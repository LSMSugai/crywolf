#!/bin/bash
# (Adapted from fpcalib/aws.) Keep a spot run alive: until s3://<bucket>/fpcalib/runs/<run>/ has a DONE or
# FAILED marker, relaunch the run (resuming from its job files) whenever its
# instance ends, trying the instance types given in order. Prints one line per event.
# Usage: dev/dyn-sim/aws/supervise.sh <bucket> <run> <script> <workers> <instance_profile> <max_hours> <instance_id> <type> [<type> ...]
#   e.g. dev/dyn-sim/aws/supervise.sh my-bucket dyn dev/dyn-sim/run_sims.R 190 pcps-sim-ec2 12 i-0abc c7a.48xlarge c7i.48xlarge m7a.48xlarge
# MAXR (environment, default 5): number of relaunches before giving up.
set -uo pipefail
BUCKET=$1; RUN=$2; SCRIPT=$3; WORKERS=$4; ROLE=$5; HOURS=$6; ID=$7; shift 7; TYPES=("$@")
MAXR=${MAXR:-5}
cd "$(git rev-parse --show-toplevel)"
S3=s3://${BUCKET}/fpcalib/runs/${RUN}
marker() { aws s3 ls "${S3}/$1" >/dev/null 2>&1; }
state() { aws ec2 describe-instances --instance-ids "$1" \
  --query 'Reservations[].Instances[].State.Name' --output text 2>/dev/null || echo gone; }
n=0
while true; do
  if marker DONE; then echo "${RUN}: finished ($(date))"; exit 0; fi
  if marker FAILED; then echo "${RUN}: script failed, see ${S3}/logs/ ($(date))"; exit 1; fi
  s=$(state "$ID")
  case "$s" in
    running|pending) sleep 300; continue;;
  esac
  [ "$n" -ge "$MAXR" ] && { echo "${RUN}: gave up after ${n} relaunches ($(date))"; exit 1; }
  n=$((n + 1))
  for T in "${TYPES[@]}"; do
    out=$(dev/dyn-sim/aws/launch_spot.sh "$BUCKET" "$RUN" "$SCRIPT" "$T" "$WORKERS" "$ROLE" "$HOURS" 2>&1 | tail -1)
    case "$out" in i-*) ID=$out; echo "${RUN}: instance ended (${s}); relaunch ${n} on ${T} ${ID} ($(date))"; break;; esac
  done
  sleep 300
done
