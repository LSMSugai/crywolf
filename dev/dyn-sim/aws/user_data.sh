#!/bin/bash
# EC2 user data (Amazon Linux 2023, x86_64). Placeholders __BUCKET__, __RUN__,
# __SCRIPT__, __JOBS__, __WORKERS__, __REPS__ and __MAX_HOURS__ are filled in
# by aws/launch_spot.sh.
# Runs dev/dyn-sim/run_sims.R inside the Docker image (adapted from fpcalib), resuming from the job files
# saved in S3; checkpoints the job folder every 10 min and on a spot
# interruption notice, copies the outputs at the end and terminates the
# instance.
set -uxo pipefail
BUCKET=__BUCKET__
RUN=__RUN__
SCRIPT=__SCRIPT__
JOBS=__JOBS__
WORKERS=__WORKERS__
REPS=__REPS__
MAX_HOURS=__MAX_HOURS__
# watchdog: terminate after MAX_HOURS whatever happens (shutdown = terminate)
shutdown -h +$((MAX_HOURS * 60))
WORK=/opt/crywolf
S3=s3://${BUCKET}/fpcalib/runs/${RUN}
TOKEN=$(curl -s -X PUT http://169.254.169.254/latest/api/token \
  -H "X-aws-ec2-metadata-token-ttl-seconds: 21600")
ID=$(curl -s -H "X-aws-ec2-metadata-token: ${TOKEN}" \
  http://169.254.169.254/latest/meta-data/instance-id)
LOG=${WORK}/run.log

dnf install -y docker
systemctl start docker
mkdir -p ${WORK}
aws s3 cp ${S3}/code.tar.gz /tmp/crywolf.tar.gz
tar xzf /tmp/crywolf.tar.gz -C ${WORK}
mkdir -p ${WORK}/${JOBS}
aws s3 sync ${S3}/jobs ${WORK}/${JOBS} --only-show-errors   # resume
touch /tmp/start_marker

checkpoint() {
  aws s3 sync ${WORK}/${JOBS} ${S3}/jobs --only-show-errors
  aws s3 cp ${LOG} ${S3}/logs/run_${ID}.log --only-show-errors || true
}

# periodic checkpoint
( while true; do sleep 600; checkpoint; done ) &

# spot interruption notice (2-minute warning): checkpoint once and stop polling
(
  while true; do
    T=$(curl -s -X PUT http://169.254.169.254/latest/api/token \
      -H "X-aws-ec2-metadata-token-ttl-seconds: 300")
    CODE=$(curl -s -o /dev/null -w "%{http_code}" -H "X-aws-ec2-metadata-token: ${T}" \
      http://169.254.169.254/latest/meta-data/spot/instance-action)
    if [ "${CODE}" = "200" ]; then echo "spot interruption notice $(date)" >> ${LOG}; checkpoint; break; fi
    sleep 5
  done
) &

cd ${WORK}
echo "instance ${ID} $(date); ${SCRIPT} WORKERS=${WORKERS} REPS=${REPS:-default}; $(ls ${JOBS} | wc -l) job files resumed" > ${LOG}
REPS_ENV=""
[ -n "${REPS}" ] && REPS_ENV="-e REPS=${REPS}"
# a job file cut short by an interruption cannot be read: drop it so the task
# runs again. R allows 128 connections by default and each future worker
# takes one, so raise the limit for the large instances.
docker build -f dev/dyn-sim/aws/Dockerfile -t crywolf-dyn . > ${WORK}/build.log 2>&1 && \
docker run --rm -v ${WORK}:/project -e WORKERS=${WORKERS} ${REPS_ENV} crywolf-dyn \
  bash -c "Rscript -e 'for (f in list.files(\"${JOBS}\", full.names = TRUE)) if (inherits(try(readRDS(f), silent = TRUE), \"try-error\")) { message(\"removed unreadable \", f); file.remove(f) }' && \
           Rscript --max-connections=512 ${SCRIPT}" >> ${LOG} 2>&1
STATUS=$?
echo "exit status ${STATUS} at $(date)" >> ${LOG}
cat ${WORK}/build.log >> ${LOG}
checkpoint
if [ "${STATUS}" = "0" ]; then MARK=DONE; else MARK=FAILED; fi
echo "${ID} exit ${STATUS} $(date)" > /tmp/${MARK}
aws s3 cp /tmp/${MARK} ${S3}/${MARK} --only-show-errors
aws s3 cp ${LOG} ${S3}/logs/run_${ID}.log --only-show-errors
shutdown -h now
