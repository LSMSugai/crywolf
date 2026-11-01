#!/bin/bash
# Upload the committed code (git HEAD) to S3 for one run (from fpcalib/aws).
# Usage: dev/dyn-sim/aws/pack_code.sh <bucket> <run>
set -euo pipefail
BUCKET=$1; RUN=$2
cd "$(git rev-parse --show-toplevel)"
if ! git diff --quiet HEAD; then echo "Uncommitted changes; commit first." >&2; exit 1; fi
TAR=$(mktemp).tar.gz
git archive --format=tar.gz -o "${TAR}" HEAD
aws s3 cp "${TAR}" "s3://${BUCKET}/fpcalib/runs/${RUN}/code.tar.gz" --only-show-errors
rm -f "${TAR}"
echo "Uploaded $(git rev-parse --short HEAD) to s3://${BUCKET}/fpcalib/runs/${RUN}/"
