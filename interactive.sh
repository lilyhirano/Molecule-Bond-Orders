#!/bin/sh
IMAGE_NAME="${IMAGE_NAME:-molecule-visualization:latest}"

docker run \
  --rm \
  -it \
  --name molecule-visualization \
  -v "$(pwd):/work" \
  -w /work \
  "${IMAGE_NAME}" \
  /bin/bash