#!/bin/sh
IMAGE_NAME="${IMAGE_NAME:-molecule-visualization:latest}"

docker run \
  --rm \
  -it \
  --name molecule-visualization \
  --user "$(id -u):$(id -g)" \
  -e HOME=/tmp \
  -v /etc/passwd:/etc/passwd:ro \
  -v /etc/group:/etc/group:ro \
  -v "$(pwd):/work" \
  -w /work \
  "${IMAGE_NAME}" \
  /bin/bash
  