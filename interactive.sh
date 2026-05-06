#!/bin/sh
IMAGE_NAME="${IMAGE_NAME:-molecule_bond_order:latest}"

docker run \
  --rm \
  -it \
  --name molecule_bond_order \
  --user "$(id -u):$(id -g)" \
  -e HOME=/tmp \
  -v /etc/passwd:/etc/passwd:ro \
  -v /etc/group:/etc/group:ro \
  -v "$(pwd):/work" \
  -w /work \
  "${IMAGE_NAME}" \
  /bin/bash
  