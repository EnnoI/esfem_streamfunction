#!/usr/bin/env bash

DUNE_DIR="/home/eioer/DUNE_DIR/2.10"
DUNE_CONTROL="${DUNE_DIR}/dune-common/bin/dunecontrol"

DUNE_CONTROL_PATH=${DUNE_DIR}:. ${DUNE_CONTROL} --opts=${DUNE_DIR}/release.opts --only=esfem_stream all
