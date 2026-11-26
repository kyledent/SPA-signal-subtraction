#!/usr/bin/env bash
# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Diamond Light Source
#
# fetch.sh -- download the EMD-32299 inputs (primary map, both half maps, the deposited mask, the
# header XML) and the fitted model PDB 7W4A into this folder. About 320 MB uncompressed.
#   bash fetch.sh
set -euo pipefail
cd "$(dirname "$0")"
B=https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-32299
curl -sSfO "$B/header/emd-32299.xml"
curl -sSf "$B/map/emd_32299.map.gz"              | gunzip > emd_32299.map
curl -sSf "$B/other/emd_32299_half_map_1.map.gz" | gunzip > emd_32299_half_map_1.map
curl -sSf "$B/other/emd_32299_half_map_2.map.gz" | gunzip > emd_32299_half_map_2.map
curl -sSfO "$B/masks/emd_32299_msk_1.map"
curl -sSf -o 7w4a.cif https://files.rcsb.org/download/7W4A.cif
ls -la emd_32299*.map 7w4a.cif
