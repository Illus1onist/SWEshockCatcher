#!/bin/bash
if [ -z "${WM_PROJECT_DIR:-}" ]; then
    source /usr/lib/openfoam/openfoam2406/etc/bashrc
fi
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR" || exit 1
mkdir -p results

for N in 128 256 512 1024 2048; do
    sed -E -i "s/hex \(0 1 2 3 4 5 6 7\) \([0-9]+ 1 1\)/hex (0 1 2 3 4 5 6 7) (${N} 1 1)/" system/blockMeshDict
    foamListTimes -rm -noZero >/dev/null 2>&1 || true
    cp 0/h.uniform 0/h
    cp 0/hU.uniform 0/hU
    blockMesh > "log.blockMesh.N${N}" 2>&1
    setFields > "log.setFields.N${N}" 2>&1
    SWEshockFoam > "log.SWEshockFoam.N${N}" 2>&1
    cp 1/h  "results/h_N${N}"
    cp 1/hU "results/hU_N${N}"
    echo "N=${N} done: $(tail -2 "log.SWEshockFoam.N${N}" | head -1)"
done
