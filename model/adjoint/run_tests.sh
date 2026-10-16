#!/bin/bash
# Run the single-face forward check (test_fwd, hord 12 / kord 8 vs GCHP) and the TL/AD tests
# (test_adv, default hord 2 / kord 17) for a list of trajectory cases.
#   run_tests.sh <devkit traj dir> <work dir> "case face step" ...
#   e.g. run_tests.sh ../../../adv_adjoint_devkit_C24/traj ../../../work "8ms 4 1" "tiltedrot 3 10"
set -uo pipefail
A=$(cd "$(dirname "$0")" && pwd)
TRAJ=$1; WORK=$2; shift 2
for c in "$@"; do
    set -- $c
    o=$WORK/${1}_f${2}_s$(printf %02d $3)
    python3 -I "$A/prep/make_face_inputs.py" --traj "$TRAJ/$1" --face $2 --step $3 --out "$o" > "$o.prep.txt" 2>&1 || { echo "$c: prep failed"; cat "$o.prep.txt"; continue; }
    f=$("$A/build/test_fwd" "$o" | grep "k < npz" | awk '{print $NF}')
    "$A/build_ad/test_adv" "$o" > "$o/test_adv_h2_k17.txt" 2>&1
    r=$(tail -1 "$o/test_adv_h2_k17.txt")
    w=$(grep -E "^  dot," "$o/test_adv_h2_k17.txt" | awk '{v=$(NF-1)+0; if (v>m) m=v} END{printf "%.1e", m}')
    printf '%-10s face %d step %2d:  fwd vs GCHP %s   TL/AD: %s (worst dot-test rel. diff %s)\n' $1 $2 $3 "$f" "$r" "$w"
done
