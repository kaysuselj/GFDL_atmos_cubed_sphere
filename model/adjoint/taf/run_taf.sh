#!/bin/bash
# Generate the TL and AD code (one TAF run, -forward -reverse) with the TAF server (staf) from the
# unmodified forward files, for the top-level top_adv (driver/adj_tops.F90).
#   taf/run_taf.sh
# The preprocessed copy of fv_mapz.F90 that goes to TAF gets the GOTO-free map1_q2 of
# taf/map1_q2_structured.F90 (see there), and ppm_limiters gets explicit dummy lengths (itot) in
# place of the assumed-size dm(*), a4(4,*) (same code otherwise); the copy is kept as taf/gen/fv_mapz_tafinput.f90 so that
# its forward can be checked against the original (test_fwd and test_adv write step_*.r8).
# Output: taf/gen/*_ad.f90 (modules with the original, TL and AD routines), taf/gen/taf.log, and
# taf/gen_gchp/*_ad_gchp.F90 (the same with renamed modules, for GCHP; tools/make_gen_gchp.py).
set -euo pipefail
STAF=${STAF:-/Users/ksuselj/TAF/staf}
D=$(cd "$(dirname "$0")" && pwd); A=$(dirname "$D"); M=$(dirname "$A")
W=$D/work; rm -rf "$W"; mkdir -p "$W"; cd "$W"
SRC="$A/stubs/stubs $M/tp_core $M/fv_fill $M/fv_mapz $M/fv_tracer2d $A/driver/adj_face_mod $A/driver/adj_tops"
FILES=""
# CPP definitions of GCHP's fvdycore build (fvdycore/CMakeLists.txt), as in the Makefile
for f in $SRC; do b=$(basename $f); gfortran -E -cpp -P -DMAPL_MODE -DSPMD -DTIMING $f.F90 > $b.f90 2>/dev/null; FILES="$FILES $b.f90"; done
gfortran -E -cpp -P "$D/map1_q2_structured.F90" > map1_q2_structured.f90
python3 -I - fv_mapz.f90 map1_q2_structured.f90 <<'PYEOF'
import re, sys
src, rep = open(sys.argv[1]).read(), open(sys.argv[2]).read()
rep = re.sub(r'(?m)^\s*!.*\n', '', rep)
pat = re.compile(r'(?ims)^[ \t]*subroutine\s+map1_q2\s*\(.*?^[ \t]*end\s+subroutine\s+map1_q2[ \t]*$')
new, nsub = pat.subn(lambda m: rep.strip('\n'), src)
assert nsub == 1, nsub
# ppm_limiters: assumed-size dummies -> explicit length itot (TAF cannot store assumed-size arrays)
new, n2 = re.subn(r'(?m)^(\s*real\s*,\s*intent\(in\)\s*::\s*dm)\(\*\)', r'\1(itot)', new)
new, n3 = re.subn(r'(?m)^(\s*real\s*,\s*intent\(inout\)\s*::\s*a4)\(4,\*\)', r'\1(4,itot)', new)
assert n2 == 1 and n3 == 1, (n2, n3)
open(sys.argv[1], 'w').write(new)
PYEOF
$STAF -forward -reverse -toplevel top_adv -input qh,qc -output qh,qc,fx,fy \
      -r8 -f08 -free -N 512 $FILES > staf.out 2>&1 || true
if grep -qE "ERROR|not present|no transformations" staf.out; then echo "TAF failed, see $W/staf.out"; exit 1; fi
if grep -E "do not depend|have no influence" taf.log; then echo "TAF dependence warnings, see $W/taf.log"; exit 1; fi
grep -E "irreducible" taf.log || true
G=$D/gen; rm -rf "$G"; mkdir -p "$G"
for f in $FILES; do cp "${f%.f90}_ad.f90" "$G/"; done
# TAF declares the adjoint of q intent(in) in start_halo_{3,4}d_ad, although the black-box adjoint
# adj_halo_q_ad overwrites it (ADKILL is honoured for adj_global_sum but not here): make it inout
python3 -I - "$G/stubs_ad.f90" <<'PYEOF'
import re, sys
t = open(sys.argv[1]).read()
def fix(m):
    return m.group(0).replace('real, intent(in) :: q_ad(', 'real, intent(inout) :: q_ad(')
t, n = re.subn(r'(?is)subroutine start_halo_[34]d_ad\(.*?end subroutine start_halo_[34]d_ad', fix, t)
assert n >= 1 and "intent(inout) :: q_ad(" in t, n
open(sys.argv[1], 'w').write(t)
PYEOF
cp taf.log "$G/taf.log"; cp fv_mapz.f90 "$G/fv_mapz_tafinput.f90"
# the back-ends (comm/*.F90) implement the black-box routines with these argument lists; the calls
# go through implicit interfaces, so check that TAF still generates them in this form
for c in "call adj_halo_q( q,help_h,help_i,help_j )" "call adj_halo_q_ad( q,q_ad,help_h,help_i,help_j )" \
         "call adj_halo_q( q_tl,help_h,help_i,help_j )" "call adj_global_sum( a_ad,n )" "call adj_global_sum( a_tl,n )"; do
    grep -qF "$c" "$G/stubs_ad.f90" || { echo "generated black-box call changed, expected: $c"; exit 1; }
done
# the copies with renamed modules that GCHP compiles (taf/gen_gchp)
python3 -I "$A/tools/make_gen_gchp.py" "$G" "$D/gen_gchp"
echo "TAF output in $G:"; wc -l "$G"/*_ad.f90 | tail -1
grep -ihE "^ *subroutine [a-z0-9_]*_(tl|ad) *\(" "$G"/*_ad.f90 | sed 's/(.*//' | sort -u | tr -s ' ' | tr '\n' ' '; echo
