#!/usr/bin/env bash
# ===========================================================================
# EXERCISE 22 — the CI script
# ===========================================================================
#
# TODO 4: complete the three tests below. Run it BEFORE fixing model.f90 --
# test 1 must FAIL. A test suite that passes against known-broken code is
# worse than no test suite, because it manufactures confidence.
#
# Usage:  ./run_tests.sh
# Exit:   0 if all tests pass, 1 otherwise (this is what CI keys on)
# ===========================================================================

set -uo pipefail          # deliberately NOT -e: a failing test must be
                          # reported, not abort the whole script

BIN=./model
MPIRUN="${MPIRUN:-mpirun --oversubscribe}"
STEPS=1000
RESTART_AT=617            # deliberately not a round number -- see TODO 5d

pass=0
fail=0

report() {
    local name="$1" ok="$2" detail="${3:-}"
    if [[ "$ok" == "yes" ]]; then
        printf '  [ PASS ]  %-28s %s\n' "$name" "$detail"
        pass=$((pass + 1))
    else
        printf '  [ FAIL ]  %-28s %s\n' "$name" "$detail"
        fail=$((fail + 1))
    fi
}

# ---------------------------------------------------------------------------
# Compare two state dumps.
#
# TODO 4a: implement the tolerance comparison.
#   - identical bytes                       -> bit-identical
#   - differs, but every value within `tol`  -> within tolerance
#   - otherwise                              -> FAIL, and print WHICH line
#
# The third case is the one that matters. "Files differ" sends someone
# hunting; "accumulated_flux differs by 3.2e-4" sends them straight to the
# line of code.
#
# awk is the shortest path: split each line into label and value, compare
# the values pairwise, print the label of any that exceed the tolerance.
# ---------------------------------------------------------------------------
compare_states() {
    local a="$1" b="$2" tol="${3:-0}"

    if [[ ! -f "$a" || ! -f "$b" ]]; then
        echo "missing file: $a or $b"
        return 2
    fi

    if cmp -s "$a" "$b"; then
        return 0                      # bit-identical
    fi

    if [[ "$tol" == "0" ]]; then
        # TODO 4a: print which labelled values differ, not just "they differ"
        echo "differs (bit-identity required)"
        return 1
    fi

    # TODO 4a: within-tolerance comparison. Suggested shape:
    #   paste "$a" "$b" | awk -v tol="$tol" '
    #     { lab=$1; x=$2; y=$4;
    #       d = (x-y); if (d<0) d=-d;
    #       s = (x<0?-x:x); if (s<1) s=1;
    #       if (d/s > tol) { print "    " lab " differs: " x " vs " y; bad=1 }
    #     } END { exit bad }'
    echo "tolerance comparison not implemented (TODO 4a)"
    return 1
}

echo "==========================================================="
echo " Exercise 22 regression suite"
echo "==========================================================="
echo

if [[ ! -x "$BIN" ]]; then
    echo "  $BIN not built -- run 'make' first."
    exit 1
fi

# ---------------------------------------------------------------------------
# TEST 1 — RESTART.
#
# A run stopped at step N and restarted must be bit-identical to an
# uninterrupted run. This is the test that catches the missing state
# variable in model.f90.
# ---------------------------------------------------------------------------
echo "--- test 1: restart reproducibility (4 ranks) ---"
$MPIRUN -n 4 $BIN run     $STEPS 0           cont.txt  > /dev/null 2>&1
$MPIRUN -n 4 $BIN restart $STEPS $RESTART_AT rst.txt   > /dev/null 2>&1

detail=$(compare_states cont.txt rst.txt 0) && ok=yes || ok=no
report "restart bit-identical" "$ok" "$detail"

# ---------------------------------------------------------------------------
# TEST 2 — REPRODUCIBILITY.
#
# The same run twice must give the same answer. This catches uninitialised
# memory, a dependence on wall-clock time, an unseeded RNG, and races.
#
# TODO 4b: run the model twice at the same rank count and compare with
# tolerance 0.
# ---------------------------------------------------------------------------
echo
echo "--- test 2: run-to-run reproducibility (4 ranks) ---"
# TODO 4b: implement. Placeholder below always fails so it cannot be
# mistaken for a passing test.
report "run-to-run identical" "no" "not implemented (TODO 4b)"

# ---------------------------------------------------------------------------
# TEST 3 — RANK INVARIANCE.
#
# The same run on 2 vs 4 vs 8 ranks. This can NOT be bit-identical while the
# global sum goes through a plain MPI_Allreduce -- see TODO 5b. So it is a
# TOLERANCE test, and choosing that tolerance honestly is the exercise.
#
# TODO 4c: run at 2, 4 and 8 ranks and compare each against the 4-rank
# reference with a justified tolerance.
# ---------------------------------------------------------------------------
echo
echo "--- test 3: rank-count invariance (2 vs 4 vs 8) ---"
# TODO 4c: implement.
report "rank invariance" "no" "not implemented (TODO 4c)"

# ---------------------------------------------------------------------------
echo
echo "==========================================================="
printf ' %d passed, %d failed\n' "$pass" "$fail"
echo "==========================================================="

if [[ $fail -gt 0 ]]; then
    echo
    echo " Expected on first run: test 1 FAILS because model.f90 has a"
    echo " genuine restart bug. Build the diagnostic (TODO 2), find the"
    echo " missing variable, fix TODO 1, and watch this go green."
    exit 1
fi
exit 0
