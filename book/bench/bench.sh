#!/bin/bash
REPO="${IMPLICIT_MANIFOLDS_ROOT:-$(cd "$(dirname "$0")/../.." && pwd)}"
# usage: bench.sh CASE [extra cmake cxx flags] [outdir]
SRC=$REPO/cpp
TC=$1; FL=$2; O=${3:-$REPO/book/runs/bench/$TC}
B=$REPO/book/runs/bench/build/$(basename $O); mkdir -p $B $O
cmake -DCMAKE_BUILD_TYPE=Release -DTEST_CASE=$TC -DCMAKE_CXX_FLAGS="$FL" -S $SRC -B $B > $O/cmake.log 2>&1
cmake --build $B --target ms ms_j1 cs pta skeleton_pta gcch reorder_hypercubes skeleton_hypercubes fch fch_to_hypercubes -- -j 2 > $O/build.log 2>&1 || { echo "$TC build FAIL"; grep -m5 error $O/build.log; exit 1; }
cd $O
T="/usr/bin/time -v"
$T $B/ms > ms.log 2>&1
$T $B/ms_j1 > ms_j1.log 2>&1
$T $B/cs - > cs.log 2>&1
$T timeout 3600 $B/pta -o pta.bin > pta.log 2>&1
$T timeout 3600 $B/gcch -o gcch.bin > gcch.log 2>&1
$T $B/fch -o fch.bin > fch.log 2>&1
$T $B/fch_to_hypercubes -i fch.bin -o fch_hc.bin > fch_hc.log 2>&1
$T $B/reorder_hypercubes -i fch_hc.bin -o fch_hc_r.bin > reorder.log 2>&1
$T $B/skeleton_hypercubes -i fch_hc_r.bin -o fch.pol -ff 15.6 > skeleton.log 2>&1
$T $B/skeleton_pta -i pta.bin -o pta.pol -ff 15.6 > skeleton_pta.log 2>&1
rm -f fch_hc.bin gcch.bin
echo "$TC done"
