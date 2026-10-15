#!/bin/bash
REPO="${IMPLICIT_MANIFOLDS_ROOT:-$(cd "$(dirname "$0")/../.." && pwd)}"
SRC=$REPO/cpp
for TC in "$@"; do
  B=$REPO/book/runs/modeling/build/$TC; mkdir -p $B
  cmake -DCMAKE_BUILD_TYPE=Release -DTEST_CASE=$TC -S $SRC -B $B > $B.cmake.log 2>&1
  cmake --build $B --target ms fch fch_to_hypercubes reorder_hypercubes skeleton_hypercubes -- -j 2 > $B.build.log 2>&1 && echo "$TC built" || { echo "$TC FAIL"; grep -m5 error $B.build.log; continue; }
  O=$REPO/book/runs/modeling/$TC; mkdir -p $O; cd $O
  ( /usr/bin/time -v $B/ms ms.pol ) > ms.log 2>&1
  ( /usr/bin/time -v timeout 3000 $B/fch -o fch.bin ) > fch.log 2>&1
  ( /usr/bin/time -v $B/fch_to_hypercubes -i fch.bin -o fch_hypercubes.bin ) > fch_to_hc.log 2>&1
  ( /usr/bin/time -v $B/reorder_hypercubes -i fch_hypercubes.bin -o fch_hypercubes_reordered.bin ) > reorder.log 2>&1
  ( /usr/bin/time -v $B/skeleton_hypercubes -i fch_hypercubes_reordered.bin -o fch.pol -ff 15.6 ) > skeleton.log 2>&1
  rm -f fch_hypercubes.bin
  echo "$TC done"; head -3 ms.log; grep -m2 "cells generated\|Time to" fch.log
done
