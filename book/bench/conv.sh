#!/bin/bash
REPO="${IMPLICIT_MANIFOLDS_ROOT:-$(cd "$(dirname "$0")/../.." && pwd)}"
SRC=$REPO/cpp
for TC in BOOK_CONV_S2 BOOK_CONV_T2; do
 for N in 8 16 32 64; do
  [ $TC = BOOK_CONV_T2 ] && [ $N = 64 ] && continue
  O=$REPO/book/runs/bench/conv/${TC}_$N; B=$REPO/book/runs/bench/build/conv_${TC}_$N; mkdir -p $O $B
  cmake -DCMAKE_BUILD_TYPE=Release -DTEST_CASE=$TC -DCMAKE_CXX_FLAGS="-DCONV_N=$N" -S $SRC -B $B > $O/cmake.log 2>&1
  cmake --build $B --target ms fch fch_to_hypercubes reorder_hypercubes skeleton_hypercubes -- -j 2 > $O/build.log 2>&1 || { echo FAIL $TC $N; continue; }
  cd $O
  $B/ms ms.pol > ms.log 2>&1
  $B/fch -o fch.bin > fch.log 2>&1
  $B/fch_to_hypercubes -i fch.bin -o hc.bin > /dev/null 2>&1
  $B/reorder_hypercubes -i hc.bin -o hcr.bin > /dev/null 2>&1
  $B/skeleton_hypercubes -i hcr.bin -o fch.pol -ff 15.6 > /dev/null 2>&1
  rm -f hc.bin hcr.bin
  echo "$TC $N done"
 done
done
