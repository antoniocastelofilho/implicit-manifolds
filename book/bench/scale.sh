#!/bin/bash
REPO="${IMPLICIT_MANIFOLDS_ROOT:-$(cd "$(dirname "$0")/../.." && pwd)}"
SRC=$REPO/cpp
for TC in BOOK_CONV_S2 BOOK_CONV_T2; do
 if [ $TC = BOOK_CONV_S2 ]; then NS="16 32 64 128"; else NS="8 16 24 32 48"; fi
 for N in $NS; do
  O=$REPO/book/runs/bench/scale/${TC}_$N; B=$REPO/book/runs/bench/build/scale_${TC}_$N; mkdir -p $O $B
  cmake -DCMAKE_BUILD_TYPE=Release -DTEST_CASE=$TC -DCMAKE_CXX_FLAGS="-DCONV_N=$N" -S $SRC -B $B > $O/cmake.log 2>&1
  cmake --build $B --target ms cs pta fch -- -j 2 > $O/build.log 2>&1 || { echo FAIL $TC $N; continue; }
  cd $O
  T="/usr/bin/time -v"
  $T $B/ms > ms.log 2>&1
  $T $B/cs - > cs.log 2>&1
  $T $B/pta -o pta.bin > pta.log 2>&1
  $T $B/fch -o fch.bin > fch.log 2>&1
  rm -f pta.bin fch.bin
  echo "$TC $N done"
 done
done
echo ALLDONE
