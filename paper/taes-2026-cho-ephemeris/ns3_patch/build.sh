#!/bin/bash
# Standalone build of the paper's Nov 2026 binary; touches nothing under ns-3-dev.
# Flags: the cmake-cache/compile_commands.json entry of the example. Links the frozen
# optimized libraries in ../bin/lib (manifest ../bin/lib_manifest.sha256).
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
B=/home/uzair/6g_ntn_ns3/ns-3-dev/build
OUT="$HERE/../bin"
python3 "$HERE/apply_edits.py"
(cd "$HERE" && diff -u orig/ntn-cho-tte-campaign.cc ntn-cho-tte-campaign.cc > ntn-cho-tte-campaign_nov.patch || true)
FLAGS="-DEIGEN_MPL2_ONLY -DHAVE_BOOST -DHAVE_BOOST_UNITS -DHAVE_EIGEN3 -DHAVE_GSL -DHAVE_LIBXML2 -DHAVE_PACKET_H -DHAVE_SQLITE3 -DHAVE_STDINT_H -DNS3_BUILD_PROFILE_OPTIMIZED -DPROJECT_SOURCE_PATH=\"/home/uzair/6g_ntn_ns3/ns-3-dev\" -D__LINUX__ -I$B/include -I/usr/include/eigen3 -I/usr/include/libxml2 -O3 -DNDEBUG -std=c++20 -fPIE -fno-semantic-interposition -march=native -mtune=native"
LIBS="$OUT/lib/libns3.43-ntn-cho-optimized.so $OUT/lib/libns3.43-ntn-constellation-optimized.so $OUT/lib/libns3.43-satellite-optimized.so $OUT/lib/libns3.43-propagation-optimized.so $OUT/lib/libns3.43-mobility-optimized.so $OUT/lib/libns3.43-network-optimized.so $OUT/lib/libns3.43-core-optimized.so"
for v in orig patched; do
  SRC="$HERE/ntn-cho-tte-campaign.cc"; NAME=p2-campaign-nov
  [ $v = orig ] && SRC="$HERE/orig/ntn-cho-tte-campaign.cc" && NAME=p2-campaign-orig
  nice -n 19 g++ $FLAGS -c "$SRC" -o "$OUT/$NAME.o"
  nice -n 19 g++ -o "$OUT/$NAME" "$OUT/$NAME.o" $LIBS -lgsl -lgslcblas -lsqlite3 -Wl,-rpath,$OUT/lib -Wl,-rpath-link,$OUT/lib
  rm -f "$OUT/$NAME.o"
  sha256sum "$OUT/$NAME"
done
