#! /bin/sh

POOMA_ROOT=$1
POOMA_ARCH=$2

grep pooma_compile_options ${POOMA_ROOT}/lib/${POOMA_ARCH}/PoomaVersions.h |
sed 's/.*"\(.*\)";$/\1/'
