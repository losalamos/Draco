#! /bin/sh

POOMA_ROOT=$1
POOMA_ARCH=$2

VERSIONFILE=${POOMA_ROOT}/lib/${POOMA_ARCH}/PoomaVersions.h
if [ ! -r ${VERSIONFILE} ]
then
   echo "Cannot find ${VERSIONFILE}" 1>&2
   exit 1
fi

raw_options=`grep pooma_compile_options \
    ${POOMA_ROOT}/lib/${POOMA_ARCH}/PoomaVersions.h |
    sed 's/.*"\(.*\)";$/\1/'`

# echo $raw_options

for option in $raw_options
do
    case $option in
    -D*=*) define=`echo $option | sed 's/-D\(.*\)=\(.*\)/#define \1 \2/'`
	;;
    -D*) define=`echo $option | sed 's/-D\(.*\)/#define \1 1/'`
	;;
    esac
    echo $define
done
