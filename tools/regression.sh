#! /bin/sh

###################################################################
#
# Script to run regression tests for draco system
#
# Usage:
#    regression.sh
#
###################################################################
# $Id$
###################################################################


#
# runregression is an sh function that will be run later in this script
#

runregression ()
{
   # directories that depend on sprnglib

   SPRNG_DEPEND_DIRS="src/imc src/rng"

   # directories that depend on pcglib

   PCG_DEPEND_DIRS="src/linalg src/P1Diffusion"

   # Remove, create the target directory, and cd into it.

   echo rm -rf $TARGETDIR
   rm -rf $TARGETDIR
   echo "Creating target directory: $TARGETDIR"
   mkdir -p $TARGETDIR || exit 1
   echo cd $TARGETDIR
   cd $TARGETDIR

   # Configure in the target directory

   echo Configuring draco in `pwd`
   rm -f config.cache
   ../../../draco/draco_config .. $CONFIGUREFLAGS

   # remove directories that depend on sprnglib if it can't be found.

   if [ ! -f ${SPRNGLIB:-"No File Exists"} ] ; then
      echo removing SPRNG dependent directories: $SPRNG_DEPEND_DIRS
      echo rm -rf $SPRNG_DEPEND_DIRS
      rm -rf $SPRNG_DEPEND_DIRS
   fi

   # remove directories that depend on pcglib if it can't be found.

   if [ ! -f ${PCGLIB:-"No File Exists"} ] ; then
      echo removing PCG dependent directories: $PCG_DEPEND_DIRS
      echo rm -rf $PCG_DEPEND_DIRS
      rm -rf $PCG_DEPEND_DIRS
   fi

   # Make the test runs

   echo gmake -k check
   gmake -k check
}

#
# Beginning of script executable
#

hostname=`hostname`
uname=`uname`

# Change to the regression directory

if [ ! -d regression ] ; then
  echo "Creating regression directory"
  mkdir regression || exit 1
fi
cd regression

# CVS checkout or update draco source tree

if [ -d draco ]; then
   echo Updating draco in `pwd`/draco
   cvs -q update -APd draco
else
   echo Checking out draco in `pwd`/draco
   cvs -q checkout draco
fi

# Set up paths to look for libraries

CONTRIB=/usr/local/contrib/radtran/$uname

case $uname in
SunOS)
    PCG_LIBPATH=${CONTRIB}/lib
    SPRNG_LIBPATH=

    BITS="0"
    C4="scalar mpi"
    ;;
IRIX64)
    PCG_LIB64PATH=${CONTRIB}/lib64
    SPRNG_LIB64PATH=

    PCG_LIBN32PATH=${CONTRIB}/libn32
    SPRNG_LIBN32PATH=

    BITS="64 N32"
#    C4="scalar mpi shmem"
    C4="scalar mpi"
    ;;
*)
    PCG_LIBPATH=
    SPRNG_LIBPATH=

    BITS="0"
    C4="scalar"
    ;;
esac

# Create a serial build

# remove the toplevel target directory

TARGETROOT=$hostname
echo rm -rf $TARGETROOT
rm -rf $TARGETROOT

for c4 in $C4
do
   for b in $BITS
   do

      CONFIGUREFLAGS="--with-c4=$c4"

      if [ "X$b" = "X0" ] ; then
         TARGETDIR=$TARGETROOT/$c4/draco
      else
         eval PCG_LIBPATH='$PCG_LIB'$b'PATH'
         TARGETDIR=$TARGETROOT/${c4}_$b/draco
      fi


      # Check if pcglib is available

      if [ -d $PCG_LIBPATH -a -f $PCG_LIBPATH/libpcg_f77.a ] ; then
         PCGLIB=$PCG_LIBPATH/libpcg_f77.a
         CONFIGUREFLAGS="--enable-pcglib --with-pcglib-lib=$PCG_LIBPATH $CONFIGUREFLAGS"
      fi

      # turn on N32 bit compilation if $b is N32

      if [ "X$b" = "XN32" ] ; then
         CONFIGUREFLAGS="--enable-32-bit $CONFIGUREFLAGS"
      fi

      # run in subshell since cd is going on.

      (runregression)

   done

done
