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

   SPRNG_DEPEND_DIRS="src/mc src/imc src/rng"

   # directories that depend on pcglib

   PCG_DEPEND_DIRS="src/linalg src/P1Diffusion"

   # directories that depend on pooma

   POOMA_DEPEND_DIRS="src/POOMA_MT"

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

   if [ "$HAS_SPRNGLIB" != "true" ] ; then
      echo removing SPRNG dependent directories: $SPRNG_DEPEND_DIRS
      echo rm -rf $SPRNG_DEPEND_DIRS
      rm -rf $SPRNG_DEPEND_DIRS
   fi

   # remove directories that depend on pcglib if it can't be found.

   if [ "$HAS_PCGLIB" != "true" ] ; then
      echo removing PCG dependent directories: $PCG_DEPEND_DIRS
      echo rm -rf $PCG_DEPEND_DIRS
      rm -rf $PCG_DEPEND_DIRS
   fi

   # remove directories that depend on pooma if it can't be found.

   if [ "$HAS_POOMA" != "true" ] ; then
      echo removing POOMA dependent directories: $POOMA_DEPEND_DIRS
      echo rm -rf $POOMA_DEPEND_DIRS
      rm -rf $POOMA_DEPEND_DIRS
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

 UNAME=`uname`
 VENDORS=/n/radtran/vendors
 SPRNG_INCPATH=${VENDORS}/sprng/include

case $uname in
SunOS)
    PCG_LIBPATH_scalar=${VENDORS}/pcglib/${UNAME}/lib/serial
    PCG_LIBPATH_mpi=${VENDORS}/pcglib/${UNAME}/lib/mpi
    SPRNG_LIBPATH=${VENDORS}/sprng/${UNAME}

    BITS="0"
    C4="scalar mpi"
    ;;
IRIX64)
    PCG_LIB64PATH_scalar=${VENDORS}/pcglib/${UNAME}/lib64/serial
    PCG_LIB64PATH_mpi=${VENDORS}/pcglib/${UNAME}/lib64/mpi
    SPRNG_LIB64PATH=${VENDORS}/sprng/${UNAME}/lib64

    PCG_LIBN32PATH_scalar=${VENDORS}/pcglib/${UNAME}/lib32/serial
    PCG_LIBN32PATH_mpi=${VENDORS}/pcglib/${UNAME}/lib32/mpi
    SPRNG_LIBN32PATH=${VENDORS}/sprng/${UNAME}/lib32

    BITS="64 N32"
    C4="scalar mpi"
#    C4="scalar mpi shmem"
    ;;
*)
    PCG_LIBPATH_scalar=${VENDORS}/pcglib/${UNAME}/lib/serial
    PCG_LIBPATH_mpi=${VENDORS}/pcglib/${UNAME}/lib/mpi
    SPRNG_LIBPATH=${VENDORS}/sprng/${UNAME}

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

      HAS_POOMA="false"
      HAS_PCGLIB="false"
      HAS_SPRNGLIB="false"

      CONFIGUREFLAGS="--with-c4=$c4"

      if [ "X$b" = "X0" ] ; then
         eval PCG_LIBPATH='$PCG_LIBPATH_'$c4
         TARGETDIR=$TARGETROOT/$c4/draco
      else
         eval PCG_LIBPATH='$PCG_LIB'$b'PATH_'$c4
         eval SPRNG_LIBPATH='$SPRNG_LIB'$b'PATH'
         TARGETDIR=$TARGETROOT/${c4}_$b/draco
      fi

      # Check if pooma is available

      # No clause exists for checking if pooma is available.
      # This will be added when POOMA_MT becomes stable.
      #
      # HAS_POOMA="true"

      # Check if pcglib is available

      if [    -d $PCG_LIBPATH \
           -a -f $PCG_LIBPATH/libpcg_f77.a ] ; then
         HAS_PCGLIB="true"
         CONFIGUREFLAGS="--enable-pcglib --with-pcglib-lib=$PCG_LIBPATH $CONFIGUREFLAGS" 
      fi

      # Check if sprnglib is available

      if [    -d $SPRNG_INCPATH            \
           -a -d $SPRNG_LIBPATH            \
           -a -f $SPRNG_LIBPATH/libcmrg.a  \
           -a -f $SPRNG_LIBPATH/liblcg.a   \
           -a -f $SPRNG_LIBPATH/liblcg64.a \
           -a -f $SPRNG_LIBPATH/liblfg.a   ] ; then
         HAS_SPRNGLIB="true"
         CONFIGUREFLAGS="--with-sprng-lib=$SPRNG_LIBPATH --with-sprng-inc=$SPRNG_INCPATH $CONFIGUREFLAGS" 
      fi

      # turn on N32 bit compilation if $b is N32

      if [ "X$b" = "XN32" ] ; then
         CONFIGUREFLAGS="--enable-32-bit $CONFIGUREFLAGS"
      fi

      # run in subshell since cd is going on.

      (runregression)

   done

done
