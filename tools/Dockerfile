#FROM kinetictheory/draco-ci-2019june

# Use ubuntu if building from scratch
FROM ubuntu:latest

# This image:
# 1. cd /D f:\work\docker (copy Dockerfile and packages.yaml to this location).
# 2. docker login -u kinetictheory (and password) # ref https://docs.docker.com/get-started/part2/
# 3. docker build --rm --pull --tag draco-ci-2019june:latest . 
#    OR: docker run -it kinetictheory/draco-ci-2019june /bin/bash -l
#        apt-get install -y --no-install-recommends [ghostview]
#        exit
# 4. docker image ls -a ==> find container name (or docker ps)
# 5. docker commit -m "added sphinx and mscgen" -a kinetictheory <hash> kinetictheory/draco-ci-2019june:latest # queues for upload
# 6. docker push kinetictheory/draco-ci-2019june:latest
# 7. docker system prune -a (remove old dangling data)

MAINTAINER KineticTheory "https://github.com/KineticTheory"

# See draco/.travis-run-tests.sh for some instructions.

## Environment needed by 'docker build' ----------------------------------------

## for apt to be noninteractive
## https://stackoverflow.com/questions/8671308/non-interactive-method-for-dpkg-reconfigure-tzdata
## https://spack.readthedocs.io/en/latest/workflows.html?highlight=docker
ENV DEBIAN_FRONTEND noninteractive
ENV DEBCONF_NONINTERACTIVE_SEEN true
ENV SPACK_ROOT=/vendors/spack
ENV DRACO_TPL="cmake@3.14.4 gsl@2.5 numdiff@5.9.0 random123@1.09 openmpi@3.1.4 netlib-lapack@3.8.0 metis@5.1.0 parmetis@4.0.3 superlu-dist trilinos mscgen"
ENV FORCE_UNSAFE_CONFIGURE=1
ENV DISTRO=bionic
ENV CLANG_FORMAT_VER=6.0
ENV APT_KEY_DONT_WARN_ON_DANGEROUS_USAGE=1

## Update packages on the raw Ubuntu image -------------------------------------------

RUN sed -i '/DISTRO/d' /etc/apt/sources.list
# try to eliminate warning about "mesg: ttyname failed: Inappropriate ioctl for device"
RUN sed -i 's/mesg/tty -s \&\& mesg/' /root/.profile
RUN cat /root/.profile

## preesed tzdata, update package index, upgrade packages and install needed software
RUN echo "tzdata tzdata/Areas select US" > /tmp/preseed.txt; \
    echo "tzdata tzdata/Zones/US select Mountain" >> /tmp/preseed.txt; \
    debconf-set-selections /tmp/preseed.txt && \
    apt-get update && \
    apt-get install -y --no-install-recommends tzdata

## Basic admin tools
RUN apt-get install -y --no-install-recommends apt-utils automake autoconf autotools-dev python3 software-properties-common flex bison ssh

## Basic developer tools
RUN apt-get install -y --no-install-recommends build-essential ca-certificates coreutils curl doxygen environment-modules gcc-8 g++-8 gfortran-8 ghostscript git grace graphviz python3-pip python3-sphinx python3-sphinx-rtd-theme tar tcl texlive tk unzip wget
# RUN apg-get upgrade
RUN if ! test -f /etc/profile.d/modules.sh; then \
      echo "source /usr/share/modules/init/bash" > /etc/profile.d/modules.sh; \
    fi

## LLVM tools like clang-format
## Note: we can't use variables in the add-apt-repository commmand as this
##       creates an invalid entry in /etc/apt/sources.list
#RUN wget -O - https://apt.llvm.org/llvm-snapshot.gpg.key | apt-key add - ; \
#    add-apt-repository 'deb http://apt.llvm.org/bionic/ llvm-toolchain-bionic-6.0 main' ; \
#    apt-get update; \
#    apt-get install -y clang-format-${CLANG_FORMAT_VER}; \
#    export PATH=$PATH:/usr/bin
RUN  apt-get install -y --no-install-recommends clang-format-${CLANG_FORMAT_VER}; \
     ln -s /usr/bin/clang-format-6.0 /usr/bin/clang-format; \
     export PATH=$PATH:/usr/bin

## code cov plugin...
RUN python3 -m pip install --upgrade pip
RUN python3 -m pip install codecov

## SPACK -----------------------------------------------------------------------------

# install/setup spack
RUN mkdir -p $SPACK_ROOT/etc/spack/linux
# Only download spack if it doesn't already exist.
RUN if ! test -d $SPACK_ROOT/opt/spack ; then \
      curl -s -L https://api.github.com/repos/spack/spack/tarball | tar xzC $SPACK_ROOT --strip 1; \
    fi
# note: if you wish to change default settings, add files alongside
#       the Dockerfile with your desired settings. Then uncomment this line
COPY packages.yaml $SPACK_ROOT/etc/spack/linux

# metis/parmetis downloads are broken right now, use a mirror.
COPY mirrors.yaml $SPACK_ROOT/etc/spack
RUN mkdir -p $SPACK_ROOT/spack.mirror/metis
RUN mkdir -p $SPACK_ROOT/spack.mirror/parmetis
COPY metis-5.1.0.tar.gz $SPACK_ROOT/spack.mirror/metis
COPY parmetis-4.0.3.tar.gz $SPACK_ROOT/spack.mirror/parmetis

RUN if ! test -f /etc/profile.d/spack.sh; then \
      echo "source $SPACK_ROOT/share/spack/setup-env.sh" > /etc/profile.d/spack.sh; \
    fi

## Provide some TPLs
RUN export PATH=$SPACK_ROOT/bin:$PATH && spack install -n ${DRACO_TPL} && spack clean -a

# image run hook: the -l will make sure /etc/profile.d/*.sh environments are loaded
CMD /bin/bash -l
