##---------------------------------------------------------------------------##
## draco/tools/configure.in
## $Id$
##---------------------------------------------------------------------------##

AC_INIT
AC_CONFIG_SRCDIR([launchtest])
AC_CONFIG_AUX_DIR(../config)
AC_DRACO_ARGS

# draco version
AC_DRACO_PREREQ

# we use the install script provided with autoconf on all machines
INSTALL='${config_dir}/install-sh -c'
INSTALL_DATA='${INSTALL} -m 444'

installdirs='${libexecdir} ${bindir}'

# required tools for draco build system testing framework
libexec_files="\${libexecdir}/launchtest"
libexec_files="\${libexecdir}/test_filter.py ${libexec_files}"
libexec_files="\${libexecdir}/apptest_modules.py ${libexec_files}"

# draco information script
bin_files='${bindir}/draco'

# useful utilities that we install
bin_files="\${bindir}/count_loc.sh ${bin_files}"
bin_files="\${bindir}/find_tags.py ${bin_files}"
bin_files="\${bindir}/find_tags_diff.py ${bin_files}"
bin_files="\${bindir}/update_release_tags.py ${bin_files}"
bin_files="\${bindir}/include_tree.py ${bin_files}"
bin_files="\${bindir}/package_depends.py ${bin_files}"
bin_files="\${bindir}/regression_filter.py ${bin_files}"
bin_files="\${bindir}/check_coverage ${bin_files}"
bin_files="\${bindir}/check_package_coverage ${bin_files}"
bin_files="\${bindir}/numdiff.py ${bin_files}"

installfiles="$libexec_files $bin_files"

configure_command="$[]0 $[]*"

AC_SUBST(configure_command)

AC_SUBST(installdirs)
AC_SUBST(installfiles)

AC_SUBST(INSTALL)
AC_SUBST(INSTALL_DATA)

AC_CONFIG_FILES([Makefile draco check_coverage check_package_coverage])
AC_OUTPUT

##---------------------------------------------------------------------------##
## end of configure.in
##---------------------------------------------------------------------------##

