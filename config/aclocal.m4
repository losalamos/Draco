dnl-------------------------------------------------------------------------dnl
dnl aclocal.m4
dnl MACROS FOR DRACO CONFIGURATION
dnl
dnl $Id$
dnl-------------------------------------------------------------------------dnl

###-------------------------------------------------------------------------###
### include the service macros used in aclocal.m4
###-------------------------------------------------------------------------###

builtin(include,ac_local.m4)dnl internal macros used in *.m4 config files
builtin(include,ac_compiler.m4)dnl C++ compilers
builtin(include, ac_f90env.m4)dnl Fortran 90 macros
builtin(include, ac_gm4.m4)dnl GNU m4 macros

###-------------------------------------------------------------------------###
### include the service macros used in configure.in scripts
###-------------------------------------------------------------------------###

builtin(include,ac_conf.m4)dnl

###-------------------------------------------------------------------------###
### include the vendor macros
###-------------------------------------------------------------------------###

builtin(include,ac_vendors.m4)dnl

###-------------------------------------------------------------------------###
### include the DRACO argument macro
###-------------------------------------------------------------------------###

builtin(include,ac_dracoarg.m4)dnl

###-------------------------------------------------------------------------###
### include the DRACO environment macro
###-------------------------------------------------------------------------###

builtin(include,ac_dracoenv.m4)dnl

###-------------------------------------------------------------------------###
### include any package-dependent macro definitions that are used 
### outside of draco, by default these must be located in the
### directory/file:
###	       <pkg>/pkg_config/ac_pkg.m4
### where <pkg> is the name of the package (solon, milagro, etc.),
### we use sinclude so that the file does not have to exist
###-------------------------------------------------------------------------###

builtin(sinclude,../pkg_config/ac_pkg.m4)dnl

###-------------------------------------------------------------------------###
### end of aclocal.m4
###-------------------------------------------------------------------------###
