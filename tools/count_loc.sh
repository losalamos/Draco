#! /bin/sh
###################################################################
#  NAME
#	count_loc.sh
#
#  SYNOPSIS
#       count_loc.sh [dir]
#
#  DESCRIPTION
#       A script to estimate lines of code.
#       The directory to be scanned is given in [dir]
#       If [dir] is absent, then the current directory is examined.
#       Attempts are made to find all files of certain types and
#       to strip out comments and blank lines, however this
#       should not be considered to be exhaustive. Ignores
#       files with names beginning with "#" or ending with "~".
#
#       Author: Mark G. Gray, Randy M. Roberts, John M. McGhee
#               Los Alamos National Laboratory
#       Date:   Wen Sep 29 10:44:35 MST 1999
#
#
###################################################################
# $Id$
###################################################################

#Count lines of code in a C++ file by counting semicolons
cpp_loc()
{
  grep ";" $* | wc -l
}

#
# awk '/<search pattern>/ {<actions>}'
#

#Count lines in text files, stripping any blank lines.
text_loc()
{
  awk '/^[ \t]*$/ {next} {print}' $* | wc -l
}

#Count lines of code in a script file, stripping any blank lines 
#and/or comments. Comment character is "#".
script_loc()
{
  awk '/^[ \t]*$/ || $1 ~ /^#/ {next} {print}' $* | wc -l
}

#Count lines of code in a html file, stripping any blank lines
#and the first line of any comments. Comment character is "<!--".
html_loc()
{
  # 1. Omit lines that start with spaces or tabs. Omit blank lines.
  # 2. OR (|| )
  # 2. omit if first record($1 ~) has an HTML comment dilemeter '<!--'

  awk 'BEGIN {skip = 0} \
    skip == 0 {if ( $1 ~ /^[ \t\n\r]*$/ )
                 {next} \
               else if ( $1 ~ /<!--/ )
                 {skip = 1} \
               else 
                 {print}; \
               next} \
    skip == 1 { if ( $1 ~ /-->/ )
                  {skip = 0};\
                next}' $* | wc -l
}

#Count lines of code in a tex file, stripping any blank lines
#and/or comments. Comment character is "%".
latex_loc()
{
  awk '/^[ \t]*$/ || $1 ~ /^%/ {next} {print}' $* | wc -l
}

#Count lines of code in a elisp file, stripping any blank lines
#and/or comments. Comment character is ";;" or ";".
elisp_loc()
{
  awk '/^[ \t]*$/ || $1 ~ /^;/ {next} {print}' $* | wc -l
}

#Count lines of code in a fortran90/95 file, stripping any blank lines
#and/or comments. Comment character is "!"
f90_loc()
{
  awk '/^[ \t]*$/ || $1 ~ /^!/ {next} {print}' $* | wc -l
}

#Count lines of code in a fortran77 file, stripping any blank lines
#and/or comments. Comment character is "C" or "c" in column 1.
f77_loc()
{
  awk '/^[ \t]*$/ || /^[Cc]/ {next} {print}' $* | wc -l
}

# Count lines of code in m4 files and configure.ac, stripping any
# blank lines and/or comments.  Comment characters are 'dnl' or '#'.
m4_loc()
{
  awk '/^[ \t]*$/ || $1 ~ /^#/ || $1 ~ /^dnl/ {next} {print}' $* | wc -l
}

#Define a fucntion which lists the names of any 
#shell script files
list_shell_scripts()
{
  for file in $*; do 
    if file $file | grep 'shell script' >/dev/null; then
      echo $file
    fi
  done
}

#Define a fucntion which  lists the names of any 
#perl script files
list_perl_scripts()
{
  for file in $*; do 
    if file $file | grep 'perl script' >/dev/null;  then
      echo $file
    fi
  done
}

#------------------------ Main Code -------------------------------------#

#Get the command line arguments (if any)
if [ $# = 0 ]; then
  TMP="."
elif [ $# = 1 ]; then
  TMP="$1"
else
  TMP="$1"
  echo "Warning -- multiple command line arguments found,"
  echo "           using first argument and ignoring others"
  echo
fi

#Check for existence of the directory to be scanned
if [ ! -d $TMP ]; then
 echo $TMP " is not a directory."
 exit 1
fi

#Get the full pathname of the directory
CODE_DIR=`(cd $TMP; pwd)`

# Print the header
echo
echo "Counting lines of source code -- " 
echo "Directory: " $CODE_DIR 
echo "Date     : " `date`
echo

# Go to the directory to be scanned.
cd $CODE_DIR

#------------------------- C++ Source Code [ .cc .hh .cpp .h ] ----------

all_cpp_fils=`find . -type f \( -name '*.cc' -o -name '*.cpp' -o -name '*.hh' -o -name '*.h' \) \
! -name '#*' -exec cat {} \;`

#Scan for C++ source code
echo -n "                 Total C++ source: "
find . -type f \
   \( -name '*.cc' -o -name '*.hh' -o -name '*.cpp' -o -name '*.h' \) \
   ! -name '#*' -exec cat {} \; |  cpp_loc 
#echo

#------------------------ C++ Test Code ----------------------------------

#Scan for C++ test code
echo -n "   C++ source in test directories: "
for i in `find . -name 'test*' -type d -print`; do
  find ${i} -type f \( -name '*.cc' -o -name '*.hh' -o -name '*.cpp' -o -name '*.h' \) ! -name '#*' \
            -exec cat {} \; 
done | cpp_loc
#echo

#Scan for C++ Design-by-Contract specifications
echo -n "      C++ contract specifications: "
for i in `find . -type f \( -name '*.cc' -o -name '*.hh' -o -name '*.cpp' -o -name '*.h' \) ! -name '#*' \
          -print`; do
#   echo ${i}
#  awk '/[Cc]heck[ ]*\(/ {print}' ${i}
#  awk '$1 ~ /[aA]ssert|[rR]equire|[eE]nsure|[Cc]heck|[iI]nsist/ {print}' ${i}
  awk '/[aA]ssert[ ]*\(|[rR]equire[ ]*\(|[eE]nsure[ ]*\(|[ ][Cc]heck[ ]*\(|[iI]nsist[ ]*\(/ {print}' ${i}
done | cpp_loc
#echo

#------------------------ C++ Documentation ----------------------------

#Scan for C++ comments (finds both C and C++ style comments)
echo -n "                     C++ comments: "
for i in `find . -type f \( -name '*.cc' -o -name '*.hh' -o -name '*.cpp' -o -name '*.h'  \) ! -name '#*' \
          -print`; do 
  awk '$1~/\/\// {print} /\/\*/, /\*\// {print}' ${i}
done | wc -l
#echo

#------------------------- f90/f95 Source Code ------------------------------

#Scan for f90/f95 source code
echo -n "             Total Fortran source: "
find . -type f \( -name '*.F' -o -name '*.f90' -o -name '*.F95' \) ! -name '#*' \
       -exec cat {} \; | f90_loc 
#echo

#------------------------ f90/f95 Documentation ----------------------------

#Scan for f90/f95 comments 
echo -n "                 Fortran comments: "
for i in `find . -type f \( -name '*.F' -o -name '*.f90' -o -name '*.f95' -o -name '*.F95' \) \
          ! -name '#*' -print`; do 
  awk '$1 ~ /!/ {print}' ${i}
done | wc -l
#echo


#------------------------ Other Documentation ----------------------------

#Scan for LaTeX source code
echo -n '              LaTeX Documentation: '
find . -type f \( -name "*.tex" -o -name '*.bib' \) ! -name '#*' -exec cat {} \; |  latex_loc
#echo

#Scan for html source
echo -n '                      HTML source: '
find . -type f \( -name '*.html' -o -name '*.htm' -o -name '*.shml'  \) \
       ! -name '#*' -exec cat {} \; |  html_loc
# echo

#Scan for plain text documents
echo -n '                        Text docs: '
find . -type f \( -name ChangeLog -o -name 'README*' -o -name '*.txt' -o -name '*.asc' \) \
  ! -name '#*' -exec cat {} \; | text_loc
#echo

#------------------------ Scripts ----------------------------------------

#Scan for Build System (M4/configure.ac/Makefile.in) scripts.
echo -n "       Build system script source: "
find . -type f ! -name '*~' ! -name '#*' \( -name '*.ac' -o -name '*.m4' -o -name Makefile -o -name '*.in' -o -name '*.cmake' -o -name CMakeLists.txt \) \
  -exec cat {} \; | m4_loc

#Scan for executable shell script source code. 
#(such as sh, bash, csh, etc.)
echo -n "   Executable shell script source: "
files=`find . -type f -perm -100  ! -name 'configure' \
      ! -name '*~' ! -name '#*' -print`
script_loc `list_shell_scripts $files` /dev/null
#echo

#Scan for Perl source code. 
echo -n "           Executable Perl source: "
files=`find . -type f -perm -100  ! -name '*~' ! -name '#*' -print`
script_loc `list_perl_scripts $files` /dev/null
#echo

#Scan for Python source code
echo -n "                    Python source: "
find . -name '*.py' ! -name '#*' -type f -exec cat {} \; | script_loc
#echo

#Scan for Expect source code
echo -n "                    Expect source: "
find . -name '*.exp' ! -name '#*' -type f -exec cat {} \; | script_loc
#echo

#------------------------ Templates --------------------------------------

#Scan for Elisp source code
echo -n "                     Elisp source: "
find . -name '*.el' ! -name '#*' -type f  -exec cat {} \; | elisp_loc
#echo

##---------------------------------------------------------------------------##

echo ""
echo "Done."


#################### end of count_loc.sh ######################################
