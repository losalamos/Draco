//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/Release.cc
 * \author Thomas Evans
 * \date   Thu Jul 15 09:31:44 1999
 * \brief  Provides the function definition for Release.
 * \note   Copyright (C) 2016 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Release.hh"
#include "ds++/config.h"
#include <sstream>
#include <cstring> // memcpy

namespace rtt_dsxx
{

//---------------------------------------------------------------------------//
// function definition for Release, define the local version number for
// this library in the form ds_#.#.# in pkg_version variable
const std::string release()
{
    std::ostringstream pkg_release;
    // Name and version
    pkg_release << "Draco-"
                << DRACO_VERSION_MAJOR << "_"
                << DRACO_VERSION_MINOR << "_"
                << DRACO_VERSION_PATCH ;

    // build date and type
    std::string const build_date( DRACO_BUILD_DATE );
    std::string const build_type( DRACO_BUILD_TYPE );
    pkg_release << ", build date " << build_date
                << "; build type: " << build_type
#ifdef DBC
                << "; DBC: " << DBC
#endif
#ifdef DRACO_DIAGNOSTICS
                << "; DRACO_DIAGNOSTICS: " << DRACO_DIAGNOSTICS
#endif
#ifdef DRACO_DIAGNOSTICS_LEVEL_3
#ifdef FPETRAP_SUPPORTED
                << "; FPE_TRAP: ON"
#endif
#endif
        ;

    return pkg_release.str();
}

//---------------------------------------------------------------------------//
/*! \brief Return a list of Draco contributing authors
 *
 * This data is assembled by hand:
 * \code
 * % files=`find . -name '*.hh' -o -name '*.cc' -o -name '*.txt' \
         -o -name '*.cmake' -o -name '*.in' -o -name '*.h'`
 * % svn annotate $files > ../file_list
 * % user_list=`cat ../file_list | awk '{print $2}' | sort -u`
 * % for name in $user_list; do numlines=`grep $name ../file_list | wc -l`;   \
        echo "$numlines: $name"; done > ../author_loc
 * % cat ../author_loc | sort -rn
 * \endcode
 *
 * Note 1: the annotate step can take a long time (do this on a local disk!)
 *
 * Note 2: I only included contributers that supplied more than 500 lines of
 * code.
 */
const std::string author_list()
{
// 77639: kellyt
// 45591: kgbudge
// 9798: tme
// 5897: bta
// 4760: lowrie
// 3954: mcghee
// 3632: mwbuksas
// 3448: gaber
// 2188: warsa
// 832: sethrj
// 799: wollaber
// 629: rsqrd
// 565: jhchang
// 434: talbotp
// 208: kwang
// 183: furnish
// 169: clevelam
// 130: bergen
// 104: phenning
// --- 100 lines ---
// 92: pahrens
// 47: batcho
// 41: jdd
// 30: tmonster
// 24: swmosher
// 8: tkelley
// 2: wawiesel

    std::stringstream alist;
//            0         1         2         3         4         5         6         7
//            01234567890123456789012345678901234567890123456789012345678901234567890123456789
    alist << "    Kelly G. Thompson, Kent G. Budge, Tom M. Evans, B. Todd Adams,\n"
          << "    Rob Lowrie, John McGhee, Mike W. Buksas, Gabriel M. Rockefeller,\n"
          << "    James S. Warsa, Seth R. Johnson, Allan B. Wollaber, Randy M. Roberts,\n"
          << "    Jae H. Chang, Paul W. Talbot, Katherine J. Wang, Jeff Furnish,\n"
          << "    Matthew A. Cleveland, Benjamin K. Bergen, Paul J. Henning.";
    return alist.str();
}

//---------------------------------------------------------------------------//
/*! \brief Print a Copyright note with an author list:
 */
const std::string copyright()
{
    std::ostringstream msg;

    msg << "Draco Contributers: \n"
        << author_list() << "\n\n"
        << "Copyright (C) 2016 Los Alamos National Security, LLC."
        << std::endl;

    return msg.str();
}

}  // end of rtt_dsxx

//---------------------------------------------------------------------------//
//! This version can be called by Fortran and wraps the C++ version.
extern "C" void ec_release( char * release_string, size_t maxlen )
{
    std::string tmp_rel = rtt_dsxx::release();
    if( tmp_rel.size() >= static_cast<size_t>(maxlen) )
    {
        tmp_rel = tmp_rel.substr(0,maxlen-1);
    }
    std::memcpy(release_string,tmp_rel.c_str(),tmp_rel.size()+1);
    return;
}

//---------------------------------------------------------------------------//
// end of Release.cc
//---------------------------------------------------------------------------//
