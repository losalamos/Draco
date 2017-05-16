//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/Release.hh
 * \author Thomas Evans
 * \date   Thu Jul 15 09:31:44 1999
 * \brief  Header file for ds++ library release function.
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#ifndef rtt_ds_Release_hh
#define rtt_ds_Release_hh

#include "ds++/config.h"
#include <algorithm>
#include <functional>
#include <map>
#include <string>

namespace rtt_dsxx {

// Typedefs
typedef std::multimap<size_t, std::string, std::greater<int>> mmdevs;
typedef std::pair<size_t, std::string> fomdev;

//! Query package for the release number.
DLL_PUBLIC_dsxx const std::string release();
//! Return a list of Draco authors
DLL_PUBLIC_dsxx const std::string author_list();
//! Return a list of Draco authors
DLL_PUBLIC_dsxx const std::string copyright();

//---------------------------------------------------------------------------//
/*!
 * \brief Format list of authors to do correct line breaks and ensures total
 *        line length is less than a specified maximum.
 *
 * \arg[in] maxlen Maximum line length
 * \arg[in] line_name Category title
 * \arg[in] list of developers
 * \return A formatted message.
 */
DLL_PUBLIC_dsxx std::string print_devs(size_t const maxlinelen,
                                       std::string const &line_name,
                                       mmdevs const &devs);

} // end of rtt_ds++

//! This version can be called by Fortran and wraps the C++ version.
extern "C" DLL_PUBLIC_dsxx void ec_release(char *release_string, size_t maxlen);

#endif // rtt_ds_Release_hh

//---------------------------------------------------------------------------//
// end of ds++/Release.hh
//---------------------------------------------------------------------------//
