//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/Shared_Lib.cc
 * \author Rob Lowrie
 * \date   Thu Apr 15 20:44:39 2004
 * \brief  Implementation of Shared_Lib.
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Shared_Lib.hh"
#include <sstream>

#include "dlfcn_support.hh"

namespace rtt_shared_lib
{

//---------------------------------------------------------------------------//
/*!
  \brief Default constructor.

  \param file_name File name of the shared lib.  If empty, open() must be
  called later in order to use a shared lib.
*/
Shared_Lib::Shared_Lib(const std::string &file_name)
    : d_handle(0)
    , d_file_name(file_name)
{
    // is_supported must be checked for all constructors.
    Insist(is_supported(), "Shared_Lib unsupported on this platform!");
    
    if ( not file_name.empty() )
    {
	open(file_name);
    }
}

//---------------------------------------------------------------------------//
/*!
  \brief Copy constructor.

  This is implemented by opening a new handle to the shared file.
*/
Shared_Lib::Shared_Lib(const Shared_Lib &from)
    : d_handle(0)
{
    open(from.d_file_name);
}

//---------------------------------------------------------------------------//
/*!
  \brief Assignment.

  This is implemented by opening a new handle to the shared file.
*/
Shared_Lib &Shared_Lib::operator=(const Shared_Lib &rhs)
{
    if ( this == &rhs )
    {
	return *this;
    }

    open(rhs.d_file_name);

    return *this;
}

//---------------------------------------------------------------------------//
/*!
  \brief Closes the shared library, if it is open.
*/
void Shared_Lib::close()
{
    if ( is_open() )
    {
	dlclose(d_handle);
    }
}

//---------------------------------------------------------------------------//
/*!
  \brief Opens a shared library.

  If a shared library is already open, that library is closed.

  \param file_name The name of the shared lib.
*/
void Shared_Lib::open(const std::string &file_name)
{
    Require(not file_name.empty());
    
    close();
    d_handle = dlopen(file_name.c_str(), RTLD_LAZY);
    d_file_name = file_name;

    if ( not is_open() )
    {
	std::ostringstream m;
	m << "Shared_Lib::open(): Error opening shared file " << file_name;
	Insist(0, m.str());
    }
}

//---------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*!
  \brief Does the dlsym() with error checking.

  The primary point of this function is to keep the dlfcn.h functions called
  within this implementation file, while this function can be called from
  within the header (i.e., within get_function()).

  \param name The name of function to load from the library.
*/
void *Shared_Lib::do_dlsym(const std::string &name)
{
    Require(is_open());
    Require(not name.empty());

    void *f = dlsym(d_handle, name.c_str());

    char *error_msg = dlerror();
    Insist(not error_msg, error_msg);

    return f;
}

} // end namespace rtt_shared_lib

//---------------------------------------------------------------------------//
//                 end of Shared_Lib.cc
//---------------------------------------------------------------------------//
