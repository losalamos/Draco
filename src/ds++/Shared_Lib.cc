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
#include "Assert.hh"
#include <string>

// Use ::dlopen, ::dlclose, ::dlerror, ::dlsym
#include <dlfcn.h>

namespace rtt_dsxx
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
    : d_file_name(from.d_file_name)
{
    open(d_file_name);
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

    d_file_name = rhs.d_file_name;
    open(d_file_name);

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

    Insist(d_handle, dlerror());
}


} // end namespace rtt_dsxx

//---------------------------------------------------------------------------//
//                 end of Shared_Lib.cc
//---------------------------------------------------------------------------//
