//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/Shared_Lib.hh
 * \author Rob Lowrie
 * \date   Thu Apr 15 20:44:39 2004
 * \brief  Header file for Shared_Lib.
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_ds_Shared_Lib_hh
#define rtt_ds_Shared_Lib_hh

#include <string>
#include "Assert.hh"

// use ::dlsym and ::dlerror
#include <dlfcn.h>

namespace rtt_dsxx
{

//===========================================================================//
/*!\class Shared_Lib
 * \brief Controls access to a shared (dynamically linked) library.
 *
 * Access to functions defined in the shared library is provided via the
 * get_function() member.
 */
/*! 
 * \example ds++/test/tstShared_Lib.cc 
 */
//===========================================================================//

class Shared_Lib 
{
    // DATA

    // The handle to the shared library.
    void *d_handle;

    // The name of the shared library.
    std::string d_file_name;
    
  public:

    // Default constructor.
    explicit Shared_Lib(const std::string &file_name = "");

    // Copy constructor.
    Shared_Lib(const Shared_Lib &from);

    //! Destructor.  Automatically closes the shared library.
    ~Shared_Lib() { close(); }

    // Assignment.
    Shared_Lib &operator=(const Shared_Lib &rhs);

    // Closes the shared library.
    void close();

    //! Returns a handle to the shared library.
    void *get_handle() const { Require(is_open()); return d_handle; }

    //! Returns the shared file name.
    std::string get_file_name() const { return d_file_name; }

    // Returns a function pointer from the shared library.
    template <class Fp_t>
    inline Fp_t get_function(const std::string &name);

    //! Returns true if library is open.
    bool is_open() const { return d_handle; }

    // Opens a shared library.
    void open(const std::string &file_name);
};

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*!
  \brief Returns a function pointer from the shared library.

  The shared library must be opened before using this function.

  \param name The name of the function in the shared lib.
  \param Fp_t The function pointer type for the function \a name.
 */
template <class Fp_t>
Fp_t Shared_Lib::get_function(const std::string &name)
{
    Require(is_open());
    Require(not name.empty());

    // This cast is so evil, one must use an old-style cast.
    Fp_t f = Fp_t(dlsym(d_handle, name.c_str()));

    char *error_msg = dlerror();
    Insist(not error_msg, error_msg);

    return f;
}

} // end namespace rtt_dsxx

#endif // rtt_ds_Shared_Lib_hh

//---------------------------------------------------------------------------//
//              end of ds++/Shared_Lib.hh
//---------------------------------------------------------------------------//
