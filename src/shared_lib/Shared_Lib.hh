//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   shared_lib/Shared_Lib.hh
 * \author Rob Lowrie
 * \date   Thu Apr 15 20:44:39 2004
 * \brief  Header file for Shared_Lib.
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef rtt_shared_lib_Shared_Lib_hh
#define rtt_shared_lib_Shared_Lib_hh

#include <shared_lib/config.h>

#include <ds++/Assert.hh>
#include <string>

namespace rtt_shared_lib {

//===========================================================================//
/*!
  \class Shared_Lib
  \brief Controls access to a shared (dynamically linked) library.

  Access to functions defined in the shared library is provided via the
  get_function() member.

  Under Draco, not all platforms support dynamic loading of shared libraries.
  Consequently, in order to write cross-platform code, one must use the
  static member function is_supported() to check whether the functionality of
  Shared_Lib is supported.  As an example,
  \code
  if ( Shared_Lib::is_supported() )
  {
     // OK, Shared_Lib is supported.
     Shared_Lib s;
     s.open("/usr/lib/libm.so");
     // ... other operations using s.
  }
  else
  {
     // Shared_Lib is unsupported!
     Shared_Lib s;               // throws an error!!!
     s.open("/usr/lib/libm.so"); // won't get this far.
  }
  \endcode
  Note that the above code should compile on all platforms, but on
  unsupported platforms, the "else" block will throw an error at run time.
 */
/*!
 * \example shared_lib/test/tstShared_Lib.cc
 *
 * This example shows how classes may be created through shared objects.
*/
//===========================================================================//

class Shared_Lib {
  // DATA

  // The handle to the shared library.
  void *d_handle;

  // The name of the shared library.
  std::string d_file_name;

public:
  // Default constructor.
  explicit Shared_Lib(const std::string &file_name = "");

  // Copy constructor.
  explicit Shared_Lib(const Shared_Lib &from);

  //! Destructor.  Automatically closes the shared library.
  ~Shared_Lib() { close(); }

  // Assignment.
  Shared_Lib &operator=(const Shared_Lib &rhs);

  // Closes the shared library.
  void close();

  //! Returns a handle to the shared library.
  void *get_handle() const {
    Require(is_open());
    return d_handle;
  }

  //! Returns the shared file name.
  std::string get_file_name() const { return d_file_name; }

  // Returns a function pointer from the shared library.
  template <class Fp_t> inline Fp_t get_function(const std::string &name);

  //! Returns true if library is open.
  bool is_open() const { return d_handle; }

  //! Returns true if platform is supported.
  static bool is_supported();

  // Opens a shared library.
  void open(const std::string &file_name);

private:
  // Does the dlsym() with error checking.
  void *do_dlsym(const std::string &name);
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
template <class Fp_t> Fp_t Shared_Lib::get_function(const std::string &name) {
  Require(is_open());

  // HACK WARNING: 5.2.10/6-7 implies that we cannot cast a
  // pointer-to-object (in this case, the void* from dlsym) to a
  // pointer-to-function.  If void* and Fp_t are different sizes, I suspect
  // the hack below may not be portable. In the end, I suspect that platforms
  // where this hack does not work don't support dlopen, anyway - lowrie

  union {
    void *vp;
    Fp_t fp;
  };
  vp = do_dlsym(name);

  return fp;
}

} // end namespace rtt_shared_lib

#endif // rtt_shared_lib_Shared_Lib_hh

//---------------------------------------------------------------------------//
//              end of shared_lib/Shared_Lib.hh
//---------------------------------------------------------------------------//
