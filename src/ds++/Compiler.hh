//----------------------------------*-C++-*----------------------------------//
/*!
  \file    ds++/Compiler.hh
  \author  Paul Henning
  \brief   Some compiler-specific definitions
  \note    Copyright 2006-2010 Los Alamos National Security, LLC.
  \version $Id$
*/
//---------------------------------------------------------------------------//

#ifndef Compiler_hh
#define Compiler_hh

/*!
 * \section compilerBackground Background
 * These are GNU C/C++ extensions that affect the visibility of functions
 * and classes in ELF objects.  Although a C++ compiler enforces the concept
 * of a "private" member functions, the generated code for those member
 * functions is still globally visible in the shared libraries.  Such
 * functions incur relocation overhead that can be supressed by making them
 * locally visible.
 *
 * To use the function visibility macros, simply append the macro to the end
 * of the function _declaration_.  For example:
 * \code
 *  int foo() const HIDE_FUNC;
 * \endcode
 * would declare the member function foo as local.

 * For classes, the syntax is slightly different.  Place the macro between the
 * "class" keyword and the name of the class.  For example:
 * \code
 *  class EXPORT_CLASS Bar { ... whatever ... };
 * \endcode
 * would give Bar global visibility.   
 *
 * NOTE: Any class/struct that is used as a type for throwing exceptions needs
 * to have an \c EXPORT_CLASS macro!
 */

/*!
 * \section compilerNote1 2010 July 1 (Kelly Thompson):
 * \sa http://gcc.gnu.org/wiki/Visibility
 * There is a more detailed discussion on the GCC wiki
 * (http://gcc.gnu.org/wiki/Visibility).  In particular, we may need to
 * treat ELF symbol visibility in a way that is portable across platforms
 * (Cygwin, Windows, etc.).
 *
 * Benefits of hidding ELF symbols:
 * - Improved load time for dynamic shared objects, DSO.
 * - The compiler can often provide better optimized code.
 * - Smaller DSO files.
 * - Reduced chance of symbol collision.
 */

#if  __GNUC__ >=4
#define HIDE_FUNC __attribute__ ((visibility ("hidden")))
#define EXPORT_FUNC __attribute__ ((visibility ("default")))
#define HIDE_CLASS __attribute ((visibility ("hidden")))
#define EXPORT_CLASS __attribute ((visibility ("default")))
#else
#define HIDE_FUNC
#define EXPORT_FUNC
#define HIDE_CLASS
#define EXPORT_CLASS
#endif

#endif
