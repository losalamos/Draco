/*!
  \file    ds++/Compiler.hh
  \author  Paul Henning
  \brief   Some compiler-specific definitions
  \note    Copyright 2006 Los Alamos National Security, LLC.
  \version $Id$
*/

#ifndef Compiler_hh
#define Compiler_hh

/* 
   These are GNU C/C++ extensions that affect the visibility of
   functions and classes in ELF objects.  Although a C++ compiler
   enforces the concept of a "private" member functions, the generated
   code for those member functions is still globally visible in the
   shared libraries.  Such functions incur relocation overhead that
   can be supressed by making them locally visible.

   To use the function visibility macros, simply append the macro to
   the end of the function _declaration_.  For example:

   int foo() const HIDE_FUNC;

   would declare the member function foo as local.

   For classes, the syntax is slightly different.  Place the macro
   between the "class" keyword and the name of the class.  For
   example:

   class EXPORT_CLASS Bar { ... whatever ... };

   would give Bar global visibility.   

   NOTE NOTE NOTE NOTE:  Any class/struct that is used as a type for
   throwing exceptions needs to have an EXPORT_CLASS macro!
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
