//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/XGetopt.hh
 * \author Katherine Wang
 * \date   Wed Nov 10 09:35:09 2010
 * \brief  Command line argument handling similar to getopt.
 * \note   Copyright (C) 2015 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_dsxx_XGetopt_hh
#define rtt_dsxx_XGetopt_hh

#include "ds++/config.h"
#include <string>
#include <map>

namespace rtt_dsxx
{

//===========================================================================//
/*!
 * \class XGetopt
 *
 * Function for parsing command line options
 * 
 * Arguments argc and argv are the argument count and array as
 * passed into the application on program invocation.  In the case
 * of Visual C++ programs, argc and argv are available via the
 * variables __argc and __argv (double underscores), respectively.
 * getopt returns the next option letter in argv that matches a
 * letter in optstring.
 * 
 * optstring is a string of recognized option letters;  if a letter
 * is followed by a colon, the option is expected to have an argument
 * that may or may not be separated from it by white space.  optarg
 * is set to point to the start of the option argument on return from
 * getopt.
 * 
 * Option letters may be combined, e.g., "-ab" is equivalent to
 * "-a -b".  Option letters are case sensitive.
 * 
 * getopt places in the external variable optind the argv index
 * of the next argument to be processed.  optind is initialized
 * to 0 before the first call to getopt.
 * 
 * When all options have been processed (i.e., up to the first
 * non-option argument), getopt returns EOF, optarg will point
 * to the argument, and optind will be set to the argv index of
 * the argument.  If there are no non-option arguments, optarg
 * will be set to NULL.
 * 
 * The special option "--" may be used to call longopts (string
 * commands, more than one character long) if defined as a mapped
 * value combined with a key value (single character command line
 * option).  "--" may be used to call the end of the options (skip
 * remaining options) when argv does not have a corresponding mapped
 * value or when it is followed by shortopts, defined or undefined.
 * 
 * For option letters contained in the string optstring, getopt
 * will return the option letter.  getopt returns a question mark (?)
 * when it encounters an option letter not included in optstring.
 * EOF is returned when processing is finished.
 * 
 * Bugs
 * 1) The GNU double-colon extension is not supported.
 * 2) TThe environment variable POSIXLY_CORRECT is not supported.
 * 3) The + syntax is not supported.
 * 4) The automatic permutation of arguments is not supported.
 * 5) This implementation of getopt() returns EOF if an error is
 *    encountered, instead of -1 as the latest standard requires.
 *
 * \code
 *	 std::map< std::string, char> long_options;
 *       long_options["add"]    = 'a';
 *       long_options["append"] = 'b';
 *   	 long_options["create"] = 'c';
 * 
 *       BOOL CMyApp::ProcessCommandLine(int argc, TCHAR *argv[])
 *       {
 *           int c;
 * 
 *           while ((c = rtt_dsxx::getopt (argc, argv, "abc:", long_options)) != -1)
 *   	     {
 *       	switch (c)
 *       	{
 *          	  case 'a':
 *                   //
 *                   // set some flag here
 *                   //
 *                   break;
 *
 *           	  case 'b':
 *                   //
 *                   // set some flag here
 *                   //
 *                   break;
 *
 *	          case 'c':
 *                   cvalue = rtt_dsxx::optarg;
 *                   break;
 *
 *           	  default:
 *                   return 0; // nothin to do.
 *       	 }
 *            }
 *
 *	    return 0;
 *
 *       }
 * \endcode
 *
 */
 /*!
 * \example ds++/test/tstXGetopt.cc
 */
//===========================================================================//
 
DLL_PUBLIC_dsxx extern int   optind;
DLL_PUBLIC_dsxx extern char *optarg;

typedef std::map< std::string, char> longopt_map;

DLL_PUBLIC_dsxx 
int getopt( int argc, char **& argv, std::string const & shortopts,
            longopt_map map = std::map<std::string, char>() );

}

#endif // rtt_dsxx_XGetopt_hh

//---------------------------------------------------------------------------//
// end of XGetOpt.hh
//---------------------------------------------------------------------------//
