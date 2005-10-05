//----------------------------------*-C++-*----------------------------------//
/*! 
 *  \file   parser/Release.hh
 *  \author Kent G. Budge
 *  \date   Thu Jul 15 09:31:44 1999
 *  \brief  Header file for parser library release function.
 */
//===========================================================================//
// $Id$
//===========================================================================//

#ifndef __parser_Release_hh__
#define __parser_Release_hh__

#include <string>

//===========================================================================//
/*!
 * \page parser-overview Overview of the Parser Package
 *
 * \version 0_0_0
 *
 * This package provides data structures and other miscellaneous support
 * for simple keyword-driven input parsing.
 *
 */
//===========================================================================//
/*!
 * \namespace rtt_parser
 *
 * \brief RTT parser namespace.
 *
 * The parser package supplies several tokenizers (scanners) as children of
 * the abstract class Token_Stream. It also supplies a Parse_Table class that
 * parses a stream of tokens based on a table of key words and associated
 * parameter parsing functions.  A number of miscellaneous useful functions
 * is supplied in utilities.hh.
 *
 * In addition, the package supplies a unit library.  This offers
 * functionality that is distinct from the Draco unit library, but perhaps it
 * should be merged with the latter at some point.
 *
 */
//===========================================================================//

namespace rtt_parser
{

//! Query package for the release number.
const std::string release();

} // end of rtt_parser

#endif                          // __parser_Release_hh__

//---------------------------------------------------------------------------//
//                              end of parser/Release.hh
//---------------------------------------------------------------------------//
