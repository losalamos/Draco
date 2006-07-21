//----------------------------------*-C++-*----------------------------------//
/*! 
 * \file Token_Equivalence.hh
 * \author Kelly G. Thompson
 * \date Thu Jul 20 9:27:29 MST 2006
 * \brief Provide services for ApplicationUnitTest framework.
 * \note Copyright © 2006 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_Token_Equivalence_hh
#define rtt_Token_Equivalence_hh

#include <string>
#include "ds++/UnitTest.hh"
#include "String_Token_Stream.hh"

namespace rtt_parser 
{
//! Search tokens for existance of keyword.
void check_token_keyword( String_Token_Stream       & tokens,
                          std::string         const & keyword,
                          rtt_dsxx::UnitTest        & ut,
                          unsigned            const & occurance=1 );

//! Search tokens for keyword and compare associated value.  Report result to UnitTest.
void check_token_keyword_value( String_Token_Stream       & tokens,
                                std::string         const & keyword,
                                int                 const   value,
                                rtt_dsxx::UnitTest        & ut,
                                unsigned            const & occurance=1 );

//! Search tokens for keyword and compare associated value.  Report result to UnitTest.
void check_token_keyword_value( String_Token_Stream       & tokens,
                                std::string         const & keyword,
                                double              const   value,
                                rtt_dsxx::UnitTest        & ut,
                                unsigned            const & occurance=1 );
}  // namespace rtt_parser

#endif  //  rtt_Token_Equivalence_hh

//--------------------------------------------------------------------//
// end of Token_Equivalence.hh
//--------------------------------------------------------------------//
