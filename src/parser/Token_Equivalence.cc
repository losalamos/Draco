//----------------------------------*-C++-*----------------------------------//
/*! 
 * \file Token_Equivalence.cc
 * \author Kelly G. Thompson
 * \date Thu Jul 20 9:27:29 MST 2006
 * \brief Provide services for ApplicationUnitTest framework.
 * \note Copyright © 2006 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <string>
#include <sstream>
#include <stdlib.h>
#include "ds++/Soft_Equivalence.hh"
#include "Token_Equivalence.hh"

namespace rtt_parser 
{

//---------------------------------------------------------------------------//
/*!
 * \brief Search token stream for keyword and compare values.
 *
 * Only works for KEYWORD Tokens.
 */
void check_token_keyword( String_Token_Stream       & tokens,
                          std::string         const & keyword,
                          rtt_dsxx::UnitTest        & ut,
                          unsigned            const & occurance )
{
    std::ostringstream msg;
    unsigned count(0);
    bool done(false);

    tokens.Rewind();
    Token token( tokens.Lookahead(1) );

    while( !done               &&
           token.Type() != END &&
           token.Type() != EXIT )
    {
        if( token.Type() == KEYWORD &&
            token.Text() == keyword &&
            ++count == occurance        )
        {
            msg << "Found the keyword \"" << keyword
                << "\"." << std::endl;
            ut.passes(msg.str());
            done = true;
        }
        
        // Get the next token in the stream.
        token = tokens.Shift();
    }

    if(!done)
    {
        msg << "Did not find the keyword \"" << keyword
            << "\" in the token stream." << std::endl;
        ut.failure(msg.str());
    }
    return;
}
//---------------------------------------------------------------------------//
/*!
 * \brief Search token stream for keyword and compare values.
 *
 * Only works for KEYWORD Tokens.
 */
void check_token_keyword_value( String_Token_Stream       & tokens,
                                std::string         const & keyword,
                                int                 const   expected_value,
                                rtt_dsxx::UnitTest        & ut,
                                unsigned            const & occurance )
{
    std::ostringstream msg;
    unsigned count(0);
    bool done(false);

    tokens.Rewind();
    Token token( tokens.Lookahead(1) );

    while( !done               &&
           token.Type() != END &&
           token.Type() != EXIT )
    {
        if( token.Type() == KEYWORD &&
            token.Text() == keyword &&
            ++count == occurance        )
        {
            // Get the value token
            token = tokens.Shift();

            // Check it's type.
            if( token.Type() != INTEGER )
            {
                msg << "Did not find the token " << keyword <<
                    " in the String_Token_Stream." << std::endl;
                ut.failure(msg.str());
                done = true;
            }

            // Get the actual value.
            int value( atoi(token.Text().c_str()));
            if( value == expected_value )
            {
                msg << "Keyword \"" << keyword
                    << "\" has the expected value of "
                    << expected_value << "." << std::endl;
                ut.passes(msg.str());
            }
            else
            {
                msg << "Keyword \"" << keyword
                    << "\" did not have the expected value of "
                    << expected_value << ".\n\t  Instead we found " << value
                    << "." << std::endl;
                ut.failure(msg.str());
            }
            done = true;
        }

        // Get the next token in the stream.
        token = tokens.Shift();
    }
    return;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Search token stream for keyword and compare values.
 *
 * Only works for KEYWORD Tokens.
 */
void check_token_keyword_value( String_Token_Stream       & tokens,
                                std::string         const & keyword,
                                double              const   expected_value,
                                rtt_dsxx::UnitTest        & ut,
                                unsigned            const & occurance )
{
    std::ostringstream msg;
    unsigned count(0);
    bool done(false);

    tokens.Rewind();
    Token token( tokens.Lookahead(1) );

    while( !done               &&
           token.Type() != END &&
           token.Type() != EXIT )
    {
        if( token.Type() == KEYWORD &&
            token.Text() == keyword &&
            ++count == occurance        )
        {
            // Get the value token
            token = tokens.Shift();

            // Check it's type.
            if( token.Type() != REAL )
            {
                msg << "Did not find the token " << keyword <<
                    " in the String_Token_Stream." << std::endl;
                ut.failure(msg.str());
                done = true;
            }

            // Get the actual value.
            double value( atof(token.Text().c_str()));
            if( rtt_dsxx::soft_equiv(value, expected_value, 1.0e-7) )
            {
                msg << "Keyword \"" << keyword
                    << "\" has the expected value of "
                    << expected_value << "." << std::endl;
                ut.passes(msg.str());
            }
            else
            {
                msg << "Keyword \"" << keyword
                    << "\" did not have the expected value of "
                    << expected_value << ".\n\t  Instead we found " << value
                    << "." << std::endl;
                ut.failure(msg.str());
            }
            done = true;
        }

        // Get the next token in the stream.
        token = tokens.Shift();
    }
    return;
}


}  // namespace rtt_parser

//--------------------------------------------------------------------//
// end of Token_Equivalence.cc
//--------------------------------------------------------------------//
