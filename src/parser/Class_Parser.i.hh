//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   utils/Class_Parser.i.hh
 * \author Kent Budge
 * \date   Mon Aug 28 07:30:16 2006
 * \brief  Member template definitions of template class Class_Parser
 * \note   Copyright (C) 2006-2014 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef utils_Class_Parser_i_hh
#define utils_Class_Parser_i_hh

// This file defines the main parse routine for a Class_Parser. This routine
// creates an object of the templated ReturnClass based on the specification
// in the input "deck".
//
// For further explanation of the three template classes, see the Class_Parser
// template definition.

//---------------------------------------------------------------------------//
/*! 
 *
 * \param tokens Token stream from which to parse the user input.
 *
 * \param allow_exit If \c true, the class parser is permitted to terminate on
 * either an END or an EXIT (end of file) token. If \c false, the class parser
 * may only terminate on an END token. This allows a host code to check for a
 * missing END in an input deck.
 *
 * \return A pointer to an object matching the user specification, or NULL if
 * the specification is not valid.
 */

template<typename Class,
         typename ReturnClass,
         bool is_reentrant,
         typename ParseTableClass>
rtt_dsxx::SP<ReturnClass>
Class_Parser<Class,
             ReturnClass,
             is_reentrant,
             ParseTableClass>::parse(rtt_parser::Token_Stream &tokens,
                                     bool const allow_exit)
{
    using rtt_dsxx::SP;
    //using namespace rtt_parser;
    using rtt_parser::Token;
    using rtt_parser::END;
    using rtt_parser::EXIT;
    
    // Is this the first call to the parser?
    static bool is_first_time = true;
    
    // The following code checks for reentrancy and treats it as a fatal
    // error, unless the is_reentrant parameter is set to show that the parse
    // has been carefully coded for safe reentrancy.
    static bool has_reentered = false;
    Insist(is_reentrant || !has_reentered,
           "Class_Parser::Parse_Specification is not reentrant");
    has_reentered = true;

    try
    {
        // Assign sentinel values to all parsed parameters if this is the
        // first time we have called this parse routine. Subsequent calls need
        // not do this since sentinel values are set at the conclusion of each
        // parse.
        if (is_first_time)
        {
            post_sentinels_();
            is_first_time = false;
        }
        
        // Save the old error count, so we can distinguish fresh errors within
        // this class keyword block from previous errors.
        unsigned const old_error_count = tokens.error_count();
        
        // Parse the class keyword block and check for completeness
        Token const terminator = parse_table_.parse(tokens);
        if (terminator.type() == END ||
            (allow_exit && terminator.type()==EXIT))
            // A class keyword block is expected to end with an END or (if
            // allow_exit is true) an EXIT.
        {
            check_completeness_(tokens);
            
            SP<ReturnClass> Result;
            if (tokens.error_count() == old_error_count)
            {
                // No fresh errors in the class keyword block.  Create the object.
                Result = create_object_();
            }
            // else there were errors in the keyword block. Don't try to
            // create a class object.  Return the null pointer.
            
            post_sentinels_();
            has_reentered = false;
            
            // In some cases we may wish to return an empty pointer.
            //Ensure(Result!=SP<ReturnClass>() ||
            //       tokens.error_count()>old_error_count);

            return Result;
        }
        else
        {
            tokens.report_syntax_error("missing 'end'?");
            return SP<ReturnClass>();
            // never reached; to eliminate spurious warning.
        }
    }
    catch (...)
    {
        // Reset the reentrancy flag so a future call will not be mistaken for
        // a reentrant call. To be safe, ensure that all parsed parameters
        // will be reset on the next call. Then rethrow the exception.
        is_first_time = true;
        has_reentered = false;
        throw;
    }
}

#endif // utils_Class_Parser_i_hh

//---------------------------------------------------------------------------//
// end of utils/Class_Parser.i.hh
//---------------------------------------------------------------------------//
