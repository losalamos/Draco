//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   parser/Class_Parser.hh
 * \author Kent Budge
 * \brief  Define template class Class_Parser
 * \note   Copyright (C) 2006-2014 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef parser_Class_Parser_hh
#define parser_Class_Parser_hh

#include "ds++/SP.hh"
#include "Token_Stream.hh"
#include "Parse_Table.hh"

namespace rtt_parser
{

//===========================================================================//
/*!
 * \class Class_Parser
 * \brief Template for parser that produces a class object.
 *
 * This template is meant to be specialized for every class for which a parser
 * is desired.  A specialization of this template for a particular class must
 * define all of the functions and the parse table described below.
 *
 * For portability reasons, we avoid direct use of this template. Instead, the
 * class being parsed is given a static method,
 * <code>Parsed_Class::parse</code>, that wraps a call to
 * <code>Class_Parser<Parsed_Class>::parse</code>. This helps ensure that
 * <code>Class_Parser<Parsed_Class></code> is instantiated exactly once.
 *
 * Note that the returned class can be different than the nominally parsed
 * class. This can be useful in inheritance hierarchies.
 *
 * A non-reentrant parser can used static storage for parsed values. Using
 * heap storage for parsed values allows a parser to be safely reentered, but
 * this requires a bit more effort, which is wasted if the parser is known to
 * be non-reentrant (and most of ours are.)
 *
 * The is_reentrant template parameter tells the boilerplate Class_Parser
 * machinery whether a parser is supposed to be reentrant or not. If not, the
 * boilerplate Class_Parser machinery will throw an exception if the parser is
 * somehow reentered anyway.
 */
//===========================================================================//

template<class Class, class ReturnClass = Class, bool is_reentrant = false>
class Class_Parser 
{
  public:

    // STATIC

    //! Parse a user specification and create the corresponding object.
    static
    rtt_dsxx::SP<ReturnClass> parse(Token_Stream &,
                                    bool allow_exit = true);

    // DATA

    //-----------------------------------------------------------------------//
    /*! Parse table for the class.
     *
     * This must be public and non-constant if we are to support various
     * mechanisms for local extensions to Capsaicin.
     */
    static Parse_Table parse_table_;

    // IMPLEMENTATION

    // The implementation must be public if we are to support derivation 

    //-----------------------------------------------------------------------//
    /*! Prepare to parse
     *
     * We assume that the parser uses static variables to store the parsed
     * parameters.  This function should reset these static variables to
     * sentinel values indicating that the corresponding parameters have not
     * been parsed yet.  This avoids weird behavior the second time we call
     * the parser, and it also allows us to detect whether a parameter has
     * been parsed yet in the current parse call.  The sentinel value should
     * be one that doesn't make sense for the parameter in question, such as 0
     * for a dimension that must be nonzero, or a negative value for a
     * parameter that must be positive.  This function should also carry out
     * any other preparations that are necessary before parsing.  It should
     * have no preconditions.
     */
    static void post_sentinels_();

    //-----------------------------------------------------------------------//
    /*! Check completeness of specification.
     *
     * This function is called to check that all required specifications were
     * found in the parsed stream.  This function must call \c
     * Token_Stream::Report_Semantic_Error at least once if it cannot
     * guarantee that the next call to create_object will be successful.
     *
     * If default values are permitted for some parameters (which is not
     * recommended) then they should be applied here.
     *
     * \param Token_Stream Stream to which to report any semantic errors in
     * the specification.
     */
    static void check_completeness_(Token_Stream &);

    //-----------------------------------------------------------------------//
    /*! Create the object after completing parsing
     *
     * Create the object from the parsed fields.  This function should have no
     * preconditions that are not guaranteed by a preceding successful call to
     * check_completeness. This function will be called only if no errors were
     * detected in the user specification.
     */
    static rtt_dsxx::SP<ReturnClass> create_object_();
};

#include "Class_Parser.i.hh"

} // end namespace rtt_parser

#endif // parser_Class_Parser_hh

//---------------------------------------------------------------------------//
// end of parser/Class_Parser.hh
//---------------------------------------------------------------------------//
