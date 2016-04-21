//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   parser/Abstract_Class_Parser.hh
 * \author Kent Budge
 * \brief  Define class Abstract_Class_Parser
 * \note   Copyright (C) 2016 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef parser_Abstract_Class_Parser_hh
#define parser_Abstract_Class_Parser_hh

#include "Parse_Table.hh"
#include "ds++/SP.hh"
#include <iostream>

namespace rtt_parser
{
using std::string;
using std::vector;
using rtt_dsxx::SP;

//===========================================================================//
/*!
 * \class Abstract_Class_Parser_Base
 * \brief Template for parser that produces a class object.
 *
 * This class exists only to serve as a base for Abstract_Class_Parser,
 * allowing all such parsers to share the same keyword table and ensuring that
 * the keyword table is properly cleaned up when the program terminates.
 */
class DLL_PUBLIC_parser Abstract_Class_Parser_Base
{
  protected:

    // TYPES

    class DLL_PUBLIC_parser c_string_vector
    {
      public:
        ~c_string_vector();
        c_string_vector(void) : data(0) {/* empty */}
        vector<char *> data;
    };

    // provide a virtual destrcutor for the base class.
    virtual ~Abstract_Class_Parser_Base() {/* empty */};

    // DATA

    //! Keywords
    static c_string_vector keys_;
};

//===========================================================================//
/*!
 * \class Abstract_Class_Parser
 * \brief Template for parser that produces a class object.
 *
 * This template is meant to be specialized for parse tables that select one
 * of a set of child classes of a single abstract class. It simplifies and
 * regularizes the task of allowing additional child classes to be added to
 * the table by a local developer working on his own version of one of the
 * Capsaicin drivers.
 *
 * \arg \a Abstract_Class The abstract class whose children are to be parsed.
 *
 * \arg \a get_parse_table A function that returns a reference to the parse
 * table for the abstract class.
 *
 * \arg \a get_parsed_object A function that returns a reference to a
 * storage location for a pointer to the abstract class.
 *
 * The key to this class is the register_child function, which is called for
 * each child class prior to attempting any parsing. It specifies a keyword
 * for selecting each child class and a function that does the actual parsing
 * of the class specification. This assumes an input grammar of the form
 *
 * <code>
 * abstract class keyword
 *   child class keyword
 *     (child class specification)
 *   end
 * end
 *
 * Note that Abstract_Class_Parser does not actually do any parsing itself. It
 * is simply a repository for keyword-parser combinations that is typically
 * used by the Class_Parser for the abstract class.
 *
 * See test/tstAbstract_Class_Parser for an example of its use.
 *
 * This template has proven useful but does not provide a fully satisfactory
 * solution to the problem of abstract class keywords other than those
 * specifying a child class.
 */
//===========================================================================//
template< typename Abstract_Class,
          Parse_Table &get_parse_table(),
          SP<Abstract_Class> &get_parsed_object() >
class Abstract_Class_Parser : private Abstract_Class_Parser_Base
{
  public:

    // TYPES

    typedef SP<Abstract_Class> Parse_Function(Token_Stream &);

    // STATIC members

    //! Register children of the abstract class
    static void register_child( string const &keyword, Parse_Function * parse_function );

    //! Check the class invariants
    static bool check_static_class_invariants();

  private:

    // IMPLEMENTATION

    //! Parse the child type
    static void parse_child_(Token_Stream &, int);

    // DATA

    //! Map of child keywords to child creation functions
    static vector<Parse_Function*> map_;
};

#include "Abstract_Class_Parser.i.hh"

} // end namespace rtt_parser

#endif // parser_Abstract_Class_Parser_hh

//---------------------------------------------------------------------------//
// end of parser/Abstract_Class_Parser.hh
//---------------------------------------------------------------------------//
