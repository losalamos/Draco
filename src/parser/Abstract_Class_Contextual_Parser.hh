//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   parser/Abstract_Class_Contextual_Parser.hh
 * \author Kent Budge
 * \brief  Define class Abstract_Class_Contextual_Parser
 * \note   Copyright (C) 2016 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef parser_Abstract_Class_Contextual_Parser_hh
#define parser_Abstract_Class_Contextual_Parser_hh

#include "Abstract_Class_Parser.hh"

namespace rtt_parser {

//===========================================================================//
/*!
 * \class Abstract_Class_Contextual_Parser
 * \brief Template for parser that produces a class object.
 *
 * This class is almost identical to Abstract_Class_Parser, except that the
 * Parse_Function typedef is to a function taking an additional context
 * argument.
 *
 * Note that this class is used only as a namespace; it should never have any
 * nonstatic members. We use a class rather than namespace because we can't
 * templatize names of namespaces.
 *
 * See test/tstAbstract_Class_Contextual_Parser for an example of its use.
 */
//===========================================================================//
template <typename Abstract_Class, Parse_Table &get_parse_table(),
          SP<Abstract_Class> &get_parsed_object(), typename Context,
          Context const &get_context()>
class Abstract_Class_Contextual_Parser {
public:
  // TYPES

  typedef SP<Abstract_Class> Parse_Function(Token_Stream &, Context const &);

  // STATIC members

  //! Register children of the abstract class
  static void register_child(string const &keyword,
                             Parse_Function *parse_function);

  //! Check the class invariants
  static bool check_static_class_invariants();

private:
  // IMPLEMENTATION

  //! Parse the child type
  static void parse_child_(Token_Stream &, int);

  // DATA

  //! Map of child keywords to child creation functions
  static vector<Parse_Function *> map_;
};

#include "Abstract_Class_Contextual_Parser.i.hh"

} // end namespace rtt_parser

#endif // parser_Abstract_Class_Contextual_Parser_hh

//---------------------------------------------------------------------------//
// end of parser/Abstract_Class_Contextual_Parser.hh
//---------------------------------------------------------------------------//
