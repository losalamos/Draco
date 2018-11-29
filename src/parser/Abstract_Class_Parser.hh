//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   parser/Abstract_Class_Parser.hh
 * \author Kent Budge
 * \brief  Define class Abstract_Class_Parser
 * \note   Copyright (C) 2016-2018 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#ifndef parser_Abstract_Class_Parser_hh
#define parser_Abstract_Class_Parser_hh

#include "Parse_Table.hh"
#include <functional>
#include <iostream>

namespace rtt_parser {
using std::pointer_to_unary_function;
using std::string;
using std::vector;

//===========================================================================//
/*!
 * \brief Closure class for wrapping context-dependent parse functions
 *
 * This template class is used to bind a get_context() function to a parse
 * function requiring a context argument, so that the function can be called
 * with the usual two parameters (Token_Stream and int) from an
 * Abstract_Class_Parse_Table.
 *
 * See test/tstAbstract_Class_Contextual_Parser.cc for an example of how it is
 * used.
 */
//===========================================================================//
template <typename Abstract_Class, typename Context,
          Context const &get_context()>
class Contextual_Parse_Functor {
public:
  Contextual_Parse_Functor(std::shared_ptr<Abstract_Class> parse_function(
      Token_Stream &, Context const &));

  std::shared_ptr<Abstract_Class> operator()(Token_Stream &) const;

private:
  std::shared_ptr<Abstract_Class> (*f_)(Token_Stream &, Context const &);
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
 * \code
 * abstract class keyword
 *   child class keyword
 *     (child class specification)
 *   end
 * end
 * \endcode
 *
 * Note that Abstract_Class_Parser does not actually do any parsing itself. It
 * is simply a repository for keyword-parser combinations that is typically used
 * by the Class_Parser for the abstract class.
 *
 * See test/tstAbstract_Class_Parser.cc for an example of its use.
 *
 * This template has proven useful but does not provide a fully satisfactory
 * solution to the problem of abstract class keywords other than those
 * specifying a child class.
 */
//===========================================================================//
template <typename Abstract_Class, Parse_Table &get_parse_table(),
          std::shared_ptr<Abstract_Class> &get_parsed_object(),
          typename Parse_Function = pointer_to_unary_function<
              Token_Stream &, std::shared_ptr<Abstract_Class>>>
class Abstract_Class_Parser {
public:
  // TYPES

  // STATIC members

  //! Register children of the abstract class
  static void register_child(string const &keyword,
                             Parse_Function parse_function);

  //! Register children of the abstract class
  static void register_child(
      string const &keyword,
      std::shared_ptr<Abstract_Class> parse_function(Token_Stream &));

  //! Check the class invariants
  static bool check_static_class_invariants();

private:
  // IMPLEMENTATION

  //! Parse the child type
  static void parse_child_(Token_Stream &, int);

  // DATA

  //! Map of child keywords to child creation functions
  static vector<Parse_Function> map_;
};

#include "Abstract_Class_Parser.i.hh"

} // end namespace rtt_parser

#endif // parser_Abstract_Class_Parser_hh

//---------------------------------------------------------------------------//
// end of parser/Abstract_Class_Parser.hh
//---------------------------------------------------------------------------//
