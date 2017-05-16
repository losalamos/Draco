//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   parser/Class_Parse_Table.hh
 * \author Kent Budge
 * \brief  Define template function parse_class
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved.
 *
 * utilities.hh declares a function template for parsers for class objects,
 * which is consistent in format with the other parse functins in utilities.hh
 * No implementation is provided. However, we recommend using the templates in
 * this header (Class_Parse_Table.hh) to provide an implementation of the form
 * \code
 * template<>
 * std::shared_ptr<Class> parse_class(Token_Stream &tokens) {
 *   return parse_class_from_table<Class_Parse_Table<Class> >(tokens);
 * }
 *
 * Why don't we define this as the default implementation? Because then any file
 * that #included Class_Parse_Table.hh would attempt to use the default
 * implementation for every call to parse_class in the file. When you are
 * parsing objects that include class subobjects, this would require including
 * the definition of the parse table class for every such subobject type. This
 * breaks encapsulation, to put it mildly. So we refrain from providing such a
 * default implementation.
 */
//----------------------------------------------------------------------------//

#ifndef parser_Class_Parse_Table_hh
#define parser_Class_Parse_Table_hh

#include "Parse_Table.hh"

namespace rtt_parser {

//----------------------------------------------------------------------------//
/*! Template for parse table classes
 *
 * No general implementation is provided. The developer wishing to make his
 * class parsable must specialize this template for his class, making sure that
 * the following members are defined:
 *
 *   Class_Parse_Table()
 *
 * Default constructor. Of course, the actual name of the class will be used
 * here. The constructor should set all data members to sentinel values
 * indicating that the corresponding parameters have not been parsed yet.  This
 * allows us to detect whether a parameter has been parsed yet.
 *
 * Where possible, the sentinel value should be one that doesn't make
 * sense for the parameter in question, such as 0 for a dimension that
 * must be nonzero, or a negative value for a parameter that must be
 * positive.  This ensures that the sentinel value can be distinguished
 * from any value a user might specify in an input "deck."
 *
 * Since Parse_Table uses static functions to parse the token stream, it is
 * generally necessary to have a static pointer somewhere that points to the
 * object being parsed and is accessible to the parser functions. By convention,
 * the pointer is declared as a private member of Class_Parse_Table with the
 * signature
 *
 *   static Class_Parse_Table *current_;
 *
 * and the parser functions are then declared as static private member functions
 * of Class_Parse_Table. However, since parse_class uses neither current_ nor
 * any of the parser functions directly, this is not strictly required -- just
 * highly recommended.
 *
 * Likewise, to ensure the parse will be reentrant, the constuctor should save
 * the previous value of current_ before changing it, and a destructor should be
 * defined that restores the previous value. Not all parsers are reentrant so
 * this is also not strictly required.
 *
 * The constructor should also carry out any other preparations that are
 * necessary before parsing.  It should have no preconditions.
 *
 *   Parse_Table const &parse_table() const;
 *
 * Returns a reference to the Parse_Table that will be used to do the actual
 * parsing. This is typically a private static member of Class_Parse_Table.
 *
 *   bool allow_exit() const;
 *
 * Is it acceptable for the input to the parser to terminate with an EXIT token,
 * or is only an END token acceptable?
 *
 *   void check_completeness(Token_Stream &);
 *
 * This function is called to check that all required specifications were found
 * in the parsed stream.  This function must call \c
 * Token_Stream::Report_Semantic_Error at least once if it cannot guarantee that
 * the subsequent call to create_object will be successful.
 *
 * If default values are permitted for some parameters (which is not
 * recommended), then they should be applied here.
 *
 *   shared_ptr<Return_Class> create_object();
 *
 * Create the object from the parsed fields.  This function should have no
 * preconditions that are not guaranteed by a preceding successful call to
 * check_completeness. create_object will be called only if no errors
 * were detected in the input "deck".
 *
 * Note that Class_Parse_Table must define the type Return_Class, which is the
 * class of the object that Class_Parse_Table is designed to parse.
 *
 * For convenience, skeletons for such a class are provided in
 * draco/environment/templates/template__parser.hh,cc
 */

template <class Class> class Class_Parse_Table;

//----------------------------------------------------------------------------//
/*! Template for helper function that produces a class object.
 *
 * \param tokens Token stream from which to parse the user input.
 *
 * \return A pointer to an object matching the user specification, or NULL if
 * the specification is not valid.
 */

template <class Class_Parse_Table>
std::shared_ptr<typename Class_Parse_Table::Return_Class>
parse_class_from_table(Token_Stream &tokens) {
  using rtt_parser::Token;
  using rtt_parser::END;
  using rtt_parser::EXIT;

  typedef typename Class_Parse_Table::Return_Class Return_Class;

  // Construct the parse object as described above.
  Class_Parse_Table parse_table;

  // Save the old error count, so we can distinguish fresh errors within
  // this class keyword block from previous errors.
  unsigned const old_error_count = tokens.error_count();

  // Parse the class keyword block and check for completeness
  Token const terminator = parse_table.parse_table().parse(tokens);
  bool allow_exit = parse_table.allow_exit(); // improve code coverage
  std::shared_ptr<Return_Class> Result;
  if (terminator.type() == END || (allow_exit && terminator.type() == EXIT))
  // A class keyword block is expected to end with an END or (if
  // allow_exit is true) an EXIT.
  {
    parse_table.check_completeness(tokens);

    if (tokens.error_count() == old_error_count) {
      // No fresh errors in the class keyword block.  Create the object.
      Result = parse_table.create_object();
    }
    // else there were errors in the keyword block. Don't try to
    // create a class object.  Return the null pointer.
  } else {
    tokens.report_syntax_error("missing 'end'?");
  }
  return Result;
}

//----------------------------------------------------------------------------//
/*! Template for helper function that produces a class object.
 *
 * \param tokens Token stream from which to parse the user input.
 *
 * \param context A context object that controls the behavior of the
 * parser. This is passed to the constructor for the Class_Parse_Table.
 *
 * \return A pointer to an object matching the user specification, or NULL if
 * the specification is not valid.
 */

template <typename Class_Parse_Table, typename Context>
std::shared_ptr<typename Class_Parse_Table::Return_Class>
parse_class_from_table(Token_Stream &tokens, Context const &context) {
  using rtt_parser::Token;
  using rtt_parser::END;
  using rtt_parser::EXIT;

  typedef typename Class_Parse_Table::Return_Class Return_Class;

  // Construct the parse object as described above.
  Class_Parse_Table parse_table(context);

  // Save the old error count, so we can distinguish fresh errors within
  // this class keyword block from previous errors.
  unsigned const old_error_count = tokens.error_count();

  // Parse the class keyword block and check for completeness
  Token const terminator = parse_table.parse_table().parse(tokens);
  bool allow_exit = parse_table.allow_exit(); // improve code coverage
  std::shared_ptr<Return_Class> Result;
  if (terminator.type() == END || (allow_exit && terminator.type() == EXIT))
  // A class keyword block is expected to end with an END or (if
  // allow_exit is true) an EXIT.
  {
    parse_table.check_completeness(tokens);

    if (tokens.error_count() == old_error_count) {
      // No fresh errors in the class keyword block.  Create the object.
      Result = parse_table.create_object();
    }
    // else there were errors in the keyword block. Don't try to
    // create a class object.  Return the null pointer.
  } else {
    tokens.report_syntax_error("missing 'end'?");
  }

  return Result;
}

} // end namespace rtt_parser

#endif // parser_Class_Parse_Table_hh

//----------------------------------------------------------------------------//
// end of parser/Class_Parse_Table.hh
//----------------------------------------------------------------------------//
