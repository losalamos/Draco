//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   parser/Abstract_Class_Parser.cc
 * \author Kent Budge
 * \brief  Define destructor for Abstract_Class_Parser_Base
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Abstract_Class_Parser.hh"

namespace rtt_parser {

c_string_vector abstract_class_parser_keys;

//---------------------------------------------------------------------------//
c_string_vector::~c_string_vector() {
  unsigned const n = data.size();
  for (unsigned i = 0; i < n; ++i) {
    delete[] data[i];
  }
}

} // end namespace rtt_parser

//---------------------------------------------------------------------------//
// end of parser/Abstract_Class_Parser.cc
//---------------------------------------------------------------------------//
