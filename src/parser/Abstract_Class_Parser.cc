//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   parser/Abstract_Class_Parser.cc
 * \author Kent Budge
 * \brief  Define destructor for Abstract_Class_Parser_Base
 * \note   Copyright (C) 2016-2020 Triad National Security, LLC.
 *         All rights reserved. */
//----------------------------------------------------------------------------//

#include "Abstract_Class_Parser.hh"

//! \bug doxygen doesn't like this file!
//! \cond doxygen_skip

namespace rtt_parser {

c_string_vector abstract_class_parser_keys;

//----------------------------------------------------------------------------//
c_string_vector::~c_string_vector() {
  Check(data.size() < UINT_MAX);
  unsigned const n = static_cast<unsigned>(data.size());
  for (unsigned i = 0; i < n; ++i) {
    delete[] data[i];
  }
}

} // end namespace rtt_parser

//! \endcond

//----------------------------------------------------------------------------//
// end of parser/Abstract_Class_Parser.cc
//----------------------------------------------------------------------------//
