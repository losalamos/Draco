//----------------------------------*-C++-*----------------------------------//
/*!
  \file   SetProps.hh
  \author lowrie
  \date   2002-04-12
  \brief  Header for SetProps.
  \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
          All rights reserved.
*/
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef INCLUDED_plot2D_SetProps_hh
#define INCLUDED_plot2D_SetProps_hh

#include "LineProps.hh"
#include "SymbolProps.hh"
#include <string>

namespace rtt_plot2D {

//===========================================================================//
/*!
  \class SetProps

  \brief Set properties for Plot2D class.

  See Grace documentation for a detailed explanation of properties.
*/
//===========================================================================//
class SetProps {
public:
  SetProps(void) : symbol(), line(), legend(){/*empty*/};

  //! The symbol properties
  SymbolProps symbol;

  //! The line properties
  LineProps line;

  //! Legend title
  std::string legend;
};

} // namespace rtt_plot2D

#endif // INCLUDED_plot2D_SetProps_hh

//---------------------------------------------------------------------------//
// end of plot2D/SetProps.hh
//---------------------------------------------------------------------------//
