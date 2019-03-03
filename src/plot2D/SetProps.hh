//----------------------------------*-C++-*----------------------------------//
/*!
  \file   SetProps.hh
  \author lowrie
  \date   2002-04-12
  \brief  Header for SetProps.
  \note   Copyright (C) 2016-2019 Triad National Security, LLC.
          All rights reserved.
*/
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//

#ifndef rtt_plot2D_SetProps_hh
#define rtt_plot2D_SetProps_hh

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

#endif // rtt_plot2D_SetProps_hh

//---------------------------------------------------------------------------//
// end of plot2D/SetProps.hh
//---------------------------------------------------------------------------//
