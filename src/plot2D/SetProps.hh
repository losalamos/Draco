//----------------------------------*-C++-*----------------------------------//
/*!
  \file   SetProps.hh
  \author lowrie
  \date   2002-04-12
  \brief  Header for SetProps.
*/
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef INCLUDED_plot2D_SetProps_hh
#define INCLUDED_plot2D_SetProps_hh

#include "LineProps.hh"
#include "SymbolProps.hh"

namespace rtt_plot2D {

//===========================================================================//
/*!
  \struct SetProps

  \brief Set properties for Plot2D class.

  See Grace documentation for a detailed explanation of properties.
*/
//===========================================================================//
struct SetProps
{
    /// The symbol properties
    SymbolProps symbol;

    /// The line properties
    LineProps line;

    /// Legend title
    std::string legend;
};

} // namespace rtt_plot2D

#endif // INCLUDED_plot2D_SetProps_hh

//---------------------------------------------------------------------------//
// end of plot2D/SetProps.hh
//---------------------------------------------------------------------------//
