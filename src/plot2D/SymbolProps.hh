//----------------------------------*-C++-*----------------------------------//
/*!
  \file   SymbolProps.hh
  \author lowrie
  \date   2002-04-12
  \brief  Header for SymbolProps.
*/
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef INCLUDED_plot2D_SymbolProps_hh
#define INCLUDED_plot2D_SymbolProps_hh

namespace rtt_plot2D {

//===========================================================================//
/*!
  \struct SymbolProps

  \brief Symbol properties for Plot2D class.

  See Grace documentation for a detailed explanation of properties.
*/
//===========================================================================//
struct SymbolProps
{
    /// Symbol type
    int type;

    /// Color of symbol border
    int color;

    /// Size of symbol
    double size;

    /// Line width for border of symbol
    double width;

    /// Fill color of symbol 
    int fillColor;

    /// Pattern for filling symbol
    int fillPattern;

    /// Constructor, uses Grace defaults for a set.
    SymbolProps()
	: type(0)
	, color(1)
	, size(1.0)
	, width(1.0)
	, fillColor(1)
	, fillPattern(0) {}
};

} // namespace rtt_plot2D

#endif // INCLUDED_plot2D_SymbolProps_hh

//---------------------------------------------------------------------------//
// end of plot2D/SymbolProps.hh
//---------------------------------------------------------------------------//
