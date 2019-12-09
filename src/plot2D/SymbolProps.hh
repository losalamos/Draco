//----------------------------------*-C++-*----------------------------------//
/*!
  \file   SymbolProps.hh
  \author lowrie
  \date   2002-04-12
  \brief  Header for SymbolProps.
  \note   Copyright (C) 2016-2019 Triad National Security, LLC.
          All rights reserved.
*/
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//

#ifndef rtt_plot2D_SymbolProps_hh
#define rtt_plot2D_SymbolProps_hh

#include "Colormap.hh"

namespace rtt_plot2D {

//===========================================================================//
/*!
  \struct SymbolProps

  \brief Symbol properties for Plot2D class.

  See Grace documentation for a detailed explanation of properties.
*/
//===========================================================================//
struct SymbolProps {
  /// Shapes available
  enum Shape {
    SHAPE_NONE,
    SHAPE_CIRCLE,
    SHAPE_SQUARE,
    SHAPE_DIAMOND,
    SHAPE_TRIANGLEUP,
    SHAPE_TRIANGLELEFT,
    SHAPE_TRIANGLEDOWN,
    SHAPE_TRIANGLERIGHT,
    SHAPE_PLUS,
    SHAPE_X,
    SHAPE_STAR,
    SHAPE_CHAR
  };

  /// Symbol shape
  Shape shape;

  /// Color of symbol border
  Colormap color;

  /// Size of symbol
  double size;

  /// Line width for border of symbol
  double width;

  /// Fill color of symbol
  Colormap fillColor;

  /// Pattern for filling symbol. 0 is none, 1 solid, ...
  int fillPattern;

  /// Constructor, uses Grace defaults for a set.
  SymbolProps()
      : shape(SHAPE_NONE), color(COLOR_BLACK), size(1.0), width(1.0),
        fillColor(COLOR_BLACK), fillPattern(0) {}
};

} // namespace rtt_plot2D

#endif // rtt_plot2D_SymbolProps_hh

//---------------------------------------------------------------------------//
// end of plot2D/SymbolProps.hh
//---------------------------------------------------------------------------//
