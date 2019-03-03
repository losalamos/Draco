//----------------------------------*-C++-*----------------------------------//
/*!
  \file   LineProps.hh
  \author lowrie
  \date   2002-04-12
  \brief  Header for LineProps.
  \note   Copyright (C) 2016-2019 Triad National Security, LLC.
          All rights reserved.
*/
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//

#ifndef rtt_plot2D_LineProps_hh
#define rtt_plot2D_LineProps_hh

#include "Colormap.hh"

namespace rtt_plot2D {

//===========================================================================//
/*!
  \struct LineProps

  \brief Line properties for Plot2D class.

  See Grace documentation for a detailed explanation of properties.
*/
//===========================================================================//
struct LineProps {
  /// Various line styles
  enum Style {
    STYLE_NONE,
    STYLE_SOLID,
    STYLE_DOT,
    STYLE_DASH,
    STYLE_LONGDASH,
    STYLE_DASHDOT,
    STYLE_LONGDASHDOT,
    STYLE_DASHDOTDOT,
    STYLE_DASHDASHDOT
  };

  /// Line style
  Style style;

  /// Line color
  Colormap color;

  /// Width of line
  double width;

  /// Constructor, uses Grace defaults for a set.
  LineProps() : style(STYLE_SOLID), color(COLOR_BLACK), width(1.0) {}
};

} // namespace rtt_plot2D

#endif // rtt_plot2D_LineProps_hh

//---------------------------------------------------------------------------//
// end of plot2D/LineProps.hh
//---------------------------------------------------------------------------//
