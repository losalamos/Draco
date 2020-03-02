//----------------------------------*-C++-*-----------------------------------//
/*!
  \file   Colormap.hh
  \author lowrie
  \date   2002-04-12
  \brief  Header for Colormap.
  \note   Copyright (C) 2016-2020 Triad National Security, LLC.
          All rights reserved.
*/
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//

#ifndef rtt_plot2D_Colormap_hh
#define rtt_plot2D_Colormap_hh

namespace rtt_plot2D {

//============================================================================//
/*!
  \enum Colormap

  \brief Defines the colormap for Grace.

  See GRACE_INSTALLATION_DIR/templates/Default.agr for default colormap.
*/
//============================================================================//
enum Colormap {
  COLOR_WHITE,
  COLOR_BLACK,
  COLOR_RED,
  COLOR_GREEN,
  COLOR_BLUE,
  COLOR_YELLOW,
  COLOR_BROWN,
  COLOR_GREY,
  COLOR_VIOLET,
  COLOR_CYAN,
  COLOR_MAGENTA,
  COLOR_ORANGE,
  COLOR_INDIGO,
  COLOR_MAROON,
  COLOR_TURQUOISE,
  COLOR_GREEN4
};

} // namespace rtt_plot2D

#endif // rtt_plot2D_Colormap_hh

//----------------------------------------------------------------------------//
// end of plot2D/Colormap.hh
//----------------------------------------------------------------------------//
