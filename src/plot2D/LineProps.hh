//----------------------------------*-C++-*----------------------------------//
/*!
  \file   LineProps.hh
  \author lowrie
  \date   2002-04-12
  \brief  Header for LineProps.
*/
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef INCLUDED_plot2D_LineProps_hh
#define INCLUDED_plot2D_LineProps_hh

namespace rtt_plot2D {

//===========================================================================//
/*!
  \struct LineProps

  \brief Line properties for Plot2D class.

  See Grace documentation for a detailed explanation of properties.
*/
//===========================================================================//
struct LineProps
{
    /// Line type
    int type;

    /// Line color
    int color;

    /// Width of line
    double width;

    /// Constructor, uses Grace defaults for a set.
    LineProps()
	: type(1)
	, color(1)
	, width(1.0) {}
    
};

} // namespace rtt_plot2D

#endif // INCLUDED_plot2D_LineProps_hh

//---------------------------------------------------------------------------//
// end of plot2D/LineProps.hh
//---------------------------------------------------------------------------//
