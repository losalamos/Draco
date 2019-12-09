//----------------------------------*-C++-*----------------------------------//
/*! \file   MathConstants.hh
 *  \author Kelly Thompson
 *  \brief  Provide a single place where physical constants (pi, speed of
 *          light, etc) are defined.
 *  \date   Fri Nov 07 10:04:52 2003
 *  \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *          All rights reserved. */
//---------------------------------------------------------------------------//

#ifndef __units_MathConstants_hh__
#define __units_MathConstants_hh__

//! \namespace rtt_units Namespace for units and physical constants
namespace rtt_units {
// Mathematical constants

// N.B. M_PI and M_E, though widely implemented as part of <cmath>, are not
// standard C++.

//! pi the ratio of a circle's circumference to its diameter (dimensionless)
static double constexpr PI = 3.141592653589793238462643383279;

// Euler's number (dimensionless)
static double constexpr N_EULER = 2.7182818284590452353602874;

} // end namespace rtt_units

#endif // __units_MathConstants_hh__

//---------------------------------------------------------------------------//
// end of units/MathConstants.hh
//---------------------------------------------------------------------------//
