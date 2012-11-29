//----------------------------------*-C++-*----------------------------------------------//
/*!
 * \file   quadrature/General_Octant_Quadrature.cc
 * \author Kelly Thompson
 * \date   Wed Sep  1 10:19:52 2004
 * \brief  
 * \note   Copyright 2004 The Regents of the University of California.
 */
//---------------------------------------------------------------------------------------//
// $Id: General_Octant_Quadrature.cc 6718 2012-08-30 20:03:01Z warsa $
//---------------------------------------------------------------------------------------//

#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>

#include "General_Octant_Quadrature.hh"

#include "ds++/Soft_Equivalence.hh"
#include "ds++/to_string.hh"
#include "units/PhysicalConstants.hh"

namespace rtt_quadrature
{
using namespace rtt_dsxx;

//---------------------------------------------------------------------------------------//
string General_Octant_Quadrature::name() const { return "General Octant Quadrature"; }

//---------------------------------------------------------------------------------------//
string General_Octant_Quadrature::parse_name()  const
{ return "general octant quadrature"; }
    
//---------------------------------------------------------------------------------------//
unsigned General_Octant_Quadrature::number_of_levels() const { return number_of_levels_; }
    
//---------------------------------------------------------------------------------------//
string General_Octant_Quadrature::as_text(string const &indent) const
{
    string Result = indent + "  type = general octant quadrature";
    Result += indent + "  number of ordinates = " + to_string(mu_.size());
    Result += indent + "  number of levels = " + to_string(number_of_levels_);

    Result += indent + "  interpolation algorithm = ";
    switch(qim())
    {
        case SN:
            Result += indent + "SN";
            break;
            
        case GQ:
            Result += indent + "GALERKIN";
            break;
            
        default:
            Insist(false, "bad case");
    }
 
    unsigned const N = mu_.size();
    for (unsigned i=0; i<N; ++i)
    {
        Result += indent + "  " + to_string(mu_[i]);
        Result += "  " + to_string(eta_[i]);
        Result += "  " + to_string(xi_[i]); 
        Result += "  " + to_string(wt_[i]);
    }

    Result += indent + "end";
    
    return Result;
}

//---------------------------------------------------------------------------------------//
/*virtual*/
void General_Octant_Quadrature::create_octant_ordinates_(vector<double> &mu,
                                                  vector<double> &eta,
                                                  vector<double> &wt) const
{
    mu = mu_;
    eta = eta_;
    wt = wt_;
}

} // end namespace rtt_quadrature

//---------------------------------------------------------------------------------------//
//                 end of General_Octant_Quadrature.cc
//---------------------------------------------------------------------------------------//
