//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_analytic/pseudo_line.cc
 * \author Kent G. Budge
 * \date   Tue Apr  5 08:42:25 MDT 2011
 * \note   Copyright (C) 2011 Los Alamos National Security, LLC.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "ode/rkqs.i.hh"
#include "ode/quad.i.hh"
#include "Pseudo_Line_Base.hh"
#include "ds++/Packing_Utils.hh"
#include "ds++/cube.hh"
#include "c4/C4_Functions.hh"
#include "cdi/CDI.hh"


namespace rtt_cdi_analytic
{
using namespace std;
using namespace rtt_ode;
using namespace rtt_dsxx;
using namespace rtt_cdi;

//---------------------------------------------------------------------------//
Pseudo_Line_Base::Pseudo_Line_Base(SP<Expression const> const &continuum,
                                   unsigned number_of_lines,
                                   double line_peak,
                                   double line_width,
                                   unsigned number_of_edges,
                                   double edge_ratio,
                                   double Tref,
                                   double Tpow,
                                   double emin,
                                   double emax,
                                   unsigned seed)
    :
    continuum_(continuum),
    seed_(seed),
    number_of_lines_(number_of_lines),
    line_peak_(line_peak),
    line_width_(line_width),
    number_of_edges_(number_of_edges),
    edge_ratio_(edge_ratio),
    Tref_(Tref),
    Tpow_(Tpow),
    center_(number_of_lines),
    edge_(number_of_edges),
    edge_factor_(number_of_edges)
{
    Require(continuum!=SP<Expression>());
    Require(line_peak>=0.0);
    Require(line_width>=0.0);
    Require(edge_ratio>=0.0);
    Require(emin>=0.0);
    Require(emax>emin);
    // Require parameter (other than emin and emax) to be same on all processors
    
    srand(seed);

    // Get global range of energy

    rtt_c4::global_min(emin);
    rtt_c4::global_max(emax);

    for (unsigned i=0; i<number_of_lines; ++i)
    {
        center_[i] = (emax-emin)*static_cast<double>(rand())/RAND_MAX + emin;
    }

    // Sort line centers
    sort(center_.begin(), center_.end());

    for (unsigned i=0; i<number_of_edges; ++i)
    {
        edge_[i] = (emax-emin)*static_cast<double>(rand())/RAND_MAX + emin;
        edge_factor_[i] = edge_ratio_*(*continuum)(vector<double>(1,edge_[i]));
    }

    // Sort edges
    sort(edge_.begin(), edge_.end());
}

//---------------------------------------------------------------------------//
// Packing function

vector<char> Pseudo_Line_Base::pack() const 
{
    throw std::range_error("sorry, pack not implemented for Pseudo_Line_Base");
    // Because we haven't implemented packing functionality for Expression
    // trees yet.
    
#if 0
// caculate the size in bytes
    unsigned const size =
        3 * sizeof(double) + 3 * sizeof(int) + continuum_->packed_size();

    vector<char> pdata(size);

    // make a packer
    rtt_dsxx::Packer packer;

    // set the packer buffer
    packer.set_buffer(size, &pdata[0]);

	
    // pack the data
    continuum_->pack(packer);
    packer << seed_;
    packer << number_of_lines_;
    packer << line_peak_;
    packer << line_width_;
    packer << number_of_edges_;
    packer << edge_ratio_;

    // Check the size
    Ensure (packer.get_ptr() == &pdata[0] + size);
	
    return pdata;
#endif
}

//---------------------------------------------------------------------------//
double Pseudo_Line_Base::monoOpacity(double const x,
                                     double const T)
    const
{
    unsigned const number_of_lines = number_of_lines_;
    double const width = line_width_;
    double const peak = line_peak_;
    
    double Result = (*continuum_)(vector<double>(1,x));
    
    for (unsigned i=0; i<number_of_lines; ++i)
    {
        double const nu0 = center_[i];
        double const d = (x - nu0)/(width*nu0);
//        Result += peak*exp(-d*d);
        Result += peak/(1+d*d);
    }
    
    unsigned const number_of_edges = number_of_edges_;
    
    for (unsigned i=0; i<number_of_edges; ++i)
    {
        double const nu0 = edge_[i];
        if (x>=nu0)
        {
            Result += edge_factor_[i]*cube(nu0/x);
        }
    }
    if (Tpow_ != 0.0)
    {
        Result = Result * pow(T/Tref_, Tpow_);
    }
    return Result;
}

} // end namespace rtt_cdi_analytic

//---------------------------------------------------------------------------//
//                              end of Pseudo_Line_Base.cc
//---------------------------------------------------------------------------//
