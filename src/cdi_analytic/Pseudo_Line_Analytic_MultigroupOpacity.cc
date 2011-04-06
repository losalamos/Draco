//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_analytic/Pseudo_Line_Analytic_MultigroupOpacity.cc
 * \author Kent G. Budge
 * \date   Tue Apr  5 08:42:25 MDT 2011
 * \note   Copyright (C) 2011 Los Alamos National Security, LLC.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "ode/rkqs.i.hh"
#include "ode/quad.i.hh"
#include "Pseudo_Line_Analytic_MultigroupOpacity.hh"
#include "ds++/Packing_Utils.hh"


namespace rtt_cdi_analytic
{
using namespace std;
using namespace rtt_ode;

typedef Analytic_MultigroupOpacity::sf_double sf_double;
typedef Analytic_MultigroupOpacity::vf_double vf_double;

//---------------------------------------------------------------------------//
class PLP_Functor
{
  public:

    typedef double return_type;
    
    PLP_Functor(Pseudo_Line_Analytic_MultigroupOpacity const *ptr)
        : ptr_(ptr)
    {}

    double operator()(double &x);

  private:
    
    Pseudo_Line_Analytic_MultigroupOpacity const *ptr_;
    
};

double PLP_Functor::operator()(double &x)
{
    unsigned const number_of_lines = ptr_->number_of_lines_;
    vector<double> const &center = ptr_->center_;
    double const width = ptr_->width_;
    double const peak = ptr_->peak_;

    double Result = ptr_->continuum_;
    
    for (unsigned i=0; i<number_of_lines; ++i)
    {
        double const nu0 = center[i];
        double const d = x - nu0;
        Result += peak*exp(-d*d/(width*width*nu0*nu0));
    }
    return Result;
}

//---------------------------------------------------------------------------//
class PLR_Functor
{
  public:

    typedef double return_type;
    
    PLR_Functor(Pseudo_Line_Analytic_MultigroupOpacity const *ptr)
        : ptr_(ptr)
    {}

    double operator()(double &x);

  private:
    
    Pseudo_Line_Analytic_MultigroupOpacity const *ptr_;
    
};

double PLR_Functor::operator()(double &x)
{
    unsigned const number_of_lines = ptr_->number_of_lines_;
    vector<double> const &center = ptr_->center_;
    double const width = ptr_->width_;
    double const peak = ptr_->peak_;

    double Result = ptr_->continuum_;
    
    for (unsigned i=0; i<number_of_lines; ++i)
    {
        double const nu0 = center[i];
        double const d = x - nu0;
        Result += peak*exp(-d*d/(width*width*nu0*nu0));
    }
    return 1.0/Result;
}

//---------------------------------------------------------------------------//
Pseudo_Line_Analytic_MultigroupOpacity::
Pseudo_Line_Analytic_MultigroupOpacity(sf_double const &group_bounds,
                                       rtt_cdi::Reaction const reaction,
                                       double continuum,
                                       unsigned number_of_lines,
                                       double peak,
                                       double width,
                                       double emin,
                                       double emax,
                                       Averaging const averaging,
                                       unsigned seed)
    :
    Analytic_MultigroupOpacity(group_bounds,
                               reaction),
    continuum_(continuum),
    seed_(seed),
    number_of_lines_(number_of_lines),
    peak_(peak),
    width_(width),
    averaging_(averaging),
    center_(number_of_lines)
{
    Require(continuum>=0.0);
    Require(peak>=0.0);
    Require(width>=0.0);
    Require(emin>=0.0);
    Require(emax>emin);
    
    srand(seed);

    for (unsigned i=0; i<number_of_lines; ++i)
    {
        center_[i] = (emax-emin)*static_cast<double>(rand())/RAND_MAX + emin;
    }

    // Precompute mass absorption coefficient
    
    unsigned const number_of_groups = group_bounds.size()-1U;
    sigma_.resize(number_of_groups, continuum_);

    switch (averaging)
    {
        case NONE:
            for (unsigned i=0; i<number_of_lines_; ++i)
            {
                double const nu0 = center_[i];
                double g1 = group_bounds[0];
                for (unsigned g=0; g<number_of_groups; ++g)
                {
                    double const g0 = g1;
                    g1 = group_bounds[g+1];
                    double const nu = 0.5*(g0 + g1);
                    double const d = nu - nu0;
                    sigma_[g] += peak_*exp(-d*d/(width_*width_*nu0*nu0));
                }
            }
            break;
            
        case ROSSELAND:
        {
            double g1 = group_bounds[0];
            for (unsigned g=0; g<number_of_groups; ++g)
            {
                double const g0 = g1;
                g1 = group_bounds[g+1];
                double eps = 1e-5;
                sigma_[g] = (g1-g0)/(
                    rtt_ode::quad(PLR_Functor(this),
                                  g0,
                                  g1,
                                  eps,
                                  rkqs<double, Quad_To_ODE<PLR_Functor> >));
            }
        }
        break;
           
        case PLANCK:
        {
            double g1 = group_bounds[0];
            for (unsigned g=0; g<number_of_groups; ++g)
            {
                double const g0 = g1;
                g1 = group_bounds[g+1];
                double eps = 1e-5;
                sigma_[g] =
                    rtt_ode::quad(PLP_Functor(this),
                                  g0,
                                  g1,
                                  eps,
                                  rkqs<double, Quad_To_ODE<PLP_Functor> >)
                    /(g1-g0);
            }
        }
        break;
        
        default:
            Insist(false, "bad case");            
    }
}

//---------------------------------------------------------------------------//
// Packing function

Analytic_MultigroupOpacity::sf_char
Pseudo_Line_Analytic_MultigroupOpacity::pack() const 
{
    sf_char pdata = Analytic_MultigroupOpacity::pack();
    
    // caculate the size in bytes
    unsigned const base_size = pdata.size();
    unsigned const size =
        3 * sizeof(double) + 2 * sizeof(int) + sizeof(Averaging);

    pdata.resize(size+base_size);

    // make a packer
    rtt_dsxx::Packer packer;

    // set the packer buffer
    packer.set_buffer(size, &pdata[base_size]);

	
    // pack the data
    packer << continuum_;
    packer << seed_;
    packer << number_of_lines_;
    packer << peak_;
    packer << width_;
    packer << averaging_;

    // Check the size
    Ensure (packer.get_ptr() == &pdata[base_size] + size);
	
    return pdata;
}

//---------------------------------------------------------------------------//
sf_double
Pseudo_Line_Analytic_MultigroupOpacity::getOpacity( double /*T*/,
                                                    double /*rho*/) const
{
    return sigma_;
}

//---------------------------------------------------------------------------//
vf_double
Pseudo_Line_Analytic_MultigroupOpacity::getOpacity(sf_double const &T,
                                                   double /*rho*/) const
{
    return vf_double(T.size(), sigma_);
}

//---------------------------------------------------------------------------//
vf_double
Pseudo_Line_Analytic_MultigroupOpacity::getOpacity(double    const /*T*/,
                                                   sf_double const &rho) const
{
    return vf_double(rho.size(), sigma_);
}

//---------------------------------------------------------------------------//
Pseudo_Line_Analytic_MultigroupOpacity::std_string
Pseudo_Line_Analytic_MultigroupOpacity::getDataDescriptor() const
{
    std_string descriptor;

    rtt_cdi::Reaction const reaction = getReactionType();
    
    if (reaction == rtt_cdi::TOTAL)
        descriptor = "Pseudo Line Multigroup Total";
    else if (reaction == rtt_cdi::ABSORPTION)
        descriptor = "Pseudo Line Multigroup Absorption";
    else if (reaction == rtt_cdi::SCATTERING)
        descriptor = "Pseudo Line Multigroup Scattering";
    else
    {
        Insist (0, "Invalid nGray multigroup model opacity!");
    }

    return descriptor;
}

} // end namespace rtt_cdi_analytic

//---------------------------------------------------------------------------//
//                              end of Pseudo_Line_Analytic_MultigroupOpacity.cc
//---------------------------------------------------------------------------//
