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
#include "ds++/cube.hh"

namespace
{

double BB(double const T, double const x)
{
    double const e = expm1(x/T);
    if (e>0.0)
        return x*x*x/e;
    else
        return x*x*T;
}

double DBB(double const T, double const x)
{
    double const e = expm1(x/T);
    if (e>0)
    {
        double const de = -exp(x/T)*x/(T*T);
        return -x*x*x*de/(e*e);
    }
    else
    {
        return x*x;
    }
}

} // end anonymous namespace


namespace rtt_cdi_analytic
{
using namespace std;
using namespace rtt_ode;
using namespace rtt_dsxx;

typedef Analytic_MultigroupOpacity::sf_double sf_double;
typedef Analytic_MultigroupOpacity::vf_double vf_double;

//---------------------------------------------------------------------------//
class PLP_Functor
{
  public:

    typedef double return_type;
    
    PLP_Functor(Pseudo_Line_Analytic_MultigroupOpacity const *ptr,
                double const T)
        : ptr_(ptr), T_(T)
    {}

    double operator()(double x);

  private:
    
    Pseudo_Line_Analytic_MultigroupOpacity const *ptr_;
    double T_;
};

double PLP_Functor::operator()(double x)
{
    return ptr_->monoOpacity(x, T_)*BB(T_, x);
}

//---------------------------------------------------------------------------//
class PLPW_Functor
{
  public:

    typedef double return_type;
    
    PLPW_Functor(double const T)
        : T_(T)
    {}

    double operator()(double x);

  private:
    
    double T_;
    
};

double PLPW_Functor::operator()(double x)
{
    return BB(T_, x);
}

//---------------------------------------------------------------------------//
class PLR_Functor
{
  public:

    typedef double return_type;
    
    PLR_Functor(Pseudo_Line_Analytic_MultigroupOpacity const *ptr,
                double const T)
        : ptr_(ptr), T_(T)
    {}

    double operator()(double x);

  private:
    
    Pseudo_Line_Analytic_MultigroupOpacity const *ptr_;
    double T_;
};

double PLR_Functor::operator()(double x)
{
    return DBB(T_,x)/ptr_->monoOpacity(x, T_);
}

//---------------------------------------------------------------------------//
class PLRW_Functor
{
  public:

    typedef double return_type;
    
    PLRW_Functor(double const T)
        : T_(T)
    {}

    double operator()(double x);

  private:
    
    double T_;
    
};

double PLRW_Functor::operator()(double x)
{
    return DBB(T_,x);
}

//---------------------------------------------------------------------------//
Pseudo_Line_Analytic_MultigroupOpacity::
Pseudo_Line_Analytic_MultigroupOpacity(sf_double const &group_bounds,
                                       rtt_cdi::Reaction const reaction,
                                       SP<Expression const> const &continuum,
                                       unsigned number_of_lines,
                                       double line_peak,
                                       double line_width,
                                       unsigned number_of_edges,
                                       double edge_ratio,
                                       double Tref,
                                       double Tpow,
                                       double emin,
                                       double emax,
                                       Averaging const averaging,
                                       unsigned const qpoints,
                                       unsigned seed)
    :
    Analytic_MultigroupOpacity(group_bounds,
                               reaction),
    continuum_(continuum),
    seed_(seed),
    number_of_lines_(number_of_lines),
    line_peak_(line_peak),
    line_width_(line_width),
    number_of_edges_(number_of_edges),
    edge_ratio_(edge_ratio),
    averaging_(averaging),
    Tref_(Tref),
    Tpow_(Tpow),
    center_(number_of_lines),
    edge_(number_of_edges),
    edge_factor_(number_of_edges),
    qpoints_(qpoints)
{
    Require(continuum>=0.0);
    Require(line_peak>=0.0);
    Require(line_width>=0.0);
    Require(edge_ratio>=0.0);
    Require(emin>=0.0);
    Require(emax>emin);
    
    srand(seed);

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
    packer << line_peak_;
    packer << line_width_;
    packer << number_of_edges_;
    packer << edge_ratio_;
    packer << averaging_;

    // Check the size
    Ensure (packer.get_ptr() == &pdata[base_size] + size);
	
    return pdata;
}

//---------------------------------------------------------------------------//
sf_double
Pseudo_Line_Analytic_MultigroupOpacity::getOpacity( double T,
                                                    double /*rho*/) const
{
    sf_double const &group_bounds = this->getGroupBoundaries();
    unsigned const number_of_groups = group_bounds.size()-1U;
    sf_double Result(number_of_groups, 0.0);

    switch (averaging_)
    {
        case NONE:
            {
                double g1 = group_bounds[0];
                for (unsigned g=0; g<number_of_groups; ++g)
                {
                    double const g0 = g1;
                    g1 = group_bounds[g+1];
                    double const nu = 0.5*(g0 + g1);
                    Result[g] = monoOpacity(nu, T);
                }
            }
            break;
            
        case ROSSELAND:
        {
            PLR_Functor rfunctor(this, T);
            PLRW_Functor wfunctor(T);
            
            double g1 = group_bounds[0];
            for (unsigned g=0; g<number_of_groups; ++g)
            {
                double const g0 = g1;
                g1 = group_bounds[g+1];
                double eps = 1e-5;
                double t = 0.0, b = 0.0;
                if (qpoints_==0)
                {
                    double x1 = g0;
                    while (x1<g1)
                    {
                        double x0 = x1;
                        x1 += 2*line_width_;
                        if (x1>g1) x1 = g1;
                        
                        t +=
                            rtt_ode::quad(rfunctor,
                                          x0,
                                          x1,
                                          eps,
                                          rkqs<double, Quad_To_ODE<PLR_Functor> >);
                        b +=
                            rtt_ode::quad(wfunctor,
                                          x0,
                                          x1,
                                          eps,
                                          rkqs<double, Quad_To_ODE<PLRW_Functor> >);
                    }
                }
                else
                {
                    for (unsigned ig=0; ig<qpoints_; ++ig)
                    {
                        double const x = (ig+0.5)*(g1-g0)/qpoints_ + g0;

                        t += rfunctor(x)*(g1-g0);
                        b += wfunctor(x)*(g1-g0);
                    }
                }

                Result[g] = b/t;
            }
        }
        break;
           
        case PLANCK:
        {
            double g1 = group_bounds[0];
            for (unsigned g=0; g<number_of_groups; ++g)
            {
                PLP_Functor pfunctor(this, T);
                PLPW_Functor wfunctor(T);
                
                double const g0 = g1;
                g1 = group_bounds[g+1];
                double eps = 1e-5;
                double t = 0.0, b = 0.0;
                if (qpoints_==0)
                {
                    double x1 = g0;
                    while (x1<g1)
                    {
                        double x0 = x1;
                        x1 += 2*line_width_;
                        if (x1>g1) x1 = g1;
                        
                        t +=
                            rtt_ode::quad(pfunctor,
                                          x0,
                                          x1,
                                          eps,
                                          rkqs<double, Quad_To_ODE<PLP_Functor> >);
                        b +=
                            rtt_ode::quad(wfunctor,
                                          x0,
                                          x1,
                                          eps,
                                          rkqs<double, Quad_To_ODE<PLPW_Functor> >);
                    }
                }
                else
                {
                    for (unsigned ig=0; ig<qpoints_; ++ig)
                    {
                        double const x = (ig+0.5)*(g1-g0)/qpoints_ + g0;
                        
                        t += pfunctor(x)*(g1-g0);
                        b += wfunctor(x)*(g1-g0);
                    }
                }
                
                Result[g] = t/b;
            }
        }
        break;
        
        default:
            Insist(false, "bad case");            
    }

    return Result;
}

//---------------------------------------------------------------------------//
vf_double
Pseudo_Line_Analytic_MultigroupOpacity::getOpacity(sf_double const &T,
                                                   double rho) const
{
    unsigned const n = T.size();
    vf_double Result(n);

    for (unsigned i=0; i<n; ++i)
    {
        double const Ti = T[i];
        Result[i] = getOpacity(Ti, rho);
    }
    return Result;
}

//---------------------------------------------------------------------------//
vf_double
Pseudo_Line_Analytic_MultigroupOpacity::getOpacity(double    const T,
                                                   sf_double const &rho) const
{
    return vf_double(rho.size(), getOpacity(T, rho[0]));
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

//---------------------------------------------------------------------------//
double Pseudo_Line_Analytic_MultigroupOpacity::monoOpacity(double const x,
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
        double const d = x - nu0;
        Result += peak*exp(-d*d/(width*width*nu0*nu0));
        //       if (d<10*width*nu0) break; // can't be significant beyond this point
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
//                              end of Pseudo_Line_Analytic_MultigroupOpacity.cc
//---------------------------------------------------------------------------//
