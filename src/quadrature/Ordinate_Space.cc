//----------------------------------*-C++-*----------------------------------------------//
/*!
 * \file   quadrature/Ordinate_Space.cc
 * \author Kent Budge
 * \date   Mon Mar 26 16:11:19 2007
 * \brief  Define methods of class Ordinate_Space
 * \note   Copyright (C) 2006-2012 Los Alamos National Security, LLC.
 */
//---------------------------------------------------------------------------------------//
// $Id: Ordinate_Space.cc 6855 2012-11-06 16:39:27Z kellyt $
//---------------------------------------------------------------------------------------//

// Vendor software
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_legendre.h>

#include "Ordinate_Space.hh"

#include "special_functions/Ylm.hh"
#include "units/PhysicalConstants.hh"

using namespace rtt_units;

namespace rtt_quadrature
{
//---------------------------------------------------------------------------------------//
/*!
 * \brief Compute the Azimuthal angle for the current quadrature direction.
 */
double Ordinate_Space::compute_azimuthalAngle( double const mu,
                                               double const eta )
{
    using rtt_units::PI;
    using rtt_dsxx::soft_equiv;

    Require( std::abs(mu)  <= 1.0 );
    Require( std::abs(eta) <= 1.0 );

    if( soft_equiv( eta, 0.0 ) ) return PI;

//-------------------------------------------------------------

    double azimuthalAngle ( std::atan2( eta, mu) );

    // For axisymmetric cooridnates only, the azimuthal angle is on [0, 2\pi]
    // It is important to remember that the positive mu axis points to the
    // left and the positive eta axis points up, when the unit sphere is
    // projected on the plane of the mu- and eta-axis. In this case, phi is
    // measured from the mu-axis counterclockwise.
    //
    // This accounts for the fact that the aziumuthal angle is discretized
    // on levels of the xi-axis, making the computation of the azimuthal angle
    // here consistent with the discretization by using the eta and mu
    // ordinates to define phi.

    if (this->geometry() == rtt_mesh_element::AXISYMMETRIC && azimuthalAngle < 0.0)
        azimuthalAngle += 2*PI;
    
    return azimuthalAngle;
}

//---------------------------------------------------------------------------------------//
/*!
 * The computation of the tau and alpha coefficients is described by Morel in
 * various technical notes on the treatment of the angle derivatives in the
 * streaming operator.
 */
/* private */
void Ordinate_Space::compute_angle_operator_coefficients_()
{
    vector<Ordinate> const &ordinates = this->ordinates();
    unsigned         const number_of_ordinates = ordinates.size();
    rtt_mesh_element::Geometry const geometry = this->geometry();

    // Compute the ordinate derivative coefficients.
    
    // Default values are for the trivial case (Cartesian geometry).
    is_dependent_.resize(number_of_ordinates, false);
    alpha_.resize(number_of_ordinates, 0.0);
    tau_.resize(number_of_ordinates, 1.0);

    // We rely on OrdinateSet to have already sorted the ordinates and
    // inserted the starting ordinates for each level. We assume that the
    // starting ordinates are distinguished by zero quadrature weight.

    levels_.resize(number_of_ordinates);
    if (geometry==rtt_mesh_element::AXISYMMETRIC)
    {
        vector<double> C;

        double Csum = 0;
        int level = -1;
        for (unsigned a=0; a<number_of_ordinates; a++)
        {
            double const mu = ordinates[a].mu();
            double const wt = ordinates[a].wt();
            if (wt!=0 || (wt == 0 && mu > 0))
                // Not a starting ordinate.  Use Morel's recurrence relations
                // to determine the next ordinate derivative coefficient.
            {
                alpha_[a] = alpha_[a-1] + mu*wt;
                Csum += wt;
                levels_[a] = level;
                is_dependent_[a] = true;
            }
            else
                // A starting ordinate. Reinitialize the recurrence relation.
            {
                Check(a==0 || std::fabs(alpha_[a-1])<1.0e-15);
                // Be sure that the previous level (if any) had a final alpha
                // of zero, to within roundoff, as expected for the Morel
                // recursion formula.

                if (a>0)
                    alpha_[a-1] = 0.0;

                alpha_[a] = 0.0;

                if (mu < 0.0)
                {
                    // Save the normalization sum for the previous level, if any. 
                    if (level>=0)
                    {
                        Check(static_cast<int>(C.size())==level);
                        C.push_back(1.0/Csum);
                    }
                    level++;
                    Csum = 0.0;
                    
                    levels_[a] = level;
                    is_dependent_[a] = false;

                    if (level> 0) 
                        first_angles_.push_back(a-1);
                }
            }
        }
        // Save the normalization sum for the final level.
        Check(static_cast<int>(C.size())==level);
        C.push_back(1.0/Csum);
        first_angles_.push_back(number_of_ordinates-1);
        number_of_levels_ = C.size();

#if DBC & 2
        unsigned const dimension = this->dimension();
        if (dimension == 2)
            // Check that the level normalizations have the expected
            // properties. 
        {
            for (unsigned n=0; n<number_of_levels_/2; n++)
            {
                Check(C[n]>0.0);
                Check(C[number_of_levels_-1-n] > 0.0);
                Check(soft_equiv(C[n], C[number_of_levels_-1-n]));
            }
        }
#endif

        double mup = -2; // sentinel
        double sinth = -2; // sentinel
        double omp(0.0);
        level = -1;

        for (unsigned a=0; a<number_of_ordinates; a++)
        {
            double const eta = ordinates[a].eta();
            double const mu = ordinates[a].mu();
            double const wt = ordinates[a].wt();
            double const omm = omp;
            double const mum = mup;
            if (wt!=0)
                // Not a new level.  Apply Morel's recurrence relation.
            {
                omp = omm - rtt_units::PI*C[level]*wt;
                if (soft_equiv(omp, 0.0))
                {
                    omp = 0.0;
                }
            }
            else if (mu < 0.0)
                // New level.  Reinitialize the recurrence relation.
            {
                    omp = rtt_units::PI;
                    Check(1-eta*eta >= 0.0);
                    sinth = std::sqrt(1-eta*eta);
                    level++;
            }
            mup = sinth*std::cos(omp);
            if (wt!=0)
            {
                tau_[a] = (mu-mum)/(mup-mum);
                //tau_[a] = 0.5;                          // old school

                Check(tau_[a] >= 0.0 && tau_[a]<1.0);
            }
                        
        }
    }
    else if (geometry == rtt_mesh_element::SPHERICAL)
    {
        number_of_levels_ = 1;
        double norm(0);
        for (unsigned a=0; a < number_of_ordinates; ++a)
        {
            levels_[a] = 0;
            double const wt(ordinates[a].wt());
            norm += wt;
        }
        double const rnorm = 1.0/norm;

        for (unsigned a=0; a < number_of_ordinates; ++a)
        {
            double const mu(ordinates[a].mu());
            double const wt(ordinates[a].wt());

            if (wt!=0)
            {
                is_dependent_[a] = true;
                alpha_[a] = alpha_[a-1] + 2*wt*mu;
            }
            else 
            {
                alpha_[a] = 0;

                if (mu < 0.0)
                {
                    is_dependent_[a] = false;
                    first_angles_.push_back(number_of_ordinates-1);
                }
                else
                {
                    is_dependent_[a] = true;
                }
                        
            }
        }

        double mup = -2; // sentinel
        for (unsigned a=0; a < number_of_ordinates; ++a)
        {
            double const mu(ordinates[a].mu());
            double const wt(ordinates[a].wt());

            double const mum = mup;

            if (wt !=0)
                mup = mum + 2*wt*rnorm;
            else 
                mup = mu;

            if (wt !=0)
            {
                tau_[a] = (mu-mum)/(2*wt*rnorm);

                Check(tau_[a]>0.0 && tau_[a]<=1.0);
            }
        }
    }
    else
    {
        Check(geometry == rtt_mesh_element::CARTESIAN);
        number_of_levels_ = 0;
    }

    Insist(first_angles_.size() == number_of_levels_, "unexpected starting direction reflection index");
}

//---------------------------------------------------------------------------------------//
vector<Moment> Ordinate_Space::compute_n2lk_(Quadrature_Class const quadrature_class,
                                             unsigned const sn_order)
{
    unsigned const dim = dimension();
    Geometry const geometry = this->geometry();
    
    if( dim == 3 )
    {
        return compute_n2lk_3D_(quadrature_class,
                                sn_order);
    }
    else if( dim == 2)
    {
        if (geometry == rtt_mesh_element::AXISYMMETRIC)
            return compute_n2lk_2Da_(quadrature_class,
                                     sn_order);
        else
            return compute_n2lk_2D_(quadrature_class,
                                    sn_order);
    }
    else
    {
        Check( dim == 1 );
        if (geometry == rtt_mesh_element::AXISYMMETRIC)
            return compute_n2lk_1Da_(quadrature_class,
                                     sn_order);
        else
            return compute_n2lk_1D_(quadrature_class,
                                    sn_order);
    }
}

//---------------------------------------------------------------------------------------//
/*! Compute the description of the moment space.
 *
 * N.B. This must not be called in the Ordinate_Space constructor, but in the
 * child class constructor, because it uses virtual functions of the child
 * class that are not set up until the child class is constructed.
 */

void Ordinate_Space::compute_moments_(Quadrature_Class const quadrature_class,
                                      unsigned const sn_order)
{
    unsigned const Lmax = expansion_order_;
    
    moments_ = compute_n2lk_(quadrature_class,
                             sn_order);
            
    moments_per_order_.resize(Lmax+1, 0U);
    number_of_moments_ = 0;
    for(unsigned n=0; n<moments_.size(); ++n)
    {
        unsigned const l = moments_[n].L();
        if (l<=Lmax)
        {
            moments_per_order_[l] += 1;
            number_of_moments_++;
        }
    }
}

//---------------------------------------------------------------------------------------//
/*!
 *
 * \param dimension Dimension of the physical problem space (1, 2, or 3)
 *
 * \param geometry Geometry of the physical problem space (spherical,
 * axisymmetric, Cartesian)
 *
 * \param ordinates Set of ordinate directions
 *
 * \param expansion_order Expansion order of the desired scattering moment
 * space.
 *
 * \param extra_starting_directions Add extra directions to each level set. In most
 * geometries, an additional ordinate is added that is opposite in direction
 * to the starting direction. This is used to implement reflection exactly in
 * curvilinear coordinates. In 1D spherical, that means an additional angle is
 * added at mu=1. In axisymmetric, that means additional angles are added that
 * are oriented opposite to the incoming starting direction on each level.
 *
 * \param ordering Ordering into which to sort the ordinates.
 */

Ordinate_Space::Ordinate_Space( unsigned const  dimension,
                                Geometry const  geometry,
                                vector<Ordinate> const &ordinates,
                                unsigned const  expansion_order,
                                bool const  extra_starting_directions,
                                Ordering const ordering)
    : Ordinate_Set(dimension,
                   geometry,
                   ordinates,
                   true, // include starting directions
                   extra_starting_directions,
                   ordering),
      expansion_order_(expansion_order),
      has_extra_starting_directions_(extra_starting_directions)
{
    Require(dimension>0 && dimension<4);
    Require(geometry!=rtt_mesh_element::END_GEOMETRY);

    compute_angle_operator_coefficients_();
    
    compute_reflection_maps_();

    Ensure(check_class_invariants());
    Ensure(has_extra_starting_directions()==extra_starting_directions);
}

//---------------------------------------------------------------------------------------//
/*!
 * The psi coefficient is used to compute the self term in the angle
 * derivative term of the streaming operator.
 */
double Ordinate_Space::psi_coefficient(unsigned const a) const
{
    Require(is_dependent(a));

    double const wt = ordinates()[a].wt();
    double const alpha_a = alpha_[a];
    double const tau_a = tau_[a];
    double const Result = alpha_a/(wt*tau_a);
    return Result;
}

//---------------------------------------------------------------------------------------//
/*!
 * The source coefficient is used to compute the previous midpoint angle term
 * in the angle derivative term of the streaming operator.
 */
double Ordinate_Space::source_coefficient(unsigned const a) const
{
    Require(is_dependent(a));

    double const wt = ordinates()[a].wt();
    double const alpha_a = alpha_[a];
    double const alpha_am1 = alpha_[a-1];
    double const tau_a = tau_[a];
    double const Result = (alpha_a*(1-tau_a)/tau_a + alpha_am1)/wt;
    return Result;
}
//---------------------------------------------------------------------------------------//
/*!
 * The bookkeeping coefficient is used to compute the next midpoint angle
 * specific intensity.
 */
double Ordinate_Space::bookkeeping_coefficient(unsigned const a) const
{
    Require(is_dependent(a));

    double const tau_a = tau_[a];
    double const Result = 1.0/tau_a;

    Require(Result>0.0);
    return Result;
}

//---------------------------------------------------------------------------------------//
bool Ordinate_Space::check_class_invariants() const
{
    if (geometry() == rtt_mesh_element::CARTESIAN)
    {
        // There are no levels in Cartesian ordinate sets.
        return number_of_levels_ == 0;
    }
    else
    {
        vector<Ordinate> const &ordinates = this->ordinates();
        unsigned const number_of_ordinates = ordinates.size();

        // Check that the number of levels is correct.
        unsigned levels = 0;
        for (unsigned a=0; a<number_of_ordinates; ++a)
        {
            if ((ordinates[a].wt() == 0) && (ordinates[a].mu() < 0.0))
            {
                ++levels;
            }
        }
        if (number_of_levels_<levels)
            return false;

        // Check that the angle derivative coefficient arrays have the correct
        // size.
        return
            is_dependent_.size()==number_of_ordinates &&
            alpha_.size()==number_of_ordinates &&
            tau_.size()==number_of_ordinates;
    }
}

//---------------------------------------------------------------------------------------//
void Ordinate_Space::compute_reflection_maps_()
{
    vector<Ordinate> const &ordinates = this->ordinates();
    unsigned const number_of_ordinates = ordinates.size();
    
    reflect_mu_.resize(number_of_ordinates);
    reflect_eta_.resize(number_of_ordinates);
    reflect_xi_.resize(number_of_ordinates);

    // Since the ordinate set will likely never number more than a few
    // hundred, we go ahead and do the simpleminded quadratic search to match
    // the ordinates up.

    for (unsigned a=0; a+1<number_of_ordinates; ++a)
    {
        for (unsigned ap=a+1; ap<number_of_ordinates; ++ap)
        {
            if (soft_equiv(ordinates[a].mu(), -ordinates[ap].mu()) &&
                soft_equiv(ordinates[a].eta(), ordinates[ap].eta()) &&
                soft_equiv(ordinates[a].xi(), ordinates[ap].xi()))
            {
                reflect_mu_[a] = ap;
                reflect_mu_[ap] = a;
            }
            else if (soft_equiv(ordinates[a].mu(), ordinates[ap].mu()) &&
                soft_equiv(ordinates[a].eta(), -ordinates[ap].eta()) &&
                soft_equiv(ordinates[a].xi(), ordinates[ap].xi()))
            {
                reflect_eta_[a] = ap;
                reflect_eta_[ap] = a;
            }
            else if (soft_equiv(ordinates[a].mu(), ordinates[ap].mu()) &&
                soft_equiv(ordinates[a].eta(), ordinates[ap].eta()) &&
                soft_equiv(ordinates[a].xi(), -ordinates[ap].xi()))
            {
                reflect_xi_[a] = ap;
                reflect_xi_[ap] = a;
            }
        }
    }
}

} // end namespace rtt_quadrature

//---------------------------------------------------------------------------------------//
//                 end of Ordinate_Space.cc
//---------------------------------------------------------------------------------------//
