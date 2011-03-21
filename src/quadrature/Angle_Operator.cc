//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/Angle_Operator.cc
 * \author Kent Budge
 * \date   Mon Mar 26 16:11:19 2007
 * \brief  Define methods of class Angle_Operator
 * \note   Copyright (C) 2006-2010 Los Alamos National Security, LLC.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include "ds++/Soft_Equivalence.hh"
#include "units/PhysicalConstants.hh"
#include "Angle_Operator.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_quadrature;
using namespace rtt_units;

namespace rtt_quadrature
{
//---------------------------------------------------------------------------//
/*!
 *
 * The computation of the tau and alpha coefficients is described by Morel in
 * various technical notes on the treatment of the angle derivatives in the
 * streaming operator.
 *
 * \param quadrature Sn quadrature set.
 *
 * \param geometry Geometry of the physical problem space.
 *
 * \param dimension Dimension of the physical problem space (1, 2, or 3)
 */

Angle_Operator::Angle_Operator( SP<Quadrature const>       const & quadrature,
                                rtt_mesh_element::Geometry const   geometry,
                                unsigned                   const   dimension)
    : OrdinateSet(quadrature, geometry, dimension),
      number_of_levels_(0),
      levels_(), is_dependent_(), alpha_(), tau_()
{
    Require(quadrature!=SP<Quadrature>());
    Require(dimension>0 && dimension<4);
    Require(geometry<rtt_mesh_element::END_GEOMETRY);

    vector<Ordinate> const &ordinates = getOrdinates();
    unsigned const number_of_ordinates = ordinates.size();

    // Compute the ordinate derivative coefficients.
    
    // Default values are for the trivial case (Cartesian geometry).
    is_dependent_.resize(number_of_ordinates, false);
    alpha_.resize(number_of_ordinates, 0.0);
    tau_.resize(number_of_ordinates, 1.0);

    // We rely on OrdinateSet to have already sorted the ordinates and
    // inserted the starting ordinates for each level. We assume that the
    // starting ordinates are distinguished by zero quadrature weight.

    levels_.resize(number_of_ordinates);
    number_of_levels_ = 0;
    if (geometry==rtt_mesh_element::AXISYMMETRIC)
    {
        number_of_levels_ = quadrature->getSnOrder();
        if (dimension == 1)
        {
            number_of_levels_ /= 2;
        }

        vector<double> C(number_of_levels_);

        double Csum = 0;
        int level = -1;
        for (unsigned a=0; a<number_of_ordinates; a++)
        {
            double const mu = ordinates[a].mu();
            double const wt = ordinates[a].wt();
            if (wt!=0)
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

                alpha_[a] = 0.0;

                if (level>=0)
                    // Save the normalization sum for the previous level, if
                    // any. 
                {
                    C[level] = 1.0/Csum;
                }
                level++;
                Csum = 0.0;

                levels_[a] = level;
                is_dependent_[a] = false;

                if (level> 0) 
                    first_angles_.push_back(a-1);
            }
        }
        // Save the normalization sum for the final level.
        C[level] = 1.0/Csum;
        first_angles_.push_back(number_of_ordinates-1);

#if DBC & 2
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
        double omp;
        level = -1;

        for (unsigned a=0; a<number_of_ordinates; a++)
        {
            double const xi = ordinates[a].xi();
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
            else
                // New level.  Reinitialize the recurrence relation.
            {
                omp = rtt_units::PI;
                Check(1-xi*xi >= 0.0);
                sinth = std::sqrt(1-xi*xi);
                level++;
            }
            mup = sinth*std::cos(omp);
            if (wt!=0)
            {
                tau_[a] = (mu-mum)/(mup-mum);
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
                is_dependent_[a] = false;
                alpha_[a] = 0;
                first_angles_.push_back(number_of_ordinates-1);
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
    }

    Insist(first_angles_.size() == number_of_levels_, "unexpected starting direction reflection index");

    //for (unsigned i=0; i < number_of_levels_; ++i)
    //{
    //    unsigned const ia=first_angles_[i];
    //    
    //    std::cout << " reflection angle for starting direction on level " << i
    //              << " out of " << first_angles_.size()
    //              << " is angle " << ia
    //              << " with mu = " << ordinates[ia].mu()
    //              << " xi = " << ordinates[ia].xi()
    //              << " wt = " << ordinates[ia].wt()
    //              << std::endl;
    //}

    Ensure(check_class_invariants());
}

//---------------------------------------------------------------------------//
/*!
 * The psi coefficient is used to compute the self term in the angle
 * derivative term of the streaming operator.
 */
double Angle_Operator::Psi_Coefficient(unsigned const a) const
{
    Require(Is_Dependent(a));

    double const wt = getOrdinates()[a].wt();
    double const alpha_a = alpha_[a];
    double const tau_a = tau_[a];
    double const Result = alpha_a/(wt*tau_a);
    return Result;
}

//---------------------------------------------------------------------------//
/*!
 * The source coefficient is used to compute the previous midpoint angle term
 * in the angle derivative term of the streaming operator.
 */
double Angle_Operator::Source_Coefficient(unsigned const a) const
{
    Require(Is_Dependent(a));

    double const wt = getOrdinates()[a].wt();
    double const alpha_a = alpha_[a];
    double const alpha_am1 = alpha_[a-1];
    double const tau_a = tau_[a];
    double const Result = (alpha_a*(1-tau_a)/tau_a + alpha_am1)/wt;
    return Result;
}
//---------------------------------------------------------------------------//
/*!
 * The bookkeeping coefficient is used to compute the next midpoint angle
 * specific intensity.
 */
double Angle_Operator::Bookkeeping_Coefficient(unsigned const a) const
{
    Require(Is_Dependent(a));

    double const tau_a = tau_[a];
    double const Result = 1.0/tau_a;

    Require(Result>0.0);
    return Result;
}

//---------------------------------------------------------------------------//
bool Angle_Operator::check_class_invariants() const
{
    if (getGeometry() == rtt_mesh_element::CARTESIAN)
    {
        return number_of_levels_ == 0;
    }
    else
    {
        vector<Ordinate> const &ordinates = getOrdinates();
        unsigned const number_of_ordinates = ordinates.size();
        
        unsigned levels = 0;
        for (unsigned a=0; a<number_of_ordinates; ++a)
        {
            if (!is_dependent_[a])
            {
                ++levels;
            }
        }
        Require(number_of_levels_>=levels);
        
        return
            is_dependent_.size()==number_of_ordinates &&
            alpha_.size()==number_of_ordinates &&
            tau_.size()==number_of_ordinates;
    }
}

//---------------------------------------------------------------------------//
/*!
 * \param quadrature Sn quadrature set.
 *
 * \param dimension Dimension of the physical problem space (1, 2, or 3)
 *
 * \param geometry Geometry of the physical problem space.
 *
 * \todo The checking done here is far from complete at present.
 */

/* static */
bool Angle_Operator::is_compatible( SP<Quadrature const> const &quadrature,
                                    rtt_mesh_element::Geometry const geometry,
                                    unsigned const dimension,
                                    ostream &cerr)
{
    Require(quadrature!=SP<Quadrature>());
    Require(dimension>0 && dimension<4);
    Require(geometry<rtt_mesh_element::END_GEOMETRY);

    bool Result = true;

    if (!soft_equiv(quadrature->getNorm(), 1.0))
    {
        cerr << "Quadrature must be normalized to unity" << endl;
        Result = false;
    }

    Ensure(check_static_class_invariants());
    return Result;
}

//---------------------------------------------------------------------------//
/*!
 *
 * In 1-D or 2-D geometry, it is useful to be able to project the
 * three-dimensional ordinate direction onto the problem geometry to produce a
 * pseudo-1D or -2D vector that can (for instance) be dotted with the
 * pseudo-1D or -2D vector representating a finite element face normal. This
 * is not entirely straightforward in 2-D, because one generally wants to map
 * the third component of the ordinate direction onto the second component of
 * the pseudo-2D vector. This function provides a handy interface for such
 * mapping that hides these details from the client.
 */
vector<double>
Angle_Operator::Projected_Ordinate(unsigned const a) const
{
    Require(a < getOrdinates().size());

    vector<Ordinate> const &ordinates = getOrdinates();
    unsigned const dimension = getDimension();
    
    Ordinate const &ordinate = ordinates[a];
    vector<double> Result(dimension);
    Result[0] = ordinate.mu();
    if (dimension==2)
    {
        Result[1] = ordinate.xi();
    }
    else if (dimension==3)
    {
        Result[1] = ordinate.eta();
        Result[2] = ordinate.xi();
    }

    Ensure(Result.size()==dimension);
    return Result;
}

} // end namespace rtt_quadrature

//---------------------------------------------------------------------------//
//                 end of Angle_Operator.cc
//---------------------------------------------------------------------------//
