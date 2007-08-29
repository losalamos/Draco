//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/Angle_Operator.cc
 * \author Kent Budge
 * \date   Mon Mar 26 16:11:19 2007
 * \brief  
 * \note   Copyright (C) 2006 Los Alamos National Security, LLC.
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
 * \param quadrature Sn quadrature set.
 *
 * \param geometry Geometry of the physical problem space.
 *
 * \param dimension Dimension of the physical problem space (1, 2, or 3)
 */

Angle_Operator::Angle_Operator( SP<Quadrature const> const &quadrature,
                                rtt_mesh_element::Geometry const geometry,
                                unsigned const dimension)
    :
    OrdinateSet(quadrature, geometry, dimension)
{
    Require(quadrature!=SP<Quadrature>());
    Require(dimension>0 && dimension<4);
    Require(geometry<rtt_mesh_element::END_GEOMETRY);

    vector<Ordinate> const &ordinates = getOrdinates();
    unsigned const number_of_ordinates = ordinates.size();

    // Compute the ordinate derivative coefficients.
    
    // Default values are for the trivial case (Cartesian geometry).
    is_dependent.resize(number_of_ordinates, false);
    alpha.resize(number_of_ordinates, 0.0);
    tau.resize(number_of_ordinates, 1.0);

    // We rely on OrdinateSet to have already sorted the ordinates and
    // inserted the starting ordinates for each level. We assume that the
    // starting ordinates are distinguished by zero quadrature weight.

    number_of_levels = 0;
    levels.resize(number_of_ordinates);
    if (geometry==rtt_mesh_element::AXISYMMETRIC)
    {
        number_of_levels = quadrature->getSnOrder();
        if (dimension == 1)
        {
            number_of_levels /= 2;
        }

        vector<double> C(number_of_levels);

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
                alpha[a] = alpha[a-1] + mu*wt;
                Csum += wt;
                levels[a] = level;
                is_dependent[a] = true;
            }
            else
                // A starting ordinate. Reinitialize the recurrence relation.
            {
                Check(a==0 || std::fabs(alpha[a-1])<1.0e-15);
                // Be sure that the previous level (if any) had a final alpha
                // of zero, to within roundoff, as expected for the Morel
                // recursion formula.

                alpha[a] = 0.0;

                if (level>=0)
                    // Save the normalization sum for the previous level, if
                    // any. 
                {
//                    Check(isFinite(1.0/Csum));
                    C[level] = 1.0/Csum;
                }
                level++;
                Csum = 0.0;

                levels[a] = level;
                is_dependent[a] = false;
            }
        }
        // Save the normalization sum for the final level.
//        Check(isFinite(1.0/Csum));
        C[level] = 1.0/Csum;

#if DBC & 2
        if (dimension == 2)
            // Check that the level normalizations have the expected
            // properties. 
        {
            for (unsigned n=0; n<number_of_levels/2; n++)
            {
                Check(C[n]>0.0);
                Check(C[number_of_levels-1-n] > 0.0);
                Check(soft_equiv(C[n], C[number_of_levels-1-n]));
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
                //Check(omp>=-1.e-15 && omp<rtt_units::PI+1.e-15);
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
//                Check(isFinite((mu-mum)/(mup-mum)));
                tau[a] = (mu-mum)/(mup-mum);
                Check(tau[a] >= 0.0 && tau[a]<1.0);
            }
        }
    }
    else if (geometry == rtt_mesh_element::SPHERICAL)
    {
        number_of_levels = 1;
        double norm(0);
        for (unsigned a=0; a < number_of_ordinates; ++a)
        {
            levels[a] = 0;
            double const wt(ordinates[a].wt());
            norm += wt;
        }

        for (unsigned a=0; a < number_of_ordinates; ++a)
        {
            double const mu(ordinates[a].mu());
            double const wt(ordinates[a].wt());

            if (wt!=0)
            {
                is_dependent[a] = true;
                alpha[a] = alpha[a-1] + 2*wt*mu;
            }
            else
            {
                is_dependent[a] = false;
                alpha[a] = 0;
            }
        }

        double mup = -2; // sentinel
        for (unsigned a=0; a < number_of_ordinates; ++a)
        {
            double const mu(ordinates[a].mu());
            double const wt(ordinates[a].wt());

            double const mum = mup;

            if (wt !=0)
                mup = mum + wt;
            else
                mup = mu;

            if (wt !=0)
            {
//                Check(isFinite((mu-mum)/wt));
                tau[a] = (mu-mum)/wt;
                Check(tau[a]>0.0 && tau[a]<=1.0);
            }
        }
    }
    else
    {
        Check(geometry == rtt_mesh_element::CARTESIAN);
    }
}

//---------------------------------------------------------------------------//
double Angle_Operator::Psi_Coefficient(unsigned const a) const
{
    Require(a!=0);

    double const wt = getOrdinates()[a].wt();
    double const alpha_a = alpha[a];
    double const tau_a = tau[a];
    double const Result = alpha_a/(wt*tau_a);
    return Result;
}

//---------------------------------------------------------------------------//
double Angle_Operator::Source_Coefficient(unsigned const a) const
{
    Require(a!=0);

    double const wt = getOrdinates()[a].wt();
    double const alpha_a = alpha[a];
    double const alpha_am1 = alpha[a-1];
    double const tau_a = tau[a];
    double const Result = (alpha_a*(1-tau_a)/tau_a + alpha_am1)/wt;
    return Result;
}
//---------------------------------------------------------------------------//
double Angle_Operator::Bookkeeping_Coefficient(unsigned const a) const
{
    Require(a!=0);

    double const tau_a = tau[a];
    double const Result = 1.0/tau_a;

    Require(Result>0.0);
    return Result;
}

//---------------------------------------------------------------------------//
bool Angle_Operator::check_class_invariants() const
{
    if (getGeometry() == rtt_mesh_element::CARTESIAN)
    {
        return number_of_levels == 0;
    }
    else
    {
        vector<Ordinate> const &ordinates = getOrdinates();
        unsigned const number_of_ordinates = ordinates.size();
        
        double levels = 0;
        for (unsigned a=0; a<number_of_ordinates; ++a)
        {
            if (!is_dependent[a])
            {
                ++levels;
            }
        }
        Require(number_of_levels>=levels);
        
        return
            is_dependent.size()==number_of_ordinates &&
            alpha.size()==number_of_ordinates &&
            tau.size()==number_of_ordinates;
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

    if ( (geometry == rtt_mesh_element::SPHERICAL) &&
         (dimension==1) &&
         (!soft_equiv(quadrature->getNorm(),2.0)))
    {
        cerr << "Quadrature must be normalized to 2.0" << endl;
        Result = false;
    }

    return Result;
}

//---------------------------------------------------------------------------//
//! Return the projection of an ordinate direction onto the mesh geometry.

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
