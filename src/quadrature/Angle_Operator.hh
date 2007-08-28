//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/Angle_Operator.hh
 * \author Kent Budge
 * \date   Mon Mar 26 16:11:19 2007
 * \brief  Definition of class Angle_Operator
 * \note   Copyright (C) 2006 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef quadrature_Angle_Operator_hh
#define quadrature_Angle_Operator_hh

#include <vector>
#include "Ordinate.hh"

namespace rtt_quadrature
{

//===========================================================================//
/*!
 * \class Angle_Operator
 * \brief Represents the angle derivative term in the sweep operator
 *
 * In curvilinear geometry, the streaming operator includes a nontrivial angle
 * operator ithat introduces dependencies between ordinates. We assume that an
 * angle operator can be cast in block bidiagonal form, so that there is not
 * more than one direct dependency per ordinate. The Angle_Operator may then
 * order the ordinates by dependency, so the first ordinate can have no
 * dependencies, the second may be directly dependent only on the first, and
 * so on. Thus a client need only check whether an ordinate is dependent on
 * the preceeding angle or not.
 *
 * Each of the blocks in the block bidiagonal form of the angle operator is
 * referred to as a "level." This is terminology held over from the particular
 * case of 2-D axisymmetric geometry, where the ordinate sets generally are
 * organized on "levels" having the same z direction cosine (xi) which are
 * coupled by the omega derivative term. It is useful to know the number of
 * such blocks to optimize storage of intermediate results.
 */
//===========================================================================//

class Angle_Operator : public rtt_quadrature::OrdinateSet
{
  public:

    // NESTED CLASSES AND TYPEDEFS

    // CREATORS

    //! Specify the ordinate quadrature.
    Angle_Operator(rtt_dsxx::SP<Quadrature const> const &quadrature,
                   rtt_mesh_element::Geometry geometry,
                   unsigned dimension);

    // MANIPULATORS

    // ACCESSORS

    unsigned Number_Of_Levels() const { return number_of_levels; }

    std::vector<unsigned> const &Levels() const
    {
        return levels;
    }

    //! Is an ordinate dependent on the preceeding ordinate?
    bool Is_Dependent(unsigned const ordinate) const
    {
        Require(ordinate<getOrdinates().size());

        return is_dependent[ordinate];
    }

    std::vector<double> const &Alpha() const { return alpha; }

    std::vector<double> const &Tau() const { return tau; }

    double Psi_Coefficient(unsigned ordinate_index) const;

    double Source_Coefficient(unsigned ordinate_index) const;

    double Bookkeeping_Coefficient(unsigned ordinate_index) const;

    std::vector<double> Projected_Ordinate(unsigned ordinate_index) const;

    bool check_class_invariants() const;

    // STATICS

    static bool is_compatible(rtt_dsxx::SP<Quadrature const> const &quadrature,
                              rtt_mesh_element::Geometry geometry,
                              unsigned dimension,
                              std::ostream &cerr);

  private:

    // NESTED CLASSES AND TYPEDEFS

    // IMPLEMENTATION 

    // DATA

    unsigned number_of_levels;
    std::vector<unsigned> levels;

    //! Is an ordinate dependent on the preceeding ordinate?
    std::vector<bool> is_dependent;
    
    /*! Coefficients for angle derivative terms.  These are defined in
     * Morel's research note of 12 May 2003 for axisymmetric geometry. 
     */
    std::vector<double> alpha;
    std::vector<double> tau;
};

} // end namespace rtt_quadrature

#endif // quadrature_Angle_Operator_hh

//---------------------------------------------------------------------------//
//              end of quadrature/Angle_Operator.hh
//---------------------------------------------------------------------------//
