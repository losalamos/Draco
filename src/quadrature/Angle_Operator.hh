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
 *
 * We illustrate with an example. In axisymmetric geometry, <a
 * href="../../pdf/transport_implementation.pdf">Morel's 
 * discretization</a> of the angle derivative term in the streaming operator is
 *
 * \f$\frac{\partial (\eta \psi)}{\partial \omega}\approx
 * \frac{1}{w_m}(\alpha_{m+1/2}\psi_{m+1/2}-\alpha_{m-1/2}\psi_{m-1/2})\f$
 *
 * where
 *
 * \f$\psi_{m+1/2} = \frac{1}{\tau_m}(\psi_m-(1-\tau_m)\psi_{m-1/2})\f$
 *
 * Thus we define
 *
 * <code>Psi_Coefficient(m)</code> = \f$P_m =
 * \frac{\alpha_{m+1/2}}{w_m\tau_m}\f$ 
 *
 * <code>Source_Coefficient(m)</code> =
 \f$S_m = \frac{\alpha_{m+1/2}\frac{1-\tau_m}{\tau_m}+\alpha_{m-1/2}}{w_m}\f$
 *
 * <code>Bookkeeping_Coefficient(m)</code> = \f$B_m = \frac{1}{\tau_m}\f$
 *
 * The angle derivative can then be coded as
 *
 * \f$\frac{\partial (\eta \psi)}{\partial \omega} =
 * P_m\psi_m-S_m\psi_{m-1/2}\f$ 
 *
 * and the next midpoint intensity as
 *
 * \f$\psi_{m+1/2} = B_m\psi_m-(1-B_m)\psi_{m-1/2})\f$
 *
 * Similar expressions can be written for spherical geometry. The
 * Angle_Operator interface hides these details, presenting only the
 * \f$P_m\f$, \f$S_m\f$, and \f$B_m\f$ coefficients required for actual
 * computation.
 *
 * Note that this discretization of the angle derivative terms must still be
 * substituted into the transport equation, which is then further discretized
 * in space. Thus the angle derivative term will generally be multiplied by an
 * additional factor arising the spatial discretization.
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

    unsigned Number_Of_Levels() const { return number_of_levels_; }

    std::vector<unsigned> const &Levels() const
    {
        return levels_;
    }

    //! Return the angle index for the most positively outward-directed angle
    //! on every level; used for reflection of starting direction boundary fluxes
    std::vector<unsigned> const &First_Angles() const
    {
        return first_angles_;
    }

    //! Is an ordinate dependent on the preceeding ordinate?
    bool Is_Dependent(unsigned const ordinate) const
    {
        Require(ordinate<getOrdinates().size());

        return is_dependent_[ordinate];
    }

    //! Return \f$\alpha_{m+1/2}\f$ for ordinate \f$m\f$
    std::vector<double> const &Alpha() const { return alpha_; }
    
    //! Return \f$P\tau_m\f$ for ordinate \f$m\f$
    std::vector<double> const &Tau() const { return tau_; }

    //! Return \f$P_m\f$ for ordinate \f$m\f$
    double Psi_Coefficient(unsigned ordinate_index) const;

    //! Return \f$S_m\f$ for ordinate \f$m\f$
    double Source_Coefficient(unsigned ordinate_index) const;

    //! Return \f$B_m\f$ for ordinate \f$m\f$
    double Bookkeeping_Coefficient(unsigned ordinate_index) const;

    //! Return the projection of an ordinate direction onto the mesh geometry.
    std::vector<double> Projected_Ordinate(unsigned ordinate_index) const;

    bool check_class_invariants() const;

    // STATICS

    //! Determine if an Angle_Operator can be generated from a specified
    //! quadrature for a specified geometry and dimension.
    static bool is_compatible(rtt_dsxx::SP<Quadrature const> const &quadrature,
                              rtt_mesh_element::Geometry geometry,
                              unsigned dimension,
                              std::ostream &cerr);

    static bool check_static_class_invariants()
    {
        return true;
    }

  private:

    // NESTED CLASSES AND TYPEDEFS

    // IMPLEMENTATION 

    // DATA

    unsigned number_of_levels_;
    std::vector<unsigned> levels_;
    std::vector<unsigned> first_angles_;

    //! Is an ordinate dependent on the preceeding ordinate?
    std::vector<bool> is_dependent_;
    
    /*! Coefficients for angle derivative terms.  These are defined in
     * Morel's research note of 12 May 2003 for axisymmetric geometry. 
     */
    std::vector<double> alpha_;
    std::vector<double> tau_;
};

} // end namespace rtt_quadrature

#endif // quadrature_Angle_Operator_hh

//---------------------------------------------------------------------------//
//              end of quadrature/Angle_Operator.hh
//---------------------------------------------------------------------------//
