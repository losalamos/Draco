//----------------------------------*-C++-*----------------------------------------------//
/*!
 * \file   quadrature/Galerkin_Ordinate_Space.hh
 * \author Kent Budge
 * \date   Mon Mar 26 16:11:19 2007
 * \brief  Definition of class Galerkin_Ordinate_Space
 * \note   Copyright (C) 2006 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------------------//
// $Id: Galerkin_Ordinate_Space.hh 6718 2012-08-30 20:03:01Z warsa $
//---------------------------------------------------------------------------------------//

#ifndef quadrature_Galerkin_Ordinate_Space_hh
#define quadrature_Galerkin_Ordinate_Space_hh

#include "Ordinate_Space.hh"

namespace rtt_quadrature
{
using std::ostream;

//=======================================================================================//
/*!
 * \class Galerkin_Ordinate_Space
 * \brief Represents ordinate operators for a Galerkin moment space. 
 *
 * The moment space contains all moments (that are not identically zero due to
 * symmetry) up to the specified scattering order, but the moment to discrete
 * operator M and discrete to moment operator D are computed as if enough additional
 * higher moments are included in the moment space to make D and M square. The
 * higher moment terms are then discarded, but the non-square D and M retain
 * the property that DM is the identity. This stabilizes the moment to
 * discrete and discrete to moment operations at high scattering orders.
 */
//=======================================================================================//

class Galerkin_Ordinate_Space : public Ordinate_Space
{
  public:

    // NESTED CLASSES AND TYPEDEFS

    // CREATORS

    //! Specify the ordinate quadrature with defaults.
    Galerkin_Ordinate_Space(unsigned dimension,
                            Geometry geometry,
                            vector<Ordinate> const &,
                            Quadrature_Class,
                            unsigned sn_order,
                            unsigned expansion_order,
                            bool extra_starting_directions=false,
                            Ordering ordering=LEVEL_ORDERED);

    // MANIPULATORS

    // ACCESSORS

    bool check_class_invariants() const;

    // SERVICES
    
    virtual QIM quadrature_interpolation_model() const;

    //! Return the discrete to moment transform matrix
    virtual vector<double> D() const;

    //! Return the moment to discrete transform matrix
    virtual vector<double> M() const;

    // STATICS

  private:

    // NESTED CLASSES AND TYPEDEFS

    // IMPLEMENTATION

    virtual void compute_M();
    virtual void compute_D();
    
    virtual vector<Moment> compute_n2lk_1D_(Quadrature_Class,
                                            unsigned sn_order);
    
    virtual vector<Moment> compute_n2lk_1Da_(Quadrature_Class,
                                             unsigned sn_order);
    
    virtual vector<Moment> compute_n2lk_2D_(Quadrature_Class,
                                            unsigned sn_order);
    
    virtual vector<Moment> compute_n2lk_3D_(Quadrature_Class,
                                            unsigned sn_order);

  private:

    // NESTED CLASSES AND TYPEDEFS

    // IMPLEMENTATION
    
    
    vector<double>  compute_M_GQ(vector<Ordinate> const &ordinates,
                                 vector< Moment > const &n2lk,
                                 unsigned const dim,
                                 double const sumwt);
    
    vector<double> compute_D_GQ(vector<Ordinate> const &ordinates,
                                vector< Moment > const &n2lk,
                                vector<double> const &mM,
                                unsigned const,
                                double const);

    // DATA

    //! Discrete to moment matrix
    vector<double> D_;
    //! Moment to discrete matrix
    vector<double> M_;
};

} // end namespace rtt_quadrature

#endif // quadrature_Galerkin_Ordinate_Space_hh

//---------------------------------------------------------------------------------------//
//              end of quadrature/Galerkin_Ordinate_Space.hh
//---------------------------------------------------------------------------------------//
