//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/OrdinateSet.hh
 * \author Kent Budge
 * \date   Tue Dec 21 14:20:03 2004
 * \brief  Declaration file for the class rtt_quadrature::Ordinate.
 * \note   Copyright (C)  2006-2012 Los Alamos National Security, LLC.
 *         All rights reserved. 
 */
//---------------------------------------------------------------------------//
// $Id: Ordinate.hh 6607 2012-06-14 22:31:45Z kellyt $
//---------------------------------------------------------------------------//

#ifndef quadrature_OrdinateSet_hh
#define quadrature_OrdinateSet_hh

#include <vector>

#include "mesh_element/Geometry.hh"
#include "ds++/SP.hh"
#include "ds++/Soft_Equivalence.hh"

#include "Quadrature.hh"
#include "QuadServices.hh"
#include "Ordinate.hh"

namespace rtt_quadrature
{

//===========================================================================//
/*!
 * \class OrdinateSet
 * \brief A collection of Ordinates that make up a complete quadrature set.
 */
//===========================================================================//

class OrdinateSet
{
  public:

    // CREATORS

    /*! Two ways to construct the object:
     * in the first one, the interpolation method is specified by the quadrature specification,
     * while in the second one, the interpolation method is specified as an argument.
     */

    OrdinateSet( rtt_dsxx::SP< Quadrature const > const quadrature,
                 rtt_mesh_element::Geometry geometry,
                 unsigned const dimension,
                 unsigned const expansion_order,
                 bool const extra_starting_directions=false,
                 comparator_t const comparator = Ordinate::SnCompare );
    
    OrdinateSet( rtt_dsxx::SP< Quadrature const > const quadrature,
                 rtt_mesh_element::Geometry geometry,
                 unsigned const dimension,
                 unsigned const expansion_order,
                 Quadrature::QIM const qim,
                 bool const extra_starting_directions=false,
                 comparator_t const comparator = Ordinate::SnCompare );
    
    //! destructor
    virtual ~OrdinateSet(){}

    // ACCESSORS

    //! Return the ordinates.
    std::vector<Ordinate> const &getOrdinates() const { return ordinates_; }

    //! Return the ordinates.
    unsigned getNumOrdinates() const { return ordinates_.size(); }

    //! Return the quadrature on which this OrdinateSet is based.
    rtt_dsxx::SP< const Quadrature > getQuadrature() const
    {
        Ensure(quadrature_ != rtt_dsxx::SP<Quadrature>());
        
        return quadrature_;
    }

    //! Return the geometry.
    rtt_mesh_element::Geometry getGeometry() const { return geometry_; }

    //! Return the dimension.
    unsigned getDimension() const { return dimension_; }

    //! Return the dimension.
    unsigned getExpansionOrder() const { return expansion_order_; }

    //! Return the norm.
    double getNorm() const { return norm_; }

    //! Return the ordering operator.
    comparator_t getComparator() const { return comparator_; }

    //! Return the ordering operator.
    bool extra_starting_directions() const
    {
        return extra_starting_directions_;
    }

    //! Return a QuadServices object.
    rtt_dsxx::SP<QuadServices const> get_qs() const { return qs_;}

    bool check_class_invariants() const;

  private:

    // Helper functions. 

    //! Initialize the ordinates, based on geometry input.
    std::vector<Ordinate> set_ordinates(rtt_mesh_element::Geometry const geometry);

    std::vector<Ordinate> create_set_from_1d_quadrature(rtt_mesh_element::Geometry const geometry);
    std::vector<Ordinate> create_set_from_2d_quadrature_for_2d_mesh(rtt_mesh_element::Geometry const geometry);
    std::vector<Ordinate> create_set_from_2d_quadrature_for_1d_mesh(rtt_mesh_element::Geometry const geometry);
    std::vector<Ordinate> create_set_from_2d_quadrature_for_3d_mesh();
    std::vector<Ordinate> create_set_from_3d_quadrature_for_3d_mesh();

    // DATA 

    // initialized
    rtt_dsxx::SP< Quadrature const > const quadrature_;
    rtt_mesh_element::Geometry const geometry_;
    unsigned const dimension_;
    unsigned const expansion_order_;
    bool const extra_starting_directions_;
    comparator_t comparator_;

    // computed during construction
    std::vector<Ordinate> const ordinates_;
    double norm_;
    rtt_dsxx::SP<QuadServices const> qs_;

};

} // end namespace rtt_quadrature

#endif // quadrature_OrdinateSet_hh

//---------------------------------------------------------------------------//
//              end of quadrature/OrdinateSet.hh
//---------------------------------------------------------------------------//
