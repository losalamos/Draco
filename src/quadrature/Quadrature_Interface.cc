//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/Quadrature_Interface.cc
 * \author Jae Chang
 * \date   Tue Jan 27 08:51:19 2004
 * \brief  Quadrature interface definitions
 * \note   Copyright (C) 2004-2015 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id: Arguments.hh 7090 2015-01-07 17:01:48Z kellyt $
//---------------------------------------------------------------------------//


#include <iostream>

#include "Quadrature_Interface.hh"
#include "Quadrature.hh"
#include "Gauss_Legendre.hh"
#include "Level_Symmetric.hh"
#include "Lobatto.hh"
#include "Square_Chebyshev_Legendre.hh"
#include "Tri_Chebyshev_Legendre.hh"
#include "Product_Chebyshev_Legendre.hh"

quadrature_data::quadrature_data()
    : dimension(0),
      type(0),
      order(0),
      azimuthal_order(0),
      geometry(0),
      mu(NULL),
      eta(NULL),
      xi(NULL),
      weights(NULL)
{

}

void init_quadrature(quadrature_data& quad)
{
    quad = quadrature_data();
}


void get_quadrature(quadrature_data& quad)
{

    using rtt_mesh_element::Geometry;
    using rtt_dsxx::SP;
    using namespace::rtt_quadrature;


    bool add_starting_directions = false;
    bool add_extra_directions = false;


    Geometry geometry;

    // Find the geometry
    switch ( quad.geometry )
    {
        case 0 :
            geometry = rtt_mesh_element::CARTESIAN;
            break;

        case 1 :
            geometry = rtt_mesh_element::AXISYMMETRIC;
            add_starting_directions = true;
            break;

        case 2:
            geometry = rtt_mesh_element::SPHERICAL;
            add_starting_directions = true;
            break;

        default :
            Insist(false,"Unrecongnized Geometry");
            geometry = rtt_mesh_element::CARTESIAN;

    }

    SP<Ordinate_Set> ordinate_set;
    unsigned size;
    vector<Ordinate> ordinates;

    if(quad.dimension==1)
    {  // 1D quadratures

        if (quad.type == 0)
        {
            Gauss_Legendre quadrature(quad.order);
            ordinate_set =
                quadrature.create_ordinate_set(1,
                                               geometry,
                                               1.0, // norm,
                                               add_starting_directions,
                                               add_extra_directions,
                                               Ordinate_Set::LEVEL_ORDERED);
        }
        else if (quad.type == 1)
        {
            Lobatto quadrature(quad.order);
            ordinate_set =
                quadrature.create_ordinate_set(1,
                                               geometry,
                                               1.0, // norm,
                                               add_starting_directions,
                                               add_extra_directions,
                                               Ordinate_Set::LEVEL_ORDERED);
        }

        ordinates = ordinate_set->ordinates();

        size = ordinates.size();

        // Copy wts and mu
        for (unsigned i = 0; i < size; i++)
        {
            Ordinate const &ordinate = ordinates[i];
            quad.weights[i] = ordinate.wt();
            quad.mu[i]      = ordinate.mu();
        }

    }
    else if( quad.dimension == 2)
    {  // 2D quadratures
        if (quad.type == 0)
        {
            Level_Symmetric quadrature(quad.order);
            ordinate_set =
                quadrature.create_ordinate_set(2,
                                               geometry,
                                               1.0, // norm,
                                               add_starting_directions,
                                               add_extra_directions,
                                               Ordinate_Set::LEVEL_ORDERED);
        }
        else if (quad.type == 1)
        {
            Tri_Chebyshev_Legendre quadrature(quad.order);
            ordinate_set =
                quadrature.create_ordinate_set(2,
                                               geometry,
                                               1.0, // norm,
                                               add_starting_directions,
                                               add_extra_directions,
                                               Ordinate_Set::LEVEL_ORDERED);
        }
        else if (quad.type == 2)
        {
            Square_Chebyshev_Legendre quadrature(quad.order);
            ordinate_set =
                quadrature.create_ordinate_set(2,
                                               geometry,
                                               1.0, // norm,
                                               add_starting_directions,
                                               add_extra_directions,
                                               Ordinate_Set::LEVEL_ORDERED);
        }
        else if (quad.type == 3)
        {
            Product_Chebyshev_Legendre quadrature(quad.order,
                                                  quad.azimuthal_order);
            ordinate_set =
                quadrature.create_ordinate_set(2,
                                               geometry,
                                               1.0, // norm,
                                               add_starting_directions,
                                               add_extra_directions,
                                               Ordinate_Set::LEVEL_ORDERED);
        }

        ordinates = ordinate_set->ordinates();

        size = ordinates.size();

        // Copy wts, mu and eta
        for (unsigned i = 0; i < size; i++)
        {
            Ordinate const &ordinate = ordinates[i];
            quad.weights[i] = ordinate.wt();
            quad.mu[i]      = ordinate.mu();
            quad.eta[i]     = ordinate.eta();
        }

    }

    return;
}
