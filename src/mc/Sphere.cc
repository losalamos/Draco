//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Sphere.cc
 * \author Mike Buksas
 * \date   Mon Jun 16 16:14:46 2003
 * \brief  Sphere implementation file
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Sphere.hh"
#include "Math.hh"
#include "Constants.hh" 

using std::vector;

namespace rtt_mc
{

Sphere::Sphere(double center_, double radius_)
    : center(center_),
      radius(radius_),
      radius_2(radius_*radius_)
{

    Check(radius > 0);

}

//---------------------------------------------------------------------------//
/*! 
 * \brief Computes the distance from a particle to the surface of the Sphere.
 *
 * \param position Position of the particle.
 * \param direction Direction of the particle.
 * \return Distance to the surface of the sphere, or huge if the surface is
 * not crossed.

 * The position and direction arguments describe a ray. This function returns
 * the positive distance along the ray where it first crosses the sphere, or
 * rtt_mc::global::huge if the ray does not cross the sphere.
 *
 * The algorithm described here uses squared quantities as much as
 * possible. This postpones the computation of a square root until, and only
 * if, it is needed.
 *
 * \image html mc/outside_crossing.png "Typical crossing configuration from outside"
 * 
 * The figure depicts a archetypal configuration when a particle is outside
 * of a spherical crossing surface. The particle (at point P) has direction
 * vector \f$\hat{\Omega}\f$ and crosses the sphere centered at C with radius
 * r. The figure is drawn in the plane containing the particle P, the sphere
 * center C and the direction vector \f$\hat{\Omega}\f$.
 *
 * Let \f$\vec{x}\f$ be the position vector of the particle and \f$\vec{c}\f$
 * be the center of the Sphere. Then define the relative particle position as
 * \f$\vec{x}'\f$. The distances in the figure (or, when more convenient,
 * their squares) can be computed as follows:
 *
 * \f[ d^2 = \vec{x}'\cdot\vec{x}' \f]
 * \f[ \beta = -\vec{\Omega}\cdot\vec{x}' \f]
 * \f[ n^2 = d^2 - \beta^2 \f]
 * \f[ \gamma^2 = r^2 - n^2 = r^2 - d^2 + \beta^2 \f]
 *
 * The two crossing distances are \f$ d_\pm=\beta\pm\gamma\f$. 
 *
 * We note that a necessary and sufficient condition for the line defined by
 * \f$(P,\vec{\Omega})\f$ to cross the sphere is \f$ n^2<r^2\f$. If the
 * particle is inside the sphere, this condition is automaticly satisfied. If
 * the particle is outside, we must ensure that the crossing happens in the
 * positive direction of motion. To do this we add the condition
 * \f$\beta>0\f$.
 *
 * For our application, we are interested in the smallest positive crossing
 * distance. For particles inside the sphere, this is always
 * \f$d_{+}=\beta+\gamma\f$. For particles outside, provided \f$ n^2<r^2\f$ and
 * \f$\beta>0\f$, this is \f$d_{-}\beta-\gamma\f$.
 *
 * The computation of \f$\gamma\f$ from \f$\gamma^2\f$ is the primary cost of
 * this algorithm. Therefore, we postpone it until the possibility that there
 * is no crossing is eliminated.
 *
 *

 */
double Sphere::distance_to(std::vector<double> position,
			   const std::vector<double>& direction) const
{

    Check(position.size()  == 3);
    Check(direction.size() == 3);

    // Shift to sphere coordinates:
    position[2] -= center;

    double center_distance_2 =  global::dot(position, position);
    double normal_projection = -global::dot(position, direction);
    double normal_distance_2 =  normal_projection * normal_projection;
    
    double normal_length_2   = center_distance_2 - normal_distance_2;
    double half_chord_2      = radius_2 - normal_length_2;

    double distance = global::huge;

    // Inside, crossing is guaranteed.
    if ( center_distance_2 < radius_2 )
    {
	Check ( half_chord_2 > 0);
	distance = normal_projection + std::sqrt(half_chord_2);
    }
    // Outside, with a crossing.
    else if ( normal_projection > 0 && normal_length_2 < radius_2 )
    {
	Check ( half_chord_2 > 0);
	distance =  normal_projection - std::sqrt(half_chord_2);
    } 
    // Else, outside, no crossing.
    
    Ensure( distance > 0 );

    return distance;

}

//---------------------------------------------------------------------------//
/*! 
 * \brief Returns the offical distance from a point to the surface of the Sphere.
 * 
 * \param position Position of the particle
 * \param direction Direction of the particle
 * \param is_inside "Offical" status of the particle with respect to the Sphere.
 * \return The "offical" distance to the surface of the sphere, or huge if
 * the sphere is not crossed.
 *
 * The algorithm in this functions differs from distance_to by the use of an
 * offical inside/outside status of the particle. This enhances the
 * robustness of detecting surface corssings by making sure that all
 * catching inconsistencies between the offical status and actual positions. 
 *
 * As an added feature, this function can be used to get both of the positive
 * crossing distances for a particle outside of the surface. Passing false as
 * the final argument gives the closer outside to inside crossing, while
 * passing true gives the more distance inside to outside crossing. This is
 * useful when detecting multiple crossings.
 *
 */
double Sphere::distance_to(std::vector<double> position, 
			   const std::vector<double>& direction,
			   bool is_inside) const
{

    Check(position.size()  == 3);
    Check(direction.size() == 3);

    // Shift to sphere coordinates:
    position[2] -= center;

    double center_distance_2 =  global::dot(position, position);
    double normal_projection = -global::dot(position, direction);
    double normal_distance_2 =  normal_projection * normal_projection;
    
    double normal_length_2   = center_distance_2 - normal_distance_2;
    
    double distance = global::huge;

    if ( is_inside ) // Offically inside
    {
	// This condition catches "leakage" of the ray to the outside when
	// it's offical status is inside.
	Require( normal_projection > 0 || center_distance_2 < radius_2 );
	
	double half_chord_2 = radius_2 - normal_length_2; 
	distance = normal_projection + std::sqrt(half_chord_2);
	
    }
    // Offically outside
    else if ( normal_projection > 0 && normal_length_2 < radius_2 ) 
    {
	// Crossing surface

	// This condition catches "leakage" of the ray to the inside when
	// it's offical status is outside.
	Require(center_distance_2 > radius_2);
	
	double half_chord_2 = radius_2 - normal_length_2; 
	distance =  normal_projection - std::sqrt(half_chord_2);
	
    }

    Ensure( distance > 0 );

    return distance;

}

//---------------------------------------------------------------------------//
bool Sphere::is_inside(std::vector<double> position) const 
{
    position[2] -= center;

    double center_distance_2 = global::dot(position, position);

    return (center_distance_2 < radius_2);
    
}

//---------------------------------------------------------------------------//
bool Sphere::is_inside(std::vector<double> position,
		       const std::vector<double> direction) const
{

    position[2] -= center;
    double normal_projection = -global::dot(position, direction);
    double center_distance_2 =  global::dot(position, position);

    if (normal_projection > 0)
	// Direction is inward, include boundary
	return ( center_distance_2 <= radius_2 );
    else
	// Direction is outward, exclude boundary
	return ( center_distance_2 < radius_2 );

    return 0;

}

//---------------------------------------------------------------------------//
double Sphere::surface_area() const
{
    return 4.0 * global::pi * radius_2; 
}

//---------------------------------------------------------------------------//
double Sphere::volume() const
{
    return 4.0/3.0 * global::pi * radius_2*radius; 
}

} // end namespace rtt_mc

//---------------------------------------------------------------------------//
//                 end of Sphere.cc
//---------------------------------------------------------------------------//
