//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Sphere.cc
 * \author Mike Buksas
 * \date   Mon Jun 16 16:14:46 2003
 * \brief  Implements a spherical surface for surface tallies
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
 * \brief Computes the distance from a particle to the surface.
 *
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
 * The position argument is modified by subtracting the vector representing
 * the center of the sphere. Hence, it becomes the relative position of the
 * particle with respect to the sphere.
 *
 * Variable center_distance_2 is the square of the distance from particle to
 * sphere center. If center_distance_2 < radius_2 the particle is inside the
 * sphere.
 *
 * Variable normal_projection is the negative dot product of the direction
 * vector and the relative position vector. If normal_projection > 0 the
 * particle is heading toward the sphere. Hence for particles outside the
 * sphere, normal_projection>0 is a necessary (but not sufficient) condition
 * for a crossing. It is also the distance along the ray to the point of
 * closest approach to the center.
 *
 * Variable normal_distance_2 is the square of normal_projection.
 *
 * Variable normal_length_2 is the square of the length of the segment
 * connecting the sphere center to the point of closest approach. It is
 * computed with the Pythagorean theorem from the distance to the center and
 * the distance to the point of closest approach.
 *
 * If the particle is outside, necessary and sufficient conditions for a
 * positive crossing distance are normal_projection>0 and normal_length_2 <
 * raidus_2
 *
 * \param position Position of the particle

 * \param direction Direction of the particle

 * \return distance to the surface of the sphere, or huge if the surface is
 * not crossed

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

double Sphere::distance_to(vector<double> position, 
			   const vector<double>& direction,
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

bool Sphere::is_inside(std::vector<double> position) const 
{
    position[2] -= center;

    double center_distance_2 = global::dot(position, position);

    return (center_distance_2 < radius_2);
    
}

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

double Sphere::surface_area() const
{
    return 4.0 * global::pi * radius_2; 
}

double Sphere::volume() const
{
    return 4.0/3.0 * global::pi * radius_2*radius; 
}

} // end namespace rtt_mc

//---------------------------------------------------------------------------//
//                 end of Sphere.cc
//---------------------------------------------------------------------------//
