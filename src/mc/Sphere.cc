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


double Sphere::distance_to(vector<double> position,
			   const vector<double>& direction) const
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

    if ( center_distance_2 < radius_2 ) // Inside
    {
	distance = normal_projection + std::sqrt(half_chord_2);
    }
    else if ( normal_projection > 0 && normal_length_2 < radius_2 )
    {
	distance =  normal_projection - std::sqrt(half_chord_2);
    } 
    
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
	Ensure( normal_projection > 0 || center_distance_2 < radius_2 );
	
	double half_chord_2 = radius_2 - normal_length_2; 
	distance = normal_projection + std::sqrt(half_chord_2);
	
    }
    // Offically outside
    else if ( normal_projection > 0 && normal_length_2 < radius_2 ) 
    {
	// Crossing surface

	// This condition catches "leakage" of the ray to the inside when
	// it's offical status is outside.
	Ensure(center_distance_2 > radius_2);
	
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
	return ( center_distance_2 <= radius_2 );
    else
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
