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
      radius(radius_)
{

    Check(radius > 0);

}


double Sphere::distance_to(vector<double> position,
			   const vector<double>& direction)
{

    Check(position.size()  == 3);
    Check(direction.size() == 3);

    // Shift to sphere coordinates:
    position[2] -= center;

    double radius_2          =  radius * radius;
    double center_distance_2 =  global::dot(position, position);
    double normal_projection = -global::dot(position, direction);
    double normal_distance_2 =  normal_projection * normal_projection;
    
    double normal_length_2   = center_distance_2 - normal_distance_2;
    double half_chord_2      = radius_2 - normal_length_2;

    if ( center_distance_2 < radius_2 ) // Inside
    {
	return normal_projection + std::sqrt(half_chord_2);
    }
    else if ( half_chord_2 > 0 && normal_projection > 0 )
    {
	return normal_projection - std::sqrt(half_chord_2);
    }

    return global::huge;

}

double Sphere::distance_to(vector<double> position, const vector<double>& direction,
			   bool is_inside)
{

    Check(position.size()  == 3);
    Check(direction.size() == 3);

    // Shift to sphere coordinates:
    position[2] -= center;

    double radius_2          =  radius * radius;
    double center_distance_2 =  global::dot(position, position);
    double normal_projection = -global::dot(position, direction);
    double normal_distance_2 =  normal_projection * normal_projection;
    
    double normal_length_2   = center_distance_2 - normal_distance_2;
    double half_chord_2      = radius_2 - normal_length_2; 
    
    double distance = global::huge;

    if ( is_inside ) // Inside
    {
	Insist(half_chord_2 >= 0, "Inconsistent data in Sphere::distance_to");
	
	distance = normal_projection + std::sqrt(half_chord_2);

	Insist(distance >= 0,  "Inconsistent data in Sphere::distance_to");
    }
    else if ( half_chord_2 > 0 && normal_projection > 0 )
    {
	distance =  normal_projection - std::sqrt(half_chord_2);

	Insist(distance >= 0,  "Inconsistent data in Sphere::distance_to");
    }

    return distance;

}

bool Sphere::is_inside(std::vector<double> position)
{
    position[2] -= center;

    double center_distance_2 = global::dot(position, position);
    double raidus_2          = radius * radius;

    return (center_distance_2 <= raidus_2);
    
}

double Sphere::surface_area()
{
    return 4.0 * global::pi * radius*radius; 
}

double Sphere::volume()
{
    return 4.0/3.0 * global::pi * radius*radius*radius; 
}

} // end namespace rtt_mc

//---------------------------------------------------------------------------//
//                 end of Sphere.cc
//---------------------------------------------------------------------------//
