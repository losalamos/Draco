//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Azimuthal_Mesh.cc
 * \author Mike Buksas
 * \date   Mon Jun 23 12:56:30 2003
 * \brief  Implements azimuthal bins centered around the z-axis.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Azimuthal_Mesh.hh" 
#include "ds++/Assert.hh"

using std::vector;

namespace rtt_imc
{

Azimuthal_Mesh::Azimuthal_Mesh(const vector<double>& cosines)
    : bins(cosines.size() + 1),
      bin_cosines(1, -1.0)
{

    Check(cosines.size() > 0);

    std::copy(cosines.begin(), cosines.end(), back_inserter(bin_cosines));
    bin_cosines.push_back(1.0);

    Require( bin_cosines.size() == bins+1 );
    Ensure( check_cosines(bin_cosines) );

}

double Azimuthal_Mesh::get_lower_cosine(int bin) const
{

    Check (bin >= 1);  Check(bin <= bins);

    return bin_cosines[bin-1];

}

double Azimuthal_Mesh::get_upper_cosine(int bin) const
{

    Check (bin >= 1);  Check(bin <= bins);

    return bin_cosines[bin];

}

bool Azimuthal_Mesh::is_in_bin(const vector<double>& direction, int bin) const
{

    Check (bin >= 1); Check(bin <= bins);

    return 
	( direction[2] >= bin_cosines[bin-1] ) &&
	( direction[2] <= bin_cosines[bin]   );

}

int Azimuthal_Mesh::find_bin(const vector<double>& direction) const
{

    for (int bin = 1; bin <= bins; ++bin)
	if ( is_in_bin(direction, bin) ) return bin;

    return 0;

}

bool Azimuthal_Mesh::check_cosines(const vector<double>& cosines) const
{

    bool is_okay = true;

    double lower_value = -1.0;
    double upper_value = cosines[0];

    is_okay = (upper_value >= lower_value) && (upper_value < lower_value + 1.0e-7);

    for (vector<double>::const_iterator i = cosines.begin() + 1;
	 i != cosines.end() && is_okay;
	 ++i)
    {

	upper_value = *i;
	is_okay = is_okay && (upper_value > lower_value);

	lower_value = upper_value;

    }

    is_okay = is_okay && (upper_value <= 1.0) && (upper_value > 1.0-1.0e-7);

    return is_okay;
}

    

} // end namespace rtt_imc

//---------------------------------------------------------------------------//
//                 end of Azimuthal_Mesh.cc
//---------------------------------------------------------------------------//
