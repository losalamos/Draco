//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Frequency.hh
 * \author Todd J. Urbatsch
 * \date   Thu Jan 10 14:05:09 2002
 * \brief  Frequency class definition.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __imc_Frequency_hh__
#define __imc_Frequency_hh__

#include "ds++/Assert.hh"
#include <vector>
#include <utility>

namespace rtt_imc
{
 
//===========================================================================//
/*!
 * \class Frequency
 *
 * \brief The Frequency class 
 *
 */
// revision history:
// -----------------
// 0) 10 Jan 2002 : original
// 
//===========================================================================//

class Frequency 
{

  public:

    // Useful typedef.
    typedef std::vector<double>        sf_double;
    typedef std::pair<double,double>   pair_doubles;

  private:

    // Boolean indicator for gray frequency treatment.
    bool gray;

    // Frequency group boundaries.
    sf_double group_boundaries;

  public:

    // Constructor.
    explicit Frequency(const sf_double & grp_bnds = sf_double())
	: group_boundaries(grp_bnds)
    {
	// Set the boolean gray indicator.
	gray = (group_boundaries.size() == 0);

	// Loose check on positivity/monotonicity of mg group boundaries. 
	Require (gray?true:
		 group_boundaries.front() < group_boundaries.back());
	Require (gray?true:group_boundaries.front() >= 0.0);
    }

    // >>> Accessors.
    
    //! Query whether frequency is gray.
    bool is_gray() const { return gray; }

    //! Query whether frequency is multigroup.
    bool is_multigroup() const { return !gray; }

    //! Get number of frequency groups.
    int get_num_groups() const { return group_boundaries.size()-1; }

    //! Get number of frequency group boundaries.
    int get_num_group_boundaries() const { return group_boundaries.size(); }

    //! Get all group boundaries.
    sf_double get_group_boundaries() const { return group_boundaries; }

    //! Get both low and high frequency group boundaries [keV].  Usage is
    //  get_group_boundaries.first and get_group_boundaries.second.
    pair_doubles get_group_boundaries(const int group) const
    {
	// Require multigroup.
	Check (!gray);

	// Require a valid group number on input.
	Check (group > 0);
	Check (group <= get_num_groups());

	// Check for monotonicity of group bounds for this group.
	Check (group_boundaries[group-1] < group_boundaries[group]); 

	// Return this group's frequency boundaries as a pair.
	return pair_doubles(group_boundaries[group-1], 
			    group_boundaries[group]); 
    }

};

} // end namespace rtt_imc

#endif                          // __imc_Frequency_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Frequency.hh
//---------------------------------------------------------------------------//
