//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Frequency.hh
 * \author Todd J. Urbatsch
 * \date   Thu Jan 10 14:05:09 2002
 * \brief  Gray_Frequency and Multigroup_Frequency class definitions.
 * \note   Copyright © 2003 The Regents of the University of California.
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
 * \class Multigroup_Frequency 
 *
 * \brief The Multigroup_Frequency class type.
 *
 */
// revision history:
// -----------------
// 0) 10 Jan 2002 : original
// 1) 16 Jan 2002 : replaced Frequency with Gray_Frequency and
//                  Multigroup_Frequency, types on which we will template.
// 
//===========================================================================//

class Multigroup_Frequency 
{

  public:

    // Useful typedef.
    typedef std::vector<double>        sf_double;
    typedef std::pair<double,double>   pair_doubles;

  private:

    // Frequency group boundaries.
    sf_double group_boundaries;

  public:

    // Constructor.
    explicit Multigroup_Frequency(const sf_double &grp_bnds)
	: group_boundaries(grp_bnds)
    {
	// Require at least one group
	Require (group_boundaries.size() > 1);

	// Loose check on positivity/monotonicity of mg group boundaries. 
	Require (group_boundaries.front() <  group_boundaries.back());
	Require (group_boundaries.front() >= 0.0);
    }

    // >>> Accessors.
    
    //! Query whether frequency is gray.
    bool is_gray() const { return false; }

    //! Query whether frequency is multigroup.
    bool is_multigroup() const { return true; }

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
	Check (group_boundaries.size() > 1);

	// Require a valid group number on input.
	Check (group > 0);
	Check (group <= get_num_groups());

	// Check for monotonicity of group bounds for this group.
	Check (group_boundaries[group-1] < group_boundaries[group]); 

	// Return this group's frequency boundaries as a pair.
	return pair_doubles(group_boundaries[group-1], 
			    group_boundaries[group]); 
    }

    //! Find a frequency group, in [1,G], given a frequency [keV]. 
    int find_group_given_a_freq(const double hnu) const
    {
	// The freqency must be multigroup.
	Require (group_boundaries.size() > 1);

	// Get the number of groups.
	int num_grps = group_boundaries.size() - 1;
	Require (num_grps > 0);

	// Return error if hnu is outside group boundaries.
	if (hnu < group_boundaries[0] || hnu > group_boundaries[num_grps])
	    return 0;

	// Binary search on groups. The search is successful when either 
	// 1) the frequency is within the group bounds or 
	// 2) the low and high search bounds are the same.
	int grp_lo  = 1;
	int grp_hi  = num_grps;
	int grp_try = static_cast<int>(0.5 * (grp_lo + grp_hi));

	while (grp_lo != grp_hi)
	{
	    // check validity of group index
	    Check (grp_try > 0 && grp_try <= num_grps);

	    // move up low search bound if frequency is above this group
	    // and split the new difference to get next estimate.
	    if (hnu > group_boundaries[grp_try])
	    {
		grp_lo  = grp_try + 1;
		grp_try = static_cast<int>(0.5 * (grp_lo + grp_hi));
	    }

	    // move down high search bound if frequency is below this group
	    // and split the new difference to get next estimate.
	    else if (hnu < group_boundaries[grp_try-1])
	    {
		grp_hi  = grp_try - 1;
		grp_try = static_cast<int>(0.5 * (grp_lo + grp_hi));
	    }

	    // otherwise, we have the correct group; return
	    else
		return grp_try;
	}

	// Check consistency of search params.
	Check (grp_lo == grp_hi && grp_try == grp_hi);
	Check (grp_try > 0      && grp_try <= num_grps);

	// check that the group actually contains the frequency.
	Check (hnu >= get_group_boundaries(grp_try).first);
	Check (hnu <= get_group_boundaries(grp_try).second);

	// return the group number (in [1,G])
	return grp_try;
    }
};
 
//===========================================================================//
/*!
 * \class Gray_Frequency 
 *
 * \brief The Gray_Frequency class type.
 *
 */
// revision history:
// -----------------
// 0) 10 Jan 2002 : original
// 1) 16 Jan 2002 : replaced Frequency with Gray_Frequency and
//                  Multigroup_Frequency, types on which we will template.
// 
//===========================================================================//

class Gray_Frequency 
{
  public:

    // Constructor.
    explicit Gray_Frequency() {/*...*/}

    // >>> Accessors.
    
    //! Query whether frequency is gray.
    bool is_gray() const { return true; }

    //! Query whether frequency is multigroup.
    bool is_multigroup() const { return false; }

    //! Get number of frequency groups.
    int get_num_groups() const { return 0; }
};

} // end namespace rtt_imc

#endif                          // __imc_Frequency_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Frequency.hh
//---------------------------------------------------------------------------//
