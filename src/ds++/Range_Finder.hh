//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/Range_finder.hh
 * \author Mike Buksas
 * \date   Thu Feb  6 12:10:56 2003
 * \brief  Header file for Range_finder
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __dsxx_Range_finder_hh__
#define __dsxx_Range_finder_hh__

#include <algorithm>
#include <iterator>

namespace rtt_dsxx
{
 
//===========================================================================//
/*!
 * \class Range_finder
 *
 * These functions locate a value in intervals described by an increasing
 * array. E.g. Let v[i] be an increasing array and r a value. These functions
 * look for i such that: v[i] < r < v[i+1]. The different versions of the
 * function have different behavior in the event that r = v[i] for some i.
 *
 * The 'left' versions assume the intervals are closed on the left, so that
 * v[i] <= r < v[i+1]. Likewise, the 'right' versions assume closure on the
 * right: v[i] < r <= v[i+1].
 *
 * The 'catch end' versions will allow equality for values at the ends of the
 * range when it would not normally be considered a part of the interval. For
 * example if r=v[0], Range_finder_left will consider r to be outside the
 * range, but Range_finder_left_catch_end will return i=0. Likewise for r=v[n-1]
 * (v's last value) Range_finder_right_catch_end will return n-1 instead of
 * considering the value to be out of range.
 *
 * The functions Range_funder and Range_finder_catch_end take an enumeration
 * parameter RANGE_DIRECTION to determine the direction: LEFT or RIGHT.
 * 
 * Values which are out of range result in an exception. 
 *
 * This function uses the equal_range algorithm, which returns an iterator
 * pair. Both iterators will point to the first value _after_ r in v, unless
 * v[i] == r for some i, then the first iterator will point to v[i] and the
 * second to v[i+1]. This enables simple verification of out of bound index
 * detection and a simple compuation for the index which works whrether or not
 * equality is obtained.
 *
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

enum RANGE_DIRECTION { LEFT = 0, RIGHT = 1};

namespace {

// Toss an exception of the iterator pair is showing out-of-range.
template <typename IT>
bool validate(const std::pair<IT,IT> &it, IT begin, IT end)
{
    return !(it.first  == end ||        // v > v[n-1]
             it.second == begin );      // v < v[0]

}

}

template <typename IT>
int Range_finder_left(
    IT begin,
    IT end, 
    typename std::iterator_traits<IT>::value_type value)
{
    const std::pair<IT,IT> it = std::equal_range(begin, end, value);
    Require(validate(it,begin,end));
    
    const int index = static_cast<int>(it.second - begin) - 1;
    Ensure(index >= 0); Ensure(begin+index < end);

    return index;
}

template <typename IT>
int Range_finder_right(
    IT begin,
    IT end,
    typename std::iterator_traits<IT>::value_type value)
{
    const std::pair<IT,IT> it = std::equal_range(begin, end, value);
    Require(validate(it,begin,end));
    
    const int index = static_cast<int>(it.first - begin) - 1;
    Ensure(index >= 0); Ensure(begin+index < end);

    return index;
}


template <typename IT>
int Range_finder_left_catch_end(
    IT begin,
    IT end,
    typename std::iterator_traits<IT>::value_type value)
{
    const std::pair<IT,IT> it = std::equal_range(begin, end, value);
    Require(validate(it,begin,end));
    
    const int index      = static_cast<int>(it.second - begin) - 1;
    const int range_size = static_cast<int>(end-begin) - 1;

    // If (index==range_size), subtract one. This indicates that the value is
    // equal to *end.
    const int fix_index = std::min(index, range_size - 1);
    Ensure(fix_index >= 0); Ensure(begin+fix_index < end);

    return fix_index;
    

}

template <typename IT>
int Range_finder_right_catch_end(
    IT begin,
    IT end,
    typename std::iterator_traits<IT>::value_type value)
{
    const std::pair<IT,IT> it = std::equal_range(begin, end, value);
    Require(validate(it,begin,end));
    
    // Extract the coordintate index (0..n-1)
    const int index = static_cast<int>(it.first - begin) - 1;

    // If we got -1 here, then v=v[0] and we want to catch this end value:
    const int fix_index = std::max(index, 0);
    Ensure(fix_index >= 0);  Ensure(begin+fix_index < end);

    return fix_index;
}


//---------------------------------------------------------------------------//
// Generic versions.
//---------------------------------------------------------------------------//

template <typename IT>
int Range_finder(
    IT begin,
    IT end,
    typename std::iterator_traits<IT>::value_type value,
    RANGE_DIRECTION direction_indicator)
{
    Check (direction_indicator == 0 || direction_indicator <= 1);

    if (direction_indicator == LEFT)
	return Range_finder_left(begin, end, value);
    else if (direction_indicator == RIGHT)
	return Range_finder_right(begin, end, value);
    else 
	Insist(0, "Invalid direction indicator in Range_finder");

    return -1;

}


template <typename IT>
int Range_finder_catch_end(
    IT begin,
    IT end,
    typename std::iterator_traits<IT>::value_type value,
    RANGE_DIRECTION direction_indicator)
{
    Check (direction_indicator == 0 || direction_indicator <= 1);

    if (direction_indicator == LEFT)
	return Range_finder_left_catch_end(begin, end, value);
    else if (direction_indicator == RIGHT)
	return Range_finder_right_catch_end(begin, end, value);
    else
	Insist(0, "Invalid direction indicator in Range_finder");

    return -1;
}


} // end namespace rtt_dsxx

#endif                          // __dsxx_Range_finder_hh__

//---------------------------------------------------------------------------//
//                              end of dsxx/Range_finder.hh
//---------------------------------------------------------------------------//
