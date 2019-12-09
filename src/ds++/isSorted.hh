//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/isSorted.hh
 * \author Randy M. Roberts
 * \date   Wed Feb 16 09:27:40 2000
 * \brief  Functions that checks to see if a container is sorted.
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC
 */
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//

#ifndef __ds_isSorted_hh__
#define __ds_isSorted_hh__

namespace rtt_dsxx {

//===========================================================================//
/*!
 * \fn    isSorted
 * \brief Checks to see if elements of [first, last) are sorted, via "<".
 */
//===========================================================================//

template <class ForwardIterator>
bool isSorted(ForwardIterator first, ForwardIterator last) {
  if (first == last)
    return true;

  bool isSorted = true;
  ForwardIterator prev = first;
  ++first;
  while (isSorted && first != last) {
    isSorted = !(*first < *prev);
    prev = first;
    ++first;
  }
  return isSorted;
}

//===========================================================================//
/*!
 * \fn isSorted
 * \brief Checks to see if elements of [first, last) are sorted, via "comp".
 */
// revision history:
// -----------------
// 0) original
//
//===========================================================================//

template <class ForwardIterator, class StrictWeakOrdering>
bool isSorted(ForwardIterator first, ForwardIterator last,
              StrictWeakOrdering comp) {
  if (first == last)
    return true;

  bool isSorted = true;
  ForwardIterator prev = first;
  ++first;
  while (isSorted && first != last) {
    isSorted = !comp(*first, *prev);
    prev = first;
    ++first;
  }
  return isSorted;
}

} // end namespace rtt_dsxx

#endif // __ds_isSorted_hh__

//---------------------------------------------------------------------------//
//                              end of ds++/isSorted.hh
//---------------------------------------------------------------------------//
