//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/Index_Set.hh
 * \author Mike Buksas
 * \date   Thu Feb  2 10:01:46 2006
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#ifndef dsxx_Index_Set_t_hh
#define dsxx_Index_Set_t_hh

#include <functional>
#include <numeric>

namespace rtt_dsxx {

//---------------------------------------------------------------------------//
/**
 * \brief Set the size of the Index_Set. Discards old size information
 * \arg sizes Pointer to unsigned integers for the index set sizes.
 */
template <unsigned D, int OFFSET>
void Index_Set<D, OFFSET>::set_size(unsigned const *const dimensions_) {
  std::copy(dimensions_, dimensions_ + D, dimensions);
  Require(sizes_okay());
  compute_size();
}

//---------------------------------------------------------------------------//
/**
 * \brief Set the size of the Index_Set to a uniform dimension. Discards old
 *        size information
 * \arg dimension The uniform dimension of the index set.
 */
template <unsigned D, int OFFSET>
void Index_Set<D, OFFSET>::set_size(const unsigned dimension) {
  for (unsigned *dim = dimensions; dim < dimensions + D; ++dim)
    *dim = dimension;
  compute_size();
}

//---------------------------------------------------------------------------//
/**
 * \brief Comparison routine
 * \arg The Index_Set object to compare to.
 */
template <unsigned D, int OFFSET>
inline bool Index_Set<D, OFFSET>::operator==(const Index_Set &rhs) const {
  if (array_size != rhs.array_size)
    return false;
  return std::equal(dimensions, dimensions + D, rhs.dimensions);
}

//---------------------------------------------------------------------------//
/**
 * \brief Make sure the indices are with the range for each dimension
 * \arg iterator An itertator to a range of indices.
 */
template <unsigned D, int OFFSET>
template <typename IT>
bool Index_Set<D, OFFSET>::indices_in_range(IT indices) const {

  int dimension = 0;
  for (IT index = indices; index != indices + D; ++index, ++dimension)
    if (!index_in_range(*index, dimension))
      return false;

  return true;
}

//---------------------------------------------------------------------------//
/**
 * \brief Return true iff the given index is within the range for the given
 *        dimension
 *
 * \arg index     The index value
 * \arg dimension The dimension of the index
 */
template <unsigned D, int OFFSET>
inline bool Index_Set<D, OFFSET>::index_in_range(int index,
                                                 unsigned dimension) const {
  Check(dimension_okay(dimension));

  return ((index >= OFFSET) &&
          (index < static_cast<int>(dimensions[dimension]) + OFFSET));
}

//---------------------------------------------------------------------------//
// IMPLEMENTAION
//---------------------------------------------------------------------------//
template <unsigned D, int OFFSET>
inline void Index_Set<D, OFFSET>::compute_size() {

  array_size = std::accumulate<unsigned *>(dimensions, dimensions + D, 1,
                                           std::multiplies<unsigned>());
  Ensure(array_size > 0);
}

} // end namespace rtt_dsxx

#endif // dsxx_Index_Set_t_hh

//---------------------------------------------------------------------------//
// end of ds++/Index_Set.hh
//---------------------------------------------------------------------------//
