//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/Index_Set.hh
 * \author Mike Buksas
 * \date   Thu Feb  2 10:01:46 2006
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#ifndef dsxx_Index_Set_hh
#define dsxx_Index_Set_hh

#include "Assert.hh"
#include <algorithm>

namespace rtt_dsxx {

//===========================================================================//
/*!
 * \class Index_Set
 * \brief Represents a D-dimensional set if indices.
 * \sa Index_Set.cc for detailed descriptions.
 */
/*!
 * \example ds++/test/tstIndex_Set.cc
 */
//===========================================================================//
template <unsigned D, int OFFSET> class Index_Set {
public:
  // CREATORS

  //! Default constructors.
  Index_Set() : array_size(0) { /* ... */
  }

  //! Construct with pointer to sizes
  Index_Set(unsigned const *const dimensions) : array_size(0) {
    set_size(dimensions);
  }

  //! Construct with all dimensions equal
  Index_Set(const unsigned dimension) : array_size(0) { set_size(dimension); }

  //! Destructor
  virtual ~Index_Set() { /* ... */
  }

  //! Comparison operator
  bool operator==(const Index_Set &rhs) const;

  //! Negative comparison operator
  bool operator!=(const Index_Set &rhs) const { return !(*this == rhs); }

  //! Re-assignment operator
  void set_size(unsigned const *const dimensions);

  //! Uniform size re-assignment operator
  void set_size(const unsigned size);

  bool index_in_range(int index) const {
    return (index >= OFFSET) && (index < static_cast<int>(array_size) + OFFSET);
  }
  bool index_in_range(int index, unsigned dimension) const;

  template <typename IT> bool indices_in_range(IT indices) const;

  int get_size() const { return array_size; }
  int min_of_index() const { return OFFSET; }
  int max_of_index() const { return OFFSET + array_size - 1; }
  int limit_of_index(const bool positive) const {
    return positive ? max_of_index() : min_of_index();
  }

  int get_size(const int d) const {
    Check(dimension_okay(d));
    return dimensions[d];
  }
  int min_of_index(const unsigned Remember(d)) const {
    Check(dimension_okay(d));
    return OFFSET;
  }
  int max_of_index(const unsigned d) const {
    Check(dimension_okay(d));
    return OFFSET + dimensions[d] - 1;
  }
  int limit_of_index(const unsigned d, const bool positive) const {
    return positive ? max_of_index(d) : min_of_index(d);
  }

  static bool direction_okay(const size_t d) { return (d > 0) && (d <= 2 * D); }
  static bool dimension_okay(const size_t d) { return d < D; }

private:
  void compute_size();

  unsigned array_size;    //!< Sizes of the whole index range
  unsigned dimensions[D]; //!< Sizes of each dimension

protected:
  // Make sure the index sizes are all positive when creating or resizing:
  bool sizes_okay() const {
    return (std::find(dimensions, dimensions + D, 0) == dimensions + D);
  }

  // Allow derived classes const access to the dimensions.
  unsigned const *get_dimensions() const { return dimensions; }
};

} // end namespace rtt_dsxx

#include "Index_Set.t.hh"

#endif // dsxx_Index_Set_hh

//---------------------------------------------------------------------------//
// end of ds++/Index_Set.hh
//---------------------------------------------------------------------------//
