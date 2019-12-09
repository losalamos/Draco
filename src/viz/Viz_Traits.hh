//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   viz/Viz_Traits.hh
 * \author Thomas M. Evans
 * \date   Fri Jan 21 17:10:54 2000
 * \brief  Viz_Traits header file.
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//

#ifndef __viz_Viz_Traits_hh__
#define __viz_Viz_Traits_hh__

#include "ds++/Assert.hh"
#include <vector>

namespace rtt_viz {

//===========================================================================//
/*!
 * \class Viz_Traits
 *
 * \brief Traits that are used by the rtt_viz package.
 *
 * These traits provide a common way to access 2D-styles arrays/fields in the
 * viz package.  Essentially, they allow vector<vector<T> > types to access data
 * using (i,j) operator overloading.  There is a general field templated type
 * class; specializations exist for vector<vector>.  Other specializations can
 * be added as needed.
 *
 * The generalized class requires the Field Type (FT) template argument to
 * have the following services:
 *
 * \arg operator()(int i, int j) where the range is [0:N-1, 0:N-1];
 * \arg nrows() returns the number of rows (i index);
 * \arg ncols(int row) returns the number of columns in row (j index);
 * \arg FT::value_type defined to the type returned by the field (int, double, etc).
 */
// revision history:
// -----------------
// 0) original
// 1) 28-JAN-00 : added explicit specializations for vector<vector<int>> and
//                vector<vector<double>> because of totalview goofiness
//
//===========================================================================//

template <typename FT> class Viz_Traits {
private:
  //! Reference to the field.
  const FT &field;

public:
  //! Constructor.
  Viz_Traits(const FT &field_in) : field(field_in) { /*...*/
  }

  //! Overloaded operator().
  typename FT::value_type operator()(size_t i, size_t j) const {
    return field(i, j);
  }

  //! Row size accessor.
  size_t nrows() const { return field.nrows(); }

  //! Column size accessor.
  size_t ncols(size_t row) const { return field.ncols(row); }
};

//---------------------------------------------------------------------------//
// Specialization for std::vector<std::vector>

template <typename T> class Viz_Traits<std::vector<std::vector<T>>> {
public:
  // Type traits
  typedef T elementType;

private:
  // Reference to vector<vector> field.
  const std::vector<std::vector<T>> &field;

public:
  // Constructor.
  Viz_Traits(const std::vector<std::vector<T>> &fin) : field(fin) {
    // Nothing to do here
  }

  // Overloaded operator().
  T operator()(size_t i, size_t j) const {
    Require(i < field.size());
    Require(j < field[i].size());
    return field[i][j];
  }

  // Row size accessor.
  size_t nrows() const { return field.size(); }

  // Column size accessor.
  size_t ncols(size_t row) const {
    Require(row < field.size());
    return field[row].size();
  }
};

//---------------------------------------------------------------------------//
// explicit specialization for vector<vector<int>> --> Need this because of
// totalviews inability to handle templated specializations

// template<>
// class Viz_Traits< std::vector<std::vector<int> > >
// {
//   public:
//     // Type traits
//     typedef int elementType;

//   private:
//     // Reference to vector<vector> field.
//     const std::vector<std::vector<int> > &field;

//   public:
//     // Constructor.
//     Viz_Traits(const std::vector<std::vector<int> > &fin) : field(fin)
//     {
//         // Nothing to do here
//     }

//     // Overloaded operator().
//     size_t operator()(size_t i, size_t j) const
//     {
//         Require(i < field.size());
//         Require(j < field[i].size());
//         return field[i][j];
//     }

//     // Row size accessor.
//     size_t nrows() const { return field.size(); }

//     // Column size accessor.
//     size_t ncols(size_t row) const
//     {
//         Require (row < field.size());
//         return field[row].size();
//     }
// };

//---------------------------------------------------------------------------//
// explicit specialization for vector<vector<double>> --> Need this because of
// totalviews inability to handle templated specializations

// template<>
// class Viz_Traits< std::vector<std::vector<double> > >
// {
//   public:
//     // Type traits
//     typedef double elementType;

//   private:
//     // Reference to vector<vector> field.
//     const std::vector<std::vector<double> > &field;

//   public:
//     // Constructor.
//     Viz_Traits(const std::vector<std::vector<double> > &fin) : field(fin)
//     {
//         // Nothing to do here
//     }

//     // Overloaded operator().
//     double operator()(size_t i, size_t j) const
//     {
//         Require(i < field.size());
//         Require(j < field[i].size());
//         return field[i][j];
//     }

//     // Row size accessor.
//     size_t nrows() const { return field.size(); }

//     // Column size accessor.
//     size_t ncols(size_t row) const
//     {
//         Require (row < field.size());
//         return field[row].size();
//     }
// };

} // namespace rtt_viz

#endif // __viz_Viz_Traits_hh__

//---------------------------------------------------------------------------//
// end of viz/Viz_Traits.hh
//---------------------------------------------------------------------------//
