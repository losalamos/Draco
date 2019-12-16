//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/VectorView.hh
 * \author Kelly G. Thompson <kgt@lanl.gov>, Steve Nolen <sdnolen@lanl.gov>
 * \date   Friday, Dec 13, 2019, 17:49 pm
 * \brief  A view class that presents contiguous 1-D as multi-D.
 * \note   Copyright (C) 2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#ifndef rtt_dsxx_VectorView_hh
#define rtt_dsxx_VectorView_hh

#include "Assert.hh"
#include <array>
#include <vector>

namespace rtt_dsxx {

//============================================================================//
/*!
 * \class VectorView
 * \brief A wrapper for std::vector<T> that presents a compatible multi-D
 *        'view'.
 *
 * \tparam T The data storage for the underlying std::vector<T>.
 * \tparam Extents A variatic template list of integral values that represent
 *         the max value of each dimension of the 'view'.
 * \todo Add operator= to allow assignment operations.
 *
 * Example:
 *
 * \code
 * vector<double> v1 = { 1,2,3,4,5,6,7,8,9 };
 * VectorView<double,3,3,3> myvecview(v1);
 * std::cout << myvecview(2,0,1) << endl;
 *
 * \endcode
 */
//============================================================================//

template <typename T, unsigned... Extents> class VectorView {
public:
  //--------------------------------------------------------------------------//
  /*! \brief Compute the 1-D index associated with the actual data storage from
   *         provided multi-D index values.
   *
   * \tparam Indices variatic template that represents a set of integral values,
   *         one for each dimension of the 'view'.
   * \param indices a list of integral types that indicate the multi-D offset.
   * \return An size_t index into the underlying std::vector<T>.
   *
   * \todo Once we adopt c++17, this logic can be simplified using 'fold
   *       expressions'.
   */
  template <typename... Indices> size_t idx(Indices... indices) {
    Insist(sizeof...(indices) == dim,
           std::string("VectorView expected ") + std::to_string(dim) +
               " arguments but only found " +
               std::to_string(sizeof...(Indices)) + ".");
    size_t idx = 0, i = 0;
    using untilFold = int[];
    // In the following line, the 'pack' is indices and the operation is the
    // comma. For an explanation see
    // https://www.codingame.com/playgrounds/2205/7-features-of-c17-that-will-simplify-your-code/fold-expressions
    // Note that we don't use the computed value 'untilFold', but its
    // construction has the side effect of setting idx.
    (void)untilFold{((idx += indices * strides[i++]), 0)...};
    Check(idx < data.size());
    return idx;
  }

  //--------------------------------------------------------------------------//
  /*! \brief Return the data associated with the provided multi-D view.
   *
   * \tparam Indices variatic template that represents a set of integral values,
   *         one for each dimension of the 'view'.
   * \param indices a list of integral types that indicate the multi-D offset.
   * \return Data from the 1-D std::vector<T> associated with the provided
   *         multi-D index set.
   */
  template <typename... Indices> T operator()(Indices... indices) {
    return data[idx(indices...)];
  }

  //--------------------------------------------------------------------------//
  /*! \brief Default constructor
   *
   * \param[in] data a reference to the underlying std::vector<T> that we are
   *               creating a 'view' for.
   *
   * When the constructor is called, we must setup the 'stride' or offset values
   * for each of indices that will be used for this view.
   */
  VectorView(std::vector<T> const &data_) : data(data_) {
    Ensure(dim > 0);
    // Member data 'strides' was initialized with Extents, but must now be
    // replaced with actual offset values.
    for (size_t i = dim - 1; i < dim; --i) {
      Insist(strides[i] > 0, "Each dimension must have a positive length.");
      strides[i] = 1;
      for (size_t j = 0; j < i; j++) {
        strides[i] *= strides[j];
      }
    }
  }

private:
  //
  // Disabled member functions to ensure correct behavior.
  //

  // Ensure that the idx member function is not defined for the template case
  // of zero arguments.
  size_t idx(void) = delete;
  // Ensure that the paren-operator member function is not defined for the
  // template case of zero arguments.
  T operator()() = delete;

  //
  // >> MEMBER DATA (state)
  //

  //! The actual data container that we are creating a view for.
  std::vector<T> const &data;
  //! Max size of each dimension for the view.
  size_t const dim = sizeof...(Extents);
  //! Offsets associated with each of the view's dimensions.
  std::array<size_t, sizeof...(Extents)> strides = {Extents...};
};

} // end namespace rtt_dsxx

#endif // rtt_dsxx_VectorView_hh

//---------------------------------------------------------------------------//
// end of ds++/VectorView.hh
//---------------------------------------------------------------------------//
