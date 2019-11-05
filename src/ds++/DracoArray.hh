//----------------------------------*-C++-*----------------------------------//
/*!                                                                              
 * \file   ds++/DracoArray.hh                                                  
 * \author Ben R. Ryan <brryan@lanl.gov> 
 * \date   2019 Nov 5                                     
 * \brief  Multidimensional container.                             
 * \note   Copyright (C) 2017-2019 Triad National Security, LLC.                 
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#ifndef rtt_dsxx_DracoArray_hh
#define rtt_dsxx_DracoArray_hh

#include "Assert.hh"
#include <vector>

namespace rtt_dsxx {

//========================================================================
/*!
 * \class DracoArray
 *
 * \brief This is a container class that provides a convenient 
 * multidimensional interpretation of an std::vector<T>.
 *
 * Based on the AthenaArray implemention in Athena++
 * https://princetonuniversity.github.io/athena/
 *
 */
//========================================================================

template <typename T> class DracoArray {
public:
  // Overload () operator for data access
  T &operator()(int i) { return data[i]; }
  T operator()(int i) const { return data[i]; }
  T &operator()(int i, int j) { return data[j + n2 * i]; }
  T operator()(int i, int j) const { return data[j + n2 * i]; }
  T &operator()(int i, int j, int k) { return data[k + n3 * (j + n2 * i)]; }
  T operator()(int i, int j, int k) const {
    return data[k + n3 * (j + n2 * i)];
  }

  // Accessor functions
  int get_dim() const { return dim; }
  int get_n1() const {
    Require(dim > 0);
    return n1;
  }
  int get_n2() const {
    Require(dim > 1);
    return n2;
  }
  int get_n3() const {
    Require(dim > 2);
    return n3;
  }
  std::vector<int> get_size() {
    std::vector<int> size_all{n1, n2, n3};
    std::vector<int> size;
    for (int n = 0; n < dim; n++) {
      size.push_back(size_all[n]);
      return size;
    }
  }
  bool is_allocated() const { return allocated; }

private:
  int dim;
  int n1, n2, n3;
  std::vector<T> data;
  bool allocated = false;

  // Constructors and memory management
public:
  DracoArray(int n1_in) {
    dim = 1;
    n1 = n1_in;
    data.resize(n1);
    allocated = true;
  }

  DracoArray(int n1_in, int n2_in) {
    dim = 2;
    n1 = n1_in;
    n2 = n2_in;
    data.resize(n1 * n2);
    allocated = true;
  }

  DracoArray(int n1_in, int n2_in, int n3_in) {
    dim = 3;
    n1 = n1_in;
    n2 = n2_in;
    n3 = n3_in;
    data.resize(n1 * n2 * n3);
    allocated = true;
  }
  void zero() { std::fill(data.begin(), data.end(), 0); }
};

} // namespace rtt_dsxx

#endif // rtt_dsxx_DracoArray_hh

//---------------------------------------------------------------------------//
// end of DracoArray.hh
//---------------------------------------------------------------------------//
