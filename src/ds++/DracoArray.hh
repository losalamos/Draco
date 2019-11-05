//----------------------------------*-C++-*----------------------------------//  
/*!                                                                              
 * \file   ds++/DracoArray.hh                                                  
 * \author Ben R. Ryan <brryan@lanl.gov                                       
 * \date   2019 Nov 5                                     
 * \brief  Multidimensional container.                             
 * \note   Copyright (C) 2017-2019 Triad National Security, LLC.                 
 *         All rights reserved. */                                               
//---------------------------------------------------------------------------//

#ifndef rtt_dsxx_DracoArray_hh
#define rtt_dsxx_DracoArray_hh

#include <vector>

namespace rtt_dsxx {

//========================================================================
/*!
 * \class DracoArray
 *
 * \brief This is a container class that provides a convenient 
 * multidimensional interpretation of an std::vector<T>.
 *
 */
//========================================================================

} // namespace rtt_dsxx

template <typename T> class DracoArray {
  public:
  MultiArray<T>(int n1_in);
  MultiArray<T>(int n1_in, int n2_in);
  MultiArray<T>(int n1_in, int n2_in, int n3_in);

  // Overload () operator for data access
  T &operator()(int i) { return data[i]; }
  T operator()(int i) const { return data[i]; }
  T &operator()(int i, int j) { return data[j + n2*i]; }
  T operator()(int i, int j) const { return data[j + n2*i]; }
  T &operator()(int i, int j, int k) { return data[k + n3*(j + n2*i)]; }
  T operator()(int i, int j, int k) const { return data[k + n3*(j + n2*i)]; }

  // Accessor functions
  int get_dim() const { return dim; }
  int get_n1() const { return n1; }
  int get_n2() const { return n2; }
  int get_n3() const { return n3; }
  std::vector<int> get_size() {
    std::vector<int> size_all{n1, n2, n3};
    std::vector<int> size;
    for (int n = 0; n < dim; n++) {
      size.push_back(size_all[n]);
      return size;
    }
  }
  void zero() {
    std::fill(data.begin(), data.end(), 0);
  }
  bool is_allocated() const { return is_allocated; }

  private:
  int dim;
  int n1, n2, n3;
  std::vector<T> data;
  bool is_allocated = false;

  template <typename T> DracoArray<T>::DracoArray(int n1_in)
  {
    dim = 1;
    n1 = n1_in;
    n2 = 1;
    n3 = 1;
    data.resize(n1);
    is_allocated = true;
  }

  template <typename T> DracoArray<T>::DracoArray(int n1_in, int n2_in)
  {
    dim = 2;
    n1 = n1_in;
    n2 = n2_in;
    n3 = 1;
    data.resize(n1*n2);
    is_allocated = true;
  }

  template <typename T> DracoArray<T>::DracoArray(int n1_in, int n2_in, int n3_in)
  {
    dim = 3;
    n1 = n1_in;
    n2 = n2_in;
    n3 = n3_in;
    data.resize(n1*n2*n3);
    is_allocated = true;
  }
};

#endif // rtt_dsxx_DracoArray_hh

//---------------------------------------------------------------------------//
// end of DracoArray.hh
//---------------------------------------------------------------------------//
