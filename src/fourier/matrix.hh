//----------------------------------*-C++-*----------------------------------//
// matrix.hh
// John Gulick
// Fri Aug  6 09:51:26 1999
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef matrix_H
#define matrix_H
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream.h>
#include <strings.h>
#include <string>
#include <list>
#include <vector.h>
#include <cstdarg>
#include <iterator>
#include <iterator.h>
#include <algobase.h>
#include <algorithm>
#include <assert.h>
#include <complex>
#include <cstring>
#include "mtl/mtl2lapack.h"
#include "mtl/dense1D.h"
#include "mtl/utils.h"
#include "mtl/mtl.h"
#include "mtl/matrix.h"
#include "gauss.hh"
#include "matrix.t.hh"

using namespace std;
using namespace mtl;
using namespace mtl2lapack;

class Matrix 
{
  public:
    //Matrix constructor.
    Matrix(int mat_order);

    //Member functions of Matrix class.
    lapack_matrix<complex<double> >::type SCB_matrix(double sigma_t, 
						     double sigma_s,
						     int quad_order, 
						     double lambda, 
						     double delta_x);
    lapack_matrix<complex<double> >::type SCB_M4S_matrix(double sigma_t, 
						         double sigma_s,
						         int quad_order, 
						         double lambda, 
						         double delta_x);
    lapack_matrix<complex<double> >::type UCB_matrix(double sigma_t, 
						     double sigma_s,
						     int quad_order, 
						     double lambda, 
						     double delta_x);
    lapack_matrix<complex<double> >::type UCB_M4S_matrix(double sigma_t, 
						         double sigma_s,
						         int quad_order, 
						         double lambda, 
						         double delta_x);
    lapack_matrix<complex<double> >::type LIN_CHAR_matrix(double sigma_t, 
							  double sigma_s,
							  int quad_order, 
							  double lambda, 
							  double delta_x);

    int             mat_order;
    double          sigma_t;
    double          sigma_s;
    double          delta_x;
    int             quad_order;
    double          *mu;
    double          *weight;
    double          gamma;
    double          lambda;
    complex<double> one_one;
    complex<double> one_two;
    complex<double> two_one;
    complex<double> two_two;
    complex<double> *T_pos;
    complex<double> *T_neg;
    complex<double> *P;
    complex<double> *H;

    double get_sigma_t()             { return sigma_t;    }
    double get_sigma_s()             { return sigma_s;    }
    double get_delta_x()             { return delta_x;    }
    int    get_quad_order()          { return quad_order; }
    double *get_mu()                 { return mu;         }
    double *get_weight()             { return weight;     }
    double get_gamma()               { return gamma;      }
    double get_lambda()              { return lambda;     }
    
  private:
    
};


#endif

//---------------------------------------------------------------------------//
//                              end of Matrix.hh
//---------------------------------------------------------------------------//
