//----------------------------------*-C++-*----------------------------------//
/*! 
 * \file fourier.hh
 * \author John Gulick
 * \date Tue Aug  3 13:44:08 1999
 * \brief Fourier package header file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef fourier_H
#define fourier_H
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
#include "gauss.hh"
#include "matrix.hh"
#include "cpgplot.h"

using namespace std;

//===========================================================================//
/*! 
 * \class Fourier
 * \brief Fourier analysis class solver.
 *
 * This class does cool stuff...
 */
//===========================================================================//


class Fourier
{
  public:
    Fourier (){}

    //! Input member function.
    void input();

    //! Solve member function.
    void solve();

    //! Plot member function.
    void plot();

    //! Print member function.
    void print();
    
  private:

    string method;
    string acceleration_type;
    string analysis_type;
    double sigma_t;
    double sigma_s;
    double delta_x;
    int    quad_order;
    double lambda_begin;
    double lambda_end;
    double lambda_step_size;
    double delta_x_begin;
    int    num_delta_x;
    int    step_counter;
    int    tot_num_delta_x;
    int    *region_number;
    int    *num_points;
    double *dx_step_size;
    double *mu,*weight;
    double *spec_rad;
    double max_lam_spec_rad;
    double max_dx_spec_rad;
    Matrix *matrix_ptr;
    double gamma;
    int    num_lambda;
    int    dx_count;
    double *eigen_test;
    double lambda;
    double max_spec_rad_lambda;
    double max_spec_rad_dx;
    double *spec_rad_lam;
    double *spec_rad_dx;

};


#endif

//---------------------------------------------------------------------------//
//                              end of fourier.hh
//---------------------------------------------------------------------------//
