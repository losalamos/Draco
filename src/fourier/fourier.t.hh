//----------------------------------*-C++-*----------------------------------//
// fourier.t.hh
// John Gulick
// Tue Aug  3 13:44:08 1999
// $Id$
//---------------------------------------------------------------------------//
// @> 
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
#include "gauss.t.hh"
#include "matrix.t.hh"
#include "cpgplot.h"

using namespace std;

class Fourier
{
  public:
    Fourier (){}

void input();
void solve();
    //void plot();
    //void print();

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
    int    tot_num_delta_x;
    int    *region_number;
    int    *num_points;
    double *dx_step_size;
    double *mu,*weight;
    double *spec_rad;
    double max_spec_rad;
    Matrix *matrix_ptr;
    double gamma;
    int num_lambda;
    int count;
    double *eigen_test;
    double lambda;
    double max_spec_rad_lambda;
    double *spec_rad_lam;
    double *spec_rad_dx;

    

  private:

};


#endif

//---------------------------------------------------------------------------//
//                              end of fourier.t.hh
//---------------------------------------------------------------------------//
