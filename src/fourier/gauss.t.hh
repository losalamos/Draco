//----------------------------------*-C++-*----------------------------------//
// gauss.t.hh
// John Gulick
// Tue Aug  3 13:28:20 1999
// $Id$
//---------------------------------------------------------------------------//
// @> Creates Gaussian Legendre Quadrature for One Dimensional Fourier 
//    Analysis.  Algorithm:  Numerical Recipes, Press.
//---------------------------------------------------------------------------//
#ifndef guass_H
#define guass_H
#include <iostream.h>
#include <math.h>


class Gauss
{
public:
    Gauss(){}

    void gauleg(int n);

    double *x;
    double *w;
    double gamma;

    double *get_mu()     { return x;     }
    double *get_weight() { return w;     }
    double get_gamma()   { return gamma; }

private:
    double x1;
    double x2;
    double z1,z,xm,xl,pp,p3,p2,p1;
    int m,j,i;
};

#endif


//---------------------------------------------------------------------------//
//                              end of gauss.t.hh
//---------------------------------------------------------------------------//
