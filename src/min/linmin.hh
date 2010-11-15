//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   utils/linmin.hh
 * \author Kent Budge
 * \date   Tue Aug 17 15:30:23 2004
 * \brief  Find minimum of a multivariate function on a specified line.
 * \note   © Copyright 2006 LANSLLC All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef utils_linmin_hh
#define utils_linmin_hh

#include <cmath>
#include <vector>

#include "ds++/Assert.hh"
#include "dbc.hh"
#include "linmin.hh"
#include "mnbrak.hh"
#include "brent.hh"

namespace rtt_utils
{

//---------------------------------------------------------------------------//
/*!
 * \class f1dim
 * \brief Helper template for line minimization
 *
 * \arg \a RandomContainer A random-access container on double.
 * \arg \a Function A function type supporting <code>double
 * operator()(RandomContainer const &)</code>.
 */
template<class RandomContainer, class Function>
class f1dim
{
  public:

    f1dim(Function func,
          RandomContainer &p,
          RandomContainer &xi)
        :
        func(func), p(p), xi(xi)
    {}

    double operator()(double const x) const
    {
        unsigned const n = p.size();
        std::vector<double> xt(n);
        for (unsigned i=0; i<n; ++i)
        {
            xt[i] = p[i] + x*xi[i];
        }
        return func(xt);
    }

  private:

    Function func;
    RandomContainer &p, &xi;
};

//---------------------------------------------------------------------------//
/*!
 * Find minimum of a multivariate function on a specified line.
 *
 * \arg \a RandomContainer A random-access container on double.
 * \arg \a Function A function type supporting <code>double
 * operator()(RandomContainer const &)</code>.
 *
 * \param[in,out] p On entry, contains a starting guess of the minimum.  On
 * exit, contains an improved estimate of the minimum.
 * \param[in,out] xi On entry, contains the search direction.  On exit,
 * contains the actual displacement to the minimum in the search direction.
 * \param[out] fret Final minimum value of the function.
 * \param[in] func Function to be minimized
 */

template <class RandomContainer, class Function>
void linmin(RandomContainer &p,
            RandomContainer &xi,
            double &fret,
            Function func)
{
    using std::numeric_limits;

    double const TOL = sqrt(numeric_limits<double>::epsilon());

    unsigned const n = p.size();

    // Initial guess for brackets
    double ax = 0.0;
    double xx = 1.0;

    double bx, fa, fx, fb, xmin;
    f1dim<RandomContainer, Function> f1(func, p, xi);
    mnbrak<f1dim<RandomContainer, Function> >(ax,
                                              xx,
                                              bx,
                                              fa,
                                              fx,
                                              fb,
                                              f1);
    
    fret = brent(ax, xx, bx, f1, TOL, xmin);
    for (unsigned j=0; j<n; ++j)
    {
        xi[j] *= xmin;
        p[j] += xi[j];
    }
}

} // end namespace rtt_utils

#endif // utils_linmin_hh

//---------------------------------------------------------------------------//
//              end of utils/linmin.hh
//---------------------------------------------------------------------------//
