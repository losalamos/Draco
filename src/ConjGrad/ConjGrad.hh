//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ConjGrad/ConjGrad.hh
 * \author Randy M. Roberts
 * \date   Mon Apr 24 08:16:51 2000
 * \brief  The basic CG algorithm.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __ConjGrad_ConjGrad_hh__
#define __ConjGrad_ConjGrad_hh__

#include "ConjGradTraits.hh"
#include "NoPrecon.hh"
#include "DefaultConvCrit.hh"

namespace rtt_ConjGrad
{
 
//===========================================================================//
/*!
 * \function ConjGrad
 *
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

template<class Field, class MatVec, class ConvCrit, class Precon,
    class ScalarMult, class DotProd, class Add, class Sub>
void conjGrad(Field &x, int &iter, const Field &b, MatVec matVec,
	      ConvCrit converged, Precon precon, ScalarMult scalarMult,
	      DotProd dotProd, Add add, Sub sub, Field &r, Field &z, Field &p,
	      Field &q)
{
    r = sub(b, matVec(x));

    typedef typename Field::value_type value_type;
    value_type rho0;
    
    for (iter=1; !converged(iter, r, x, b); ++iter)
    {
	z = precon(r);
	value_type rho1 = dotProd(r, z);
	if (iter == 1)
	{
	    p = z;
	}
	else
	{
	    value_type beta = rho1/rho0;
	    p = add(z, scalarMult(beta, p));
	}
	q = matVec(p);
	value_type alpha = rho1/dotProd(p, q);
	x = add(x, scalarMult(alpha, p));
	r = sub(r, scalarMult(alpha, q));

	rho0 = rho1;
    }
}

template<class Field, class MatVec, class ConvCrit, class Precon>
inline void conjGrad(Field &x, int &iter, const Field &b, MatVec matVec,
		     ConvCrit converged, Precon precon, Field &r)
{
    Field z = ConjGradTraits<Field>::create(b);
    Field p = ConjGradTraits<Field>::create(b);
    Field q = ConjGradTraits<Field>::create(b);

    conjGrad(x, iter, b, matVec, converged, precon,
	     ConjGradTraits<Field>::ScalarMultiplies(),
	     ConjGradTraits<Field>::Dot(),
	     ConjGradTraits<Field>::Add(),
	     ConjGradTraits<Field>::Sub(), r, z, p, q);
}

// High level interface to conjGrad with no preconditioner

template<class Field, class MatVec, class ConvCrit>
inline void conjGrad(Field &x, int &iter, const Field &b, MatVec matVec,
		     ConvCrit converged, Field &r)
{
    conjGrad(x, iter, b, matVec, converged, NoPrecon<Field>(), r);
}

// High level interface to conjGrad with default convergence criteria

template<class Field, class MatVec, class Precon>
inline void conjGrad(Field &x, int &iter, const Field &b, MatVec matVec,
		     int maxIters, typename Field::value_type eps,
		     Precon precon, Field &r)
{
    typedef ConjGradTraits<Field>::Norm Norm;
    conjGrad(x, iter, b, matVec,
	     DefaultConvCrit<Field, Norm>(Norm(), maxIters, eps), precon, r);
}

// High level interface to conjGrad with no preconditioner and
// default convergence criteria

template<class Field, class MatVec>
inline void conjGrad(Field &x, int &iter, const Field &b, MatVec matVec,
		     int maxIters, typename Field::value_type eps,
		     Field &r)
{
    typedef ConjGradTraits<Field>::Norm Norm;
    conjGrad(x, iter, b, matVec,
	     DefaultConvCrit<Field, Norm>(Norm(), maxIters, eps), r);
}

} // end namespace rtt_ConjGrad

#endif                          // __ConjGrad_ConjGrad_hh__

//---------------------------------------------------------------------------//
//                              end of ConjGrad/ConjGrad.hh
//---------------------------------------------------------------------------//
