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
 * \brief A function to perform conjugate gradient calculations.
 *
 * There are several templated conjGrad functions of varying degrees of
 * generality. Users should normally use one of the high level interfaces,
 * which will make use of various default parameters and implementations. The 
 * high level functions use a cascade of calls to lower level wrappers until
 * the lowest level function actually implements the CG functionality in a
 * very general manner. Because the cascade of calls passes certain
 * user-defined objects by value (e.g. MatVec), these probably should be
 * lightweight wrappers to the actual implementation classes in order to
 * avoid excessive copying costs.
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

// Medium level interface that uses traits for Add, Dot, etc.

template<class Field, class MatVec, class ConvCrit, class Precon, class Reduce>
inline void conjGrad(Field &x, int &iter, const Field &b, MatVec matVec,
		     ConvCrit converged, Precon precon, Field &r,
		     Reduce reduce)
{
    Field z = ConjGradTraits<Field>::create(b);
    Field p = ConjGradTraits<Field>::create(b);
    Field q = ConjGradTraits<Field>::create(b);

    conjGrad(x, iter, b, matVec, converged, precon,
	     ConjGradTraits<Field>::ScalarMultiplies(),
	     ConjGradTraits<Field>::Dot<Reduce>(reduce),
	     ConjGradTraits<Field>::Add(),
	     ConjGradTraits<Field>::Sub(), r, z, p, q);
}

// High level interface to conjGrad with no preconditioner

template<class Field, class MatVec, class ConvCrit, class Reduce>
inline void conjGrad(Field &x, int &iter, const Field &b, MatVec matVec,
		     ConvCrit converged, Field &r, Reduce reduce)
{
    // Calls Medium level interface that uses traits for Add, Dot, etc.
    conjGrad(x, iter, b, matVec, converged, NoPrecon<Field>(), r, reduce);
}

// High level interface to conjGrad with default convergence criteria

template<class Field, class MatVec, class Precon, class Reduce>
inline void conjGrad(Field &x, int &iter, const Field &b, MatVec matVec,
		     int maxIters, typename Field::value_type eps,
		     Precon precon, Field &r, Reduce reduce)
{
    typedef ConjGradTraits<Field>::Norm<Reduce> Norm;

    // Calls Medium level interface that uses traits for Add, Dot, etc.
    conjGrad(x, iter, b, matVec,
	     DefaultConvCrit<Field, Norm>(Norm(reduce), maxIters, eps),
	     precon, r, reduce);
}

// High level interface to conjGrad with no preconditioner and
// default convergence criteria

template<class Field, class MatVec, class Reduce>
inline void conjGrad(Field &x, int &iter, const Field &b, MatVec matVec,
		     int maxIters, typename Field::value_type eps,
		     Field &r, Reduce reduce)
{
    typedef ConjGradTraits<Field>::Norm<Reduce> Norm;

    // Calls High level interface to conjGrad with no preconditioner
    conjGrad(x, iter, b, matVec,
	     DefaultConvCrit<Field, Norm>(Norm(reduce), maxIters, eps),
	     r, reduce);
}

} // end namespace rtt_ConjGrad

#endif                          // __ConjGrad_ConjGrad_hh__

//---------------------------------------------------------------------------//
//                              end of ConjGrad/ConjGrad.hh
//---------------------------------------------------------------------------//
