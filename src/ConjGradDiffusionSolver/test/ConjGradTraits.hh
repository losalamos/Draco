//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ConjGradDiffusionSolver/ConjGradTraits.hh
 * \author Randy M. Roberts
 * \date   Mon Apr 24 16:58:46 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __ConjGradDiffusionSolver_ConjGradTraits_hh__
#define __ConjGradDiffusionSolver_ConjGradTraits_hh__

#include "ConjGrad/ConjGradTraits.hh"
#include "mesh/Mesh_XYZ.hh"

namespace rtt_ConjGrad
{
 
// Specialization for SolverP1Diff::ccsf

template<class T>
class ConjGradTraits<Mesh_XYZ:: template cctf<T> >
{

    // NESTED CLASSES AND TYPEDEFS

  public:

    typedef typename Mesh_XYZ:: template cctf<T> cctf;
    typedef typename cctf::value_type value_type;

    struct ScalarMultiplies : public std::binary_function<value_type,cctf,cctf>
    {
	cctf operator()(const value_type v, const cctf &f) const
	{
	    cctf ret(f);
	    ret *= v;
	    return ret;
	}
    };

    struct Add : public std::binary_function<cctf,cctf,cctf>
    {
	cctf operator()(const cctf &f1, const cctf &f2) const
	{
	    cctf ret(f1);
	    ret += f2;
	    return ret;
	}
    };

    struct Sub : public std::binary_function<cctf,cctf,cctf>
    {
	cctf operator()(const cctf &f1, const cctf &f2) const
	{
	    cctf ret(f1);
	    ret -= f2;
	    return ret;
	}
    };

    struct Dot : public std::binary_function<cctf,cctf,value_type>
    {
	value_type operator()(const cctf &f1, const cctf &f2) const
	{
	    value_type ret = value_type();
	    typename cctf::const_iterator it1 = f1.begin();
	    typename cctf::const_iterator it2 = f2.begin();
	    while (it1 != f1.end())
	    {
		ret = ret + (*it1)*(*it2);
		++it1;
		++it2;
	    }

	    C4::gsum(ret);
	    return ret;
	}
    };

    struct Norm : public std::unary_function<cctf,value_type>
    {
	value_type operator()(const cctf &f) const
	{
	    return std::sqrt(Dot()(f,f));
	}
    };

  private:
    
    // DATA
    
  public:

    // CREATORS
    
    //Defaulted: ConjGradTraits();
    //Defaulted: ConjGradTraits(const ConjGradTraits &rhs);
    //Defaulted: ~ConjGradTraits();

    // MANIPULATORS
    
    //Defaulted: ConjGradTraits& operator=(const ConjGradTraits &rhs);

    // ACCESSORS

    static cctf create(const cctf &x) { return x; }

  private:
    
    // IMPLEMENTATION
};
} // end namespace rtt_ConjGrad

#endif // __ConjGradDiffusionSolver_ConjGradTraits_hh__

//---------------------------------------------------------------------------//
// end of ConjGradDiffusionSolver/ConjGradTraits.hh
//---------------------------------------------------------------------------//
