//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ConjGrad/ConjGradTraits.hh
 * \author Randy M. Roberts
 * \date   Mon Apr 24 09:14:20 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __ConjGrad_ConjGradTraits_hh__
#define __ConjGrad_ConjGradTraits_hh__

#include "DefaultScalarMultiplies.hh"

#include "traits/MT_traits.hh"
#include "ds++/Mat.hh"
#include "c4/global.hh"
#include <functional>
#include <cmath>

namespace rtt_ConjGrad
{
 
//===========================================================================//
/*!
 * \class ConjGradTraits
 *
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

template<class Field>
class ConjGradTraits 
{

    // NESTED CLASSES AND TYPEDEFS

  public:
    
    typedef typename Field::value_type value_type;
    
    typedef DefaultScalarMultiplies<value_type, Field> ScalarMultiplies;
    typedef std::plus<Field> Add;
    typedef std::minus<Field> Sub;

    struct Dot : public std::binary_function<Field,Field,value_type>
    {
	value_type operator()(const Field &f1, const Field &f2) const
	{
	    return rtt_traits::vector_traits<Field>::dot(f1, f2);
	}
    };

    struct Norm : public std::unary_function<Field,value_type>
    {
	value_type operator()(const Field &f) const
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

    static Field create(const Field &x) { return x; }

  private:
    
    // IMPLEMENTATION
};

// Specialization for Mat1<T>

template<class T>
class ConjGradTraits<rtt_dsxx::Mat1<T> >
{

    // NESTED CLASSES AND TYPEDEFS

  public:

    typedef rtt_dsxx::Mat1<T> Mat;
    typedef typename Mat::value_type value_type;

    struct ScalarMultiplies : public std::binary_function<value_type,Mat,Mat>
    {
	Mat operator()(const value_type v, const Mat &f) const
	{
	    Mat ret(f);
	    ret *= v;
	    return ret;
	}
    };

    struct Add : public std::binary_function<Mat,Mat,Mat>
    {
	Mat operator()(const Mat &f1, const Mat &f2) const
	{
	    Mat ret(f1);
	    ret += f2;
	    return ret;
	}
    };

    struct Sub : public std::binary_function<Mat,Mat,Mat>
    {
	Mat operator()(const Mat &f1, const Mat &f2) const
	{
	    Mat ret(f1);
	    ret -= f2;
	    return ret;
	}
    };

    struct Dot : public std::binary_function<Mat,Mat,value_type>
    {
	value_type operator()(const Mat &f1, const Mat &f2) const
	{
	    f1.assert_conformality(f2);
	
	    value_type ret = value_type();
	    typename Mat::const_iterator it1 = f1.begin();
	    typename Mat::const_iterator it2 = f2.begin();
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

    struct Norm : public std::unary_function<Mat,value_type>
    {
	value_type operator()(const Mat &f) const
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

    static Mat create(const Mat &x) { return x; }

  private:
    
    // IMPLEMENTATION
};

} // end namespace rtt_ConjGrad

#endif                          // __ConjGrad_ConjGradTraits_hh__

//---------------------------------------------------------------------------//
//                              end of ConjGrad/ConjGradTraits.hh
//---------------------------------------------------------------------------//
