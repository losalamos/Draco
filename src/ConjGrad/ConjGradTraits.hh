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

#include "traits/ContainerTraits.hh"
#include "ds++/Mat.hh"
#include "ds++/Assert.hh"
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
    typedef rtt_traits::ContainerTraits<Field> ContainerTraits;
    
    struct ScalarMultiplies
	: public std::binary_function<value_type,Field,Field>
    {
	Field operator()(const value_type v, const Field &f) const
	{
	    Field ret(f);
	    for (typename Field::iterator rit = ret.begin();
		 rit != ret.end(); ++rit)
	    {
		*rit = v * (*rit);
	    }
	    return ret;
	}
    };

    struct Add : public std::binary_function<Field,Field,Field>
    {
	Field operator()(const Field &f1, const Field &f2) const
	{
	    Assert(ContainerTraits::conformal(f1, f2));
	    
	    Field ret(f1);

	    typename Field::iterator rit = ret.begin();
	    typename Field::const_iterator f2it = f2.begin();
	    while (rit != ret.end())
	    {
		*rit = *rit + *f2it;
		++rit;
		++f2it;
	    }
	    return ret;
	}
    };

    struct Sub : public std::binary_function<Field,Field,Field>
    {
	Field operator()(const Field &f1, const Field &f2) const
	{
	    Assert(ContainerTraits::conformal(f1, f2));
	    
	    Field ret(f1);

	    typename Field::iterator rit = ret.begin();
	    typename Field::const_iterator f2it = f2.begin();
	    while (rit != ret.end())
	    {
		*rit = *rit - *f2it;
		++rit;
		++f2it;
	    }
	    return ret;
	}
    };

    template <class Reduce>
    class Dot : public std::binary_function<Field,Field,value_type>
    {
      private:

	Reduce reducer;

      public:

	Dot(const Reduce reducer_in) : reducer(reducer_in) { /* empty */ }

	value_type operator()(const Field &f1, const Field &f2) const
	{
	    Assert(ContainerTraits::conformal(f1, f2));
	    
	    value_type ret = value_type();
	    typename Field::const_iterator it1 = f1.begin();
	    typename Field::const_iterator it2 = f2.begin();
	    while (it1 != f1.end())
	    {
		ret = ret + (*it1)*(*it2);
		++it1;
		++it2;
	    }

	    ret = reducer.sum(ret);
	    return ret;
	}
    };

    template <class Reduce>
    class Norm : public std::unary_function<Field,value_type>
    {
      private:

	Reduce reducer;

      public:

	Norm(const Reduce reducer_in) : reducer(reducer_in) { /* empty */ }

	value_type operator()(const Field &f) const
	{
	    return std::sqrt(Dot<Reduce>(reducer)(f,f));
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

    template <class Reduce>
    class Dot : public std::binary_function<Mat,Mat,value_type>
    {
      private:

	Reduce reducer;

      public:

	Dot(const Reduce reducer_in) : reducer(reducer_in) { /* empty */ }

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

	    ret = reducer.sum(ret);
	    return ret;
	}
    };

    template <class Reduce>
    class Norm : public std::unary_function<Mat,value_type>
    {
      private:

	Reduce reducer;

      public:

	Norm(const Reduce reducer_in) : reducer(reducer_in) { /* empty */ }

	value_type operator()(const Mat &f) const
	{
	    return std::sqrt(Dot<Reduce>(reducer)(f,f));
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
