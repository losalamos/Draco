//----------------------------------*-C++-*----------------------------------//
// BilinearInterpTable.hh
// Randy M. Roberts
// Tue Apr  7 12:59:40 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __matprops_BilinearInterpTable_hh__
#define __matprops_BilinearInterpTable_hh__

#include "ds++/Mat.hh"
#include "ds++/Assert.hh"
#include <vector>
#include <algorithm>
#include <functional>

#ifndef BEGIN_NS_XTM
#define BEGIN_NS_XTM namespace XTM  {
#define END_NS_XTM }
#endif

BEGIN_NS_XTM
    
//===========================================================================//
// class BilinearInterpTable - 
//
// Date created : Tue Apr  7 12:59:40 1998
// Purpose      : Binary Interpolation Table
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class BilinearInterpTable
{

    //=======================================================================//
    // NESTED CLASSES AND TYPEDEFS
    //=======================================================================//

  public:
    
    // This class is part of the Memento pattern
    // from p. 283, "Design Patterns", E. Gamma, et. al., 1995
    //
    // This class is meant to be "the internal state" of the
    // BilinearInterpTable (to be held externally),
    // and therefore, must be seen only by that class.
    // I am using it to contain the internal state for a given
    // (x1, x2) pair's interpolation.
    //
    // For faster access to the interpolation table, a user can create
    // a field of Memento to store the intermediate interpolation table
    // calculations.

    class Memento;

  private:

    // Internal class needed for STL algorithms.
    
    template<class T>
    class UnaryInterpolator;

    // Internal class needed for STL algorithms.
    
    template<class T1, class T2>
    class BinaryInterpolator;

    //=======================================================================//
    // DATA
    //=======================================================================//

  private:
    
    // x1vals is the x1 axis of the 2-dimensional grid of tabulated y values.
    std::vector<double> x1vals;

    // x2vals is the x2 axis of the 2-dimensional grid of tabulated y values.
    std::vector<double> x2vals;

    // yvals is the 2-dimensional grid of tabulated y values.
    dsxx::Mat2<double> yvals;
    
    //=======================================================================//
    // CREATORS
    //=======================================================================//
    
  public:

    //------------------------------------------------------------------------//
    // BilinearInterpTable:
    //    Constructor supplying the two axes grids and the two-dimensional
    //    table of evaluations.
    //------------------------------------------------------------------------//

    BilinearInterpTable(const std::vector<double> &x1vals_,
			const std::vector<double> &x2vals_,
			const dsxx::Mat2<double> &yvals_)
	: x1vals(x1vals_), x2vals(x2vals_), yvals(yvals_)
    {
	Require(x1vals.size() >= 2);
	Require(x2vals.size() >= 2);
	Require(yvals.size() == x1vals.size()*x2vals.size());
	Require(axesAreOrdered());
    }
    
    //------------------------------------------------------------------------//
    // BilinearInterpTable:
    //    Constructor supplying the two axes grids and a binary functional
    //    to perform the calculation.
    //------------------------------------------------------------------------//

    template<class BinaryOperation>
    BilinearInterpTable(const std::vector<double> &x1vals_,
			const std::vector<double> &x2vals_,
			BinaryOperation binary_op)
	: x1vals(x1vals_), x2vals(x2vals_), yvals(x1vals.size(),x2vals.size())
    {
	Require(x1vals.size() >= 2);
	Require(x2vals.size() >= 2);
	Require(yvals.size() == x1vals.size()*x2vals.size());
	Require(axesAreOrdered());

	for (int i=0; i<x1vals.size(); i++)
	    for (int j=0; j<x2vals.size(); j++)
		yvals(i,j) = binary_op(x1vals[i], x2vals[j]);
    }
    
    //=======================================================================//
    // MANIPULATORS
    //=======================================================================//
    
    //------------------------------------------------------------------------//
    // setTable:
    //    Manipulator supplying a binary functional to modify the table
    //    y values.
    //------------------------------------------------------------------------//

    template<class BinaryOperation>
    void setTable(BinaryOperation binary_op)
    {
	Require(x1vals.size() >= 2);
	Require(x2vals.size() >= 2);
	Require(yvals.size() == x1vals.size()*x2vals.size());
	Require(axesAreOrdered());

	for (int i=0; i<x1vals.size(); i++)
	    for (int j=0; j<x2vals.size(); j++)
		yvals(i,j) = binary_op(x1vals[i], x2vals[j]);
    }
    
    //=======================================================================//
    // ACCESSORS
    //=======================================================================//

    //------------------------------------------------------------------------//
    // getMemento:
    //    Accessor to obtain the Memento from an x1, x2 interpolation.
    //    This allows the user to avoid the repeated calculation of
    //    interpolation quantities.
    //------------------------------------------------------------------------//

    Memento getMemento(double x1, double x2) const;
    
    inline double interpolate(const Memento &memento) const;

    inline double interpolate(std::pair<double,double> p) const;

    inline double interpolate(double x1, double x2) const;

    template<class InputIterator, class OutputIterator>
    inline void interpolate(InputIterator first,
			    InputIterator last,
			    OutputIterator result) const;
    
    template<class FTIN, class FTOUT>
    inline void interpolate(const FTIN &args, FTOUT &ret_vals) const;

    template<class InputIterator1, class InputIterator2, class OutputIterator>
    inline void interpolate(InputIterator1 first1,
			    InputIterator1 last1,
			    InputIterator2 first2,
			    OutputIterator result) const;
    
    template<class FTIN1, class FTIN2, class FTOUT>
    inline void interpolate(const FTIN1 &arg1s,
			    const FTIN2 &arg2s,
			    FTOUT &ret_vals) const;

  private:
    
    //=======================================================================//
    // IMPLEMENTATION
    //=======================================================================//

    bool axesAreOrdered() const;

};

//===========================================================================//
// class BilinearInterpTable::Memento
//     This class is part of the Memento pattern
//     from "Design Patterns", E. Gamma, et. al., 1995
//
//     This class is meant to be "the internal state" of the
//     BilinearInterpTable (to be held externally),
//     and therefore, must be seen only by that class.
//     I am using it to contain the internal state for a given
//     (x1, x2) pair's interpolation.
//
//     For faster access to the interpolation table, a user can create
//     a field of Memento to store the intermediate interpolation table
//     calculations.
//===========================================================================//

class BilinearInterpTable::Memento
{
    friend class BilinearInterpTable;

    //=======================================================================//
    // DATA
    //=======================================================================//
    
  private:

    // The j index into the x1 axis that is just less than the x1 value.
    int j;

    // The k index into the x2 axis that is just less than the x2 value.
    int k;

    // t and u are defined so that the interpolation to a given (x1,x2)
    // is given by:
    //    y = (1-t)*(1-u)*y1 + t*(1-u)*y2 + t*u*y3 + (1-t)*u*y4
    // where y1 ... y4 are numbered counterclockwise, with y1 in the
    // lower left corner.
    
    double t;
    double u;

    
    //=======================================================================//
    // CREATORS
    //=======================================================================//

    // all Creators are private, so only usable by BilinearInterpTable

    Memento(int j_, int k_, double t_, double u_)
	: j(j_), k(k_), t(t_), u(u_)
    {
	// empty
    }

  public:
    Memento()
    {
	Assert(0);
    }
    
  private:

    //=======================================================================//
    // MANIPULATORS
    //=======================================================================//
    
    // *** none ***
    
    //=======================================================================//
    // ACCESSORS
    //=======================================================================//

    // *** none ***
    
  private:
    
    //=======================================================================//
    // IMPLEMENTATION
    //=======================================================================//

    // *** none ***
};

//=======================================================================//
// INLINE METHODS
//=======================================================================//

inline double BilinearInterpTable::interpolate(std::pair<double,double> p) const
{
    return interpolate(getMemento(p.first, p.second));
}
    
inline double BilinearInterpTable::interpolate(double x1, double x2) const
{
    return interpolate(getMemento(x1, x2));
}

inline double BilinearInterpTable::interpolate(const Memento &memento) const
{
    const double t = memento.t;
    const double u = memento.u;
    const int j = memento.j;
    const int k = memento.k;

    Require(j >= 0 && j < x1vals.size());
    Require(k >= 0 && k < x2vals.size());
    Require(t >= 0.0 && t <= 1.0);
    Require(u >= 0.0 && u <= 1.0);
    
    return (1-t)*(1-u) * yvals(j  ,k)
	+  t*(1-u)     * yvals(j+1,k)
	+  t*u         * yvals(j+1,k+1)
	+  (1-t)*u     * yvals(j  ,k+1);
}

template<class T>
class BilinearInterpTable::UnaryInterpolator
    : public std::unary_function<T, double>
{
  private:
    const BilinearInterpTable &table;
  public:
    UnaryInterpolator(const BilinearInterpTable &table_)
	: table(table_)
    {
	// empty
    }
    double operator()(const T &arg) const
    {
	return table.interpolate(arg);
    }
};
	
template<class InputIterator, class OutputIterator>
inline void BilinearInterpTable::interpolate(InputIterator first,
					     InputIterator last,
					     OutputIterator result) const
{
    typedef typename std::iterator_traits<InputIterator>::value_type value_type;

    // This class turns the interploation table in to a unary functor.

    std::transform(first, last, result, UnaryInterpolator<value_type>(*this));
}

template<class FTIN, class FTOUT>
inline void BilinearInterpTable::interpolate(const FTIN &args,
					     FTOUT &ret_vals) const
{
    Require(args.size() == ret_vals.size());
    interpolate(args.begin(), args.end(), ret_vals.begin());
}

// This class turns the interploation table in to a binary functor.

template<class T1, class T2>
class BilinearInterpTable::BinaryInterpolator
    : public std::binary_function<T1, T2, double>
{
  private:
    const BilinearInterpTable &table;
  public:
    BinaryInterpolator(const BilinearInterpTable &table_)
	: table(table_)
    {
	// empty
    }
    double operator()(const T1 &arg1,
		      const T2 &arg2) const
    {
	return table.interpolate(arg1, arg2);
    }
};
	
template<class InputIterator1, class InputIterator2, class OutputIterator>
inline void BilinearInterpTable::interpolate(InputIterator1 first1,
					     InputIterator1 last1,
					     InputIterator2 first2,
					     OutputIterator result) const
{
    typedef typename std::iterator_traits<InputIterator1>::value_type
	value_type1;
    typedef typename std::iterator_traits<InputIterator2>::value_type
	value_type2;
    
    std::transform(first1, last1, first2, result,
		   BinaryInterpolator<value_type1, value_type2>(*this));
}

template<class FTIN1, class FTIN2, class FTOUT>
inline void BilinearInterpTable::interpolate(const FTIN1 &arg1s,
					     const FTIN2 &arg2s,
					     FTOUT &ret_vals) const
{
    Require(arg1s.size() == arg2s.size());
    Require(arg1s.size() == ret_vals.size());
    interpolate(arg1s.begin(), arg1s.end(), arg2s.begin(), ret_vals.begin());
}

END_NS_XTM  // namespace XTM

#endif                          // __matprops_BilinearInterpTable_hh__

//---------------------------------------------------------------------------//
//                              end of matprops/BilinearInterpTable.hh
//---------------------------------------------------------------------------//
