//----------------------------------*-C++-*----------------------------------//
// BilinearInterpTable.hh
// Randy M. Roberts
// Tue Apr  7 12:59:40 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __matprops_BilinearInterpTable_hh__
#define __matprops_BilinearInterpTable_hh__

#include "BilinearInterpGrid.hh"
#include "ds++/Mat.hh"
#include "ds++/SP.hh"
#include "ds++/Assert.hh"
#include <algorithm>
#include <functional>

namespace rtt_matprops {
    
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
    
    typedef BilinearInterpGrid::Memento Memento;

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

    // The grid axis values on which to do the interpolation.
    
    dsxx::SP<BilinearInterpGrid> grid;
    
    // yvals is the 2-dimensional grid of tabulated y values.

    dsxx::Mat2<double> yvals;

    // Is this an empty table?

    bool isEmpty;
    
    //=======================================================================//
    // CREATORS
    //=======================================================================//
    
  public:

    //------------------------------------------------------------------------//
    // BilinearInterpTable:
    //    Creates an empty table.
    //------------------------------------------------------------------------//

    BilinearInterpTable()
	: isEmpty(true)
    {
	// *** empty ***
    }
    
    //------------------------------------------------------------------------//
    // BilinearInterpTable:
    //    Constructor supplying the two axes grids and the two-dimensional
    //    table of evaluations.
    //------------------------------------------------------------------------//

    BilinearInterpTable(const std::vector<double> &x1vals_,
			const std::vector<double> &x2vals_,
			const dsxx::Mat2<double> &yvals_);
    
    //------------------------------------------------------------------------//
    // BilinearInterpTable:
    //    Constructor supplying the two axes grids and a binary functional
    //    to perform the calculation.
    //------------------------------------------------------------------------//

    template<class BinaryOperation>
    inline BilinearInterpTable(const std::vector<double> &x1vals_,
			       const std::vector<double> &x2vals_,
			       BinaryOperation binary_op);
    
    //------------------------------------------------------------------------//
    // BilinearInterpTable:
    //    Constructor supplying the grid and the two-dimensional
    //    table of evaluations.
    //------------------------------------------------------------------------//

    BilinearInterpTable(const dsxx::SP<BilinearInterpGrid> &grid_,
			const dsxx::Mat2<double> &yvals_);
    
    //------------------------------------------------------------------------//
    // BilinearInterpTable:
    //    Constructor supplying the grid and a binary functional
    //    to perform the calculation.
    //------------------------------------------------------------------------//

    template<class BinaryOperation>
    inline BilinearInterpTable(const dsxx::SP<BilinearInterpGrid> &grid_,
			       BinaryOperation binary_op);
    
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
	Require(yvals.size() == grid->size());

	for (int i=0; i<grid->size(1); i++)
	    for (int j=0; j<grid->size(2); j++)
		yvals(i,j) = binary_op(grid->x1(i), grid->x2(j));
    }
    
    //=======================================================================//
    // ACCESSORS
    //=======================================================================//

    //------------------------------------------------------------------------//
    // hasData:
    //    Return true if this table has data.
    //------------------------------------------------------------------------//

    bool hasData() const { return !isEmpty; }
    
    //------------------------------------------------------------------------//
    // getGrid:
    //    Return a const reference to the grid.
    //------------------------------------------------------------------------//

    const dsxx::SP<BilinearInterpGrid> getGrid() const { return grid; }
    
    //------------------------------------------------------------------------//
    // interpolate:
    //    Given a Memento, return an interpolated value.
    //------------------------------------------------------------------------//

    inline double interpolate(const Memento &memento) const;

    //------------------------------------------------------------------------//
    // interpolate:
    //    Given a pair of doubles, return an interpolated value.
    //------------------------------------------------------------------------//

    inline double interpolate(std::pair<double,double> p) const;

    //------------------------------------------------------------------------//
    // interpolate:
    //    Given a pair of doubles, return an interpolated value.
    //------------------------------------------------------------------------//

    inline double interpolate(double x1, double x2) const;

    //------------------------------------------------------------------------//
    // interpolate:
    //    Given begin and end iterators of "things" (probably Memento's
    //    or std::pair<double,double>'s), return interpolations to the
    //    result iterator.
    //------------------------------------------------------------------------//

    template<class InputIterator, class OutputIterator>
    inline void interpolate(InputIterator first,
			    InputIterator last,
			    OutputIterator result) const;
    
    //------------------------------------------------------------------------//
    // interpolate:
    //    Given a field of "things" (probably Memento's
    //    or std::pair<double,double>'s), return interpolations to the
    //    ret_vals field.
    //------------------------------------------------------------------------//

    template<class FTIN, class FTOUT>
    inline void interpolate(const FTIN &args, FTOUT &ret_vals) const;

    //------------------------------------------------------------------------//
    // interpolate:
    //    Given begin and end iterators of something (probably double's)
    //    and a begin iterator of something else (also probably double;s),
    //    return interpolations to the result iterator.
    //------------------------------------------------------------------------//

    template<class InputIterator1, class InputIterator2, class OutputIterator>
    inline void interpolate(InputIterator1 first1,
			    InputIterator1 last1,
			    InputIterator2 first2,
			    OutputIterator result) const;
    
    //------------------------------------------------------------------------//
    // interpolate:
    //    Given a field of something (probably double's)
    //    and another field of something else (also probably double;s),
    //    return interpolations to the ret_vals field.
    //------------------------------------------------------------------------//

    template<class FTIN1, class FTIN2, class FTOUT>
    inline void interpolate(const FTIN1 &arg1s,
			    const FTIN2 &arg2s,
			    FTOUT &ret_vals) const;

  private:
    
    //=======================================================================//
    // IMPLEMENTATION
    //=======================================================================//

};

//=======================================================================//
// INLINE METHODS
//=======================================================================//

//------------------------------------------------------------------------//
// BilinearInterpTable:
//    Constructor supplying the two axes grids and a binary functional
//    to perform the calculation.
//------------------------------------------------------------------------//

template<class BinaryOperation>
inline
BilinearInterpTable::BilinearInterpTable(const std::vector<double> &x1vals_,
					 const std::vector<double> &x2vals_,
					 BinaryOperation binary_op)
    : grid(new BilinearInterpGrid(x1vals_, x2vals_)),
      yvals(x1vals_.size(),x2vals_.size()), isEmpty(false)
{
    Require(yvals.size() == grid->size());

    setTable(binary_op);
}
    
//------------------------------------------------------------------------//
// BilinearInterpTable:
//    Constructor supplying the grid and a binary functional
//    to perform the calculation.
//------------------------------------------------------------------------//

template<class BinaryOperation>
inline
BilinearInterpTable::BilinearInterpTable(const dsxx::SP<BilinearInterpGrid> &grid_,
					 BinaryOperation binary_op)
    : grid(grid_), yvals(grid_->size(1),grid_->size(2)),
      isEmpty(false)
{
    Require(yvals.size() == grid->size());

    setTable(binary_op);
}
    
//------------------------------------------------------------------------//
// interpolate:
//    Given a pair of doubles, return an interpolated value.
//------------------------------------------------------------------------//

inline double BilinearInterpTable::interpolate(std::pair<double,double> p) const
{
    return interpolate(grid->getMemento(p.first, p.second));
}
    
//------------------------------------------------------------------------//
// interpolate:
//    Given a pair of doubles, return an interpolated value.
//------------------------------------------------------------------------//

inline double BilinearInterpTable::interpolate(double x1, double x2) const
{
    return interpolate(grid->getMemento(x1, x2));
}

//------------------------------------------------------------------------//
// interpolate:
//    Given a Memento, return an interpolated value.
//------------------------------------------------------------------------//

inline double BilinearInterpTable::interpolate(const Memento &memento) const
{
    int j, k;
    grid->getIndices(memento, j, k);

    double t, u;
    grid->getCoefficients(memento, t, u);
    
    return (1-t)*(1-u) * yvals(j  ,k)
	+  t*(1-u)     * yvals(j+1,k)
	+  t*u         * yvals(j+1,k+1)
	+  (1-t)*u     * yvals(j  ,k+1);
}

//===========================================================================//
// class BilinearInterpTable::UnaryInterpolator
//     A nested class used with the std::transform to repeatedly call
//     the interpolation over a container of "things" (probably Memento's
//     or std::pair<double,double>'s).
//===========================================================================//

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
	
//------------------------------------------------------------------------//
// interpolate:
//    Given begin and end iterators of "things" (probably Memento's
//    or std::pair<double,double>'s, return interpolations to the
//    result iterator.
//------------------------------------------------------------------------//

template<class InputIterator, class OutputIterator>
inline void BilinearInterpTable::interpolate(InputIterator first,
					     InputIterator last,
					     OutputIterator result) const
{
    typedef typename std::iterator_traits<InputIterator>::value_type value_type;

    // This class turns the interploation table in to a unary functor.

    std::transform(first, last, result, UnaryInterpolator<value_type>(*this));
}

//------------------------------------------------------------------------//
// interpolate:
//    Given a field of "things" (probably Memento's
//    or std::pair<double,double>'s), return interpolations to the
//    ret_vals field.
//------------------------------------------------------------------------//

template<class FTIN, class FTOUT>
inline void BilinearInterpTable::interpolate(const FTIN &args,
					     FTOUT &ret_vals) const
{
    Require(args.size() == ret_vals.size());
    interpolate(args.begin(), args.end(), ret_vals.begin());
}

//===========================================================================//
// class BilinearInterpTable::BinaryInterpolator
//     A nested class used with the std::transform to repeatedly call
//     the interpolation over two containers (probably both of double's).
//===========================================================================//

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
	
//------------------------------------------------------------------------//
// interpolate:
//    Given begin and end iterators of something (probably double's)
//    and a begin iterator of something else (also probably double;s),
//    return interpolations to the result iterator.
//------------------------------------------------------------------------//

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

//------------------------------------------------------------------------//
// interpolate:
//    Given a field of something (probably double's)
//    and another field of something else (also probably double;s),
//    return interpolations to the ret_vals field.
//------------------------------------------------------------------------//

template<class FTIN1, class FTIN2, class FTOUT>
inline void BilinearInterpTable::interpolate(const FTIN1 &arg1s,
					     const FTIN2 &arg2s,
					     FTOUT &ret_vals) const
{
    Require(arg1s.size() == arg2s.size());
    Require(arg1s.size() == ret_vals.size());
    interpolate(arg1s.begin(), arg1s.end(), arg2s.begin(), ret_vals.begin());
}

} // end of rtt_matprops namespace

#endif                          // __matprops_BilinearInterpTable_hh__

//---------------------------------------------------------------------------//
//                              end of matprops/BilinearInterpTable.hh
//---------------------------------------------------------------------------//
