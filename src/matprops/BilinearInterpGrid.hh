//----------------------------------*-C++-*----------------------------------//
// BilinearInterpGrid.hh
// Randy M. Roberts
// Tue Apr 14 15:28:53 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __matprops_BilinearInterpGrid_hh__
#define __matprops_BilinearInterpGrid_hh__

#include "ds++/Assert.hh"
#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>

#ifndef BEGIN_NS_XTM
#define BEGIN_NS_XTM namespace XTM  {
#define END_NS_XTM }
#endif

BEGIN_NS_XTM

// Forward Reference

class BilinearInterpTable;

//===========================================================================//
// class BilinearInterpGrid - 
//
// Date created :
// Purpose      :
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class BilinearInterpGrid
{

    friend class BilinearInterpTable;

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


    // DATA

  private:

    bool errorOnOutOfBounds;
    
    // x1vals is the x1 axis of the 2-dimensional grid of tabulated y values.

    std::vector<double> x1vals;

    // x2vals is the x2 axis of the 2-dimensional grid of tabulated y values.

    std::vector<double> x2vals;

  public:

    // CREATORS
    
    //------------------------------------------------------------------------//
    // BilinearInterpGrid:
    //    Creates an empty grid.
    //------------------------------------------------------------------------//

    BilinearInterpGrid(bool errorOnOutOfBounds_ = true)
	: errorOnOutOfBounds(errorOnOutOfBounds_)
    {
	// *** empty ***
    }

    //------------------------------------------------------------------------//
    // BilinearInterpGrid:
    //    Constructor supplying the two axes grids and the two-dimensional
    //    table of evaluations.
    //------------------------------------------------------------------------//

    BilinearInterpGrid(const std::vector<double> &x1vals_,
		       const std::vector<double> &x2vals_,
		       bool errorOnOutOfBounds_ = true)
	: x1vals(x1vals_), x2vals(x2vals_),
	  errorOnOutOfBounds(errorOnOutOfBounds_)
    {
	Require(x1vals.size() >= 2);
	Require(x2vals.size() >= 2);
	Require(axesAreOrdered());
    }
    
    //=======================================================================//
    // MANIPULATORS
    //=======================================================================//

    // *** none ***
    
    //=======================================================================//
    // ACCESSORS
    //=======================================================================//

    //------------------------------------------------------------------------//
    // getGrid:
    //   Return the grid to the user.
    //------------------------------------------------------------------------//

    const std::vector<double> &getGrid(int i) const
    {
	Require(i == 1 || i==2);
	return (i==1) ? x1vals : x2vals;
    }
    
    //------------------------------------------------------------------------//
    // getMemento:
    //    Accessor to obtain the Memento from an x1, x2 interpolation.
    //    This allows the user to avoid the repeated calculation of
    //    interpolation quantities.
    //------------------------------------------------------------------------//

    Memento getMemento(double x1, double x2) const;
    

    //------------------------------------------------------------------------//
    // size:
    //    Return the size of the grid's axes.
    //    With size(1) return size of x1 axis;
    //    with size(2) return size of x2 axis.
    //------------------------------------------------------------------------//

    int size(int dim) const
    {
	return ( dim == 1 ? x1vals.size() : ( dim == 2 ? x2vals.size() : 0 ) );
    }
    
    //------------------------------------------------------------------------//
    // size:
    //    Return the total size of the grid's axes.
    //    With size() return size(1) * size(2) (See above).
    //------------------------------------------------------------------------//

    int size() const { return size(1)*size(2); }

    //------------------------------------------------------------------------//
    // x1:
    //    Return the j'th value on the x1 axis.
    //------------------------------------------------------------------------//

    double x1(int j) const
    {
	Require(j >= 0 && j < size(1));
	return x1vals[j];
    }
    
    //------------------------------------------------------------------------//
    // x2:
    //    Return the k'th value on the x2 axis.
    //------------------------------------------------------------------------//

    double x2(int k) const
    {
	Require(k >= 0 && k < size(2));
	return x2vals[k];
    }

  private:
    
    //=======================================================================//
    // IMPLEMENTATION
    //=======================================================================//

    //------------------------------------------------------------------------//
    // getIndices:
    //    Return the j and k indices into the x1 and x2 axes from
    //    the memento.
    //------------------------------------------------------------------------//

    inline void getIndices(const Memento &memento, int &j, int &k) const;

    //------------------------------------------------------------------------//
    // getCoefficients:
    //    Return the bilinear interpolation coefficients from the memento.
    //------------------------------------------------------------------------//

    inline void getCoefficients(const Memento &memento, double &t,
				double &u) const;
    
    //------------------------------------------------------------------------//
    // axesAreOrdered:
    //    Return true if both axes are monotomically increasingly ordered.
    //------------------------------------------------------------------------//

    bool axesAreOrdered() const;
};

//===========================================================================//
// class BilinearInterpGrid::Memento
//     This class is part of the Memento pattern
//     from "Design Patterns", E. Gamma, et. al., 1995
//
//     This class is meant to be "the internal state" of the
//     BilinearInterpGrid (to be held externally),
//     and therefore, must be seen only by that class.
//     I am using it to contain the internal state for a given
//     (x1, x2) pair's interpolation.
//
//     For faster access to the interpolation table, a user can create
//     a field of Memento to store the intermediate interpolation table
//     calculations.
//===========================================================================//

class BilinearInterpGrid::Memento
{
    friend class BilinearInterpGrid;

    friend std::ostream &operator<<(std::ostream &os, const Memento &rhs);
    
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

    // all Creators are private, so only usable by BilinearInterpGrid

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

//------------------------------------------------------------------------//
// operator<<
//    Use this to emit debug output.
//------------------------------------------------------------------------//

inline std::ostream &operator<<(std::ostream &os,
				const BilinearInterpGrid::Memento &rhs)
{
    os << "(" << rhs.j << ","
       << rhs.t << ","
       << rhs.k << ","
       << rhs.u << ")";
    return os;
}

//------------------------------------------------------------------------//
// getIndices:
//    Return the j and k indices into the x1 and x2 axes from
//    the memento.
//------------------------------------------------------------------------//

inline void BilinearInterpGrid::getIndices(const Memento &memento,
					   int &j, int &k) const
{
    j = memento.j;
    k = memento.k;

    Assert(j >= 0 && j < size(1));
    Assert(k >= 0 && k < size(2));
}

//------------------------------------------------------------------------//
// getCoefficients:
//    Return the bilinear interpolation coefficients from the memento.
//------------------------------------------------------------------------//

inline void BilinearInterpGrid::getCoefficients(const Memento &memento,
						double &t, double &u) const
{
    t = memento.t;
    u = memento.u;

    Assert(t >= 0.0 && t <= 1.0);
    Assert(u >= 0.0 && u <= 1.0);
}

END_NS_XTM  // namespace XTM

#endif                          // __matprops_BilinearInterpGrid_hh__

//---------------------------------------------------------------------------//
//                              end of matprops/BilinearInterpGrid.hh
//---------------------------------------------------------------------------//
