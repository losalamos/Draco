//----------------------------------*-C++-*----------------------------------//
// BilinearInterpGrid.cc
// Randy M. Roberts
// Tue Apr 14 15:28:53 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "matprops/BilinearInterpGrid.hh"

#include <algorithm>
#include <numeric>

using rtt_matprops::BilinearInterpGrid;
using std::lower_bound;

typedef BilinearInterpGrid::Memento Memento;

//------------------------------------------------------------------------//
// getMemento:
//    Accessor to obtain the Memento from an x1, x2 interpolation.
//    This allows the user to avoid the repeated calculation of
//    interpolation quantities.
//------------------------------------------------------------------------//

Memento BilinearInterpGrid::getMemento(double x1, double x2) const
{
    Require(x1vals.size() >= 2);
    Require(x2vals.size() >= 2);
    Require(axesAreOrdered());

    // cerr << endl << "getMemento: x1=" << x1 << " x2=" << x2 << endl;
    
    // find the index into the x1 axis before x1

    // Because of the STL definition of lower_bound, we must manually
    // check for the equality with the first element to ensure
    // consistent indexing to element smaller or equal to the value.
    
    int j = 0;
    if (x1vals[0] != x1)
	j = lower_bound(x1vals.begin(), x1vals.end(), x1) - x1vals.begin() - 1;
    
    // find the index into the x2 axis before x2

    // Because of the STL definition of lower_bound, we must manually
    // check for the equality with the first element to ensure
    // consistent indexing to element smaller or equal to the value.
    
    int k = 0;
    if (x2vals[0] != x2)
	k = lower_bound(x2vals.begin(), x2vals.end(), x2) - x2vals.begin() - 1;

    double t;
    double u;

    // With errorOnOutOfBounds == false we will flatten
    // the table as if there were duplicate first
    // and last y1 and y2 values at +/- infinity.
	
    if (!errorOnOutOfBounds && j < 0)
    {
	j = 0;
	t = 0.0;
    }
    else if (!errorOnOutOfBounds && (j >= x1vals.size()-1))
    {
	j = x1vals.size()-2;
	t = 1.0;
    }
    else
    {
	// Assert that we are not doing an extrapolation.

	Assert(j >= 0);
	Assert(j < x1vals.size()-1);
	Assert(x1vals[j] <= x1);
	Assert(x1vals[j+1] >= x1);

	t = (x1 - x1vals[j]) / (x1vals[j+1] - x1vals[j]);
    }
	
    if (!errorOnOutOfBounds && k < 0)
    {
	k = 0;
	u = 0.0;
    }
    else if (!errorOnOutOfBounds && (k >= x2vals.size()-1))
    {
	k = x2vals.size()-2;
	u = 1.0;
    }
    else
    {
	// Assert that we are not doing an extrapolation.

	Assert(k >= 0);
	Assert(k < x2vals.size()-1);
	Assert(x2vals[k] <= x2);
	Assert(x2vals[k+1] >= x2);

	u = (x2 - x2vals[k]) / (x2vals[k+1] - x2vals[k]);
    }
    
    // These extra assertions also check that the x1vals and x2vals
    // vectors are ordered.
    
    Assert(t >= 0.0 && t <= 1.0);
    Assert(u >= 0.0 && u <= 1.0);
    
    return Memento(j, k, t, u);    
}

//------------------------------------------------------------------------//
// axesAreOrdered:
//    Return true if both axes are monotomically increasingly ordered.
//------------------------------------------------------------------------//

bool BilinearInterpGrid::axesAreOrdered() const
{
    Require(x1vals.size() >= 2);
    Require(x2vals.size() >= 2);

    using std::inner_product;
    using std::logical_and;
    using std::less;

    // I know....
    // This is a bit obfusticated, but it is rather simple to explain.
    // The std::inner_product has been co-opted to perform the following:
    //    compare, using "less than", each element with the one at
    //    one index greater, and logically "and" all of those comparisons.
    // The results is true if the vector is monotomically increasing.
    
    bool ordered = inner_product(x1vals.begin(), x1vals.end()-1,
				 x1vals.begin()+1,
				 true,
				 logical_and<bool>(),
				 less<double>());

    ordered = ordered && inner_product(x2vals.begin(), x2vals.end()-1,
				       x2vals.begin()+1,
				       true,
				       logical_and<bool>(),
				       less<double>());
    
    return ordered;
}

//---------------------------------------------------------------------------//
//                              end of BilinearInterpGrid.cc
//---------------------------------------------------------------------------//
