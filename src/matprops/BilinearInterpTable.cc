//----------------------------------*-C++-*----------------------------------//
// BilinearInterpTable.cc
// Randy M. Roberts
// Tue Apr  7 12:59:40 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "matprops/BilinearInterpTable.hh"
#include <functional>
#include <algorithm>
#include <numeric>

// #include <iostream>
// using std::cerr;
// using std::endl;

using XTM::BilinearInterpTable;
using std::lower_bound;

typedef BilinearInterpTable::Memento Memento;

Memento BilinearInterpTable::getMemento(double x1, double x2) const
{
    // cerr << endl << "getMemento: x1=" << x1 << " x2=" << x2 << endl;
    
    // find the index into the x1 axis before x1
    
    int j;
    if (x1vals[0] == x1)
	j = 0;
    else
	j = lower_bound(x1vals.begin(), x1vals.end(), x1) - x1vals.begin() - 1;
    
    // find the index into the x2 axis before x2
    
    int k;
    if (x2vals[0] == x2)
	k = 0;
    else
	k = lower_bound(x2vals.begin(), x2vals.end(), x2) - x2vals.begin() - 1;

    // cerr << "j=" << j << " k=" << k << endl;
    // cerr << "x1vals[j]=" << x1vals[j] << " x2vals[k]=" << x2vals[k] << endl;
    // cerr << "x1vals[j+1]=" << x1vals[j+1]
    //      << " x2vals[k+1]=" << x2vals[k+1] << endl;

    // Assert that we are not doing an extrapolation

    Assert(j >= 0);
    Assert(j < x1vals.size()-1);
    Assert(x1vals[j] <= x1);
    Assert(x1vals[j+1] >= x1);
    
    Assert(k >= 0);
    Assert(k < x2vals.size()-1);
    Assert(x2vals[k] <= x2);
    Assert(x2vals[k+1] >= x2);
    
    double t = (x1 - x1vals[j]) / (x1vals[j+1] - x1vals[j]);
    double u = (x2 - x2vals[k]) / (x2vals[k+1] - x2vals[k]);

    // These extra assertions also check that the x1vals and x2vals
    // vectors are ordered.
    
    Assert(t >= 0.0 && t <= 1.0);
    Assert(u >= 0.0 && u <= 1.0);
    
    return Memento(j, k, t, u);    
}

bool BilinearInterpTable::axesAreOrdered() const
{
    using std::inner_product;
    using std::logical_and;
    using std::less;
    
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
//                              end of BilinearInterpTable.cc
//---------------------------------------------------------------------------//
