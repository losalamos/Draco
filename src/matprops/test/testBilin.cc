#include "matprops/BilinearInterpTable.hh"
#include "ds++/Assert.hh"
#include "ds++/Mat.hh"

#include <iostream>
#include <vector>
#include <list>

#ifndef BEGIN_NS_XTM
#define BEGIN_NS_XTM namespace XTM  {
#define END_NS_XTM }
#endif

using namespace XTM;

using dsxx::Mat2;
using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::back_inserter;
using std::list;
using std::pair;

template<class C>
void showSequence(const char *str, const C &c,
		  std::ostream &os = cout, const char *sep=" ")
{
    os << str;
    C::const_iterator iter = c.begin();
    while (iter != c.end())
	os << *iter++ << sep;
    os << endl;
}

inline double func(double x1, double x2)
{
#if 0
    const double pi = 3.1415;
    return std::sin(x1*x2*pi/100.0);
#else
    // The following function is bi-linear, and should
    // be exactly interpolated.
    
    return  10.2*x1*x2 + 3.4*x1 + 9.2*x2 + 4.0;
#endif
}

BilinearInterpTable getTable()
{
    const double x1[] = {0., 2., 4., 6., 8., 9.};
    const double x2[] = {1., 3., 5., 7., 9., 10., 11.};

    const int sz1 = sizeof(x1)/sizeof(double);
    const int sz2 = sizeof(x2)/sizeof(double);
    
    // cout << sizeof(x1) << " " << sizeof(x1)/sizeof(double) << endl;
    
    vector<double> vx1(x1, x1+sz1);
    vector<double> vx2(x2, x2+sz2);

    BilinearInterpGrid grid(vx1, vx2);

#if 0
    
    Mat2<double> mat(sz1, sz2);

    for (int i=0; i<sz1; i++)
    {
	for (int j=0; j<sz2; j++)
	{
	    mat(i, j) = func(vx1[i], vx2[j]);
	}
    }

    return BilinearInterpTable(grid, mat);

#else

    dsxx::SP<BilinearInterpGrid> spgrid(new BilinearInterpGrid(grid));

    return BilinearInterpTable(spgrid, func);
    
#endif
}

void testBilinDoit()
{
    BilinearInterpTable biLinTable = getTable();
    
    const double testvals1[] = {1.1, 5.5, 3.3, 6.6, 1.0, 9.0, 9.0};
    const double testvals2[] = {1.7, 4.0, 5.5, 3.3, 6.6, 1.0, 9.0};

    const int sztest = sizeof(testvals1)/sizeof(double);
    Assert(sizeof(testvals2)/sizeof(double) == sztest);

    vector<BilinearInterpTable::Memento> mementos;

    // Load up the mementos

    const dsxx::SP<BilinearInterpGrid> grid = biLinTable.getGrid();

    for (int i=0; i<sztest; i++)
	mementos.push_back(grid->getMemento(testvals1[i], testvals2[i]));
    
    vector<double> results;
    
    biLinTable.interpolate(testvals1, testvals1+sztest, testvals2,
			   back_inserter(results));

    showSequence(" results: ", results, cerr);

    // Uncomment the next line if you want to trigger an assertion failure.
    // results.push_back(999.9);
    
    for (int i=0; i<results.size(); i++)
	results[i] = 0.0;

    biLinTable.interpolate(mementos, results);

    showSequence(" results: ", results, cerr);

    list<pair<double,double> > testlist;
    for (int i=0; i<sztest; i++)
	testlist.push_back(pair<double,double>(testvals1[i], testvals2[i]));

    for (int i=0; i<results.size(); i++)
	results[i] = 0.0;

    biLinTable.interpolate(testlist, results);

    showSequence(" results: ", results, cerr);

    vector<double> expected;
    for (int i=0; i<sztest; i++)
	expected.push_back(func(testvals1[i], testvals2[i]));

    showSequence("expected: ", expected, cerr);
}
