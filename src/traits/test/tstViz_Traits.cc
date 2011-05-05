//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   traits/test/tstViz_Traits.cc
 * \author Thomas M. Evans
 * \date   Fri Jan 21 17:51:52 2000
 * \brief  Viz_Traits test.
 * \note   Copyright (C) 2000-2010 Los Alamos National Security, LLC.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Traits_Test.hh"
#include "ds++/Release.hh"
#include "../Viz_Traits.hh"

#include <iostream>
#include <string>
#include <vector>

using namespace std;

using rtt_traits::Viz_Traits;

// passing condition
bool passed = true;
#define ITFAILS passed = rtt_traits_test::fail(__LINE__);

//---------------------------------------------------------------------------//
// simple test field class for checking viz traits

template<class T>
class Test_Field
{
  public:
    typedef T value_type;

  private:
    vector<vector<T> > data;

  public:
    Test_Field(const vector<vector<T> > &data_in) : data(data_in) {}

    T operator()(size_t i, size_t j) const { return data[i][j]; }
    size_t nrows() const { return data.size(); }
    size_t ncols(size_t r) const { return data[r].size(); }
};

//---------------------------------------------------------------------------//
// test vector traits specialization

template<class VVF>
void test_vector()
{
    VVF field(3);
    for (size_t i = 0; i < field.size(); i++)
    {
	field[i].resize(i+2);
	for (size_t j = 0; j < field[i].size(); j++)
	    field[i][j] = 2 * i + 4 * j;
    }

    Viz_Traits<VVF> vdf(field);

    if (vdf.nrows() != field.size())         ITFAILS; 
    for (size_t i = 0; i < vdf.nrows(); i++)
    {
	if (vdf.ncols(i) != field[i].size()) ITFAILS;
	for (size_t j = 0; j < vdf.ncols(i); j++)
	{
	    if (static_cast<int>(vdf(i, j)) != field[i][j])       ITFAILS;
	    if (vdf(i, j) != 2*i + 4*j)                           ITFAILS;
	}
    }				   
}

//---------------------------------------------------------------------------//
// standard Viz_Traits field test

template<class T>
void test_FT()
{
    vector<vector<T> > field(3);
    for (size_t i = 0; i < field.size(); i++)
    {
	field[i].resize(i+2);
	for (size_t j = 0; j < field[i].size(); j++)
	    field[i][j] = 2 * i + 4 * j;
    }

    Test_Field<T> test_field(field);
    
    Viz_Traits<Test_Field<T> > vt(test_field);

    if (vt.nrows() != 3)                                         ITFAILS;
    for (size_t i = 0; i < vt.nrows(); i++)
    {
	if (vt.ncols(i) != field[i].size())                          ITFAILS;
	for (size_t j = 0; j < vt.ncols(i); j++)
	    if (vt(i, j) != field[i][j])                             ITFAILS;
    }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    // version tag
    cout << argv[0] << ": version " << rtt_dsxx::release() << endl; 
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	    return 0;
    
    // tests
    test_vector<vector<vector<int> > >   ();
    test_vector<vector<vector<double> > >();
    test_vector<vector<vector<float> > > ();

    test_FT<int>   ();
    test_FT<double>();

    // status of test
    cout <<     "\n*************************************";
    if (passed) 
        cout << "\n****Viz_Traits Self Test: PASSED ****";
    cout <<     "\n*************************************\n\n"
         << "Done testing Viz_Traits." << endl;
}

//---------------------------------------------------------------------------//
//                              end of tstViz_Traits.cc
//---------------------------------------------------------------------------//
