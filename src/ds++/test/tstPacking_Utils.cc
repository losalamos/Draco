//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstPacking_Utils.cc
 * \author Thomas M. Evans
 * \date   Fri Jul 20 17:22:36 2001
 * \brief  Packing Utils test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "ds++_Test.hh"
#include "../Release.hh"
#include "../Packing_Utils.hh"

#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <cstdlib>

using namespace std;

using rtt_dsxx::Packer;
using rtt_dsxx::Unpacker;

// passing condition
bool passed = true;
#define ITFAILS passed = rtt_dsxx_test::fail(__LINE__);

//---------------------------------------------------------------------------//
 
void packing_test()
{
    // make some data
    double x = 102.45;
    double y = 203.89;
    double z = 203.88;

    int ix   = 10;
    int iy   = 11; 
    int iz   = 12;

    // make 2 buffers for data
    int s1   = 2 * sizeof(double) + 2 * sizeof(int);
    char *b1 = new char[s1];
    int s2   = sizeof(double) + sizeof(int);
    char *b2 = new char[s2];

    // pack the data
    {
	Packer p;

	p.set_buffer(s1, b1);
	p << x << ix;
	p.pack(y);
	p.pack(iy);

	if (p.get_ptr() != b1+s1) ITFAILS;

	// try catching a failure
	bool caught = false;
	try
	{
	    p << iz;
	}
	catch (const rtt_dsxx::assertion &a)
	{
	    cout << "Should catch this: " << a.what() << endl;
	    caught = true;
	}
	if (!caught) ITFAILS;

	p.set_buffer(s2, b2);
	p << iz << z;
	
	if (p.get_ptr() != b2+s2) ITFAILS;
    }
    
    // unpack the data
    {
	Unpacker u;

	double   d = 0;
	int      i = 0;

	u.set_buffer(s1, b1);
	u >> d >> i;
	if (d != 102.45)          ITFAILS;
	if (i != 10)              ITFAILS;

	u.unpack(d);
	u.unpack(i); 
	if (d != 203.89)          ITFAILS;
	if (i != 11)              ITFAILS;

	if (u.get_ptr() != s1+b1) ITFAILS;

	// try catching a failure
	bool caught = false;
	try
	{
	     u.unpack(i);
	}
	catch (const rtt_dsxx::assertion &a)
	{
	    cout << "Should catch this: " << a.what() << endl;
	    caught = true;
	}
	if (!caught) ITFAILS;

	u.set_buffer(s2, b2);
	u >> i >> d;
	if (i != 12)              ITFAILS;
	if (d != 203.88)          ITFAILS;
	
	if (u.get_ptr() != s2+b2) ITFAILS;
    }

    delete [] b1;
    delete [] b2;

    // try packing a vector and char array
    double r = 0;
    srand(125);

    vector<double> vx(100, 0.0);
    vector<double> ref(100, 0.0);

    char c[4] = {'c','h','a','r'};

    for (int i = 0; i < vx.size(); i++)
    {
	r      = rand();
	vx[i]  = r / RAND_MAX;
	ref[i] = vx[i];
    }

    int size     = 100 * sizeof(double) + 4;
    char *buffer = new char[size];
    
    // pack
    {
	Packer p;
	p.set_buffer(size, buffer);

	for (int i = 0; i < vx.size(); i++)
	    p << vx[i];

	for (int i = 0; i < 4; i++)
	    p << c[i];

	if (p.get_ptr() != buffer+size) ITFAILS;
    }
    
    // unpack
    {
	char cc[4];
	vector<double> x(100, 0.0);

	Unpacker u;
	u.set_buffer(size, buffer);

	for (int i = 0; i < x.size(); i++)
	    u >> x[i];

	for (int i = 0; i < 4; i++)
	    u >> cc[i];

	if (u.get_ptr() != buffer+size) ITFAILS;
	
	for (int i = 0; i < x.size(); i++)
 	    if (x[i] != ref[i]) ITFAILS;

	if (c[0] != 'c') ITFAILS;
	if (c[1] != 'h') ITFAILS;
	if (c[2] != 'a') ITFAILS;
	if (c[3] != 'r') ITFAILS;
    }

    delete [] buffer;
}

//---------------------------------------------------------------------------//
// main

int main(int argc, char *argv[])
{
    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	{
	    cout << argv[0] << ": version " << rtt_dsxx::release() 
		 << endl; 
	    return 0;
	}

    // tests
    try
    {
	packing_test();
    }
    catch (const rtt_dsxx::assertion &ass)
    {
	cout << "Test: assertion failure at line " 
	     << ass.what() << endl;
	return 1;
    }
    catch(...)
    {
	cout << "HELP ME" << endl;
	return 1;
    }

    // status of test
    cout << endl;
    cout <<     "************************************" << endl; 
    if (passed) 
    {
        cout << "**** Packing Utils Self Test: PASSED" << endl;
    }
    cout <<     "************************************" << endl;
    cout << endl;

    cout << "Done testing Packing_Utils." << endl;
}

//---------------------------------------------------------------------------//
//                              end of tstPacking_Utils.cc
//---------------------------------------------------------------------------//
