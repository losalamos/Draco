//----------------------------------*-C++-*----------------------------------//
// Copyright 1996-2006 The Regents of the University of California.
// Copyright 2006-2010 LANS, LLC
// All rights reserved.
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// Created on: Thu Sep 12 10:56:25 1996
// Created by: Geoffrey Furnish
// Also maintained by:
//
//---------------------------------------------------------------------------//

#include "UserVec.hh"

#include <iostream>
#include <cmath>
#include <string>

using std::cout;
using std::endl;
using std::exp;

int sz = 5;
bool passed = true;

//---------------------------------------------------------------------------//
// Function to print out the contents of a UserVec.
//---------------------------------------------------------------------------//

template<class T> void dump( UserVec<T>& u )
{
    for( int i=0; i < u.size(); i++ )
	cout << "u[" << i << "]=" << u[i] << endl;
}

//---------------------------------------------------------------------------//
// Check that we can handle UserVec the usual way.
//---------------------------------------------------------------------------//

void t1()
{
    cout << "t1: start\n";

    UserVec<float> u(sz), v(sz);

    for( int i=0; i < u.size(); i++ )
	u[i] = 2.*i;
    passed &= (u[1] == 2.);

    u = 1.;
    passed &= (u[0] == 1.);

    v = u;
    passed &= (v[0] == u[0]);

    cout << "t1: end\n";
}

//---------------------------------------------------------------------------//
// Check the simple binary operations with assignments.
//---------------------------------------------------------------------------//

void t2()
{
    cout << "t2: start\n";

    UserVec<float> u(sz), v(sz), w(sz);

    float a = 2.0;

    for( int i=0; i < v.size(); i++ ) {
	v[i] = 2.*i;
	w[i] = 10.-i;
    }

    u = 1.;
    u += a;
    passed &= (u[0] == 1. + a);

    u -= a;
    passed &= (u[0] == 1.);

    u *= a;
    passed &= (u[0] == a);

    u /= a;
    passed &= (u[0] == 1.);


    u += v;
    passed &= (u[2] == 5.);

    u -= v;
    passed &= (u[2] == 1.);

    u *= v;
    passed &= (u[2] == 4.);

    u /= v;
    passed &= (u[2] == 1.);

    u = v;
    passed &= (u[2] == 4.);


    u = v + w;
    passed &= (u[1] == 11.);

    u = v - w;
    passed &= (u[1] == -7.);

    u = v * w;
    passed &= (u[1] == 18.);

    u = w / v;
    passed &= (u[2] == 2.);


    u = 1.;
    u += v + w;
    passed &= (u[1] == 12.);

    u -= v + w;
    passed &= (u[1] == 1.);

    u *= v + w;
    passed &= (u[1] == 11.);

    u /= v + w;
    passed &= (u[1] == 1.);


    cout << "t2: end\n";
}

//---------------------------------------------------------------------------//
// Now a more involved test of binary operations.
//---------------------------------------------------------------------------//

void t3()
{
    cout << "t3: start\n";

    UserVec<float> a(sz), b(sz), c(sz), d(sz), e(sz), f(sz);

    a = 4.;
    b = 2.;
    c = 2.;
    d = 5.;
    e = 2.;

    f = (a+b)*c/(d-e);
    passed &= (f[0] == 4.);

    //dump(f);

    cout << "t3: end\n";
}

//---------------------------------------------------------------------------//
// Check the simple unary operations with assignments.
//---------------------------------------------------------------------------//

void t4()
{
    cout << "t4: start\n";

    UserVec<float> u(sz), v(sz);

    for( int i=0; i < v.size(); i++ ) {
	v[i] = 2.*i;
    }

    u = +v;
    passed &= (u[1] == v[1]);

    u = -v;
    passed &= (u[1] == -v[1]);

    cout << "t4: end\n";
}

//---------------------------------------------------------------------------//
// Check the other binary operations with assignments.
//---------------------------------------------------------------------------//

void t5()
{
    cout << "t5: start\n";

    UserVec<float> a(sz), b(sz), c(sz), d(sz);

    a = 4.;
    b = 3.;

    c = pow( a, 3.f );
    passed &= (c[0] == 64.);

    c = 2.;
    d = pow( c, 3 );
    passed &= (d[0] == 8.);

    d = 1.f + pow( c, 3 );
    passed &= (d[0] == 9.);

    c = pow(a,b);
    passed &= (c[0] == 64.);

    a = 0.;
    b = 1.;
    c = atan2(a,b);
    passed &= (c[0] == 0.);
    a = 4.;
    b = 3.;

    c = min(a,b);
    passed &= (c[0] == 3.);

    c = max(a,b);
    passed &= (c[0] == 4.);

    b = 7.;
    c = fmod(b,a);
    passed &= (c[0] == 3.);
    b = 3.;

    d = pow(a,min(a,b));
    passed &= (d[0] == 64.);

    cout << "t5: end\n";
}

//---------------------------------------------------------------------------//
// Check the other unary operations with assignments.
//---------------------------------------------------------------------------//

void t6()
{
    cout << "t6: start\n";

    UserVec<float> a(sz), b(sz);
    UserVec<int> ai(sz), bi(sz);
    UserVec<long> al(sz), bl(sz);

    a = 0.;

    b = sin(a);
    passed &= (b[0] == 0.);

    b = cos(a);
    passed &= (b[0] == 1.);

    b = tan(a);
    passed &= (b[0] == 0.);

    b = asin(a);
    passed &= (b[0] == 0.);

    a = 1.;
    b = acos(a);
    passed &= (b[0] == 0.);
    a = 0.;

    b = atan(a);
    passed &= (b[0] == 0.);

    b = sinh(a);
    passed &= (b[0] == 0.);

    b = cosh(a);
    passed &= (b[0] == 1.);

    b = tanh(a);
    passed &= (b[0] == 0.);

    b = exp(a);
    passed &= (b[0] == 1.);

    a = exp(1.);
    b = log(a);
    passed &= (std::fabs(b[0] - 1.) < 0.00001);

    a = 10.;
    b = log10(a);
    passed &= (b[0] == 1.);

    a = 9.;
    b = sqrt(a);
    passed &= (b[0] == 3.);

    a = 3.4;
    b = ceil(a);
    passed &= (b[0] == 4.);

    ai = -3;
    bi = abs(ai);
    passed &= (bi[0] == 3);

    al = -3;
    bl = labs(al);
    passed &= (bl[0] == 3);

    a = -3.4;
    b = fabs(a);
    passed &= (std::fabs(b[0] - 3.4) < 0.00001);

    a = 3.4;
    b = floor(a);
    passed &= (b[0] == 3.);

    cout << "t6: end\n";
}

//---------------------------------------------------------------------------//
// Check incompatible participation.
//---------------------------------------------------------------------------//

void t7()
{
    cout << "t7: start\n";

    UserVec<float> a(sz), b(sz);
    FooBar<float>  c(sz);

    a = -99.;

    b = 2.; c = 3.;

// This statement should be illegal!
//     a = b + c;

    cout << "t7: end\n";
}

//---------------------------------------------------------------------------//
// Check operator+= with various arguments.
//---------------------------------------------------------------------------//

void t8()
{
    cout << "t8: start\n";

    UserVec<float> u(sz), v(sz), w(sz);

    float a = 2.0;

    for( int i=0; i < v.size(); i++ ) {
	v[i] = 2.*i;
	w[i] = 5.-i;
    }
    u = 0;

    cout << "-----------------------------------" << endl;
    cout << "Output of u"                         << endl;
    cout << "-----------------------------------" << endl;
    dump(u);

    cout << "-----------------------------------" << endl;
    cout << "Output of v"                         << endl;
    cout << "-----------------------------------" << endl;
    dump(v);

    cout << "-----------------------------------" << endl;
    cout << "Output of w"                         << endl;
    cout << "-----------------------------------" << endl;
    dump(w);

    cout << "-----------------------------------" << endl;
    cout << "Output of u += 2.0"                  << endl;
    cout << "-----------------------------------" << endl;
    u += a;
    dump(u);

    cout << "-----------------------------------" << endl;
    cout << "Output of u += v"                    << endl;
    cout << "-----------------------------------" << endl;
    u += v;
    dump(u);

    cout << "-----------------------------------" << endl;
    cout << "Output of u += v + w"                << endl;
    cout << "-----------------------------------" << endl;
    u += v + w;
    dump(u);

    cout << "t8: end\n";
}

void version(const std::string &progname)
{
    std::string version = "1.0.0";
    cout << progname << ": version " << version << endl;
}

//---------------------------------------------------------------------------//
// Main program, just run through each test in turn.
//---------------------------------------------------------------------------//

int main( int argc, char *argv[] )
{
    version(argv[0]);
    for (int arg=1; arg < argc; arg++)
	if (std::string(argv[arg]) == "--version")
	    return 0;
    
    t1();
    t2();
    t3();
    t4();
    t5();
    t6();
    t7();
    //t8();

    // Print the status of the test.
    cout <<     "\n***********************************************";
    if (passed) 
        cout << "\n**** Expression Template Self Test: PASSED ****";
    else
        cout << "\n**** Expression Template Self Test: FAILED ****";
    cout <<     "\n***********************************************\n" << endl;
    return 0;
}

//---------------------------------------------------------------------------//
//                              end of tstUV.cc
//---------------------------------------------------------------------------//
