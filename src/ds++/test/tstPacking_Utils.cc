//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstPacking_Utils.cc
 * \author Thomas M. Evans
 * \date   Wed Nov  7 15:58:08 2001
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "ds_test.hh"
#include "../Release.hh"
#include "../Packing_Utils.hh"
#include "../Assert.hh"
#include "../Soft_Equivalence.hh"

#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <sstream>
#include <typeinfo>

using namespace std;

using rtt_dsxx::Packer;
using rtt_dsxx::Unpacker;
using rtt_dsxx::pack_data;
using rtt_dsxx::unpack_data;
using rtt_dsxx::soft_equiv;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void do_some_packing(Packer &p,
		     const vector<double> &vd,
		     const vector<int> &vi)
{
    for ( int i = 0; i < vd.size(); ++i )
	p << vd[i];

    for ( int i = 0; i < vi.size(); ++i )
	p << vi[i];
}

void compute_buffer_size_test()
{
    // make data

    int num_vd = 5;
    vector<double> vd(num_vd, 2.3432);
    vd[3] = 22.4;

    int num_vi = 3;
    vector<int> vi(num_vi, 6);
    vi[0] = 7;
    vi[1] = 22;

    unsigned int total_size = num_vi * sizeof(int) + num_vd * sizeof(double);

    Packer p;

    // Compute the required buffer size.

    p.compute_buffer_size_mode();
    do_some_packing(p, vd, vi); // computes the size

    if ( total_size != p.size() ) ITFAILS;

    vector<char> buffer(p.size());

    // Pack into buffer.

    p.set_buffer(p.size(), &buffer[0]);
    do_some_packing(p, vd, vi); // actually does the packing

    // Unpack

    Unpacker u;

    u.set_buffer(p.size(), &buffer[0]);

    for ( int i = 0; i < vd.size(); ++i )
    {
	double d;
	u >> d;
	if ( ! soft_equiv(d, vd[i]) ) ITFAILS;
    }

    for ( int i = 0; i < vi.size(); ++i )
    {
	int j;
	u >> j;
	if ( j != vi[i] ) ITFAILS;
    }

    if (rtt_ds_test::passed)
	PASSMSG("compute_buffer_size_test() worked fine.");
}

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
#if DBC
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
#endif
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

#if DBC
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
#endif 

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

void std_string_test()
{
    vector<char> pack_string;

    {
	// make a string
	string hw("Hello World"); 

	// make a packer 
	Packer packer;

	// make a char to write the string into
	pack_string.resize(hw.size() + 1 * sizeof(int));

	packer.set_buffer(pack_string.size(), &pack_string[0]);

	packer << static_cast<int>(hw.size());

	// pack it
	for (string::const_iterator it = hw.begin(); it != hw.end(); it++)
	    packer << *it;
	
	if (packer.get_ptr() != &pack_string[0] + pack_string.size()) ITFAILS;
	if (packer.get_ptr() != packer.begin()  + pack_string.size()) ITFAILS;
    }

    // now unpack it
    Unpacker unpacker;
    unpacker.set_buffer(pack_string.size(), &pack_string[0]);
    
    // unpack the size of the string
    int size;
    unpacker >> size;

    string nhw;
    nhw.resize(size);

    // unpack the string
    for (string::iterator it = nhw.begin(); it != nhw.end(); it++)
	unpacker >> *it;

    if (unpacker.get_ptr() != &pack_string[0] + pack_string.size()) ITFAILS;

    // test the unpacked string
    // make a string
    string hw("Hello World"); 

    if (hw == nhw)
    {
	ostringstream message;
	message << "Unpacked string " << nhw << " that matches original "
		<< "string " << hw;
	PASSMSG(message.str());
    }
    else
    {
	ostringstream message;
	message << "Failed to unpack string " << hw << " correctly. Instead "
		<< "unpacked " << nhw;
	FAILMSG(message.str());
    }
}

//---------------------------------------------------------------------------//
 
void packing_functions_test()
{
    // make a packing container
    vector<char> total_packed;

    vector<double> x_ref;
    string         y_ref("Todd Urbatsch is an ass!");
    
    vector<double> x_new;
    string         y_new;

    // pack some data
    {
	vector<char> packed_int(sizeof(int));
	vector<char> packed_vector;
	vector<char> packed_string;

	vector<double> x(5);
	string         y("Todd Urbatsch is an ass!");

	x_ref.resize(5);

	for (int i = 0; i < 5; i++)
	{
	    x[i]     = 100.0 * (static_cast<double>(i) + x[i]) + 2.5;
	    x_ref[i] = x[i];
	}

	pack_data(x, packed_vector);
	pack_data(y, packed_string);

	if (packed_vector.size() != 5 * sizeof(double) + sizeof(int)) ITFAILS;
	if (packed_string.size() != y_ref.size() + sizeof(int))       ITFAILS;

	Packer p;
	p.set_buffer(sizeof(int), &packed_int[0]);

	// pack the size of the vector
	p << static_cast<int>(packed_vector.size());

	// push the vector onto the total packed
	total_packed.insert(total_packed.end(),
			    packed_int.begin(), packed_int.end());
	total_packed.insert(total_packed.end(),
			    packed_vector.begin(), packed_vector.end());

	// reset the packer
	p.set_buffer(sizeof(int), &packed_int[0]);
	
	// pack the size of the string
	p << static_cast<int>(packed_string.size());

	// push the string onto the total packed
	total_packed.insert(total_packed.end(),
			    packed_int.begin(), packed_int.end());
	total_packed.insert(total_packed.end(),
			    packed_string.begin(), packed_string.end());

	if (total_packed.size() != 2*packed_int.size() + packed_vector.size()
	    + packed_string.size()) ITFAILS;
    }

    // unpack the data
    {
	int size;
	Unpacker u;
	u.set_buffer(total_packed.size(), &total_packed[0]);
	
	// unpack the packed vector
	u >> size;
	vector<char> packed_vector(size);

	for (int i = 0; i < size; i++)
	    u >> packed_vector[i];

	// unpack the packed string
	u >> size;
	vector<char> packed_string(size);

	for (int i = 0; i < size; i++)
	    u >> packed_string[i];

	if (u.get_ptr() != &total_packed[0] + total_packed.size()) ITFAILS;

	unpack_data(x_new, packed_vector);
	unpack_data(y_new, packed_string);
    }

    if (!soft_equiv(x_new.begin(), x_new.end(), x_ref.begin(), x_ref.end()))
	ITFAILS;

    if (y_new == y_ref)
    {
	ostringstream message;
	message << "Correctly unpacked string " << y_new << " from "
		<< "packed " << y_ref;
	PASSMSG(message.str());
    }
    else
    {
	ostringstream message;
	message << "Failed to unpack string " << y_ref << " correctly. Instead "
		<< "unpacked " << y_new;
	FAILMSG(message.str());
    }

    if (rtt_ds_test::passed)
	PASSMSG("pack_data and unpack_data work fine.");
}

//---------------------------------------------------------------------------//

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

    try
    {
	// >>> UNIT TESTS
	packing_test();
	
	std_string_test();

	packing_functions_test();

	compute_buffer_size_test();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tstPacking_Utils, " << ass.what()
	     << endl;
	return 1;
    }

    // status of test
    cout << endl;
    cout <<     "*********************************************" << endl;
    if (rtt_ds_test::passed) 
    {
        cout << "**** tstPacking_Utils Test: PASSED" 
	     << endl;
    }
    cout <<     "*********************************************" << endl;
    cout << endl;

    cout << "Done testing tstPacking_Utils." << endl;
}   

//---------------------------------------------------------------------------//
//                        end of tstPacking_Utils.cc
//---------------------------------------------------------------------------//
