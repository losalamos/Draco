//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstRCF.cc
 * \author Thomas M. Evans
 * \date   Wed Jan 28 10:53:26 2004
 * \brief  Test of RCF (reference counted field) class.
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

#include "../Release.hh"
#include "../Soft_Equivalence.hh"
#include "../RCF.hh"
#include "ds_test.hh"

using namespace std;

using rtt_dsxx::RCF;
using rtt_dsxx::soft_equiv;

typedef vector<double> dbl_field;

//---------------------------------------------------------------------------//
// TESTING INFRASTRUCTURE
//---------------------------------------------------------------------------//

int nfields = 0;

class Field
{
  public:
    typedef dbl_field::value_type     value_type;
    typedef dbl_field::size_type      size_type;
    typedef dbl_field::iterator       iterator;
    typedef dbl_field::const_iterator const_iterator;

  private:
    dbl_field d;

  public:
    Field() : d(5, 1.0) { nfields++; }
    Field(int n, value_type v = value_type()) : d(n,v) { nfields++; }
    ~Field() { nfields--; }

    value_type& operator[](int i) { return d[i]; }
    const value_type& operator[](int i) const { return d[i]; }
    size_t size() const { return d.size(); }
    bool empty() const { return d.empty(); }

    const_iterator begin() const { return d.begin(); }
    iterator begin() { return d.begin(); }

    const_iterator end() const { return d.end(); }
    iterator end() { return d.end(); }
};

//---------------------------------------------------------------------------//

RCF<Field> get_field()
{
    RCF<Field> f(new Field);
    if (nfields != 1) ITFAILS;

    return f;
}

//---------------------------------------------------------------------------//

void use_const_field(const RCF<dbl_field> &f, const dbl_field &ref)
{
    // test const_iterator access
    if (!soft_equiv(f.begin(), f.end(), ref.begin(), ref.end())) ITFAILS;

    // get the field to test const get_field
    const dbl_field &field = f.get_field();
    if (!soft_equiv(field.begin(), field.end(), ref.begin(), ref.end())) 
        ITFAILS;

    // check constant operator[] access
    for (int i = 0; i < f.size(); i++)
    {
        if (!soft_equiv(f[i], ref[i])) ITFAILS;
    }
}


//---------------------------------------------------------------------------//

void use_const_field(const Field &f, const dbl_field &ref)
{
    if (!soft_equiv(f.begin(), f.end(), ref.begin(), ref.end())) ITFAILS;
}

//---------------------------------------------------------------------------//

void use_non_const_field(Field &f)
{
    // change element 2
    f[1] = 13.231;
} 

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
// Test the following:
//
//   RCF()
//   RCF(Field_t *)
//   RCF(n, v)
//   RCF(begin, end)
//   operator=(Field_t *)
//   begin() const
//   end() const
//   begin()
//   end()
//   size()
//   operator[]()
//   operator[]() const
//   empty()
//   get_field()
//   get_field() const
//   assigned()
// 

void test_simple_construction_copy()
{
    // make a smart field on a vector of doubles
    RCF<dbl_field> sf;
    RCF<dbl_field> y;
    if (sf.assigned()) ITFAILS;

    {
        sf = new dbl_field(10, 5.2);
        if (!sf.assigned()) ITFAILS;

        dbl_field ref(10, 5.2);

        if (!soft_equiv(sf.begin(), sf.end(), ref.begin(), ref.end())) ITFAILS;

        // fill in 2.2 (tests non-const begin and end)
        fill(sf.begin(), sf.end(), 2.2);
        fill(ref.begin(), ref.end(), 2.2);

        if (!soft_equiv(sf.begin(), sf.end(), ref.begin(), ref.end())) ITFAILS;
        
        // check size
        if (sf.size() != 10) ITFAILS;
        if (sf.empty())      ITFAILS;

        // check with subscript access
        for (int i = 0; i < sf.size(); i++)
        {
            if (!soft_equiv(sf[i], 2.2)) ITFAILS;
            
            // reassign and check
            sf[i] = 12.46;
            if (!soft_equiv(sf[i], 12.46)) ITFAILS;
        }

        fill(ref.begin(), ref.end(), 12.46);
        if (!soft_equiv(sf.begin(), sf.end(), ref.begin(), ref.end())) ITFAILS;

        // check const functions
        use_const_field(sf, ref);

        // get field and empty it
        sf.get_field().resize(0);
        if (!sf.empty()) ITFAILS;

        // make a field, using alternative ctor.
        RCF<dbl_field> x(10, 12.46);
        if (!x.assigned()) ITFAILS;
        if (!soft_equiv(x.begin(), x.end(), ref.begin(), ref.end())) ITFAILS;

        // assign it to x
        y = x;

        // change x (which also changes y)
        x.get_field().resize(2);
        x[0] = 1.1;
        y[1] = 1.2;

        if (y.size() != 2) ITFAILS;

        if (y[0] != 1.1) ITFAILS;
        if (x[1] != 1.2) ITFAILS;

	// check range constructor
	RCF<dbl_field> z(x.begin(), x.end());
        if (!soft_equiv(x.begin(), x.end(), z.begin(), z.end())) ITFAILS;
    }

    if (!sf.assigned()) ITFAILS;
    if (!y.assigned())  ITFAILS;

    if (y.size() != 2) ITFAILS;
    if (y[0] != 1.1) ITFAILS;
    if (y[1] != 1.2) ITFAILS;

    if (rtt_ds_test::passed)
        PASSMSG("Simple construction and copy ok.");
}

//---------------------------------------------------------------------------//

void test_counting()
{
    if (nfields != 0) ITFAILS;

    RCF<Field> f = get_field();
    if (!f.assigned()) ITFAILS;

    if (nfields != 1) ITFAILS;

    {
        RCF<Field> g = f;
        
        if (nfields != 1) ITFAILS;
    }
    
    if (nfields != 1) ITFAILS;

    dbl_field ref(5, 1.0);
    
    // check const field access
    use_const_field(f.get_field(), ref);
    if (!soft_equiv(f.begin(), f.end(), ref.begin(), ref.end())) ITFAILS;
    
    // check non-const field access
    use_non_const_field(f.get_field());
    ref[1] = 13.231;
    if (!soft_equiv(f.begin(), f.end(), ref.begin(), ref.end())) ITFAILS;

    RCF<Field> g;
    if (nfields != 1) ITFAILS;

    // test copying and assignment
    {
        g = f;
        if (nfields != 1) ITFAILS;
        f = new Field();
        if (nfields != 2) ITFAILS;
    }

    if (nfields != 2) ITFAILS;

    g = RCF<Field>();
    if (g.assigned()) ITFAILS;

    if (nfields != 1) ITFAILS;

    if (rtt_ds_test::passed)
        PASSMSG("Reference counting and copy construction ok.");
}

//---------------------------------------------------------------------------//

void test_constness()
{
    if (nfields != 0) ITFAILS;

    RCF<const Field> f = get_field();
    if (!f.assigned()) ITFAILS;

    if (nfields != 1) ITFAILS;

    {
        RCF<const Field> g = f;
        
        if (nfields != 1) ITFAILS;
    }
    
    if (nfields != 1) ITFAILS;

    dbl_field ref(5, 1.0);
    
    // check const field access
    use_const_field(f.get_field(), ref);
    if (!soft_equiv(f.begin(), f.end(), ref.begin(), ref.end())) ITFAILS;

    RCF<const Field> g;
    if (nfields != 1) ITFAILS;

    // test copying and assignment
    {
        g = f;
        if (nfields != 1) ITFAILS;
        f = new Field();
        if (nfields != 2) ITFAILS;
    }

    if (nfields != 2) ITFAILS;

    g = RCF<const Field>();
    if (g.assigned()) ITFAILS;

    if (nfields != 1) ITFAILS;

    if (rtt_ds_test::passed)
        PASSMSG("Constness tests ok.");
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (std::string(argv[arg]) == "--version")
	{
	    std::cout << argv[0] << ": version " 
		      << rtt_dsxx::release() 
		      << std::endl;
	    return 0;
	}

    try
    {
	// >>> UNIT TESTS

        test_simple_construction_copy();
        test_counting();
        test_constness();

        // make sure that the field number is zero
        if (nfields == 0)
        {
            PASSMSG("All fields destroyed.");
        }
        else
        {
            FAILMSG("Error in reference counting of fields.");
        }
    }
    catch (rtt_dsxx::assertion &ass)
    {
	std::cout << "While testing tstRCF, " << ass.what()
		  << std::endl;
	return 1;
    }

    // status of test
    std::cout << std::endl;
    std::cout <<     "*********************************************" 
	      << std::endl;
    if (rtt_ds_test::passed) 
    {
        std::cout << "**** tstRCF Test: PASSED" 
		  << std::endl;
    }
    std::cout <<     "*********************************************" 
	      << std::endl;
    std::cout << std::endl;
    
    std::cout << "Done testing tstRCF." << std::endl;
    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of tstRCF.cc
//---------------------------------------------------------------------------//
