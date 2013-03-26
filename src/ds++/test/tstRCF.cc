//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstRCF.cc
 * \author Thomas M. Evans
 * \date   Wed Jan 28 10:53:26 2004
 * \brief  Test of RCF (reference counted field) class.
 * \note   Copyright (C) 2003-2013 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "../Release.hh"
#include "../Soft_Equivalence.hh"
#include "../RCF.hh"
#include "../ScalarUnitTest.hh"
#include <vector>

using namespace std;
using namespace rtt_dsxx;

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
    Field(void) : d(5, 1.0) { nfields++; }
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

RCF<Field> get_field(UnitTest &ut)
{
    RCF<Field> f(new Field);
    if (nfields != 1) ut.failure("test fails");

    return f;
}

//---------------------------------------------------------------------------//

void use_const_field(RCF<dbl_field> const & f,
                     dbl_field      const & ref,
                     UnitTest &ut)
{
    // test const_iterator access
    if (!soft_equiv(f.begin(), f.end(), ref.begin(), ref.end()))
        ut.failure("test fails");

    // get the field to test const get_field
    const dbl_field &field = f.get_field();
    if (!soft_equiv(field.begin(), field.end(), ref.begin(), ref.end())) 
        ut.failure("test fails");

    // check constant operator[] access
    for (size_t i = 0; i < f.size(); i++)
        if (!soft_equiv(f[i], ref[i])) ut.failure("test fails");
}


//---------------------------------------------------------------------------//

void use_const_field(const Field &f,
                     const dbl_field &ref,
                     UnitTest &ut)
{
    if (!soft_equiv(f.begin(), f.end(), ref.begin(), ref.end()))
        ut.failure("test fails");
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

void test_simple_construction_copy(UnitTest &ut)
{
    // make a smart field on a vector of doubles
    RCF<dbl_field> sf;
    RCF<dbl_field> y;
    if (sf.assigned()) ut.failure("test fails");
    {
        sf = new dbl_field(10, 5.2);
        if (!sf.assigned()) ut.failure("test fails");

        dbl_field ref(10, 5.2);

        if (!soft_equiv(sf.begin(), sf.end(), ref.begin(), ref.end()))
            ut.failure("test fails");

        // fill in 2.2 (tests non-const begin and end)
        fill(sf.begin(), sf.end(), 2.2);
        fill(ref.begin(), ref.end(), 2.2);

        if (!soft_equiv(sf.begin(), sf.end(), ref.begin(), ref.end()))
            ut.failure("test fails");
        
        // check size
        if (sf.size() != 10) ut.failure("test fails");
        if (sf.empty())      ut.failure("test fails");

        // check with subscript access
        for (size_t i = 0; i < sf.size(); i++)
        {
            if (!soft_equiv(sf[i], 2.2)) ut.failure("test fails");
            
            // reassign and check
            sf[i] = 12.46;
            if (!soft_equiv(sf[i], 12.46)) ut.failure("test fails");
        }

        fill(ref.begin(), ref.end(), 12.46);
        if (!soft_equiv(sf.begin(), sf.end(), ref.begin(), ref.end()))
            ut.failure("test fails");

        // check const functions
        use_const_field(sf, ref, ut);

        // get field and empty it
        sf.get_field().resize(0);
        if (!sf.empty()) ut.failure("test fails");

        // make a field, using alternative ctor.
        RCF<dbl_field> x(10, 12.46);
        if (!x.assigned()) ut.failure("test fails");
        if (!soft_equiv(x.begin(), x.end(), ref.begin(), ref.end()))
            ut.failure("test fails");

        // assign it to x
        y = x;

        // change x (which also changes y)
        x.get_field().resize(2);
        x[0] = 1.1;
        y[1] = 1.2;

        if (y.size() != 2) ut.failure("test fails");

        if (y[0] != 1.1) ut.failure("test fails");
        if (x[1] != 1.2) ut.failure("test fails");

	// check range constructor
	RCF<dbl_field> z(x.begin(), x.end());
        if (!soft_equiv(x.begin(), x.end(), z.begin(), z.end()))
            ut.failure("test fails");
    }

    if (!sf.assigned()) ut.failure("test fails");
    if (!y.assigned())  ut.failure("test fails");

    if (y.size() != 2) ut.failure("test fails");
    if (y[0] != 1.1) ut.failure("test fails");
    if (y[1] != 1.2) ut.failure("test fails");

    // test some RCF< const Field > functions.
    {
        // reference values
        int const len(5);
        dbl_field ref( len, 3.1415 );

        // test copy ctor
        {
            RCF< const dbl_field > cf( y );
            // test empty()
            if( cf.empty() ) ut.failure("test fails");
            // test size() -- same size as y!
            if( cf.size() != 2 ) ut.failure("test fails");
            // test bracket operator
            if( ! soft_equiv( cf[0], y[0] ) ) ut.failure("test fails");
        }

        // create a new RCF using alternate ctor
        {
            RCF< const dbl_field > cf2( 5, 3.1415 );
            if( ! cf2.assigned() ) ut.failure("test fails");
            if( ! soft_equiv( cf2.begin(), cf2.end(), 
                              ref.begin(), ref.end() )) ut.failure("test fails");
        }
        
        // check range constructor
        {
            RCF< const dbl_field > cf3( ref.begin(), ref.end() );
            if( ! cf3.assigned() ) ut.failure("test fails");
            if( ! soft_equiv( cf3.begin(), cf3.end(), 
                              ref.begin(), ref.end() )) ut.failure("test fails");
        }

        // check constructor from ptr to field
        {
            RCF< const Field > cf( new Field );
            if (!cf.assigned()) ut.failure("test fails");
        }
    }
    
    if (ut.numFails==0)
        ut.passes("Simple construction and copy ok.");
}

//---------------------------------------------------------------------------//

void test_counting(UnitTest &ut)
{
    if (nfields != 0) ut.failure("test fails");

    RCF<Field> f = get_field(ut);
    if (!f.assigned()) ut.failure("test fails");

    if (nfields != 1) ut.failure("test fails");

    {
        RCF<Field> g = f;
        
        if (nfields != 1) ut.failure("test fails");
    }
    
    if (nfields != 1) ut.failure("test fails");

    dbl_field ref(5, 1.0);
    
    // check const field access
    use_const_field(f.get_field(), ref, ut);
    if (!soft_equiv(f.begin(), f.end(), ref.begin(), ref.end()))
        ut.failure("test fails");
    
    // check non-const field access
    use_non_const_field(f.get_field());
    ref[1] = 13.231;
    if (!soft_equiv(f.begin(), f.end(), ref.begin(), ref.end()))
        ut.failure("test fails");

    RCF<Field> g;
    if (nfields != 1) ut.failure("test fails");

    // test copying and assignment
    {
        g = f;
        if (nfields != 1) ut.failure("test fails");
        Field *ptr = new Field();
        f = ptr;
        if (nfields != 2) ut.failure("test fails");
        f = ptr;
        if (nfields != 2) ut.failure("test fails");
    }

    if (nfields != 2) ut.failure("test fails");

    g = RCF<Field>();
    if (g.assigned()) ut.failure("test fails");

    if (nfields != 1) ut.failure("test fails");

    if (ut.numFails==0)
        ut.passes("Reference counting and copy construction ok.");
}

//---------------------------------------------------------------------------//

void test_constness(UnitTest &ut)
{
    if (nfields != 0) ut.failure("test fails");

    RCF<const Field> f = get_field(ut);
    if (!f.assigned()) ut.failure("test fails");

    if (nfields != 1) ut.failure("test fails");

    {
        RCF<const Field> g = f;
        
        if (nfields != 1) ut.failure("test fails");
    }
    
    if (nfields != 1) ut.failure("test fails");

    dbl_field ref(5, 1.0);
    
    // check const field access
    use_const_field(f.get_field(), ref, ut);
    if (!soft_equiv(f.begin(), f.end(), ref.begin(), ref.end()))
        ut.failure("test fails");

    RCF<const Field> g;
    if (nfields != 1) ut.failure("test fails");

    // test copying and assignment
    {
        g = f;
        if (nfields != 1) ut.failure("test fails");
        f = new Field();
        if (nfields != 2) ut.failure("test fails");
        g = const_cast<Field *>(&g.get_field());
        if (nfields != 2) ut.failure("test fails");
    }

    if (nfields != 2) ut.failure("test fails");

    g = RCF<const Field>();
    if (g.assigned()) ut.failure("test fails");

    if (nfields != 1) ut.failure("test fails");

    if (ut.numFails==0)
        ut.passes("Constness tests ok.");
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    ScalarUnitTest ut( argc, argv, release );
    try
    {
        test_simple_construction_copy(ut);
        test_counting(ut);
        test_constness(ut);

        // make sure that the field number is zero
        if (nfields == 0)
            ut.passes("All fields destroyed.");
        else
            ut.failure("Error in reference counting of fields.");
    }
    UT_EPILOG(ut);
}   

//---------------------------------------------------------------------------//
// end of tstRCF.cc
//---------------------------------------------------------------------------//
