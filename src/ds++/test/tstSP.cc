//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstSP.cc
 * \author Thomas M. Evans
 * \date   Wed Feb  5 17:29:59 2003
 * \brief  SP test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "ds_test.hh"
#include "../Release.hh"
#include "../SP.hh"
#include "../Assert.hh"

#include <iostream>
#include <vector>
#include <cmath>
#include <typeinfo>
#include <sstream>

using namespace std;

using rtt_dsxx::SP;

//---------------------------------------------------------------------------//
// TEST HELPERS
//---------------------------------------------------------------------------//

int nfoos = 0;
int nbars = 0;
int nbazs = 0;
int nbats = 0;

#define CHECK_0_OBJECTS if (nfoos != 0) ITFAILS; if (nbars != 0) ITFAILS; if (nbazs != 0) ITFAILS; if (nbats != 0) ITFAILS;

#define CHECK_N_OBJECTS(nf, nb, nbz, nbt) if (nfoos != nf) ITFAILS; if (nbars != nb) ITFAILS; if (nbazs != nbz) ITFAILS; if (nbats != nbt) ITFAILS;

//---------------------------------------------------------------------------//

class Foo
{
  private:
    int v;

  public:
    Foo() 
	: v(0)
    {
	nfoos++;
    }

    explicit Foo(int i)
	: v(i)
    {
	nfoos++;
    }

    Foo(const Foo &f)
	: v(f.v)
    {
	nfoos++;
    }

    virtual ~Foo()
    {
	nfoos--;
    }

    virtual int vf() { return v; }

    int f() { return v+1; }
};

//---------------------------------------------------------------------------//

class Bar : public Foo
{
  private:
    Bar(const Bar &);

  public:
    explicit Bar(int i) 
	: Foo(i)
    {
	nbars++;
    }

    virtual ~Bar()
    {
	nbars--;
    }

    virtual int vf() { return Foo::f() + 1; }

    int f() { return Foo::f() + 2; }
};

//---------------------------------------------------------------------------//

class Baz : public Bar
{
  private:
    Baz(const Baz &);

  public:
    explicit Baz(int i) 
	: Bar(i)
    {
	nbazs++;
    }

    virtual ~Baz()
    {
	nbazs--;
    }

    virtual int vf() { return Bar::f() + 1; }

    int f() { return Bar::f() + 2; }
};

//---------------------------------------------------------------------------//

class Wombat 
{
  private:
    Wombat(const Wombat &);

  public:
    Wombat() { nbats++; }
    virtual ~Wombat() { nbats--; }
};

//---------------------------------------------------------------------------//
 
SP<Foo> get_foo()
{
    SP<Foo> f(new Foo(10));
    return f;
}

//---------------------------------------------------------------------------//

SP<Bar> get_bar()
{
    SP<Bar> b(new Bar(20));
    return b;
}

//---------------------------------------------------------------------------//

void test_foobar(SP<Foo> f, int v)
{
    if (f->vf() != v) ITFAILS;
}

//---------------------------------------------------------------------------//

void kill_SPBar(SP<Bar> &b)
{
    b = SP<Bar>();
}

//---------------------------------------------------------------------------//

void temp_change_SP(SP<Foo> f)
{
    CHECK_N_OBJECTS(1, 1, 0, 0);

    // this is a temporary change
    f = new Foo(100);

    CHECK_N_OBJECTS(2, 1, 0, 0);

    if (f->vf() != 100) ITFAILS;

    if (rtt_ds_test::passed)
	PASSMSG("SP<Bar> successfully (temporarily) reassigned to SP<Foo>.");
}

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
// here we test the following SP members:
//
//    SP();
//    SP(T *);
//    SP(const SP<T> &);
//    SP<T>& operator=(T *);
//    SP<T>& operator=(const SP<T> &);
//    T* operator->() const;
//    bool operator==(const T *) const;
//    bool operator!=(const T *) const;
//    bool operator==(const SP<T> &) const;
//    bool operator!=(const SP<T> &) const;
// 
// plus we test 
//
//    bool operator==(const T *, const SP<T> &);
//    bool operator!=(const T *, const SP<T> &);
//
void type_T_test()
{
    CHECK_0_OBJECTS;

    // test explicit constructor for type T *
    {
	// make a Foo, Bar, and Baz
	SP<Foo> spfoo(new Foo(1));
	SP<Bar> spbar(new Bar(2));
	SP<Baz> spbaz(new Baz(3));
	
	// there should be 3 Foos, 2 Bars and 1 Baz
	CHECK_N_OBJECTS(3, 2, 1, 0);
}

    // now all should be destroyed
    CHECK_0_OBJECTS;

    if (rtt_ds_test::passed)
	PASSMSG("Explicit constructor for type T * ok.");

    // test copy constructor for type T *
    {
	SP<Foo> rspfoo;
	SP<Bar> rspbar;
	SP<Baz> rspbaz;
	{
	    // no objects yet
	    CHECK_0_OBJECTS;

	    Foo *f  = new Foo(1);
	    Bar *b  = new Bar(2);
	    Baz *bz = new Baz(3);
	    
	    SP<Foo> spfoo(f);
	    SP<Bar> spbar(b);
	    SP<Baz> spbaz(bz);
	    
	    // there should be 3 Foos, 2 Bars and 1 Baz
	    CHECK_N_OBJECTS(3, 2, 1, 0);

	    // now assign
	    rspfoo = spfoo;
	    rspbar = spbar;
	    rspbaz = spbaz;
	    
	    // there are no additional objects made because the SP will make
	    // additional references
	    CHECK_N_OBJECTS(3, 2, 1, 0);

	    if (rtt_ds_test::passed)
		PASSMSG("Assignment of SP<T> ok.");

	    // now copy construct
	    SP<Foo> ispfoo = rspfoo;
	    SP<Bar> ispbar = rspbar;

	    // still no new foos created
	    CHECK_N_OBJECTS(3, 2, 1, 0);
	    if (ispfoo->f() != 2) ITFAILS;
	    if (spfoo->f()  != 2) ITFAILS;
	    if (rspfoo->f() != 2) ITFAILS;

	    if (rtt_ds_test::passed)
		PASSMSG("Copy construct of SP<T> ok.");

	    // now make a foo pointer and assign
	    Foo *ff = new Foo(10);
	    ispfoo  = ff;

	    // still no new foos created
	    CHECK_N_OBJECTS(4, 2, 1, 0);
	    if (ispfoo->f() != 11) ITFAILS;
	    if (spfoo->f()  != 2)  ITFAILS;
	    if (rspfoo->f() != 2)  ITFAILS;

	    if (rtt_ds_test::passed)
		PASSMSG("Assignment of T* ok.");

	    // now we can check equality 
	    if (rspfoo == spfoo)
	    {
		PASSMSG("Equality operation ok.");
	    }
	    else
	    {
		FAILMSG("Equality operation failed.");
	    }
	    
	    // now check inequality
	    if (rspfoo != spfoo)
	    {
		FAILMSG("Equality operation failed.");
	    }
	    else
	    {
		PASSMSG("Equality operation ok.");
	    }

	    if (rspbar != spbar) ITFAILS;
	    if (rspbaz != spbaz) ITFAILS;

	    if (spfoo != f)      ITFAILS;
	    if (spbar != b)      ITFAILS;
	    if (spbaz != bz)     ITFAILS;

	    if (spfoo == b)      ITFAILS; // this is ok because a Bar * can
					  // be passed to Foo *

	    if (spbar == dynamic_cast<Bar *>(f)) ITFAILS;

	    if (f  != spfoo)     ITFAILS;
	    if (b  != spbar)     ITFAILS;
	    if (bz != spbaz)     ITFAILS;

	    if (f == spfoo)
	    {
		PASSMSG("Overloaded equality operators ok.");
	    }
	    else
	    {
		FAILMSG("Overloaded equality operators failed.");
	    }

	    if (rtt_ds_test::passed)
		PASSMSG("Equality/Inequality operations ok.");
	}
	
	// we should still have objects left even because we still have
	// viable SPs in scope
	CHECK_N_OBJECTS(3, 2, 1, 0);
    }

    // now all should be destroyed
    CHECK_0_OBJECTS;

    if (rtt_ds_test::passed)
	PASSMSG("Operations on type T ok");
}

//---------------------------------------------------------------------------//
// here we test the following SP members:
//
//    SP();
//    SP(X *);
//    SP<const SP<X> &);
//    SP<T>& operator=(T *);
//    SP<T>& operator=(X *);
//    SP<T>& operator=(const SP<X> &);
//    T* operator->() const;
//    T& operator*() const;
//    T* bp() const;
//    operator bool() const;
//    bool operator!() const;
// 
void type_X_test()
{
    CHECK_0_OBJECTS;

    // check explicit constructor
    {
	// make a foo pointer
	SP<Foo> spfoo(new Bar(10));
	CHECK_N_OBJECTS(1, 1, 0, 0);

	if (spfoo->vf() != 12) ITFAILS;
	if (spfoo->f() != 11)  ITFAILS;

	Foo &f = *spfoo;
	if (f.f() != 11)  ITFAILS;
	if (f.vf() != 12) ITFAILS;

	Foo ff = *spfoo;
	if (ff.vf() != 10) ITFAILS;

	Bar *b = dynamic_cast<Bar *>(spfoo.bp());
	if (b->vf() != 12) ITFAILS;
	if (b->f() != 13)  ITFAILS;

	if (typeid(spfoo.bp()) != typeid(Foo *)) ITFAILS;
	if (typeid(*spfoo.bp()) != typeid(Bar))  ITFAILS;

	CHECK_N_OBJECTS(2, 1, 0, 0);
    }

    CHECK_0_OBJECTS;

    if (rtt_ds_test::passed)
	PASSMSG("Explicit constructor for type X * ok.");

    // check SP<X> constructor and assignment
    {
	// make some objects
	SP<Foo> spfoo;
	SP<Bar> spbar;
	SP<Foo> spfoo2;

	if (spfoo)  ITFAILS;
	if (spbar)  ITFAILS;
	if (spfoo2) ITFAILS;
	{
	    spbar = new Bar(50);
	    CHECK_N_OBJECTS(1, 1, 0, 0);

	    if (spbar->f() != 53)  ITFAILS;
	    if (spbar->vf() != 52) ITFAILS;

	    // now assign to base class SP
	    spfoo = spbar;
	    CHECK_N_OBJECTS(1, 1, 0, 0);

	    if (spfoo->f() != 51)  ITFAILS;
	    if (spfoo->vf() != 52) ITFAILS;

	    if (typeid(spfoo.bp()) != typeid(Foo *)) ITFAILS;
	    if (typeid(*spfoo.bp()) != typeid(Bar))  ITFAILS;
	    if (typeid(spbar.bp()) != typeid(Bar *)) ITFAILS;

	    if (rtt_ds_test::passed)
		PASSMSG("Assignment with SP<X> ok.");

	    // now do copy construction
	    SP<Foo> rspfoo(spbar);
	    CHECK_N_OBJECTS(1, 1, 0, 0);

	    if (rspfoo->f() != 51)  ITFAILS;
	    if (rspfoo->vf() != 52) ITFAILS;

	    if (typeid(rspfoo.bp()) != typeid(Foo *)) ITFAILS;
	    if (typeid(*rspfoo.bp()) != typeid(Bar))  ITFAILS;
	    
	    if (rtt_ds_test::passed)
		PASSMSG("Copy constructor with SP<X> ok.");

	    // now check assignment with X *
	    rspfoo = new Bar(12);
	    CHECK_N_OBJECTS(2, 2, 0, 0);

	    if (rspfoo->f() != 13)  ITFAILS;
	    if (rspfoo->vf() != 14) ITFAILS;

	    if (typeid(rspfoo.bp()) != typeid(Foo *)) ITFAILS;
	    if (typeid(*rspfoo.bp()) != typeid(Bar))  ITFAILS;

	    if (rtt_ds_test::passed)
		PASSMSG("Assignment with X * ok.");

	    // assign SPfoo2 to a bar
	    spfoo2 = new Bar(20);
	    CHECK_N_OBJECTS(3, 3, 0, 0);
	    
	}
	// still have 2 object
	CHECK_N_OBJECTS(2, 2, 0, 0);

	// assign spfoo to a baz
	spfoo2 = new Baz(45);
	CHECK_N_OBJECTS(2, 2, 1, 0);

	if (spfoo2->f() != 46)  ITFAILS;
	if (spfoo2->vf() != 49) ITFAILS;

	if (typeid(*spfoo2.bp()) != typeid(Baz)) ITFAILS;

	// assign spbar to NULL
	spbar = SP<Bar>();
	CHECK_N_OBJECTS(2, 2, 1, 0);

	// spfoo should still point to the same bar
	if (spfoo->f() != 51)  ITFAILS;
	if (spfoo->vf() != 52) ITFAILS;
	
	if (typeid(spfoo.bp()) != typeid(Foo *)) ITFAILS;
	if (typeid(*spfoo.bp()) != typeid(Bar))  ITFAILS;
	if (typeid(spbar.bp()) != typeid(Bar *)) ITFAILS;
	
	if (rtt_ds_test::passed)
	    PASSMSG("Set to SP<>() releases pointer.");

	// assign spfoo to NULL
	spfoo = SP<Foo>();
	CHECK_N_OBJECTS(1, 1, 1, 0);

	if (spfoo) ITFAILS;
	if (spbar) ITFAILS;

	if (rtt_ds_test::passed)
	    PASSMSG("Overloaded bool ok.");

	if (!spfoo2) ITFAILS;

	if (rtt_ds_test::passed)
	    PASSMSG("Overloaded ! (not) ok."); 
    }

    CHECK_0_OBJECTS;

    if (rtt_ds_test::passed)
	PASSMSG("Operations on type X ok");
}

//---------------------------------------------------------------------------//

void fail_modes_test()
{
    // make an object and try to reference it
    SP<Foo>    spfoo;
    SP<Bar>    spbar;
    SP<Baz>    spbaz;
    SP<Wombat> spbat;

    if (spfoo) ITFAILS;
    if (spfoo) ITFAILS;
    if (spfoo) ITFAILS;

    // try to reference a function
    bool caught = false;
#if DBC
    try
    {
	spfoo->f();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	caught = true;
	ostringstream m;
	m << "Good, caught the following exception, " << endl
	  << ass.what();
	PASSMSG(m.str());
    }
    if (!caught)
	FAILMSG("Failed to catch illegal access exception.");

    // try assigning a derived NULL to a base; the spfoo base pointer is
    // still a foo in the case 
    spfoo = spbar;
    if (typeid(spfoo.bp()) != typeid(Foo *)) ITFAILS;

    CHECK_0_OBJECTS;

    // now try assigning to a non-derived class of Foo that is NULL,
    // unfortunately, even though this shouldn't be allowed we get away with
    // it because wombat has some virtual functions and is NULL; however,
    // this isn't really dangerous because spfoo still doesn't point to
    // anything 
    spfoo  = spbat;
    caught = false;
    try
    {
	spfoo->f();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	caught = true;
	ostringstream m;
	m << "Good, caught the following exception, " << endl
	  << ass.what();
	PASSMSG(m.str());
    }
    if (!caught)
	FAILMSG("Failed to catch illegal access exception.");
#endif
    // now make a wombat and try
    spbat  = new Wombat;
    CHECK_N_OBJECTS(0, 0, 0, 1);

    caught = false;
    try
    {
	spfoo = spbat;
    }
    catch (rtt_dsxx::assertion &ass)
    {
	caught = true;
	ostringstream m;
	m << "Good, caught the following exception, " << endl
	  << ass.what();
	PASSMSG(m.str());
    }
    if (!caught)
	FAILMSG("Failed to catch illegal assignment.");

    // now try copy construction
    caught = false;
    try
    {
	SP<Foo> s(spbat);
    }
    catch (rtt_dsxx::assertion &ass)
    {
	caught = true;
	ostringstream m;
	m << "Good, caught the following exception, " << endl
	  << ass.what();
	PASSMSG(m.str());
    }
    if (!caught)
	FAILMSG("Failed to catch illegal assignment.");

    // now try copy and assignment on X *
    Wombat *bat = new Wombat();
    CHECK_N_OBJECTS(0, 0, 0, 2);

    caught = false;
    try
    {
	spfoo = bat;
    }
    catch (rtt_dsxx::assertion &ass)
    {
	caught = true;
	ostringstream m;
	m << "Good, caught the following exception, " << endl
	  << ass.what();
	PASSMSG(m.str());
    }
    if (!caught)
	FAILMSG("Failed to catch illegal assignment.");

    // now try copy construction
    caught = false;
    try
    {
	SP<Foo> s(bat);
    }
    catch (rtt_dsxx::assertion &ass)
    {
	caught = true;
	ostringstream m;
	m << "Good, caught the following exception, " << endl
	  << ass.what();
	PASSMSG(m.str());
    }
    if (!caught)
	FAILMSG("Failed to catch illegal assignment.");
    
    // assign wombat to a pointer to clean it up
    spbat = bat;
    CHECK_N_OBJECTS(0, 0, 0, 1);
    
    if (rtt_ds_test::passed)
	PASSMSG("Failure modes work ok.");
}

//---------------------------------------------------------------------------//

void equality_test()
{
    CHECK_0_OBJECTS;

    // try some more advanced stuff
    SP<Foo> f1;
    SP<Foo> f2;

    Foo *f  = new Foo(5);
    Foo *ff = new Foo(5);

    f1 = f;
    f2 = f1;

    if (f2 != f1) ITFAILS;

    // now f and ff are equivalent, but the smart pointers won't be because
    // they don't point to the same instance of Foo *
    f2 = ff;

    if (f2 == f1) ITFAILS;

    CHECK_N_OBJECTS(2, 0, 0, 0);

    if (rtt_ds_test::passed)
	PASSMSG("Equality tests work ok.");
}

//---------------------------------------------------------------------------//

void get_test()
{
    CHECK_0_OBJECTS;

    // get a foo and bar
    {
	
	SP<Foo> f  = get_foo();
	SP<Foo> fb = get_bar();
	SP<Bar> b  = get_bar();
	
	CHECK_N_OBJECTS(3, 2, 0, 0);
	
	if (fb == b) ITFAILS;
	
	if (f->f() != 11)   ITFAILS;
	if (fb->vf() != 22) ITFAILS;
	if (b->vf() != 22)  ITFAILS;
	
	if (fb->f() != 21)  ITFAILS;
	if (b->f() != 23)   ITFAILS;
	
    }

    if (rtt_ds_test::passed)
	PASSMSG("Get/factory tests work ok.");
}

//---------------------------------------------------------------------------//

void access_test()
{
    CHECK_0_OBJECTS;

    SP<Bar> b(new Bar(10));
    CHECK_N_OBJECTS(1, 1, 0, 0);

    test_foobar(b, 12);
    CHECK_N_OBJECTS(1, 1, 0, 0);

    kill_SPBar(b);
    CHECK_0_OBJECTS;

    b = new Bar(12);
    temp_change_SP(b); // this temporarily changes to a Foo

    if (b->vf() != 14)                   ITFAILS;
    if (typeid(b.bp()) != typeid(Bar *)) ITFAILS;

    CHECK_N_OBJECTS(1, 1, 0, 0);

    if (rtt_ds_test::passed)
	PASSMSG("Accessor/set-style tests work ok.");
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

	CHECK_0_OBJECTS;

	type_T_test();
	cout << endl;

	type_X_test();
	cout << endl;
	
	fail_modes_test();
	cout << endl;

	equality_test();
	cout << endl;

	get_test();
	cout << endl;

	access_test();

	CHECK_0_OBJECTS;
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tstSP, " << ass.what()
	     << endl;
	return 1;
    }

    // status of test
    cout << endl;
    cout <<     "*********************************************" << endl;
    if (rtt_ds_test::passed) 
    {
        cout << "**** tstSP Test: PASSED" 
	     << endl;
    }
    cout <<     "*********************************************" << endl;
    cout << endl;

    cout << "Done testing tstSP." << endl;
    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of tstSP.cc
//---------------------------------------------------------------------------//
