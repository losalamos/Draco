//----------------------------------*-C++-*----------------------------------//
// tstSP.cc
// Geoffrey Furnish
// 3 December 1994
//---------------------------------------------------------------------------//
// @> Test program for the SP<T> class
//---------------------------------------------------------------------------//

#include <iostream.h>

#include "SP.hh"
// #include "List.cc"
// #include "Map.cc"

// autodoc: noprint foo

int nfoos = 0;

class foo {

    int v;

    foo( const foo& );
    
  public:
    foo() { v=0; nfoos++; cout << "foo constructed.\n"; }
    foo( int i ) { v=i; nfoos++; cout << "foo constructed.\n";; }
    virtual ~foo() { nfoos--; cout << "foo destroyed.\n"; }
    void method() { cout << "foo::method invoked.\n"; }
    int val() const { return v; }
};

class bar : public foo
{
};

class baz : public bar
{
};

class wombat {};

#ifdef __GNUC__
INSTANTIATE_SP(foo);
INSTANTIATE_SPList(foo);
INSTANTIATE_Map(int,SP<foo>);
#endif

#ifdef __DECCXX
#pragma define_template SP<foo>
#pragma define_template SPrep<foo>
#pragma define_template SPList<foo>
#pragma define_template SPList_iter<foo>
#pragma define_template SPLink<foo>
#pragma define_template Map<int, SP<foo> >
#pragma define_template Link<int, SP<foo> >
#pragma define_template Map_iter<int, SP<foo> >
#endif

//---------------------------------------------------------------------------//
// Function to help us keep track of the foo pool.
//---------------------------------------------------------------------------//

int expect( int n ) {
    cout << "Expecting " << n << " foo's, ";
    if (nfoos == n) {
	cout << "good.\n";
	return 1;
    } else {
	cout << "but there are " << nfoos << ". <<<<<<<<<<<<<\n";
	return 0;
    }
}

SP<foo> getafoo()
{
    SP<foo> p = new foo;
    return p;
}

SP<foo> useafoo( SP<foo> s )
{
    s->method();

    SP<foo> p = new foo;

    p = s;

    cout << "a foo should've just been destroyed.\n";

    return p;
}

#ifdef __GNUC__
INSTANTIATE_Slist(SP<foo>);
#endif

void tst_map();

void x1()
{
    cout << "\n\n test x1 -- mostly same-pointer-type ops.\n";

    {
	cout << "\n Check foo ctor/dtor." << endl;
	expect(0);
	foo *pf = new foo;
	expect(1);
	delete pf;
	expect(0);
    }

    {
	cout << "\n Check SP def ctor/dtor." << endl;
	expect(0);
	{
	    SP<foo> spf;
	    expect(0);
	}
	expect(0);
    }

    {
	cout << "\n Check SP(T*) ctor/dtor." << endl;
	expect(0);
	{
	    foo *pf = new foo;
	    SP<foo> spf( pf );
	    expect(1);
	}
	cout << "A foo should've just been destroyed." << endl;
	expect(0);
    }

    {
	cout << "\n Check SP<foo> s = new foo, which ops are those?" << endl;
	expect(0);
	{
	    SP<foo> spf = new foo;
	    expect(1);
	}
	expect(0);
    }

    {
	cout << "\n Check assign from T *." << endl;
	expect(0);
	{
	    SP<foo> spf;
	    foo *pf = new foo;
	    expect(1);
	    spf = pf;
	    expect(1);
	}
	expect(0);
    }
    
    {
	cout << "\n Check SP assign from like type." << endl;
	expect(0);
	{
	    SP<foo> sp1 = new foo;
	    expect(1);
	    SP<foo> sp2 = new foo;
	    expect(2);
	    sp1 = sp2;
	    expect(1);
	}
	expect(0);
    }

    {
	cout << " \n Check SP def ctor + copy ctor," << endl;
	expect(0);
	{
	    SP<foo> sp1;
	    expect(0);
	    SP<foo> sp2 = new foo;
	    expect(1);
	    sp1 = sp2;
	    expect(1);
	}
	expect(0);
    }

    {
	cout << "\n Check SP copy ctor from like type." << endl;
	expect(0);
	{
	    SP<foo> sp1 = new foo;
	    expect(1);
	    SP<foo> sp2( sp1 );
	    expect(1);
	}
	expect(0);
    }
    
    {
	cout << "\n Check getafoo." << endl;
	expect(0);
	{
	    SP<foo> spf = getafoo();
	    expect(1);
	}
	expect(0);
    }
}

void x2()
{
    cout << "\n\n test x2 -- test compatible-pointer-type ops.\n";

    {
	cout << "\n Check derived class ctor/dtor." << endl;
	expect(0);
	{
	    bar *b = new bar;
	    expect(1);
	    delete b;
	}
	expect(0);
    }

    {
	cout << "\n Check X* ctor." << endl;
	expect(0);
	{
	    SP<foo> spf = new bar;
	    expect(1);
	}
	expect(0);
    }

    {
	cout << "\n Check SP<X> copy ctor." << endl;
	expect(0);
	{
	    SP<foo> spf = new bar;
	    expect(1);
	    SP<bar> spb( spf );
	    expect(1);
	}
	expect(0);
    }

    {
	cout << "\n Check assign from SP<X>." << endl;
	expect(0);
	{
	    SP<foo> spf;
	    SP<bar> spb = new bar;
	    expect(1);
	    spf = spb;
	    expect(1);
	}
	expect(0);
    }

    {
	cout << "\n Check downcast copy ctor." << endl;
	expect(0);
	{
	    SP<foo> spf = new bar;
	    SP<bar> spb( spf );
	    expect(1);
	}
	expect(0);
    }

    {
	cout << "\n Check downcast assingment." << endl;
	expect(0);
	{
	    SP<foo> spf = new bar;
	    SP<bar> spb = new bar;
	    expect(2);
	    spb = spf;
	    expect(1);
	}
	expect(0);
    }
}

void x3()
{
    cout << "\n\n test x3 -- tricky copy/assign cases.\n";

    {
	cout << "\n Check assign to self." << endl;
	expect(0);
	{
	    SP<foo> sp1 = new foo;
	    expect(1);
	    sp1 = sp1;
	    expect(1);
	}
	expect(0);
    }

    {
	cout << "\n Check assign to self, T *." << endl;
	expect(0);
	{
	    SP<foo> sp1 = new foo;
	    SP<foo> sp2( sp1 );
	    expect(1);
	    sp1 = sp2;
	    expect(1);
	}
	expect(0);
    }

    {
	cout << "\n Check assign to self, X * (upcast)." << endl;
	expect(0);
	{
	    SP<foo> spf = new bar;
	    SP<bar> spb( spf );
	    spf = spb;
	    expect(1);
	}
	expect(0);
    }

    {
	cout << "\n Check assign to self, X * (downcast)." << endl;
	expect(0);
	{
	    SP<foo> spf = new bar;
	    SP<bar> spb( spf );
	    spb = spf;
	    expect(1);
	}
	expect(0);
    }
}

void x4()
{
    cout << "\n\n test x4 -- comparison operators, and suchlike.\n";

    {
	cout << "\n Check operator bool." << endl;
	expect(0);
	{
	    SP<foo> spf;
	    expect(0);
	    if (spf) cout << "bgous\n"; else cout << "good\n";
	    spf = new foo;
	    expect(1);
	    if (spf) cout << "good\n"; else cout << "bogus\n";
	}
	expect(0);
    }

    {
	cout << "\n Check operator!." << endl;
	expect(0);
	{
	    SP<foo> spf;
	    if (!spf) cout << "good\n"; else cout << "bogus\n";
	    spf = new foo;
	    if (!spf) cout << "bogus\n"; else cout << "good\n";
	}
	expect(0);
    }

    {
	cout << "\n Check for nullness." << endl;
	expect(0);
	{
	    SP<foo> spf;
// 	    if (spf == 0) cout << "good\n"; else cout << "bogus\n";
	}
	expect(0);
    }

    {
	cout << "\n Check comparison to dumb pointer." << endl;
	expect(0);
	{
	    SP<foo> spf;
	    foo *pf = new foo;
	    spf = pf;

	    if (spf == pf) cout << "good\n"; else cout << "bogus\n";
	    if (spf != pf) cout << "bogus\n"; else cout << "good\n";

	    if (pf == spf) cout << "good\n"; else cout << "bogus\n";
	    if (pf != spf) cout << "bogus\n"; else cout << "good\n";
	}
	expect(0);
    }

    {
	cout << "\n Check comparison to same type." << endl;
	expect(0);
	{
	    SP<foo> sp1 = new foo;
	    SP<foo> sp2 = new foo;
	    expect(2);
	    if (sp1 == sp2) cout << "bogus\n"; else cout << "good\n";
	    if (sp1 != sp2) cout << "good\n"; else cout << "bogus\n";
	}
	expect(0);
    }

    {
	cout << "\n Check SP comparison between compatible types." << endl;
	expect(0);
	{
	    SP<foo> spf = new bar;
	    SP<bar> spb = spf;

	    cout << "checking same object comparison." << endl;
	    if (spf == spb) cout << "good\n"; else cout << "bogus\n";
	    if (spf != spb) cout << "bogus\n"; else cout << "good\n";
	}
	expect(0);
	{
	    SP<foo> spf = new bar;
	    SP<bar> spb = new bar;
	    expect(2);

	    cout << "checking different object comparison." << endl;
	    if (spf == spb) cout << "bogus\n"; else cout << "good\n";
	    if (spf != spb) cout << "good\n"; else cout << "bogus\n";
	}
	expect(0);
    }
}

main()
{
#ifndef __KCC
    ios::sync_with_stdio();
#endif

    cout << "\n\n\n    tstSP starting.\n";

    x1();
    x2();
    x3();
    x4();
    
    tst_map();
}

//---------------------------------------------------------------------------//
// Test interaction of smart pointers with Maps.
//---------------------------------------------------------------------------//

void u1()
{
#if 0
    int i;
    expect(0);
    {
	Map< int, SP<foo> > mif;

	cout << "Putting some foo's in a Map.\n";

	for( i=0; i < 5; i++ ) {
	    SP<foo> spf = new foo(i);
	    mif[i] = spf;
	}

	expect(5);

	cout << "Now the Map will go out of scope.\n";
    }
    expect(0);
    {
	Map< int, SP<foo> > mif;

	cout << "Putting some foo's in a new Map.\n";

	for( i=0; i < 5; i++ ) {
	    SP<foo> spf = new foo(i);
	    mif[i] = spf;
	}

	expect(5);

	cout << "Walking a Map.\n";
	Map_iter<int,SP<foo> > mi( mif );
	for( ; mi; mi++ ) {
	    SP<foo> spf = mi.value();
	    cout << spf->val() << ' ';
	}
	cout << endl;

	{
	    cout << "Using copy ctor to make a new Map.\n";
	    Map<int,SP<foo> > mif2( mif );
	    expect(5);
	    cout << "Ready to let aliasing Map go out of scope.\n";
	}
	expect(5);

	SPList<foo> s;
	for( mi.reset(); mi; mi++ ) {
	    SP<foo> spf = mi.value();
	    s.insert( spf );
	}
	expect(5);

	cout << "Now the Map and the SPList will go out of scope.\n";
    }

    expect(0);
#endif
}

void tst_map()
{
    cout << "\n\n Testing SP<T> with Map<K,V>.\n\n";

    u1();
}

//---------------------------------------------------------------------------//
//                              end of tstSP.cc
//---------------------------------------------------------------------------//
