//----------------------------------*-C++-*----------------------------------//
// tstSP.cc
// Geoffrey Furnish
// 3 December 1994
//---------------------------------------------------------------------------//
// @> Test program for the SP<T> class
//---------------------------------------------------------------------------//

#include <iostream>
#include <string>

using std::cout;
using std::endl;

#include "../SP.hh"

using dsxx::SP;

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

//---------------------------------------------------------------------------//
// Function to help us keep track of the foo pool.
//---------------------------------------------------------------------------//

int expect( int n ) {
    cout << "Expecting " << n << " foo's, ";
    if (nfoos == n) {
	cout << "and found them. --> test: passed \n";
	return 1;
    } else {
	cout << "but there are " << nfoos << ". --> test: failed \n";
	return 0;
    }
}

SP<foo> getafoo()
{
    SP<foo> p(new foo);
    return p;
}

SP<foo> useafoo( SP<foo> s )
{
    s->method();

    SP<foo> p(new foo);

    p = s;

    cout << "a foo should've just been destroyed.\n";

    return p;
}

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
	    SP<foo> spf(new foo);
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
	    SP<foo> sp1(new foo);
	    expect(1);
	    SP<foo> sp2(new foo);
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
	    SP<foo> sp2(new foo);
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
	    SP<foo> sp1(new foo);
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
	    SP<foo> spf(new bar);
	    expect(1);
	}
	expect(0);
    }

    {
	cout << "\n Check SP<X> copy ctor." << endl;
	expect(0);
	{
	    SP<foo> spf(new bar);
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
	    SP<bar> spb(new bar);
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
	    SP<foo> spf(new bar);
	    SP<bar> spb( spf );
	    expect(1);
	}
	expect(0);
    }

    {
	cout << "\n Check downcast assingment." << endl;
	expect(0);
	{
	    SP<foo> spf(new bar);
	    SP<bar> spb(new bar);
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
	    SP<foo> sp1(new foo);
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
	    SP<foo> sp1(new foo);
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
	    SP<foo> spf(new bar);
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
	    SP<foo> spf(new bar);
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
	    if (spf) cout << "test: failed\n"; else cout << "test: passed\n";
	    spf = new foo;
	    expect(1);
	    if (spf) cout << "test: passed\n"; else cout << "test: failed\n";
	}
	expect(0);
    }

    {
	cout << "\n Check operator!." << endl;
	expect(0);
	{
	    SP<foo> spf;
	    if (!spf) cout << "test: passed\n"; else cout << "test: failed\n";
	    spf = new foo;
	    if (!spf) cout << "test: failed\n"; else cout << "test: passed\n";
	}
	expect(0);
    }

    {
	cout << "\n Check for nullness." << endl;
	expect(0);
	{
	    SP<foo> spf;
// 	    if (spf == 0) cout << "test: passed\n"; else cout << "test: failed\n";
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

	    if (spf == pf) cout << "test: passed\n"; else cout << "test: failed\n";
	    if (spf != pf) cout << "test: failed\n"; else cout << "test: passed\n";

	    if (pf == spf) cout << "test: passed\n"; else cout << "test: failed\n";
	    if (pf != spf) cout << "test: failed\n"; else cout << "test: passed\n";
	}
	expect(0);
    }

    {
	cout << "\n Check comparison to same type." << endl;
	expect(0);
	{
	    SP<foo> sp1(new foo);
	    SP<foo> sp2(new foo);
	    expect(2);
	    if (sp1 == sp2) cout << "test: failed\n"; else cout << "test: passed\n";
	    if (sp1 != sp2) cout << "test: passed\n"; else cout << "test: failed\n";
	}
	expect(0);
    }

    {
	cout << "\n Check SP comparison between compatible types." << endl;
	expect(0);
	{
	    SP<foo> spf(new bar);
	    SP<bar> spb = spf;

	    cout << "checking same object comparison." << endl;
	    if (spf == spb) cout << "test: passed\n"; else cout << "test: failed\n";
	    if (spf != spb) cout << "test: failed\n"; else cout << "test: passed\n";
	}
	expect(0);
	{
	    SP<foo> spf(new bar);
	    SP<bar> spb(new bar);
	    expect(2);

	    cout << "checking different object comparison." << endl;
	    if (spf == spb) cout << "test: failed\n"; else cout << "test: passed\n";
	    if (spf != spb) cout << "test: passed\n"; else cout << "test: failed\n";
	}
	expect(0);
    }
}

void version(const std::string &progname)
{
    std::string version = "1.0.0";
    cout << progname << ": version " << version << endl;
}

int main(int argc, char *argv[])
{

    for (int arg=1; arg < argc; arg++)
	{
	    if (std::string(argv[arg]) == "--version")
		{
		    version(argv[0]);
		    return 0;
		}
	}

    cout << "\n\n\n    tstSP starting.\n";

    x1();
    x2();
    x3();
    x4();
    return 0;
}

//---------------------------------------------------------------------------//
//                              end of tstSP.cc
//---------------------------------------------------------------------------//
