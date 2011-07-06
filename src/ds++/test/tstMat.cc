//----------------------------------*-C++-*----------------------------------//
/*! \file  tstMat.cc
 * \author Geoffrey Furnish
 * \date   Wed Apr  2 12:48:48 1997
 * \note   Copyright (c) 1997-2010 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "../Assert.hh"
#include "../Release.hh"
#include "../Allocators.hh"
#include "ds_test.hh"

#include <iostream>
#include <string>
#include <sstream>
#include <vector>

// This is heinous beyond description, but I /really/ want to test all the
// parts of the Matrix classes here.  Wantom abuse of CPP like that which
// follows, is just to prove that I mean business!

#define private public
#define protected public

#include "../Mat.hh"

using namespace rtt_dsxx;
using namespace std;


//---------------------------------------------------------------------------//

template<class T, class A>
void ikv1( Mat1<T,A>& /* x */ )
{
    PASSMSG("Successfully passed a reference to a Mat1<T,A> to the function ikv1.");
}

//---------------------------------------------------------------------------//

void t1()
{
    cout << "t1: beginning.\n";

    // Default construction
    {
        int newsize( 5 );
        
        // Test Mat1<int>;

	Mat1<int> x;
        PASSMSG("Successfully contructed a Mat1<int>.");
        
	x.redim( newsize-1 );
	x.redim( newsize );
        PASSMSG("Successfully resized a Mat1<int>.");

	ikv1( x );

        // test conformality function
        x.assert_conformality( x );
        if( x.conformal( newsize ) )
            PASSMSG("Mat1<int> correctly reported size.");
        else
            FAILMSG("Mat1<int> incorrectly reported size.");

        int value(-1);
        x.elements( value );
        if( value == newsize )
            PASSMSG("Mat1<int> correctly reported size.");
        else
            FAILMSG("Mat1<int> incorrectly reported size.");

        if( x.get_xmin() == 0 )
            PASSMSG("Mat1.get_xmin() returned 0.");
        else
            FAILMSG("Mat1.get_xmin() did not return 0.");

        if( x.get_xlen() == newsize )
            PASSMSG("Mat1.get_xlen() returned 5.");
        else
            FAILMSG("Mat1.get_xlen() did not return 5.");
    }

    // Construction using Guarded Allocator
    //  This does not work because the Mat class has been designed to only use
    //  the simple allocator.  assert_conformality takes a Mat1<T> instead of
    //  Mat1<T,Allocator>. 
    // {
    //     int newsize( 5 );
        
    //     // Test Mat1<int>;

    //     Mat1< int, Guarded_Allocator<int> > x;
    //     PASSMSG("Successfully contructed a Mat1< int,Guarded_Allocator<int> >.");
        
    //     x.redim( newsize-1 );
    //     x.redim( newsize );
    //     PASSMSG("Successfully resized a Mat1< int,Guarded_Allocator<int> >.");

    //     ikv1( x );

    //     // test conformality function
    //     x.assert_conformality( x );
    //     if( x.conformal( newsize ) )
    //         PASSMSG("Mat1< int,Guarded_Allocator<int> > correctly reported size.");
    //     else
    //         FAILMSG("Mat1< int,Guarded_Allocator<int> > incorrectly reported size.");

    //     int value(-1);
    //     x.elements( value );
    //     if( value == newsize )
    //         PASSMSG("Mat1< int,Guarded_Allocator<int> > correctly reported size.");
    //     else
    //         FAILMSG("Mat1< int,Guarded_Allocator<int> > incorrectly reported size.");

    //     if( x.get_xmin() == 0 )
    //         PASSMSG("Mat1.get_xmin() returned 0.");
    //     else
    //         FAILMSG("Mat1.get_xmin() did not return 0.");

    //     if( x.get_xlen() == newsize )
    //         PASSMSG("Mat1.get_xlen() returned 5.");
    //     else
    //         FAILMSG("Mat1.get_xlen() did not return 5.");
    // }
    
    // Construction using Bounds object.
    {
        int const bmin(0);
        int const bmax(9);
        
        Mat1<int> x( Bounds(bmin,bmax) );
        PASSMSG("Successfully contructed a Mat1<int> using the Bounds constructor.");
        Mat1<int> const y( Bounds(bmin,bmax) );
        PASSMSG("Successfully contructed a const Mat1<int> using the Bounds constructor.");

        // Test paren operator
        if( x(0) == 0 )
            PASSMSG("Successful test of paren operator.");
        else
            FAILMSG("Unsuccessful test of paren operator.");
        
        if( y(0) == 0 )
            PASSMSG("Successful test of const paren operator.");
        else
            FAILMSG("Unsuccessful test of const paren operator.");

        // Test *= operator
        x*=5;
        bool mybool( true );
        for( int i=bmin; i<bmax; ++i )
            if( x(i) != 0 )
            {
                mybool = false; break;
            }
        if( mybool )
            PASSMSG("Successful test of *= operator (scalar).");
        else
            FAILMSG("Unsuccessful test of *= operator (scalar).");

        x*=y;
        mybool = true;
        for( int i=bmin; i<bmax; ++i )
            if( x(i) != 0 )
            {
                mybool = false; break;
            }
        if( mybool )
            PASSMSG("Successful test of *= operator (Mat1).");
        else
            FAILMSG("Unsuccessful test of *= operator (Mat1).");

        // Test redim(Bounds)
        int const newbmax( 5 );
        x.redim( Bounds( bmin, newbmax ) );
        if (x.get_xlen() == newbmax-bmin+1)
            PASSMSG("Successful test redim(Bounds).");
        else
            FAILMSG("Unsuccessful test redim(Bounds).");
        
        // What does this do?
        // Mat1<int> const z( Bounds(0,1) );
        // x*=z;
    }

    // Construction from C array and other constructor corner cases
    {
        int raw_x[3] = {0, 1, 2};
        Mat1<int> x(raw_x, 3);
        if (x[0]!=0)
            FAILMSG("construction from C array FAILS");

        x.redim(2);
        if (x.size()!=2)
            FAILMSG("redim from C array FAILS");

        Mat1<int> y;
        y.redim(Bounds(1,2));
        if (y.size()!=2)
            FAILMSG("redim from C array FAILS");
        
        Mat1<int> z(raw_x, 3);
        z.redim(Bounds(2,2));
        if (z.size()!=1)
            FAILMSG("redim from C array FAILS");
    }

    // Assignment
    {
        int raw_x[3] = {0, 1, 2};
        Mat1<int> x(raw_x, 3);
        x = x;
        if (x[0]!=0)
            FAILMSG("assignment from C array FAILS");

        if (x.conformal(4))
            FAILMSG("conformal FAILS");

        x = Mat1<int>(Bounds(1,2));
        if (x.size() != 2)
            FAILMSG("assignment from variable bounds FAILS");

        if (x.conformal(4))
            FAILMSG("conformal FAILS");
    }

    // Test equivalence
    {
        int raw_x[3] = {0, 1, 2};
        Mat1<int> x(raw_x, 3);
        Mat1<int> y(2);
        if (y==x)
            FAILMSG("equivalence FAILS");

        y = x;
        y[0] = 1;
        if (y==x)
            FAILMSG("equivalence FAILS");        
    }

    // Test += operator
    {
        int const bmax( 5 );
        Mat1<int> x( bmax );
        Mat1<int> y( bmax );

        y+=1;
        bool mybool( true );
        for( int i=0; i<bmax; ++i )
            if( y(i) != 1 )
            {
                mybool = false; break;
            }
        if( mybool )
            PASSMSG("Successful test of += operator (scalar).");
        else
            FAILMSG("Unsuccessful test of += operator (scalar).");

        x+=y;
        mybool = true;
        for( int i=0; i<bmax; ++i )
            if( x(i) != 1 )
            {
                mybool = false; break;
            }
        if( mybool )
            PASSMSG("Successful test of += operator (Mat1).");
        else
            FAILMSG("Unsuccessful test of += operator (Mat1).");
    }

    // Test -= operator
    {
        int const bmax( 5 );
        Mat1<int> x( bmax );
        Mat1<int> y( bmax );

        y-=-1;
        bool mybool( true );
        for( int i=0; i<bmax; ++i )
            if( y(i) != 1 )
            {
                mybool = false; break;
            }
        if( mybool )
            PASSMSG("Successful test of -= operator (scalar).");
        else
            FAILMSG("Unsuccessful test of -= operator (scalar).");

        x-=y;
        mybool = true;
        for( int i=0; i<bmax; ++i )
            if( x(i) != -1 )
            {
                mybool = false; break;
            }
        if( mybool )
            PASSMSG("Successful test of -= operator (Mat1).");
        else
            FAILMSG("Unsuccessful test of -= operator (Mat1).");
    }

    // Test /= operator
    {
        int const bmax( 5 );
        Mat1<int> x( bmax );
        Mat1<int> y( bmax );

        y+=10;
        y/=2;
        bool mybool( true );
        for( int i=0; i<bmax; ++i )
            if( y(i) != 5 )
            {
                mybool = false; break;
            }
        if( mybool )
            PASSMSG("Successful test of /= operator (scalar).");
        else
            FAILMSG("Unsuccessful test of /= operator (scalar).");

        x=5;
        x/=y;
        mybool = true;
        for( int i=0; i<bmax; ++i )
            if( x(i) != 1 )
            {
                mybool = false; break;
            }
        if( mybool )
            PASSMSG("Successful test of /= operator (Mat1).");
        else
            FAILMSG("Unsuccessful test of /= operator (Mat1).");
    }

    PASSMSG("Done with test t1.");
}

//---------------------------------------------------------------------------//
// Test the Mat2, every way we can think of.
//---------------------------------------------------------------------------//

void t2()
{
    cout << "t2: beginning.\n";
    
    {
	Mat2<int> x;
        PASSMSG("Successfully contructed a Mat2<int> using default contructor.");
    }
    
    {
        int const bmax(3);
           
        // Check various fundamental computations.
	Mat2<int> x(bmax,bmax);
        
	Assert( x.nx() == bmax );
	Assert( x.ny() == bmax );
	Assert( x.index(0,0) == 0 );
	Assert( x.size() == bmax*bmax );

        if (x.conformal(bmax+1, bmax)) FAILMSG("Mat2<int> conformal check FAILS");
        if (x.conformal(bmax, bmax+2)) FAILMSG("Mat2<int> conformal check FAILS");

	int k=0;
	for( int j=0; j < bmax; j++ )
	    for( int i=0; i < bmax; i++ )
            {
		x(i,j) = k++;
                cout << x(i,j) << endl;
            }

	k = 0;
	for( Mat2<int>::iterator xi = x.begin(); xi != x.end(); )
	    if (*xi++ != k++)
		throw "bogus";

        ostringstream msg;
        msg << "Successfully contructed a Mat2<int> using "
            << "specified size contstructor.";
        PASSMSG(msg.str());

        // Test conformality
        
        x.assert_conformality( x );
        if( x.conformal( bmax, bmax ) )
            PASSMSG("Mat2<int> correctly reported size.");
        else
            FAILMSG("Mat2<int> incorrectly reported size.");

        // test elements, get_xlen and get_xmin member functions.
        vector<int> velem;
        velem.push_back(0);
        velem.push_back(0);
        
        x.elements( velem[0], velem[1] );
        if( velem[0] == bmax )
            PASSMSG("Mat2<int> correctly reported size for 1st dim.");
        else
            FAILMSG("Mat2<int> incorrectly reported size for 1st dim.");
        if( velem[1] == bmax )
            PASSMSG("Mat2<int> correctly reported size for 2nd dim.");
        else
            FAILMSG("Mat2<int> incorrectly reported size for 2nd dim.");
        
        if( x.get_xmin() == 0 )
            PASSMSG("Mat2.get_xmin() returned 0.");
        else
            FAILMSG("Mat2.get_xmin() did not return 0.");
        if( x.get_ymin() == 0 )
            PASSMSG("Mat2.get_ymin() returned 0.");
        else
            FAILMSG("Mat2.get_ymin() did not return 0.");
        if( x.get_xlen() == bmax )
            PASSMSG("Mat2.get_xlen() returned 3.");
        else
            FAILMSG("Mat2.get_xlen() did not return 3.");
        if( x.get_ylen() == bmax )
            PASSMSG("Mat2.get_ylen() returned 3.");
        else
            FAILMSG("Mat2.get_ylen() did not return 3.");        
    }

    // Test construction from Bounds object.
    {
        int const bmax(3);

        Mat2<int> x( Bounds(1,bmax), Bounds(1,bmax) );

	Assert( x.nx() == bmax );
	Assert( x.ny() == bmax );
        Assert( x.index(1,1) == 4 );
	Assert( x.size() == 9 );

	int k=0;
	for( int j=1; j < 4; j++ )
	    for( int i=1; i < 4; i++ )
		x(i,j) = k++;

	k = 0;
	for( Mat2<int>::iterator xi = x.begin(); xi != x.end(); )
	    if (*xi++ != k++)
		throw "bogus";

        PASSMSG("Successfully constructed a Mat2<int> using Bounds constructor.");
    }

    // Test construction from C array
    {
        int const bmax(4);
        
        Mat2<int> x;
        x.redim( Bounds(1,bmax), Bounds(1,bmax+1));
        if (x.size() != bmax*(bmax+1))
            FAILMSG("Mat2 redim from Bounds FAILS");

        int raw_y[16] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
        
        Mat2<int> y(raw_y, 4, 4);
        x.redim( Bounds(1,4), Bounds(1,4+1));
        if (x.size() != 4*(4+1))
            FAILMSG("Mat2 redim from Bounds FAILS");
        
        Mat2<int> z(raw_y, 4, 4);
        z.redim(4,5);
        if (z.size() != 4*5)
            FAILMSG("Mat2 redim from Bounds FAILS");
        
        Mat2<int> w;
        w.redim(4,5);
        if (w.size() != 4*5)
            FAILMSG("Mat2 redim from Bounds FAILS");

        w.redim(2,3);
        if (x==y)
            FAILMSG("Mat2 operator== FAILS");

        y(1,1) = 3;
        z(1,1) = 4;
        if (y==z)
            FAILMSG("Mat2 operator== FAILS");
    }

    // Test paren operator
    {
        int const bmin(0);
        int const bmax(3);
        Mat2<int> x( Bounds(bmin,bmax), Bounds(bmin,bmax) );
        Mat2<int> const y( Bounds(bmin,bmax), Bounds(bmin,bmax) );

        if( x(0,0) == y(0,0) )
            PASSMSG("Successfull test of paren operator.");
        else
            FAILMSG("Unsuccessfull test of paren operator.");
    }

    // Test = operator
    {
        int const bmin(0);
        int const bmax(3);
        int value(5);
        Mat2<int> x( bmax, bmax );

        x=value;
        bool mybool(true);
        for( int i=bmin; i< bmax; ++i )
            for( int j=bmin; j< bmax; ++j )
                if( x(i,j) != 5 )
                {
                    mybool = false;
                    break;
                }
        if( mybool )
            PASSMSG("Successfull test of = operator.");
        else
            FAILMSG("Unsuccessfull test of = operator.");
    }

    // Test *= operator
    {
        int const bmin(0);
        int const bmax(3);
        int value(1);
        Mat2<int> x( bmax, bmax );
        Mat2<int> y( bmax, bmax );

        x=value;
        x*=5;
        bool mybool(true);
        for( int i=bmin; i< bmax; ++i )
            for( int j=bmin; j< bmax; ++j )
                if( x(i,j) != 5 )
                {
                    mybool = false;
                    break;
                }
        if( mybool )
            PASSMSG("Successfull test of *= operator (scalar).");
        else
            FAILMSG("Unsuccessfull test of *= operator (scalar).");

        y=2;
        x*=y;
        mybool = true;
        for( int i=bmin; i< bmax; ++i )
            for( int j=bmin; j< bmax; ++j )
                if( x(i,j) != 10 )
                {
                    mybool = false;
                    break;
                }
        if( mybool )
            PASSMSG("Successfull test of *= operator (Mat2).");
        else
            FAILMSG("Unsuccessfull test of *= operator (Mat2).");
    }
    
    // Test += operator
    {
        int const bmin(0);
        int const bmax(3);
        int value(1);
        Mat2<int> x( bmax, bmax );
        Mat2<int> y( bmax, bmax );

        x=value;
        x+=5;
        bool mybool(true);
        for( int i=bmin; i< bmax; ++i )
            for( int j=bmin; j< bmax; ++j )
                if( x(i,j) != 6 )
                {
                    mybool = false;
                    break;
                }
        if( mybool )
            PASSMSG("Successfull test of += operator (scalar).");
        else
            FAILMSG("Unsuccessfull test of += operator (scalar).");

        y=2;
        x+=y;
        mybool = true;
        for( int i=bmin; i< bmax; ++i )
            for( int j=bmin; j< bmax; ++j )
                if( x(i,j) != 8 )
                {
                    mybool = false;
                    break;
                }
        if( mybool )
            PASSMSG("Successfull test of += operator (Mat2).");
        else
            FAILMSG("Unsuccessfull test of += operator (Mat2).");
    }

    // Test -= operator
    {
        int const bmin(0);
        int const bmax(3);
        int value(5);
        Mat2<int> x( bmax, bmax );
        Mat2<int> y( bmax, bmax );

        x=value;
        x-=4;
        bool mybool(true);
        for( int i=bmin; i< bmax; ++i )
            for( int j=bmin; j< bmax; ++j )
                if( x(i,j) != 1 )
                {
                    mybool = false;
                    break;
                }
        if( mybool )
            PASSMSG("Successfull test of -= operator (scalar).");
        else
            FAILMSG("Unsuccessfull test of -= operator (scalar).");

        y=1;
        x-=y;
        mybool = true;
        for( int i=bmin; i< bmax; ++i )
            for( int j=bmin; j< bmax; ++j )
                if( x(i,j) != 0 )
                {
                    mybool = false;
                    break;
                }
        if( mybool )
            PASSMSG("Successfull test of -= operator (Mat2).");
        else
            FAILMSG("Unsuccessfull test of -= operator (Mat2).");
    }

    // Test /= operator
    {
        int const bmin(0);
        int const bmax(3);
        int value(10);
        Mat2<int> x( bmax, bmax );
        Mat2<int> y( bmax, bmax );

        x=value;
        x/=2;
        bool mybool(true);
        for( int i=bmin; i< bmax; ++i )
            for( int j=bmin; j< bmax; ++j )
                if( x(i,j) != 5 )
                {
                    mybool = false;
                    break;
                }
        if( mybool )
            PASSMSG("Successfull test of /= operator (scalar).");
        else
            FAILMSG("Unsuccessfull test of /= operator (scalar).");

        y=5;
        x/=y;
        mybool = true;
        for( int i=bmin; i< bmax; ++i )
            for( int j=bmin; j< bmax; ++j )
                if( x(i,j) != 1 )
                {
                    mybool = false;
                    break;
                }
        if( mybool )
            PASSMSG("Successfull test of /= operator (Mat2).");
        else
            FAILMSG("Unsuccessfull test of /= operator (Mat2).");
    }

    // test redim operators
    {
        int const bmin(0);
        int const bmax(3);
        Mat2<int> x( bmax, bmax );

        int const b2min(0);
        int const b2max(4);
        x.redim(b2max,b2max);

        if( x.get_xmin() == b2min         &&
            x.get_xlen() == b2max-b2min   &&
            x.get_ymin() == b2min         &&
            x.get_ylen() == x.get_xlen() )
        {
            PASSMSG("Successfull test of redim operation.");
        }
        else
        {
            ostringstream msg;
            msg << "Unsuccessfull test of redim operation.\n"
                << "\t x.get_xmin() = " << x.get_xmin()
                << ", expected " << b2min << endl
                << "\t x.get_xlen() = " << x.get_xlen()
                << ", expected " << b2max-b2min << endl
                << "\t x.get_ymin() = " << x.get_ymin()
                << ", expected " << b2min << endl
                << "\t x.get_ylen() = " << x.get_ylen()
                << ", expected " << x.get_xlen() << endl;
            FAILMSG(msg.str());
        }

        int const b3min(-3);
        int const b3max(1);
        Bounds b(b3min,b3max);
        Mat2<int> y(Bounds(bmin,bmax),Bounds(bmin,bmax));
        y.redim(b,b);
        if( y.get_xmin() == b3min         &&
            y.get_xlen() == b3max-b3min+1 &&
            y.get_ymin() == b3min         &&
            y.get_ylen() == y.get_xlen() )
        {
            PASSMSG("Successfull test of redim operation (Bounds).");
        }
        else
        {
            ostringstream msg;
            msg << "Unsuccessfull test of redim operation (Bounds).\n"
                << "\t y.get_xmin() = " << y.get_xmin()
                << ", expected " << b3min << endl
                << "\t y.get_xlen() = " << y.get_xlen()
                << ", expected " << b3max-b3min+1 << endl
                << "\t y.get_ymin() = " << y.get_ymin()
                << ", expected " << b3min << endl
                << "\t y.get_ylen() = " << y.get_ylen()
                << ", expected " << y.get_xlen() << endl;
            FAILMSG(msg.str());
        }
    }
    
    PASSMSG("Successful completion of t2 tests.");
}


void version(const std::string &progname)
{
    std::string version = rtt_dsxx::release();
    cout << progname << ": version " << version << endl;
}

//---------------------------------------------------------------------------//

int main( int argc, char *argv[] )
{
    // version tag

    version( argv[0] );
    for (int arg=1; arg < argc; arg++)
        if (std::string(argv[arg]) == "--version")
            return 0;
    
    try {
	cout << "Initiating test of the Mat family.\n";
	t1();
	t2();
    }
    catch( assertion& a )
    {
	cout << "Test: Failed assertion: " << a.what() << endl;
    }
    catch( ... )
    {
        cout << "FATAL error occured in tstMat.\n";
    }
    
    // status of test
    cout << endl;
    cout <<     "*********************************************" << endl;
    if( rtt_ds_test::passed ) 
    {
        cout << "**** tstMat Test: PASSED" << endl;
    }
    cout <<     "*********************************************" << endl;
    cout << endl;

    cout << "Done testing Mat family.\n";

    return 0;
}

//---------------------------------------------------------------------------//
//                              end of tstMat.cc
//---------------------------------------------------------------------------//
