//----------------------------------*-C++-*----------------------------------//
// nmtest.cc
// Geoffrey Furnish
// 25 January 1994
//---------------------------------------------------------------------------//
// @> A simple program to test the namelist facilities
//---------------------------------------------------------------------------//

#include "ds++/String.hh"
#include "Items.hh"
#include "Group.hh"

int a1, a2, a3, a4, a5;

int bool1, flog1;
float flt1;
double db1;
String s1;

int main( int argc, char *argv[] )
{
    NML_Group nmg( "test group" );

    nmg.addblock( "aaa",
		  NMI_INT, "a1", 1, &a1,
		  NMI_INT, "a2", 2, &a2,
		  NMI_INT, "a3", 3, &a3,
		  NMI_END );

    nmg.addblock( "aaa",
		  NMI_INT, "a4", 4, &a4,
		  NMI_INT, "a5", 5, &a5,
//		  NMI_INT, "a3", 3, &a3,     // test duplicate rejection.
		  NMI_END );

    nmg.addblock( "bbb",
		  N_float( flt1, 3.14 ),
		  N_double( db1, 49.00000000005 ),
		  N_bool( bool1, on ),
		  N_flog( flog1, off ),
		  N_String( s1, "tokamak" ),
		  NMI_END );

    nmg.writegroup( "x1.dat" );

    cout << "Done writing x1.dat\n" << flush;

    nmg.readgroup( "y.in" );

    nmg.writegroup( "x2.dat" );

    nmg.set_defaults();
    nmg.modify_block_defaults( "aaa",
			       "a1", 21,
			       NMI_MOD_END );

    nmg.writegroup( "x3.dat" );

    return 0;
}

extern "C" int Tcl_AppInit() { return 0; }

//---------------------------------------------------------------------------//
//                              end of nmtest.cc
//---------------------------------------------------------------------------//
