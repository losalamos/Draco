//---------------------------------*-C++-*---------------------------------//
// oldtest.cc
// Geoffrey Furnish
// 27 February 1992
//-------------------------------------------------------------------------//
// @> This is the oritinal nmtest program.  
//-------------------------------------------------------------------------//

// #ifndef OLD_GNU
#include <stdio.h>
// #else
// #if defined(LINUX) || defined(HPUX) || defined(UNICOS) || \
//     defined(AIX) || defined(SUNOS)
// #define ATT_io
// #endif
// #endif
// 
// #ifdef ATT_io
// #include <iostream.h>
// #include <fstream.h>
// #else
// #include <stream.h>
// #endif

#include "Items.hh"
#include "Group.hh"
//#include <namelist.h>

// block aaa

int a1;
float a2;

int a3;

#define MAX_STRING 40

//char a4[ MAX_STRING ];
String a4;

// block bbb

int b1;
float b2;
int b3;
int b4;
int b5;
float b6;
float b7;
float b8;
float b9;
int b10, b11, b12, b13, b14, b15, b16;
int b20, b21, b22, b23, b24, b25, b26;

// block ccc

int color;
#define RED 1
#define GREEN 2
#define BLUE 3
#define MAX_NSELECT 10
int nselect[ MAX_NSELECT ];
int n_nselect;
int fortran_logical;

NML_Group *mygroup;

void setup_namelist(void);
void print_program_values(void);
int main( int argc, char *argv[] );
void callback_a(void *d);
void callback_b(void *d);
void myquit(void *d);

int cba = 'a';
int cbb = 'b';
int cbq = 'q';

main( int argc, char *argv[] )
{
    cout << "\n Namelist demonstrator program.\n\n";

    setup_namelist();

    argc--; argv++;
    
    if (argc && !strcmp(*argv,"-tk") ) {
	cout << "Initiating Tk interface.\n\n";
	mygroup->popup_Tk(argc, argv);
    } else {
	cout << "argc = " << argc << '\n';
	cout << "*argv = " << *argv << '\n';
    }

    mygroup->writegroup( "x.out" );

    print_program_values();
    mygroup->readgroup( "x.in" );
    print_program_values();

    mygroup->writegroup( "y.out" );

    mygroup->set_defaults();
    print_program_values();
    mygroup->writegroup( "z.out" );
}

void setup_namelist()
{
    mygroup = new NML_Group( "My only group" );

    mygroup->addblock( "aaa",
			 NMI_INT, "A1", 1, &a1,
			 NMI_FLOAT, "a2", 2., &a2,
			 NMI_END );

    mygroup->addblock( "bbb",
			 NMI_INT, "b1", 11, &b1,
			 NMI_FLOAT, "b2", 12., &b2,
		         NMI_PAGE_BREAK,
			 NMI_INT, "b3", 13, &b3,
			 NMI_INT, "b4", 14, &b4,
			 NMI_INT, "b5", 15, &b5,
			 NMI_FLOAT, "b6", 16., &b6,
			 NMI_FLOAT, "b7", 17., &b7,
			 NMI_FLOAT, "b8", 18., &b8,
			 NMI_FLOAT, "b9", 19., &b9,
			 NMI_INT, "b10", 10, &b10,
			 NMI_INT, "b11", 10, &b11,
			 NMI_INT, "b12", 10, &b12,
			 NMI_INT, "b13", 10, &b13,
			 NMI_INT, "b14", 10, &b14,
			 NMI_INT, "b15", 10, &b15,
			 NMI_INT, "b16", 10, &b16,
			 NMI_INT, "b20", 10, &b20,
			 NMI_INT, "b21", 10, &b21,
			 NMI_INT, "b22", 10, &b22,
			 NMI_INT, "b23", 10, &b23,
			 NMI_INT, "b24", 10, &b24,
			 NMI_INT, "b25", 10, &b25,
			 NMI_INT, "b26", 10, &b26,
			 NMI_END );

    mygroup->addblock( "ccc",
		       NMI_INT_SET, "color", "red", 0, 3, &color,
				"red", RED,
				"green", GREEN,
				"blue", BLUE,
// 			 NMI_INT_SET, "color", 1, 3, &color,
// 				"red", RED,
// 				"green", GREEN,
// 				"blue", BLUE,
		       NMI_INT_ARRAY, "nselect", "1,2,3,4,7", nselect, 
		                MAX_NSELECT, &n_nselect,
		       NMI_FLOG, "fortran_logical", "on", &fortran_logical,
		       NMI_END );

    mygroup->addblock( "aaa",
		       NMI_INT, "a3", 3, &a3,
// 		       NMI_STRING, "a4", "default_string", a4, MAX_STRING,
 		       NMI_STRING, "a4", "default_string", &a4,
		       NMI_END );

    mygroup->add_callback( "Execute funciton A", callback_a, &cba );
    mygroup->add_callback( "Execute funciton B", callback_b, &cbb );
    mygroup->add_callback( "QUIT", myquit, &cbq );
    
    cout << "Namelist blocks have been defined.\n";
}

void print_program_values()
{
    cout << "\n\nThe values of the program variables are:\n";

    cout << "For block aaa, a1, a2 and a3\n";
    cout << a1 << " " << a2 << " " << a3 << '\n';;

    cout << "For block bbb, b1 to b9\n";
    cout << b1 << " "
	 << b2 << " "
	 << b3 << " "
	 << b4 << " "
	 << b5 << " "
	 << b6 << " "
	 << b7 << " "
	 << b8 << " "
	 << b9 << '\n';
    
    cout << "For block ccc, color\n";
    cout << color << '\n';;
    if (color==RED) cout << "color is RED!\n";
    if (color==GREEN) cout << "color is GREEN!\n";
    if (color==BLUE) cout << "color is BLUE!\n";

    cout << "ccc:nselect.  #= " << n_nselect << '\n';
    for( int i=0; i < n_nselect; i++ )
      cout << nselect[i] << " ";

    cout << "\n fortran_logical=" << fortran_logical << '\n';

    cout << '\n';
}
void callback_a( void *d )
{
    cout << "\n\nReached callback_a.  Do anything you like!\n";
    cout << "Data passed was :" << *(char *)d << ":.\n" << flush;
}

void callback_b( void *d )
{
    cout << "\n\nReached callback_b.  Jump in the lake.\n";
    cout << "Data passed was :" << *(char *)d << ":.\n" << flush;
}

void myquit( void *d )
{
    cout << "Reached user defined quit action.\n";
    cout << "Your favorite stuff here.\n";
    cout << "Data passed was :" << *(char *)d << ":.\n" << flush;
}
