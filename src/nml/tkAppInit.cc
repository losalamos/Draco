/* 
 * tkAppInit.c --
 *
 *	Provides a default version of the Tcl_AppInit procedure for
 *	use in wish and similar Tk-based applications.
 *
 * Copyright (c) 1993 The Regents of the University of California.
 * All rights reserved.
 *
 * Permission is hereby granted, without written agreement and without
 * license or royalty fees, to use, copy, modify, and distribute this
 * software and its documentation for any purpose, provided that the
 * above copyright notice and the following two paragraphs appear in
 * all copies of this software.
 * 
 * IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY FOR
 * DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES ARISING OUT
 * OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE UNIVERSITY OF
 * CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * THE UNIVERSITY OF CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES,
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
 * AND FITNESS FOR A PARTICULAR PURPOSE.  THE SOFTWARE PROVIDED HEREUNDER IS
 * ON AN "AS IS" BASIS, AND THE UNIVERSITY OF CALIFORNIA HAS NO OBLIGATION TO
 * PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
 */

#ifndef lint
static char rcsid[] = "$Header$ SPRITE (Berkeley)";
#endif /* not lint */

#include "tk.h"
#include "itcl.h"
#include "itk.h"

#include "nml/nmitcl.hh"

#include "ds++/String.hh"
#include "nml/Items.hh"
#include "nml/Group.hh"

/*
 * The following variable is a special hack that allows applications
 * to be linked using the procedure "main" from the Tk library.  The
 * variable generates a reference to "main", which causes main to
 * be brought in from the library (and all of Tk and Tcl with it).
 */

// Not needed, since I ...  hmm, why don't I need this ???
// extern int main();
// int *tclDummyMainPtr = (int *) main;

#if TK_MAJOR_VERSION >= 4

int Tcl_AppInit( Tcl_Interp *interp );

main( int argc, char *argv[] )
{
    Tk_Main( argc, argv, Tcl_AppInit );
    return 0;
}

#endif

/*
 *----------------------------------------------------------------------
 *
 * Tcl_AppInit --
 *
 *	This procedure performs application-specific initialization.
 *	Most applications, especially those that incorporate additional
 *	packages, will have their own version of this procedure.
 *
 * Results:
 *	Returns a standard Tcl completion code, and leaves an error
 *	message in interp->result if an error occurs.
 *
 * Side effects:
 *	Depends on the startup script.
 *
 *----------------------------------------------------------------------
 */

extern "C" int Blt_Init( Tcl_Interp *interp );

int nmltest();
void s1_test_cb( const String& v, void *data );

int Tcl_AppInit( Tcl_Interp *interp )
{
    Tk_Window main;

    main = Tk_MainWindow(interp);

    /*
     * Call the init procedures for included packages.  Each call should
     * look like this:
     *
     * if (Mod_Init(interp) == TCL_ERROR) {
     *     return TCL_ERROR;
     * }
     *
     * where "Mod" is the name of the module.
     */

    if (Tcl_Init(interp) == TCL_ERROR) {
	return TCL_ERROR;
    }
    if (Tk_Init(interp) == TCL_ERROR) {
	return TCL_ERROR;
    }
    if (Itcl_Init(interp) == TCL_ERROR) {
	return TCL_ERROR;
    }
    if (Itk_Init(interp) == TCL_ERROR) {
	return TCL_ERROR;
    }
    if (Blt_Init(interp) == TCL_ERROR) {
	return TCL_ERROR;
    }
    if (nmltk_Init(interp) == TCL_ERROR) {
	return TCL_ERROR;
    }

    nmltest();

    /*
     * Call Tcl_CreateCommand for application-specific commands, if
     * they weren't already created by the init procedures called above.
     */

    /*
     * Specify a user-specific startup file to invoke if the application
     * is run interactively.  Typically the startup file is "~/.apprc"
     * where "app" is the name of the application.  If this line is deleted
     * then no user-specific startup file will be run under any conditions.
     */

    tcl_RcFileName = "~/.wishrc";
    return TCL_OK;
}

int a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14;

int bool1, flog1;
float flt1;
double db1;
String s1;
int shape;

int nmltest()
{
    NML_Group *nmg = new NML_Group( "test" );

    nmg->addblock( "aaa",
		   NMI_INT, "a1", 1, &a1,
		   NMI_INT, "a2", 2, &a2,
		   NMI_INT, "a3", 3, &a3,
		   NMI_INT, "a4", 3, &a4,
		   NMI_PAGE_BREAK, "second screen",
		   NMI_INT, "a5", 3, &a5,
		   NMI_INT, "a6", 3, &a6,
		   NMI_INT, "a7", 3, &a7,
		   NMI_INT, "a8", 3, &a8,
		   NMI_INT, "a9", 3, &a9,
		   NMI_PAGE_BREAK, "third screen",
		   NMI_INT, "a10", 3, &a10,
		   NMI_INT, "a11", 3, &a11,
		   NMI_INT, "a12", 3, &a12,
		   NMI_INT, "a13", 3, &a13,
		   NMI_INT, "a14", 3, &a14,
		   NMI_END );

//     nmg->addblock( "aaa",
// 		   NMI_INT, "a4", 4, &a4,
// 		   NMI_INT, "a5", 5, &a5,
// //		  NMI_INT, "a3", 3, &a3,     // test duplicate rejection.
// 		   NMI_END );

    nmg->addblock( "bbb",
		   NMI_PAGE_BREAK, "A one-page block",
		   N_float( flt1, 3.14 ),
		   N_double( db1, 49.00000000005 ),
		   N_bool( bool1, on ),
		   N_flog( flog1, off ),
		   N_String( s1, "tokamak" ),
		   NMI_INT_SET, "shape", "circle", 0, 4, &shape,
		   "circle", 0,
		   "square", 1,
		   "triangle", 2,
		   "diamond", 3,
		   NMI_END );

// Now let's install a callback for an item, switch to defered
// callback, and install a binding which will invoke it on demand.

    NML_Item *ps1 = nmg->Get_item( "bbb", "s1" );
    ps1->defered_callback = 1;
    ps1->callback_binding = "<Return>";

    NML_Callback cb( (void *) 7, s1_test_cb );
    nmg->add_callback( "bbb", "s1", cb );

    return 1;
}

void s1_test_cb( const String& v, void *data )
{
    cout << "s1_test_cb: callback activated\n";
    cout << "v = " << v << ", data = " << (int) data 
	 << endl << flush;
}
