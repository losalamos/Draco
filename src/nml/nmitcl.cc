//----------------------------------*-C++-*----------------------------------//
// nmitcl.cc
// Geoffrey Furnish
// 18 April 1994
//---------------------------------------------------------------------------//
// @> Extended wish interface for the namelist utility.  Uses itcl.
//---------------------------------------------------------------------------//

#include "tk.h"
#include "itcl.h"
#include "nml/nmitcl.hh"
#include "nml/Group.hh"
#include "nml/TraceInfo.hh"

#include "Assert.hh"

extern map<String,NML_Group *> NML_grplist;

#include <stdio.h>
#include <string.h>
#include <ctype.h>

// Evals the specified command, aborting on an error.

static int tcl_cmd( Tcl_Interp *interp, char *cmd );

// Evals the specified string, returning the result.

static int tcl_eval( Tcl_Interp *interp, char *cmd );

int nml_testCmd( ClientData cd, Tcl_Interp *interp,
		 int argc, char **argv );

//---------------------------------------------------------------------------//
// nmltk_Init()

// This is the Physcial Dynamics Namelist WISH extension
// initialization routine.  Nothing particularly spectacular happens
// here.  We set up a couple of extension commands and try to patch
// the auto_path to be able to find all our cool [incr Tcl] classes.
//---------------------------------------------------------------------------//

int nmltk_Init( Tcl_Interp *interp )
{
// Now add in the namelist api to the tcl interpreter.

    Tk_Window main;

    main = Tk_MainWindow(interp);

// Set up access to the Physical Dynamics C++ Namelist Tcl library.

    char *libdir = Tcl_GetVar2( interp, "env", "NML_LIBRARY",
				TCL_GLOBAL_ONLY);
    if (libdir == NULL)
	libdir = NML_LIBRARY;

    Tcl_SetVar( interp, "nml_library", libdir, TCL_GLOBAL_ONLY );

// Now add our extension commands.  Need to "itclize" these someday.

    Tcl_CreateCommand( interp, "nml_get_blocks", nml_get_blocks,
		       (ClientData) main, (void (*)(ClientData)) NULL);

//     Tcl_CreateCommand( interp, "nml_get_items", nml_get_items,
// 		       (ClientData) main, (void (*)(ClientData)) NULL);

// New extension facilities for [incr Tcl] 2.0.

    if ( Itcl_RegisterC( interp, "nml_test", nml_testCmd) != TCL_OK ) {
	return TCL_ERROR;
    }

    if ( Itcl_RegisterC( interp, "nml:Block:get_items",
			 nml_get_items ) != TCL_OK ) {
	return TCL_ERROR;
    }

// Now run the NML library Tcl initiailzation code.

    static char initCmd[] =
        "if [file exists $nml_library/init.nml] {\n\
	     source $nml_library/init.nml\n\
         } else {\n\
             set msg \"can't find $nml_library/init.nml\\n\"\n\
             append msg \"Perhaps your GTS build environment is not \\n\"\n\
             append msg \"configured correctly or you need to \\n\"\n\
             append msg \"set your NML_LIBRARY environment variable?\"\n\
             error $msg\n\
         }";

    return Tcl_Eval( interp, initCmd );
}

int nml_testCmd( ClientData cd, Tcl_Interp *interp,
		 int argc, char **argv )
{
    Tcl_SetVar2( interp, "itmdata", "foo", "wow", 0 );

    Tcl_SetResult( interp, "Made it into the extension command", TCL_STATIC );
    return TCL_OK;    
}

//---------------------------------------------------------------------------//
// nml_get_blocks()

// Return a list of block names for the indicated group.  This is
// needed in order for the Group class to set up the list of Blocks.
//---------------------------------------------------------------------------//

int nml_get_blocks( ClientData cd, Tcl_Interp *interp,
		    int argc, char **argv )
{
    if (argc != 2) {
	Tcl_SetResult( interp, "wrong # of args.", TCL_STATIC );
	return TCL_ERROR;
    }

    NML_Group *pg = NML_grplist[ argv[1] ];
    if ( pg == 0 ) {
	Tcl_AppendResult( interp, "No such group :", *argv[1],
			  ":.", 0 );
	return TCL_ERROR;
    }

    list<NML_Block *> pblist = pg->Blocklist();
    for( list<NML_Block *>::iterator pbi = pblist.begin();
	 pbi != pblist.end(); pbi++ ) {
	String name = (*pbi)->Name();
	Tcl_AppendElement( interp, &name[0] );
    }

    return TCL_OK;
}

//---------------------------------------------------------------------------//
// nml_get_items()

// This is one of the main portions of the extension facility.  Here
// we locate all the information which the Tcl side is going to need,
// and stuff it into a big associative array which the various Tcl
// sections can examine and utilize.  We also set up traces on those
// members of the array which hold the Tcl representation of the item
// values.  This allows us to perform various actions such as updating
// the compiled side, and performing cleanup when the Block window
// goes away.
//
// The format of the information communicated back to the Tcl side is:
//        ...  <fill in later>  ...
//---------------------------------------------------------------------------//

// @symbol = nml:Block:get_items.  This is now a method in Block.itk
// Specifically, it implements the method get_items in nml::Block

int nml_get_items( ClientData cd, Tcl_Interp *interp,
		   int argc, char **argv )
{
    if (argc != 1) {
	Tcl_SetResult( interp, "wrong # of args.", TCL_STATIC );
	return TCL_ERROR;
    }
    for( int i=0; i < argc; i++ )
	cout << "arg#" << i << " = " << argv[i] << endl << flush;

    char *group = Tcl_GetVar( interp, "group", 0 );
    char *block = Tcl_GetVar( interp, "name", 0 );

    cout << "Want block " << block << " in group " << group << endl << flush;

    NML_Group *pg = NML_grplist[ group ];
    if ( pg == 0 ) {
	Tcl_AppendResult( interp, "No such group :", *argv[1],
			  ":.", 0 );
	return TCL_ERROR;
    }

    int block_found = 0;

    NML_Block *pb;
    list<NML_Block *> pblist = pg->Blocklist();
    for( list<NML_Block *>::iterator pbi = pblist.begin();
	 pbi != pblist.end(); pbi++ ) {
	pb = *pbi;
	if ( pb->Name() == block ) {
	    block_found = 1;
	    break;
	}
    }
    if (!block_found) {
	Tcl_AppendResult( interp, "no such block :", *argv[2], ":.", 0 );
	return TCL_ERROR;
    }

    String allnames;

// Keep track of pages and number of items processed.  That way if
// they put a PageBreak first, we know it's only for the purpose of
// supplying a description tag, and don't make a page with no entries.

    int page=1, processed=0;

// Okay, thanks to the braindead tracing in itcl 1.3, we have to use a
// "unique global array" for storing the info.  Will base this on the
// $this for the active object.

    String This = Tcl_GetVar( interp, "oldthis", 0 );

    String items = String( &This[1] ) + "_items";

//     String items = "itmdata";
    cout << "array base name >" << items << "<.\n";

    char buf[ 200 ];
    sprintf( buf, "global ::%s\n", &items[0] );
    cout << "buf=" << buf << flush;
    if ( Tcl_Eval( interp, buf ) != TCL_OK ) {
	interp->result = "Not able to set global access.";
	return TCL_ERROR;
    }

// Set a default description for the first page.

    if (Tcl_SetVar2( interp, &items[0], "pagedescrip1", "", TCL_GLOBAL_ONLY )
	== NULL) {
	interp->result = "Unable to set global array element.";
	return TCL_ERROR;
    }

// Now walk the list.

    list<NML_Item *> itmlist = pb->Itemlist();
    for( list<NML_Item *>::iterator ili = itmlist.begin();
	 ili != itmlist.end(); ili++ ) {
	NML_Item *pi = *ili;
	String name = pi->Name();
	String type = pi->Type();
	String value = pi->Value();
	String def_val = pi->Default();
	String choices = pi->Choices();
	String widget = pi->Widget();
	String active = pi->Active() ? "yes" : "no";
	String bind = pi->callback_binding;

// Some things require special handling...

	if (type == "bool" || type == "flog") {
	    // For these classes, force to "on" or "off", since the
	    // checkbutton widget is kind of limited.

	    if (value == "n" || value == "no" || value == "f" ||
		value == ".f" || value == "false" ||
		value == ".false")
		value = "off";
	    if (value == "y" || value == "yes" || value == "t" ||
		value == ".t" || value == "true" ||
		value == ".true")
		value = "on";
	}

// Now handle the case of page breaks.

	if (name == "PageBreak") {
	    char buf[40];
  	    if (!processed) {
// No items processed yet, so this is just to tag the first page.

		sprintf( buf, "pagedescrip%d", page );

		Tcl_SetVar2( interp, &items[0], buf, &value[0],
			     TCL_GLOBAL_ONLY );
		page++;
	    } else {
// Start a new page and tag it.
		page++;
		sprintf( buf, "pagedescrip%d", page );
		
		Tcl_SetVar2( interp, &items[0], buf, &value[0],
			     TCL_GLOBAL_ONLY );
		allnames += name + ' ';
	    }
	    continue;
	}

	allnames += name + ' ';

	for( char *p = &name[0]; *p; p++ )
	    *p = tolower( *p );

	String kval = name + "_value";
	Tcl_SetVar2( interp, &items[0], &kval[0], &value[0], TCL_GLOBAL_ONLY );

	String ktype = name + "_type";
	Tcl_SetVar2( interp, &items[0], &ktype[0], &type[0], TCL_GLOBAL_ONLY );

	String kdv = name + "_default";
	Tcl_SetVar2( interp, &items[0], &kdv[0], &def_val[0],
		     TCL_GLOBAL_ONLY );

	String kchoice = name + "_choices";
	Tcl_SetVar2( interp, &items[0], &kchoice[0], &choices[0],
		     TCL_GLOBAL_ONLY );

	String kwid = name + "_widget";
	Tcl_SetVar2( interp, &items[0], &kwid[0], &widget[0],
		     TCL_GLOBAL_ONLY );

	String kactive = name + "_active";
	Tcl_SetVar2( interp, &items[0], &kactive[0], &active[0],
		     TCL_GLOBAL_ONLY );

// Now let's set a trace on writes to the value.

	NML_TraceInfo *pti = new NML_TraceInfo( pb, pi );

	if ( Tcl_TraceVar2( interp, &items[0], &kval[0],
			    TCL_TRACE_WRITES | TCL_TRACE_UNSETS |
			    TCL_GLOBAL_ONLY,
			    nml_catch_write,
			    (ClientData) pti ) != TCL_OK ) {
	    interp->result = "Can't trace a var in itcl.cc";
	    return TCL_ERROR;
	}

// If callback_binding is set for this item, then set bind string so
// Tcl can see it, and put a trace on writes to the variable which
// will be written by the bind.  That way we can trap on the write,
// and effect the callback.  The issue here is that we want to have
// some ClientData associated with the binding, and that's not
// possible from Tcl.  So we have to install this handler.

	if (pi->callback_binding != "") {
	    NML_TraceInfo *pti = new NML_TraceInfo( pb, pi );

	    String kbind = name + "_bind";
	    Tcl_SetVar2( interp, &items[0], &kbind[0], &bind[0], 0 );
	    String kbindx = kbind + 'x';
	    if ( Tcl_TraceVar2( interp, &items[0], &kbindx[0],
				TCL_TRACE_WRITES | TCL_TRACE_UNSETS |
				TCL_GLOBAL_ONLY,
				nml_catch_write,
				(ClientData) pti ) != TCL_OK ) {
		interp->result = "Can't trace a var in itcl.cc";
		return TCL_ERROR;
	    }
	}

	processed++;
    }

    Tcl_SetVar2( interp, &items[0], "items", &allnames[0], TCL_GLOBAL_ONLY );

    return TCL_OK;
}

//---------------------------------------------------------------------------//
// nml_catch_writes()

// Catch writes to a variable from the GUI side.  Should propagate
// these back to the compiled side.  But that doesn't seem to be being
// done just now.  Have to look into it.
//---------------------------------------------------------------------------//

char *nml_catch_write( ClientData cd, Tcl_Interp *interp,
		       char *name1, char *name2, int flags )
{
    NML_TraceInfo *pti = (NML_TraceInfo *) cd;

    NML_Block *pb = pti->Get_pb();
    NML_Item  *pi = pti->Get_pi();

// Now set the new value of the variable on the compiled side, make
// sure callbacks are run, etc.

    if (flags & TCL_TRACE_UNSETS) {
	cout << name1 << '(' << name2 << ") has been unset.\n";

	delete pti;
    }
    else if (flags & TCL_TRACE_WRITES) {

	cout << name1 << '(' << name2 << ") has been written.\n";

// Now lets figure out if this was a write to the value, or if it was
// a trap on the binding.

	int len2 = strlen(name2);
	char *s = name2 + len2 - 6;
	if (len2 > 6 && !strcmp(s,"_bindx")) {
	    // We've got a callback_binding trap

	    String nval = Tcl_GetVar2( interp, name1, name2, flags );
	    pi->do_callbacks();
	}
	else {
	    // We've got a normal value write trap.

	    flags &= TCL_GLOBAL_ONLY;
	    String nval = Tcl_GetVar2( interp, name1, name2, flags );

	    cout << "The new value is " << nval << endl;
	    pi->set_value( nval );
	    pi->update_value();
	    if (!pi->defered_callback)
		pi->do_callbacks();
	}
    }
    else
	throw( "Unrecognized trace directive." );

    return NULL;
}

// Copied from PLplot's tkshell.c, author mjl.

/*----------------------------------------------------------------------*\
* tcl_cmd
*
* Evals the specified command, aborting on an error.
\*----------------------------------------------------------------------*/

static int tcl_cmd( Tcl_Interp *interp, char *cmd )
{
#ifdef DEBUG_ENTER
    fprintf(stderr, "evaluating command %s\n", cmd);
#endif

    if (tcl_eval(interp, cmd)) {
	fprintf(stderr, "TCL command \"%s\" failed:\n\t %s\n",
		cmd, interp->result);
	return TCL_ERROR;
    }
    return TCL_OK;
}

/*----------------------------------------------------------------------*\
* tcl_eval
*
* Evals the specified string, returning the result.
* Use a static string buffer to hold the command, to ensure it's in
* writable memory (grrr...).
\*----------------------------------------------------------------------*/

static String cmdbuf;

static int tcl_eval( Tcl_Interp *interp, char *cmd )
{
    cmdbuf = cmd;

    return Tcl_VarEval(interp, &cmdbuf[0], (char **) NULL);
}

//---------------------------------------------------------------------------//
//                              end of nmitcl.cc
//---------------------------------------------------------------------------//
