//----------------------------------*-C++-*----------------------------------//
// nmitcl.hh
// Geoffrey Furnish
// 18 April 1994
//---------------------------------------------------------------------------//
// @> Extended wish interface for the namelist utility.  Uses itcl.
//---------------------------------------------------------------------------//

#ifndef __nml_nmitcl_hh__
#define __nml_nmitcl_hh__

int nmltk_Init( Tcl_Interp *interp );

// The rest of these comprise the extended wish api for namelist tk.

int nml_get_blocks( ClientData cd, Tcl_Interp *interp,
		    int argc, char **argv );
int nml_get_items( ClientData cd, Tcl_Interp *interp,
		   int argc, char **argv );

char *nml_catch_write( ClientData cd, Tcl_Interp *interp,
		       char *name1, char *name2, int flags );

#endif                          // __nml_nmitcl_hh__

//---------------------------------------------------------------------------//
//                              end of nmitcl.hh
//---------------------------------------------------------------------------//
