//----------------------------------*-C++-*----------------------------------//
// dump.hh
// Maurice LeBrun
// Thu Feb  9 15:53:47 1995
//---------------------------------------------------------------------------//
// @> Functions for dumping various things in various ways.
//---------------------------------------------------------------------------//

#ifndef __base_dump_hh__
#define __base_dump_hh__

#include "c4/global.hh"
#include "c4/SpinLock.hh"

#define NODE0 if (!node)

#include "ds++/String.hh"
#include "ds++/Mat.hh"

template<class T>
void print_Mat( const dsxx::Mat1<T>& a, dsxx::String name );

template<class T>
void print_Mat( const dsxx::Mat2<T>& m, dsxx::String name );

template<class T>
void print2_Mat( const dsxx::Mat2<T>& m, dsxx::String name );

#include "util/ADFile.hh"

template<class T>
void ADdump_Mat( ADFile *f, dsxx::Mat3<T>& m, dsxx::String id );

template<class T>
void ADdump_Mat( ADFile *f, dsxx::Mat2<T>& m, dsxx::String id );

void dump( ADFile *f, char *buf, dsxx::String id );

//---------------------------------------------------------------------------//
// Used to dump an Mat1<T> which is understood to be one field extent
// wide on each processor.  For example, ak in phitor.
//---------------------------------------------------------------------------//

template<class T, class A>
void dump_fe( ADFile *f, dsxx::Mat1<T,A>& a, dsxx::String id )
{
    int nfu = a.nx();

    ADKey adk;
    int node = C4_node(), nodes = C4_nodes(), group = C4_group();
    int ncu = nfu * nodes;

    int i;
    dsxx::Mat2<T> mm( ncu, 1 );

    if (!node) {
	sprintf( adk.s, "Mat1f(%d) fe, %s", ncu, &id[0] );
	f->Start_new_entry(adk);
    }

    if (node == 0) {

	for( i=0; i < nfu; i++ )
	    mm(i,0) = a[i];

	for( int n=1; n < nodes; n++ )
	    C4_Recv( &mm(nfu*n,0), nfu*sizeof(T), n, 202, group );

	f->write( (void *) &mm(0,0), sizeof(T)*ncu );
	f->End_of_entry();
    }
    else
	C4_Send( &a[0], nfu*sizeof(T), 0, 202, group );
}

#endif                          // __base_dump_hh__

//---------------------------------------------------------------------------//
//                              end of base/dump.hh
//---------------------------------------------------------------------------//
