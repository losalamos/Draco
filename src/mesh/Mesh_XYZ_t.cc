//----------------------------------*-C++-*----------------------------------//
// Mesh_XYZ_t.cc
// Geoffrey M. Furnish
// Wed May 27 11:02:08 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "c4/global.hh"

template<class T>
void dump( const Mesh_XYZ::cell_array<T>& data, char *name )
{
    cout << "dumping a Mesh_XYZ::cell_array: " << name << endl;
    {
	HTSyncSpinLock h;
	char buf[80];
	for( int i=0; i < data.size(); i++ ) {
	    sprintf( buf, "node %d, cell %d, value=%lf \n",
		     C4::node(), i, data(i) );
	    cout << buf;
	}
    }
}

template<class T>
Mesh_XYZ::guarded_cell_array<T>& 
Mesh_XYZ::guarded_cell_array<T>::operator=( const Mesh_XYZ::cell_array<T>& c )
{
    for( int k=zoff; k < zoff + nczp; k++ )
        for( int j=0; j < ncy; j++ )
            for( int i=0; i < ncx; i++ )
                data(i,j,k) = c(i,j,k);

    update_guard_cells();
    
    return *this;
}

template<class T>
void Mesh_XYZ::guarded_cell_array<T>::update_guard_cells()
{
// Not implemented yet...
}

//---------------------------------------------------------------------------//
//                              end of Mesh_XYZ_t.cc
//---------------------------------------------------------------------------//
