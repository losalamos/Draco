//----------------------------------*-C++-*----------------------------------//
// Mesh_XYZ.t.cc
// Geoffrey M. Furnish
// Wed May 27 11:02:08 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "c4/global.hh"
#include "c4/C4_Req.hh"
#include "ds++/Assert.hh"

template<class T>
void dump( const Mesh_XYZ::cctf<T>& data, char *name )
{
    cout << "dumping a Mesh_XYZ::cctf: " << name << endl;
    {
    //	HTSyncSpinLock h;
	char buf[80];
	for( int i=0; i < data.size(); i++ ) {
        //	    sprintf( buf, "node %d, cell %d, value=%lf \n",
        //		     C4::node(), i, data(i) );
	    sprintf( buf, "cell %d, value=%lf \n",
		     i, data(i) );
	    cout << buf;
	}
    }
}

template<class T>
void dump( const Mesh_XYZ::fcdtf<T>& data, char *name )
{
    cout << "dumping a Mesh_XYZ::fcdtf: " << name << endl;
    {
    //	HTSyncSpinLock h;
	char buf[80];
	for( int i=0; i < data.size()/6; i++ ) {
            for( int j=0; j < 6; j++ ) {
	        sprintf( buf, "cell %d, face %d, value=%lf \n",
		         i, j, data(i,j) );
                cout << buf;
            }
	}
    }
}

template<class T>
Mesh_XYZ::gcctf<T>& 
Mesh_XYZ::gcctf<T>::operator=( const Mesh_XYZ::cctf<T>& c )
{
    for( int k=zoff; k < zoff + nczp; k++ )
        for( int j=0; j < ncy; j++ )
            for( int i=0; i < ncx; i++ )
                data(i,j,k) = c(i,j,k);

    update_guard_cells();
    
    return *this;
}

template<class T>
void Mesh_XYZ::gcctf<T>::update_guard_cells()
{
    using namespace C4;
    C4_Req lrcv, rrcv;

    dsxx::Mat2<T> lrbuf( &data(0,0,zoff-1), ncx, ncy );
    dsxx::Mat2<T> rrbuf( &data(0,0,zoff+nczp), ncx, ncy );

    if (node > 0)
        RecvAsync( lrcv, &lrbuf(0,0), ncx*ncy, node-1 );

    if (node < lastnode)
        RecvAsync( rrcv, &rrbuf(0,0), ncx*ncy, node+1 );

    if (node > 0)
        Send( &data(0,0,zoff), ncx*ncy, node-1 );

    if (node < lastnode)
        Send( &data(0,0,zoff+nczp-1), ncx*ncy, node+1 );
}

template<class T>
Mesh_XYZ::gfcdtf<T>& 
Mesh_XYZ::gfcdtf<T>::operator=( const Mesh_XYZ::fcdtf<T>& c )
{
    for( int k=zoff; k < zoff + nczp; k++ )
        for( int j=0; j < ncy; j++ )
            for( int i=0; i < ncx; i++ )
                for ( int f=0; f < 6; f++ )
                    data(f,i,j,k) = c(i,j,k,f);

    update_gfcdtf();
    
    return *this;
}

template<class T>
void Mesh_XYZ::gfcdtf<T>::update_gfcdtf()
{
    using namespace C4;
    C4_Req lrcv, rrcv;

    dsxx::Mat3<T> lrbuf( &data(0,0,0,zoff-1), 6, ncx, ncy );
    dsxx::Mat3<T> rrbuf( &data(0,0,0,zoff+nczp), 6, ncx, ncy );

    if (node > 0)
        RecvAsync( lrcv, &lrbuf(0,0,0), 6*ncx*ncy, node-1 );

    if (node < lastnode)
        RecvAsync( rrcv, &rrbuf(0,0,0), 6*ncx*ncy, node+1 );

    if (node > 0)
        Send( &data(0,0,0,zoff), 6*ncx*ncy, node-1 );

    if (node < lastnode)
        Send( &data(0,0,0,zoff+nczp-1), 6*ncx*ncy, node+1 );
}

template<class T>
Mesh_XYZ::bstf<T>&
Mesh_XYZ::bstf<T>::operator=( T x )
{
    // left face
    for ( int j = 0; j < ncy; ++j )
        for ( int k = zoff; k < zoff + nczp; ++k )
            data(0,j,k,0) = x;

    // right face
    for ( int j = 0; j < ncy; ++j )
        for ( int k = zoff; k < zoff + nczp; ++k )
            data(ncx-1,j,k,1) = x;

    // front face
    for ( int i = 0; i < ncx; ++i )
        for ( int k = zoff; k < zoff + nczp; ++k )
            data(i,0,k,2) = x;

    // back face
    for ( int i = 0; i < ncx; ++i )
        for ( int k = zoff; k < zoff + nczp; ++k )
            data(i,ncy-1,k,3) = x;

    // bottom face
    if (node == 0)
        for ( int i = 0; i < ncx; ++i )
            for ( int j = 0; j < ncy; ++j )
                data(i,j,0,4) = x;

    // top face
    if (node == lastnode)
        for ( int i = 0; i < ncx; ++i )
            for ( int j = 0; j < ncy; ++j )
                data(i,j,ncz-1,5) = x;

    return *this;
}

template<class T>
T& Mesh_XYZ::bstf<T>::operator()( int i, int j, int k, int f )
{
    Assert (   (i == 0 && f == 0) || (i == ncx - 1 && f == 1)
            || (j == 0 && f == 2) || (j == ncy - 1 && f == 3)
            || (k == 0 && f == 4) || (k == ncz - 1 && f == 5));

    return data(i,j,k,f);
}

template<class T>
T Mesh_XYZ::bstf<T>::operator()( int i, int j, int k, int f ) const
{
    Assert (   (i == 0 && f == 0) || (i == ncx - 1 && f == 1)
            || (j == 0 && f == 2) || (j == ncy - 1 && f == 3)
            || (k == 0 && f == 4) || (k == ncz - 1 && f == 5));

    return data(i,j,k,f);
}

template <class T1, class T2, class Op>
void Mesh_XYZ::scatter
( Mesh_XYZ::fcdtf<T1>& to, const Mesh_XYZ::cctf<T2>& from, const Op& op )
{
    Mesh_XYZ::gcctf<T2> gfrom(from);
    for ( int i = 0; i < to.ncx; ++i )
      for ( int j = 0; j < to.ncy; ++j )
        for ( int k = to.zoff; k < to.zoff + to.nczp; ++k )
	{
          for ( int f = 0; f < 6; ++f )
	     op(to(i,j,k,f), gfrom(i,j,k));
          if (i != 0) op(to(i,j,k,0), gfrom(i-1,j,k));
          if (i != to.ncx - 1) op(to(i,j,k,1), gfrom(i+1,j,k));
          if (j != 0) op(to(i,j,k,2), gfrom(i,j-1,k));
          if (j != to.ncy - 1) op(to(i,j,k,3), gfrom(i,j+1,k));
          if (k != 0) op(to(i,j,k,4), gfrom(i,j,k-1));
          if (k != to.ncz - 1) op(to(i,j,k,5), gfrom(i,j,k+1));
        }
}

template <class T1, class T2, class Op>
void Mesh_XYZ::scatter
( Mesh_XYZ::cctf<T1>& to, const Mesh_XYZ::fcdtf<T2>& from, const Op& op )
{
    for ( int i = 0; i < to.ncx; ++i )
      for ( int j = 0; j < to.ncy; ++j )
        for ( int k = to.zoff; k < to.zoff + to.nczp; ++k )
          for ( int f = 0; f < 6; ++f )
	    op(to(i,j,k), from(i,j,k,f));
}

template <class T1, class T2, class Op> 
void Mesh_XYZ::scatter
( Mesh_XYZ::fcdtf<T1>& to, const Mesh_XYZ::vctf<T2>& from, const Op& op )
{
    for ( int i = 0; i < to.ncx; ++i )
      for ( int j = 0; j < to.ncy; ++j )
        for ( int k = to.zoff; k < to.zoff + to.nczp; ++k )
        {
          op(to(i,j,k,0), from(i,j,k,0));
          op(to(i,j,k,0), from(i,j,k,2));
          op(to(i,j,k,0), from(i,j,k,4));
          op(to(i,j,k,0), from(i,j,k,6));
          op(to(i,j,k,1), from(i,j,k,1));
          op(to(i,j,k,1), from(i,j,k,3));
          op(to(i,j,k,1), from(i,j,k,5));
          op(to(i,j,k,1), from(i,j,k,7));
          op(to(i,j,k,2), from(i,j,k,0));
          op(to(i,j,k,2), from(i,j,k,1));
          op(to(i,j,k,2), from(i,j,k,4));
          op(to(i,j,k,2), from(i,j,k,5));
          op(to(i,j,k,3), from(i,j,k,2));
          op(to(i,j,k,3), from(i,j,k,3));
          op(to(i,j,k,3), from(i,j,k,6));
          op(to(i,j,k,3), from(i,j,k,7));
          op(to(i,j,k,4), from(i,j,k,0));
          op(to(i,j,k,4), from(i,j,k,1));
          op(to(i,j,k,4), from(i,j,k,2));
          op(to(i,j,k,4), from(i,j,k,3));
          op(to(i,j,k,5), from(i,j,k,4));
          op(to(i,j,k,5), from(i,j,k,5));
          op(to(i,j,k,5), from(i,j,k,6));
          op(to(i,j,k,5), from(i,j,k,7));
        }  
}

template <class T1, class T2, class Op>
void Mesh_XYZ::gather
( Mesh_XYZ::fcdtf<T1>& to, const Mesh_XYZ::cctf<T2>& from, const Op& op )
{
    for ( int i = 0; i < to.ncx; ++i )
      for ( int j = 0; j < to.ncy; ++j )
        for ( int k = to.zoff; k < to.zoff + to.nczp; ++k )
          for ( int f = 0; f < 6; ++f )
	    op(to(i,j,k,f), from(i,j,k));
}

template <class T>
void Mesh_XYZ::swap
( Mesh_XYZ::fcdtf<T>& to, const Mesh_XYZ::fcdtf<T>& from )
{
    Mesh_XYZ::gfcdtf<T> gfrom(from);
    for ( int i = 0; i < to.ncx; ++i )
      for ( int j = 0; j < to.ncy; ++j )
        for ( int k = to.zoff; k < to.zoff + to.nczp; ++k )
	{
          if (i != 0) to(i,j,k,0) = gfrom(i-1,j,k,1);
          else to(i,j,k,0) = 0;
          if (i != to.ncx - 1) to(i,j,k,1) = gfrom(i+1,j,k,0);
          else to(i,j,k,1) = 0;
          if (j != 0) to(i,j,k,2) = gfrom(i,j-1,k,3);
          else to(i,j,k,2) = 0;
          if (j != to.ncy - 1) to(i,j,k,3) = gfrom(i,j+1,k,2);
          else to(i,j,k,3) = 0;
          if (k != 0) to(i,j,k,4) = gfrom(i,j,k-1,5);
          else to(i,j,k,4) = 0;
          if (k != to.ncz - 1) to(i,j,k,5) = gfrom(i,j,k+1,4);
          else to(i,j,k,5) = 0;
	}
}

template <class T1, class T2, class Op>
void Mesh_XYZ::strip
( Mesh_XYZ::bstf<T1>& to, const Mesh_XYZ::fcdtf<T2>& from, const Op& op )
{
    // left face
    for ( int j = 0; j < to.ncy; ++j )
        for ( int k = to.zoff; k < to.zoff + to.nczp; ++k )
            op(to(0,j,k,0), from(0,j,k,0));

    // right face
    for ( int j = 0; j < to.ncy; ++j )
        for ( int k = to.zoff; k < to.zoff + to.nczp; ++k )
            op(to(to.ncx-1,j,k,1), from(to.ncx-1,j,k,1));

    // front face
    for ( int i = 0; i < to.ncx; ++i )
        for ( int k = to.zoff; k < to.zoff + to.nczp; ++k )
            op(to(i,0,k,2), from(i,0,k,2));

    // back face
    for ( int i = 0; i < to.ncx; ++i )
        for ( int k = to.zoff; k < to.zoff + to.nczp; ++k )
            op(to(i,to.ncy-1,k,3), from(i,to.ncy-1,k,3));

    // bottom face
    if (to.node == 0)
        for ( int i = 0; i < to.ncx; ++i )
            for ( int j = 0; j < to.ncy; ++j )
                op(to(i,j,0,4), from(i,j,0,4));

    // top face
    if (to.node == to.lastnode)
        for ( int i = 0; i < to.ncx; ++i )
            for ( int j = 0; j < to.ncy; ++j )
                op(to(i,j,to.ncz-1,5), from(i,j,to.ncz-1,5));
}

template <class T1, class T2, class Op>
void Mesh_XYZ::coat
( Mesh_XYZ::fcdtf<T1>& to, const Mesh_XYZ::bstf<T2>& from, const Op& op )
{
    // left face
    for ( int j = 0; j < to.ncy; ++j )
        for ( int k = to.zoff; k < to.zoff + to.nczp; ++k )
            op(to(0,j,k,0), from(0,j,k,0));

    // right face
    for ( int j = 0; j < to.ncy; ++j )
        for ( int k = to.zoff; k < to.zoff + to.nczp; ++k )
            op(to(to.ncx-1,j,k,1), from(to.ncx-1,j,k,1));

    // front face
    for ( int i = 0; i < to.ncx; ++i )
        for ( int k = to.zoff; k < to.zoff + to.nczp; ++k )
            op(to(i,0,k,2), from(i,0,k,2));

    // back face
    for ( int i = 0; i < to.ncx; ++i )
        for ( int k = to.zoff; k < to.zoff + to.nczp; ++k )
            op(to(i,to.ncy-1,k,3), from(i,to.ncy-1,k,3));

    // bottom face
    if (to.node == 0)
        for ( int i = 0; i < to.ncx; ++i )
            for ( int j = 0; j < to.ncy; ++j )
                op(to(i,j,0,4), from(i,j,0,4));

    // top face
    if (to.node == to.lastnode)
        for ( int i = 0; i < to.ncx; ++i )
            for ( int j = 0; j < to.ncy; ++j )
                op(to(i,j,to.ncz-1,5), from(i,j,to.ncz-1,5));
}

//---------------------------------------------------------------------------//
//                              end of Mesh_XYZ.t.cc
//---------------------------------------------------------------------------//
