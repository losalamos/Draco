//----------------------------------*-C++-*----------------------------------//
// Mesh_XYZ.t.cc
// Geoffrey M. Furnish
// Wed May 27 11:02:08 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "Mesh_XYZ.hh"
#include "c4/global.hh"
#include "c4/C4_Req.hh"
#include "c4/SpinLock.hh"
#include "ds++/Assert.hh"

template<class T>
void dump( const Mesh_XYZ::cctf<T>& data, char *name )
{
    cout << "dumping a Mesh_XYZ::cctf: " << name << endl;
    {
        C4::HTSyncSpinLock h;
        char buf[80];
        for( int i=0; i < data.size(); i++ ) {
            sprintf( buf, "node %d, cell %d, value=%lf \n",
                     C4::node(), i, data(i) );
            cout << buf;
	}
    }
}

template<class T>
void dump( const Mesh_XYZ::fcdtf<T>& data, char *name )
{
    cout << "dumping a Mesh_XYZ::fcdtf: " << name << endl;
    {
        C4::HTSyncSpinLock h;
        char buf[80];
        for( int i=0; i < data.size()/6; i++ ) {
            for( int j=0; j < 6; j++ ) {
                sprintf( buf, "node %d, cell %d, face %d, value=%lf \n",
                         C4::node(), i, j, data(i,j) );
                cout << buf;
            }
        }
    }
}

template<class T>
void dump( const Mesh_XYZ::nctf<T>& data, char *name )
{
    cout << "dumping a Mesh_XYZ::nctf: " << name << endl;
    {
        C4::HTSyncSpinLock h;
        char buf[80];
        for( int i=0; i < data.size(); i++ ) {
            sprintf( buf, "processor %d, node %d, value=%lf \n",
                     C4::node(), i, data(i) );
            cout << buf;
	}
    }
}

template<class T>
void dump( const Mesh_XYZ::vctf<T>& data, char *name )
{
    cout << "dumping a Mesh_XYZ::vctf: " << name << endl;
    {
        C4::HTSyncSpinLock h;
        char buf[80];
        for( int i=0; i < data.size()/8; i++ ) {
            for( int j=0; j < 8; j++ ) {
                sprintf( buf, "node %d, cell %d, vertex %d, value=%lf \n",
                         C4::node(), i, j, data(i,j) );
                cout << buf;
            }
        }
    }
}

template<class T>
bool Mesh_XYZ::fcdtf<T>::operator==( const Mesh_XYZ::fcdtf<T>& x ) const
{
    if (this == &x)
        return true;

    if ( x.size() != this->size() )
        return false;

    if (x.mesh != this->mesh)
        return false;

    if (x.data != this->data)
        return false;

    return true;
}

template<class T>
bool Mesh_XYZ::fcdtf<T>::operator!=( const Mesh_XYZ::fcdtf<T>& x ) const
{
    return !(*this == x);
}

template<class T>
bool Mesh_XYZ::fcdtf<T>::operator<( const Mesh_XYZ::fcdtf<T>& x ) const
{
    return this->data < x.data;
}

template<class T>
bool Mesh_XYZ::fcdtf<T>::operator>( const Mesh_XYZ::fcdtf<T>& x ) const
{
    return (x < *this);
}

template<class T>
bool Mesh_XYZ::fcdtf<T>::operator<=( const Mesh_XYZ::fcdtf<T>& x ) const
{
    return !(x < *this);
}

template<class T>
bool Mesh_XYZ::fcdtf<T>::operator>=( const Mesh_XYZ::fcdtf<T>& x ) const
{
    return !(*this < x);
}

template<class T>
bool Mesh_XYZ::cctf<T>::operator==( const Mesh_XYZ::cctf<T>& x ) const
{
    if (this == &x)
        return true;

    if ( x.size() != this->size() )
        return false;

    if (x.mesh != this->mesh)
        return false;

    if (x.data != this->data)
        return false;

    return true;
}

template<class T>
bool Mesh_XYZ::cctf<T>::operator!=( const Mesh_XYZ::cctf<T>& x ) const
{
    return !(*this == x);
}

template<class T>
bool Mesh_XYZ::cctf<T>::operator<( const Mesh_XYZ::cctf<T>& x ) const
{
    return this->data < x.data;
}

template<class T>
bool Mesh_XYZ::cctf<T>::operator>( const Mesh_XYZ::cctf<T>& x ) const
{
    return (x < *this);
}

template<class T>
bool Mesh_XYZ::cctf<T>::operator<=( const Mesh_XYZ::cctf<T>& x ) const
{
    return !(x < *this);
}

template<class T>
bool Mesh_XYZ::cctf<T>::operator>=( const Mesh_XYZ::cctf<T>& x ) const
{
    return !(*this < x);
}

template<class T>
bool Mesh_XYZ::nctf<T>::operator==( const Mesh_XYZ::nctf<T>& x ) const
{
    if (this == &x)
        return true;

    if ( x.size() != this->size() )
        return false;

    if (x.mesh != this->mesh)
        return false;

    if (x.data != this->data)
        return false;

    return true;
}

template<class T>
bool Mesh_XYZ::nctf<T>::operator!=( const Mesh_XYZ::nctf<T>& x ) const
{
    return !(*this == x);
}

template<class T>
bool Mesh_XYZ::nctf<T>::operator<( const Mesh_XYZ::nctf<T>& x ) const
{
    return this->data < x.data;
}

template<class T>
bool Mesh_XYZ::nctf<T>::operator>( const Mesh_XYZ::nctf<T>& x ) const
{
    return (x < *this);
}

template<class T>
bool Mesh_XYZ::nctf<T>::operator<=( const Mesh_XYZ::nctf<T>& x ) const
{
    return !(x < *this);
}

template<class T>
bool Mesh_XYZ::nctf<T>::operator>=( const Mesh_XYZ::nctf<T>& x ) const
{
    return !(*this < x);
}

template<class T>
bool Mesh_XYZ::vctf<T>::operator==( const Mesh_XYZ::vctf<T>& x ) const
{
    if (this == &x)
        return true;

    if ( x.size() != this->size() )
        return false;

    if (x.mesh != this->mesh)
        return false;

    if (x.data != this->data)
        return false;

    return true;
}

template<class T>
bool Mesh_XYZ::vctf<T>::operator!=( const Mesh_XYZ::vctf<T>& x ) const
{
    return !(*this == x);
}

template<class T>
bool Mesh_XYZ::vctf<T>::operator<( const Mesh_XYZ::vctf<T>& x ) const
{
    return this->data < x.data;
}

template<class T>
bool Mesh_XYZ::vctf<T>::operator>( const Mesh_XYZ::vctf<T>& x ) const
{
    return (x < *this);
}

template<class T>
bool Mesh_XYZ::vctf<T>::operator<=( const Mesh_XYZ::vctf<T>& x ) const
{
    return !(x < *this);
}

template<class T>
bool Mesh_XYZ::vctf<T>::operator>=( const Mesh_XYZ::vctf<T>& x ) const
{
    return !(*this < x);
}

template<class T>
bool Mesh_XYZ::bstf<T>::operator==( const Mesh_XYZ::bstf<T>& x ) const
{
    if (this == &x)
        return true;

    if ( x.size() != this->size() )
        return false;

    if (x.mesh != this->mesh)
        return false;

    if (x.data != this->data)
        return false;

    return true;
}

template<class T>
bool Mesh_XYZ::bstf<T>::operator!=( const Mesh_XYZ::bstf<T>& x ) const
{
    return !(*this == x);
}

template<class T>
bool Mesh_XYZ::bstf<T>::operator<( const Mesh_XYZ::bstf<T>& x ) const
{
    return this->data < x.data;
}

template<class T>
bool Mesh_XYZ::bstf<T>::operator>( const Mesh_XYZ::bstf<T>& x ) const
{
    return (x < *this);
}

template<class T>
bool Mesh_XYZ::bstf<T>::operator<=( const Mesh_XYZ::bstf<T>& x ) const
{
    return !(x < *this);
}

template<class T>
bool Mesh_XYZ::bstf<T>::operator>=( const Mesh_XYZ::bstf<T>& x ) const
{
    return !(*this < x);
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
Mesh_XYZ::gvctf<T>& 
Mesh_XYZ::gvctf<T>::operator=( const Mesh_XYZ::vctf<T>& c )
{
    for( int k=zoff; k < zoff + nczp; k++ )
        for( int j=0; j < ncy; j++ )
            for( int i=0; i < ncx; i++ )
                for ( int v=0; v < 8; v++ )
                    data(v,i,j,k) = c(i,j,k,v);

    update_gvctf();
    
    return *this;
}

template<class T>
void Mesh_XYZ::gvctf<T>::update_gvctf()
{
    using namespace C4;
    C4_Req lrcv, rrcv;

    dsxx::Mat3<T> lrbuf( &data(0,0,0,zoff-1), 8, ncx, ncy );
    dsxx::Mat3<T> rrbuf( &data(0,0,0,zoff+nczp), 8, ncx, ncy );

    if (node > 0)
        RecvAsync( lrcv, &lrbuf(0,0,0), 8*ncx*ncy, node-1 );

    if (node < lastnode)
        RecvAsync( rrcv, &rrbuf(0,0,0), 8*ncx*ncy, node+1 );

    if (node > 0)
        Send( &data(0,0,0,zoff), 8*ncx*ncy, node-1 );

    if (node < lastnode)
        Send( &data(0,0,0,zoff+nczp-1), 8*ncx*ncy, node+1 );
}

template<class T>
void Mesh_XYZ::bstf<T>::next_element( int& i, int& j, int& k, int& f) const
{
    if (f == 0)
    // on left face
    {
        if (j != ncy - 1)
        {
            ++j;
        }
        else if (k != zoff + nczp - 1)
        {
            j = 0;
            ++k;
        }
        else
        // go to right face
        {
            i = ncx - 1;
            j = 0;
            k = zoff;
            f = 1;
        }
    }
    else if (f == 1)
    // on right face
    {
        if (j != ncy - 1)
        {
            ++j;
        }
        else if (k != zoff + nczp - 1)
        {
            j = 0;
            ++k;
        }
        else
        // go to front face
        {
            i = 0;
            j = 0;
            k = zoff;
            f = 2;
        }
    }
    else if (f == 2)
    // on front face
    {
        if (i != ncx - 1)
        {
            ++i;
        }
        else if (k != zoff + nczp - 1)
        {
            i = 0;
            ++k;
        }
        else
        {
        // go to back face
            i = 0;
            j = ncy - 1;
            k = zoff;
            f = 3;
        }
    }
    else if (f == 3)
    // on back face
    {
        if (i != ncx - 1)
        {
            ++i;
        }
        else if (k != zoff + nczp - 1)
        {
            i = 0;
            ++k;
        }
        else
        {
            if (node == 0)
            // go to bottom face
            {
                i = 0;
                j = 0;
                k = 0;
                f = 4;
            }
            else if (node == lastnode)
            // go to top face
            {
                i = 0;
                j = 0;
                k = ncz - 1;
                f = 5;
            }
            else
            // at end; use dummy face
            {
                f = -1;
            }
        }
    }
    else if (f == 4)
    // on bottom face
    {
        if (i != ncx - 1)
        {
            ++i;
        }
        else if (j != ncy - 1)
        {
            i = 0;
            ++j;
        }
        else
        {
            if (node == lastnode)
            // go to top face
            {
                i = 0;
                j = 0;
                k = ncz - 1;
                f = 5;
            }
            else
            // at end; use dummy face
            {
                f = -1;
            }
        }
    }
    else
    // on top face
    {
        if (i != ncx - 1)
        {
            ++i;
        }
        else if (j != ncy - 1)
        {
            i = 0;
            ++j;
        }
        else
        // at end; use dummy face
        {
            f = -1;
        }
    }
}

template<class T>
void Mesh_XYZ::bstf<T>::get_indexes
( int& i, int& j, int& k, int& f, const int index) const
{
    int locindex;

    if (index < ncy*nczp)
    // on left face
    {
        locindex = index;
        f = 0;
        k = locindex/ncy + zoff;
        j = locindex - (k-zoff)*ncy;
        i = 0;
    }
    else if (index < 2*ncy*nczp)
    // on right face
    {
        locindex = index - ncy*nczp;
        f = 1;
        k = locindex/ncy + zoff;
        j = locindex - (k-zoff)*ncy;
        i = ncx - 1;
    }
    else if (index < (2*ncy+ncx)*nczp)
    // on front face
    {
        locindex = index - 2*ncy*nczp;
        f = 2;
        k = locindex/ncx + zoff;
        j = 0;
        i = locindex - (k-zoff)*ncx;
    }
    else if (index < 2*(ncy+ncx)*nczp)
    // on back face
    {
        locindex = index - (2*ncy+ncx)*nczp;
        f = 3;
        k = locindex/ncx + zoff;
        j = ncy - 1;
        i = locindex - (k-zoff)*ncx;
    }
    else if (node == 0 && index < 2*(ncy+ncx)*nczp + ncx*ncy)
    // on bottom face
    {
        locindex = index - 2*(ncy+ncx)*nczp;
        f = 4;
        k = 0;
        j = locindex/ncx;
        i = locindex - j*ncx;
    }
    else
    // on top face
    {
        locindex = index - 2*(ncy+ncx)*nczp;
        if (node == 0)
            locindex -= ncx*ncy;
        f = 5;
        k = ncz - 1;
        j = locindex/ncx;
        i = locindex - j*ncx;
    }
}

template<class T>
Mesh_XYZ::bstf<T>::iterator& Mesh_XYZ::bstf<T>::iterator::operator++()
{
    bfield->next_element(i,j,k,f);
    if (f == -1)
        *this = bfield->end();
    else
        p = &(*bfield)(i,j,k,f);

    return *this;
}

template<class T>
Mesh_XYZ::bstf<T>::iterator
Mesh_XYZ::bstf<T>::iterator::operator++( int dummy )
{
    Mesh_XYZ::bstf<T>::iterator temp = *this;

    bfield->next_element(i,j,k,f);
    if (f == -1)
        *this = bfield->end();
    else
        p = &(*bfield)(i,j,k,f);

    return temp;
}

template<class T>
Mesh_XYZ::bstf<T>::iterator&
Mesh_XYZ::bstf<T>::iterator::operator=( const iterator& iter )
{
    p = iter.p;
    i = iter.i;
    j = iter.j;
    k = iter.k;
    f = iter.f;
    bfield = iter.bfield;

    return *this;
}

template<class T>
Mesh_XYZ::bstf<T>::const_iterator::const_iterator( const iterator& iter )
{
    p = iter.p;
    i = iter.i;
    j = iter.j;
    k = iter.k;
    f = iter.f;
    bfield = iter.bfield;
}

template<class T>
Mesh_XYZ::bstf<T>::const_iterator&
Mesh_XYZ::bstf<T>::const_iterator::operator++()
{
    bfield->next_element(i,j,k,f);
    if (f == -1)
        *this = bfield->end();
    else
        p = &(*bfield)(i,j,k,f);

    return *this;
}

template<class T>
Mesh_XYZ::bstf<T>::const_iterator
Mesh_XYZ::bstf<T>::const_iterator::operator++( int dummy )
{
    Mesh_XYZ::bstf<T>::const_iterator temp = *this;

    bfield->next_element(i,j,k,f);
    if (f == -1)
        *this = bfield->end();
    else
        p = &(*bfield)(i,j,k,f);

    return temp;
}

template<class T>
Mesh_XYZ::bstf<T>::const_iterator&
Mesh_XYZ::bstf<T>::const_iterator::operator=( const const_iterator& iter )
{
    p = iter.p;
    i = iter.i;
    j = iter.j;
    k = iter.k;
    f = iter.f;
    bfield = iter.bfield;

    return *this;
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
const T& Mesh_XYZ::bstf<T>::operator()( int i, int j, int k, int f ) const
{
    Assert (   (i == 0 && f == 0) || (i == ncx - 1 && f == 1)
            || (j == 0 && f == 2) || (j == ncy - 1 && f == 3)
            || (k == 0 && f == 4) || (k == ncz - 1 && f == 5));

    return data(i,j,k,f);
}

template<class T>
const T& Mesh_XYZ::bstf<T>::operator[]( int index ) const
{
    int i, j, k, f;
    get_indexes(i,j,k,f,index);
    return data(i,j,k,f);
}

template<class T>
T& Mesh_XYZ::bstf<T>::operator[]( int index )
{
    int i, j, k, f;
    get_indexes(i,j,k,f,index);
    return data(i,j,k,f);
}

template<class T>
Mesh_XYZ::bstf<T>::iterator Mesh_XYZ::bstf<T>::begin()
{
    return Mesh_XYZ::bstf<T>::iterator(&data(0,0,zoff,0),0,0,zoff,0,this);
}

template<class T>
Mesh_XYZ::bstf<T>::iterator Mesh_XYZ::bstf<T>::end()
{
    int i, j, k, f;
    T* p;

    i = ncx - 1;
    j = ncy - 1;
    k = zoff + nczp - 1;
    f = 5;
    p = &data(i,j,k,f);
    ++p;

    return Mesh_XYZ::bstf<T>::iterator(p,i,j,k,f,this);
}

template<class T>
Mesh_XYZ::bstf<T>::const_iterator Mesh_XYZ::bstf<T>::begin() const
{
    return Mesh_XYZ::bstf<T>::const_iterator
                              (&data(0,0,zoff,0),0,0,zoff,0,this);
}

template<class T>
Mesh_XYZ::bstf<T>::const_iterator Mesh_XYZ::bstf<T>::end() const
{
    int i, j, k, f;
    const T* p;

    i = ncx - 1;
    j = ncy - 1;
    k = zoff + nczp - 1;
    f = 5;
    p = &data(i,j,k,f);
    ++p;

    return Mesh_XYZ::bstf<T>::const_iterator(p,i,j,k,f,this);
}

template<class T>
Mesh_XYZ::bstf<T>::size_type Mesh_XYZ::bstf<T>::size() const
{
    int tempsize = 2*(ncx+ncy)*nczp;
    if (node == 0)
        tempsize += ncx*ncy;
    if (node == lastnode)
        tempsize += ncx*ncy;

    return tempsize;
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
          if (i != 0)
              op(to(i,j,k,0), gfrom(i-1,j,k));
          if (i != to.ncx - 1)
              op(to(i,j,k,1), gfrom(i+1,j,k));
          if (j != 0)
              op(to(i,j,k,2), gfrom(i,j-1,k));
          if (j != to.ncy - 1)
              op(to(i,j,k,3), gfrom(i,j+1,k));
          if (k != 0)
              op(to(i,j,k,4), gfrom(i,j,k-1));
          if (k != to.ncz - 1)
              op(to(i,j,k,5), gfrom(i,j,k+1));
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
void Mesh_XYZ::scatter
( Mesh_XYZ::vctf<T1>& to, const Mesh_XYZ::fcdtf<T2>& from, const Op& op )
{
    for ( int i = 0; i < to.ncx; ++i )
      for ( int j = 0; j < to.ncy; ++j )
        for ( int k = to.zoff; k < to.zoff + to.nczp; ++k )
        {
          op(to(i,j,k,0), from(i,j,k,0));
          op(to(i,j,k,0), from(i,j,k,2));
          op(to(i,j,k,0), from(i,j,k,4));

          op(to(i,j,k,1), from(i,j,k,1));
          op(to(i,j,k,1), from(i,j,k,2));
          op(to(i,j,k,1), from(i,j,k,4));

          op(to(i,j,k,2), from(i,j,k,0));
          op(to(i,j,k,2), from(i,j,k,3));
          op(to(i,j,k,2), from(i,j,k,4));

          op(to(i,j,k,3), from(i,j,k,1));
          op(to(i,j,k,3), from(i,j,k,3));
          op(to(i,j,k,3), from(i,j,k,4));

          op(to(i,j,k,4), from(i,j,k,0));
          op(to(i,j,k,4), from(i,j,k,2));
          op(to(i,j,k,4), from(i,j,k,5));

          op(to(i,j,k,5), from(i,j,k,1));
          op(to(i,j,k,5), from(i,j,k,2));
          op(to(i,j,k,5), from(i,j,k,5));

          op(to(i,j,k,6), from(i,j,k,0));
          op(to(i,j,k,6), from(i,j,k,3));
          op(to(i,j,k,6), from(i,j,k,5));

          op(to(i,j,k,7), from(i,j,k,1));
          op(to(i,j,k,7), from(i,j,k,3));
          op(to(i,j,k,7), from(i,j,k,5));

        }  
}

template <class T1, class T2, class Op>
void Mesh_XYZ::scatter
( Mesh_XYZ::nctf<T1>& to, const Mesh_XYZ::vctf<T2>& from, const Op& op )
{
    Mesh_XYZ::gvctf<T2> gfrom(from);
    for ( int i = 0; i < from.ncx; ++i )
      for ( int j = 0; j < from.ncy; ++j )
      {
          if (from.zoff != 0)
          {
              op(to(i,j,from.zoff), gfrom(i,j,from.zoff-1,4));
              op(to(i+1,j,from.zoff), gfrom(i,j,from.zoff-1,5));
              op(to(i,j+1,from.zoff), gfrom(i,j,from.zoff-1,6));
              op(to(i+1,j+1,from.zoff), gfrom(i,j,from.zoff-1,7));
          }
          for ( int k = from.zoff; k < from.zoff + from.nczp; ++k )
          {
              op(to(i,j,k), gfrom(i,j,k,0));
              op(to(i+1,j,k), gfrom(i,j,k,1));
              op(to(i,j+1,k), gfrom(i,j,k,2));
              op(to(i+1,j+1,k), gfrom(i,j,k,3));
              op(to(i,j,k+1), gfrom(i,j,k,4));
              op(to(i+1,j,k+1), gfrom(i,j,k,5));
              op(to(i,j+1,k+1), gfrom(i,j,k,6));
              op(to(i+1,j+1,k+1), gfrom(i,j,k,7));
          }
          if (from.zoff + from.nczp != from.ncz)
          {
              op(to(i,j,from.zoff + from.nczp),
                 gfrom(i,j,from.zoff + from.nczp,0));
              op(to(i+1,j,from.zoff + from.nczp),
                 gfrom(i,j,from.zoff + from.nczp,1));
              op(to(i,j+1,from.zoff + from.nczp),
                 gfrom(i,j,from.zoff + from.nczp,2));
              op(to(i+1,j+1,from.zoff + from.nczp),
                 gfrom(i,j,from.zoff + from.nczp,3));
          }
      }
}

template <class T1, class T2, class Op>
void Mesh_XYZ::scatter
( Mesh_XYZ::cctf<T1>& to, const Mesh_XYZ::vctf<T2>& from, const Op& op )
{
    for ( int i = 0; i < to.ncx; ++i )
      for ( int j = 0; j < to.ncy; ++j )
        for ( int k = to.zoff; k < to.zoff + to.nczp; ++k )
          for ( int v = 0; v < 8; ++v )
	    op(to(i,j,k), from(i,j,k,v));
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

template <class T1, class T2, class Op>
void Mesh_XYZ::gather
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
void Mesh_XYZ::gather
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

template <class T1, class T2, class Op>
void Mesh_XYZ::gather
( Mesh_XYZ::vctf<T1>& to, const Mesh_XYZ::nctf<T2>& from, const Op& op )
{
    for ( int i = 0; i < to.ncx; ++i )
      for ( int j = 0; j < to.ncy; ++j )
        for ( int k = to.zoff; k < to.zoff + to.nczp; ++k )
        {
            op(to(i,j,k,0), from(i,j,k));
            op(to(i,j,k,1), from(i+1,j,k));
            op(to(i,j,k,2), from(i,j+1,k));
            op(to(i,j,k,3), from(i+1,j+1,k));
            op(to(i,j,k,4), from(i,j,k+1));
            op(to(i,j,k,5), from(i+1,j,k+1));
            op(to(i,j,k,6), from(i,j+1,k+1));
            op(to(i,j,k,7), from(i+1,j+1,k+1));
        }
}

template <class T1, class T2, class Op>
void Mesh_XYZ::gather
( Mesh_XYZ::vctf<T1>& to, const Mesh_XYZ::cctf<T2>& from, const Op& op )
{
    for ( int i = 0; i < to.ncx; ++i )
      for ( int j = 0; j < to.ncy; ++j )
        for ( int k = to.zoff; k < to.zoff + to.nczp; ++k )
          for ( int v = 0; v < 8; ++v )
            op(to(i,j,k,v), from(i,j,k));
}

template <class T>
void Mesh_XYZ::swap_faces
( Mesh_XYZ::fcdtf<T>& to, const Mesh_XYZ::fcdtf<T>& from )
{
    Mesh_XYZ::gfcdtf<T> gfrom(from);
    for ( int i = 0; i < to.ncx; ++i )
      for ( int j = 0; j < to.ncy; ++j )
        for ( int k = to.zoff; k < to.zoff + to.nczp; ++k )
        {
          if (i != 0)
              to(i,j,k,0) = gfrom(i-1,j,k,1);
          else
              to(i,j,k,0) = 0;
          if (i != to.ncx - 1)
              to(i,j,k,1) = gfrom(i+1,j,k,0);
          else
              to(i,j,k,1) = 0;
          if (j != 0)
              to(i,j,k,2) = gfrom(i,j-1,k,3);
          else
              to(i,j,k,2) = 0;
          if (j != to.ncy - 1)
              to(i,j,k,3) = gfrom(i,j+1,k,2);
          else
              to(i,j,k,3) = 0;
          if (k != 0)
              to(i,j,k,4) = gfrom(i,j,k-1,5);
          else
              to(i,j,k,4) = 0;
          if (k != to.ncz - 1)
              to(i,j,k,5) = gfrom(i,j,k+1,4);
          else
              to(i,j,k,5) = 0;
	}
}

template <class T>
T Mesh_XYZ::sum( const Mesh_XYZ::cctf<T>& from )
{
    T sum = 0;

    for (Mesh_XYZ::cctf<T>::const_iterator iter = from.begin();
         iter != from.end(); ++iter)
        sum += *iter;

    C4::gsum<T>(sum);

    return sum;
}

template <class T>
T Mesh_XYZ::sum( const Mesh_XYZ::fcdtf<T>& from )
{
    T sum = 0;

    for (Mesh_XYZ::fcdtf<T>::const_iterator iter = from.begin();
         iter != from.end(); ++iter)
        sum += *iter;

    C4::gsum<T>(sum);

    return sum;
}

template <class T>
T Mesh_XYZ::sum( const Mesh_XYZ::nctf<T>& from )
{
    T sum = 0;
    int bound;

    if (C4::node() == (C4::nodes()-1))
        bound = 1;
    else
        bound = 0;

    for ( int i = 0; i <= from.ncx; ++i )
      for ( int j = 0; j <= from.ncy; ++j )
        for ( int k = from.zoff; k < from.zoff + from.nczp + bound; ++k )
          sum += from(i,j,k);

    C4::gsum<T>(sum);

    return sum;
}

template <class T>
T Mesh_XYZ::sum( const Mesh_XYZ::vctf<T>& from )
{
    T sum = 0;

    for (Mesh_XYZ::vctf<T>::const_iterator iter = from.begin();
         iter != from.end(); ++iter)
        sum += *iter;

    C4::gsum<T>(sum);

    return sum;
}

template <class T>
T Mesh_XYZ::sum( const Mesh_XYZ::bstf<T>& from )
{
    T sum = 0;

    for (Mesh_XYZ::bstf<T>::const_iterator iter = from.begin();
         iter != from.end(); ++iter)
        sum += *iter;

    C4::gsum<T>(sum);

    return sum;
}

template <class T>
T Mesh_XYZ::min( const Mesh_XYZ::cctf<T>& from )
{
    Mesh_XYZ::cctf<T>::const_iterator iter = from.begin();
    T minimum = *iter;

    for (Mesh_XYZ::cctf<T>::const_iterator iter = from.begin();
         iter != from.end(); ++iter)
        minimum = (minimum < *iter) ? minimum : *iter;

    C4::gmin<T>(minimum);

    return minimum;
}

template <class T>
T Mesh_XYZ::min( const Mesh_XYZ::fcdtf<T>& from )
{
    Mesh_XYZ::fcdtf<T>::const_iterator iter = from.begin();
    T minimum = *iter;

    for (Mesh_XYZ::fcdtf<T>::const_iterator iter = from.begin();
         iter != from.end(); ++iter)
        minimum = (minimum < *iter) ? minimum : *iter;

    C4::gmin<T>(minimum);

    return minimum;
}

template <class T>
T Mesh_XYZ::min( const Mesh_XYZ::nctf<T>& from )
{
    Mesh_XYZ::nctf<T>::const_iterator iter = from.begin();
    T minimum = *iter;

    for (Mesh_XYZ::nctf<T>::const_iterator iter = from.begin();
         iter != from.end(); ++iter)
        minimum = (minimum < *iter) ? minimum : *iter;

    C4::gmin<T>(minimum);

    return minimum;
}

template <class T>
T Mesh_XYZ::min( const Mesh_XYZ::vctf<T>& from )
{
    Mesh_XYZ::vctf<T>::const_iterator iter = from.begin();
    T minimum = *iter;

    for (Mesh_XYZ::vctf<T>::const_iterator iter = from.begin();
         iter != from.end(); ++iter)
        minimum = (minimum < *iter) ? minimum : *iter;

    C4::gmin<T>(minimum);

    return minimum;
}

template <class T>
T Mesh_XYZ::min( const Mesh_XYZ::bstf<T>& from )
{
    Mesh_XYZ::bstf<T>::const_iterator iter = from.begin();
    T minimum = *iter;

    for (Mesh_XYZ::bstf<T>::const_iterator iter = from.begin();
         iter != from.end(); ++iter)
        minimum = (minimum < *iter) ? minimum : *iter;

    C4::gmin<T>(minimum);

    return minimum;
}

template <class T>
T Mesh_XYZ::max( const Mesh_XYZ::cctf<T>& from )
{
    Mesh_XYZ::cctf<T>::const_iterator iter = from.begin();
    T maximum = *iter;

    for (Mesh_XYZ::cctf<T>::const_iterator iter = from.begin();
         iter != from.end(); ++iter)
        maximum = (maximum > *iter) ? maximum : *iter;

    C4::gmax<T>(maximum);

    return maximum;
}

template <class T>
T Mesh_XYZ::max( const Mesh_XYZ::fcdtf<T>& from )
{
    Mesh_XYZ::fcdtf<T>::const_iterator iter = from.begin();
    T maximum = *iter;

    for (Mesh_XYZ::fcdtf<T>::const_iterator iter = from.begin();
         iter != from.end(); ++iter)
        maximum = (maximum > *iter) ? maximum : *iter;

    C4::gmax<T>(maximum);

    return maximum;
}

template <class T>
T Mesh_XYZ::max( const Mesh_XYZ::nctf<T>& from )
{
    Mesh_XYZ::nctf<T>::const_iterator iter = from.begin();
    T maximum = *iter;

    for (Mesh_XYZ::nctf<T>::const_iterator iter = from.begin();
         iter != from.end(); ++iter)
        maximum = (maximum > *iter) ? maximum : *iter;

    C4::gmax<T>(maximum);

    return maximum;
}

template <class T>
T Mesh_XYZ::max( const Mesh_XYZ::vctf<T>& from )
{
    Mesh_XYZ::vctf<T>::const_iterator iter = from.begin();
    T maximum = *iter;

    for (Mesh_XYZ::vctf<T>::const_iterator iter = from.begin();
         iter != from.end(); ++iter)
        maximum = (maximum > *iter) ? maximum : *iter;

    C4::gmax<T>(maximum);

    return maximum;
}

template <class T>
T Mesh_XYZ::max( const Mesh_XYZ::bstf<T>& from )
{
    Mesh_XYZ::bstf<T>::const_iterator iter = from.begin();
    T maximum = *iter;

    for (Mesh_XYZ::bstf<T>::const_iterator iter = from.begin();
         iter != from.end(); ++iter)
        maximum = (maximum > *iter) ? maximum : *iter;

    C4::gmax<T>(maximum);

    return maximum;
}

//---------------------------------------------------------------------------//
//                              end of Mesh_XYZ.t.cc
//---------------------------------------------------------------------------//
