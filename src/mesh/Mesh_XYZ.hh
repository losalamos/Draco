//----------------------------------*-C++-*----------------------------------//
// Mesh_XYZ.hh
// RHS Linux User
// Thu Nov 20 22:08:34 1997
//---------------------------------------------------------------------------//
// @> A 3-d cartesian structured mesh facility.
//---------------------------------------------------------------------------//

#ifndef __mesh_Mesh_XYZ_hh__
#define __mesh_Mesh_XYZ_hh__

#include "Mesh_DB.hh"

#include "traits/MT_traits.hh"
#include "ds++/Mat.hh"
#include "ds++/SP.hh"

// using rtt_dsxx::Mat1;
// using rtt_dsxx::Mat2;
// using rtt_dsxx::Mat3;

#include "c4/NodeInfo.hh"

#include "xm/xm.hh"

#include <iterator>

#include "./MeshXYZConnFacesAroundVertices.hh"

struct XYZ_Mapper : public C4::NodeInfo
{
    typedef rtt_dsxx::Mat1<int>::size_type size_type;

    typedef rtt_mesh::Mesh_DB Mesh_DB;

    int nct;			// # of total cells in problem.
    int nxy;                    // # of cells in an x-y plane.
    int nczp;                   // # of cells in z this processor.
    int ncp;			// # of cells on this processor.

    int ncx;
    int ncy;
    int ncz;

    int zoff;                   // z index offset of this processor's cells.
    int goff;			// global offset of this processor's cells.

    Mesh_DB mdb;

// Methods

    XYZ_Mapper( const Mesh_DB& mdb_in );

    const Mesh_DB &get_Mesh_DB() const { return mdb; }

// Convert a local cell index to its (i,j,k) indexes in whole domain.

    int I(int x) const { x += goff; return (x % nxy) % ncx; }
    int J(int x) const { x += goff; return (x % nxy) / ncx; }
    int K(int x) const { x += goff; return  x / nxy; }

    int goffset( int i, int j, int k ) const { return ncx*(k*ncy + j) +i; }
    int goffset( int nc ) const
    {
	nc += goff;
	return goffset( I(nc), J(nc), K(nc) );
    }

    int local_cell_index( int i, int j, int k ) const
    {
	return goffset(i,j,k) - goff;
    }

// Stuff needed to instantiate a Banded_Matrix

    int row_offset() const { return goff; }
    int nrows_this_processor() const { return ncp; }
    int nrows_total() const { return nct; }
    bool verbose() const { return false; }
};

//===========================================================================//
// class Mesh_XYZ - A 3-d cartesian mesh class

// This is a 3-d cartesisan structured mesh.  The main purpose of having this
// class is in order for it to be instantiated by various test articles
// throughout the Draco system.  It could also be useful to Draco clients as
// a point of reference during the construction of their own mesh classes.
//===========================================================================//

class Mesh_XYZ : private XYZ_Mapper
{

  public:

    template<class T> class cctf;
    template<class T> class gcctf;
    template<class T> class fcdtf;
    template<class T> class gfcdtf;
    template<class T> class nctf;
    template<class T> class gnctf;
    template<class T> class vctf;
    template<class T> class gvctf;
    template<class T> class bstf;
    template<class T, int N> class tiny_vec;

    typedef rtt_mesh::Mesh_DB Mesh_DB;

    typedef XYZ_Mapper::size_type size_type;
    typedef rtt_dsxx::SP<const Mesh_XYZ> FieldConstructor;

    typedef cctf<double> ccsf;
    typedef fcdtf<double> fcdsf;
    typedef nctf<double> ncsf;
    typedef vctf<double> vcsf;
    typedef bstf<double> bssf;
    typedef cctf<int> ccif;
    typedef fcdtf<int> fcdif;
    typedef nctf<int> ncif;
    typedef vctf<int> vcif;
    typedef bstf<int> bsif;
    typedef cctf<tiny_vec<double, 3> > ccvsf;
    typedef fcdtf<tiny_vec<double, 3> > fcdvsf;
    typedef nctf<tiny_vec<double, 3> > ncvsf;
    typedef vctf<tiny_vec<double, 3> > vcvsf;
    typedef bstf<tiny_vec<double, 3> > bsvsf;
    typedef gcctf<double> gccsf;
    typedef gvctf<double> gvcsf;
    typedef tiny_vec<double, 3> vec3;

    // ConnFacesAroundVertices is a class that defines objects
    // that will iterate through a face-centered field around each vertex,
    // before going onto the next vertex's faces.

    template<class FaceField>
    class ConnFacesAroundVertices
	: public rtt_mesh::MeshXYZConnFacesAroundVertices<FaceField>
    {
	friend class Mesh_XYZ;

      public:

	ConnFacesAroundVertices(FaceField& field_):
	    rtt_mesh::MeshXYZConnFacesAroundVertices<FaceField>(field_)
	{
	    // empty
	}
    };

// Face centered discontinuous field
// Has a value on each face in each cell.

    template<class T>
    class fcdtf : private XYZ_Mapper,
                  public xm::Indexable< T, fcdtf<T> > {
	friend class Mesh_XYZ;
	friend class gfcdtf<T>;

	rtt_dsxx::Mat2<T> data;
        const Mesh_XYZ* mesh;
        FieldConstructor spm;

	fcdtf( const Mesh_XYZ *m )
	    : XYZ_Mapper( m->get_Mesh_DB() ),
	      data( m->get_ncells(), 6 ), mesh( m ), spm( 0 )
	{}

      public:
	typedef T value_type;
        typedef T& reference;
        typedef const T& const_reference;
        typedef typename rtt_dsxx::Mat2<T>::pointer pointer;
        typedef typename rtt_dsxx::Mat2<T>::const_pointer const_pointer;
        typedef typename rtt_dsxx::Mat2<T>::iterator iterator;
        typedef typename rtt_dsxx::Mat2<T>::const_iterator const_iterator;
        typedef typename rtt_dsxx::Mat2<T>::difference_type difference_type;
        typedef typename rtt_dsxx::Mat2<T>::size_type size_type;

	fcdtf( const FieldConstructor& spm_ )
	    : XYZ_Mapper( spm_->get_Mesh_DB() ),
	      data( spm_->get_ncells(), 6 ), mesh( spm_.bp() ), spm( spm_ )
	{}

	fcdtf( const fcdtf& f ) 
          : XYZ_Mapper( f.get_Mesh_DB() ),
            data( f.data ), mesh( f.mesh ), spm( f.spm ) {}

        ~fcdtf() {}

        const Mesh_XYZ& get_Mesh() const { return *mesh; }

        const FieldConstructor& get_FieldConstructor() const { return spm; }

	fcdtf& operator=( T x ) { data = x; return *this; }

        fcdtf& operator=( const fcdtf& f )
	{
            data = f.data;
            return *this;
	}

        template<class X>
        fcdtf& operator=( const xm::Xpr< T, X, fcdtf >& x )
        {
            return assign_from( x );
        }

    // i, j, k == global xyz cell indicies
    // f == face index
    // c == local cell index

	T&       operator()( int c, int f )       { return data(c,f); }
	const T& operator()( int c, int f ) const { return data(c,f); }

	T& operator()( int i, int j, int k, int f )
	{
	    return data( local_cell_index(i,j,k), f );
	}
	const T& operator()( int i, int j, int k, int f ) const 
	{
	    return data( local_cell_index(i,j,k), f );
	}

        const T& operator[]( int i ) const { return data[i]; }
        T& operator[]( int i ) { return data[i]; }

        iterator begin() { return data.begin(); }
        iterator end() { return data.end(); }

        const_iterator begin() const { return data.begin(); }
        const_iterator end() const { return data.end(); }

        size_type size() const { return data.size(); }
        size_type max_size() const { return data.size(); }
        bool empty() const { return data.empty(); }

        void swap ( fcdtf& x ) { data.swap(x.data); }

        bool operator==( const fcdtf& x ) const;
        bool operator!=( const fcdtf& x ) const;
        bool operator<( const fcdtf& x ) const;
        bool operator>( const fcdtf& x ) const;
        bool operator<=( const fcdtf& x ) const;
        bool operator>=( const fcdtf& x ) const;
    };

// Guarded face centered discontinuous field
// Has a value on each face in each cell.

    template<class T>
    class gfcdtf : private XYZ_Mapper,
                   public xm::Indexable< T, gfcdtf<T> > {
	rtt_dsxx::Mat4<T> data;
        const Mesh_XYZ* mesh;
        FieldConstructor spm;

      public:
        typedef typename rtt_dsxx::Mat4<T>::iterator iterator;
        typedef typename rtt_dsxx::Mat4<T>::const_iterator const_iterator;

	gfcdtf( const FieldConstructor& spm_ )
	    : XYZ_Mapper( spm_->get_Mesh_DB() ),
	      data( rtt_dsxx::Bounds( 0, 5 ),
                    rtt_dsxx::Bounds( 0, ncx - 1 ),
                    rtt_dsxx::Bounds( 0, ncy - 1 ),
                    rtt_dsxx::Bounds( zoff - 1, zoff + nczp ) ),
              mesh( spm_.bp() ), spm( spm_ )
	{}

        const Mesh_XYZ& get_Mesh() const { return *mesh; }

        gfcdtf<T>& operator=( const fcdtf<T>& c );
        void update_gfcdtf();

        gfcdtf<T>( const fcdtf<T>& f )
            : XYZ_Mapper( f.get_Mesh_DB() ),
              data( rtt_dsxx::Bounds( 0, 5 ),
                    rtt_dsxx::Bounds( 0, ncx - 1 ),
                    rtt_dsxx::Bounds( 0, ncy - 1 ),
                    rtt_dsxx::Bounds( zoff - 1, zoff + nczp ) ),
              mesh( f.mesh ), spm( f.spm )
        { *this = f; }

    // i, j, k == global xyz cell indicies
    // f == face index

    // Note that the order of the indexing is different than the
    // order of the data layout.

        const T& operator()( int i, int j, int k, int f ) const
        { return data(f,i,j,k); }
        T& operator()( int i, int j, int k, int f )
        { return data(f,i,j,k); }

        const T& operator[]( int i ) const { return data[i]; }
        T& operator[]( int i ) { return data[i]; }

        iterator begin() { return data.begin(); }
        iterator end() { return data.end(); }

        const_iterator begin() const { return data.begin(); }
        const_iterator end() const { return data.end(); }

        int size() const { return data.size(); }
    };

    typedef gfcdtf<double> gfcdsf;

// Cell centered field
// Has a value in each cell.

    template<class T>
    class cctf : private XYZ_Mapper,
                       public xm::Indexable< T, cctf<T> >
    {
        rtt_dsxx::Mat1<T> data;
        const Mesh_XYZ* mesh;
        FieldConstructor spm;

	cctf( const Mesh_XYZ *m ) 
          : XYZ_Mapper( m->get_Mesh_DB() ),
            data( m->get_ncells() ), mesh( m ), spm( 0 )
        {}

      public:
	typedef T value_type;
        typedef T& reference;
        typedef const T& const_reference;
        typedef typename rtt_dsxx::Mat1<T>::pointer pointer;
        typedef typename rtt_dsxx::Mat1<T>::const_pointer const_pointer;
        typedef typename rtt_dsxx::Mat1<T>::iterator iterator;
        typedef typename rtt_dsxx::Mat1<T>::const_iterator const_iterator;
        typedef typename rtt_dsxx::Mat1<T>::difference_type difference_type;
        typedef typename rtt_dsxx::Mat1<T>::size_type size_type;

	cctf( const FieldConstructor& spm_ ) 
          : XYZ_Mapper( spm_->get_Mesh_DB() ),
            data( spm_->get_ncells() ), mesh( spm_.bp() ), spm( spm_ ) {}

	cctf( const cctf& c ) 
          : XYZ_Mapper( c.get_Mesh_DB() ),
            data( c.data ), mesh( c.mesh ), spm( c.spm ) {}

        ~cctf() {}

        const Mesh_XYZ& get_Mesh() const { return *mesh; }

        const FieldConstructor& get_FieldConstructor() const { return spm; }

        cctf& operator=( T x )
        {
            data = x;
            return *this;
        }

        cctf& operator=( const cctf& c )
	{
            data = c.data;
            return *this;
	}

        template<class X>
        cctf& operator=( const xm::Xpr< T, X, cctf >& x )
        {
            return assign_from( x );
        }

        const T& operator()( int i ) const { return data(i); }
        T& operator()( int i ) { return data(i); }

	T& operator()( int i, int j, int k )
	{
	    return data( local_cell_index(i,j,k) );
	}
	const T& operator()( int i, int j, int k ) const 
	{
	    return data( local_cell_index(i,j,k) );
	}

        const T& operator[]( int i ) const { return data(i); }
        T& operator[]( int i ) { return data(i); }

        iterator begin() { return data.begin(); }
        iterator end() { return data.end(); }

        const_iterator begin() const { return data.begin(); }
        const_iterator end() const { return data.end(); }

        size_type size() const { return data.size(); }
        size_type max_size() const { return data.size(); }
        bool empty() const { return data.empty(); }

        void swap ( cctf& x ) { data.swap(x.data); }

        bool operator==( const cctf& x ) const;
        bool operator!=( const cctf& x ) const;
        bool operator<( const cctf& x ) const;
        bool operator>( const cctf& x ) const;
        bool operator<=( const cctf& x ) const;
        bool operator>=( const cctf& x ) const;

	friend class Mesh_XYZ;
	friend class gcctf<T>;
    };

// Guarded cell centered field
// Has a value in each cell.

    template<class T>
    class gcctf
        : private XYZ_Mapper,
          public xm::Indexable< T, gcctf<T> >
    {
        rtt_dsxx::Mat3<T> data;
        const Mesh_XYZ* mesh;
        FieldConstructor spm;

      public:
        typedef typename rtt_dsxx::Mat3<T>::iterator iterator;
        typedef typename rtt_dsxx::Mat3<T>::const_iterator const_iterator;

        gcctf( const FieldConstructor& spm_ )
            : XYZ_Mapper( spm_->get_Mesh_DB() ),
              data( rtt_dsxx::Bounds( 0, ncx - 1 ),
                    rtt_dsxx::Bounds( 0, ncy - 1 ),
                    rtt_dsxx::Bounds( zoff - 1, zoff + nczp ) ),
              mesh( spm_.bp() ), spm( spm_ )
        {}

        const Mesh_XYZ& get_Mesh() const { return *mesh; }

        gcctf<T>& operator=( const cctf<T>& c );
        void update_guard_cells();

        gcctf<T>( const cctf<T>& c )
            : XYZ_Mapper( c.get_Mesh_DB() ),
              data( rtt_dsxx::Bounds( 0, ncx - 1 ),
                    rtt_dsxx::Bounds( 0, ncy - 1 ),
                    rtt_dsxx::Bounds( zoff - 1, zoff + nczp ) ),
              mesh( c.mesh ), spm( c.spm )
        { *this = c; }

        const T& operator()( int i, int j, int k ) const
        { return data(i,j,k); }
        T& operator()( int i, int j, int k )
        { return data(i,j,k); }

    // Operators needed for glommable expression templates.
        const T& operator[]( int i ) const { return data[i]; }
        T& operator[]( int i ) { return data[i]; }

        iterator begin() { return data.begin(); }
        iterator end()   { return data.end(); }

        const_iterator begin() const { return data.begin(); }
        const_iterator end()   const { return data.end(); }

        int size() const { return data.size(); }
    };

// Node centered field
// Has a value at each node.

    template<class T>
    class nctf : private XYZ_Mapper,
                       public xm::Indexable< T, nctf<T> >
    {
        rtt_dsxx::Mat1<T> data;
        const Mesh_XYZ* mesh;
        FieldConstructor spm;

	static int calcSize(const Mesh_XYZ &m)
	{
	    const int bound = (C4::node() == (C4::nodes()-1)) ? 1 : 0;
	    return (m.get_ncx()+1)*(m.get_ncy()+1)*(m.get_nczp()+bound);
	}

	bool inRange(int i, int j, int k) const
	{
	    const int bound = (C4::node() == (C4::nodes()-1)) ? 1 : 0;
	    bool inRng = i >= 0 && i < (ncx + 1);
	    inRng = inRng && j >= 0 && j < (ncy + 1);
	    inRng = inRng && k >= 0 && k < (nczp + bound);
	    return inRng;
	}

	nctf( const Mesh_XYZ *m ) 
          : XYZ_Mapper( m->get_Mesh_DB() ),
            data( calcSize(*m) ),
            mesh( m ), spm( 0 )
        {}

      public:
	typedef T value_type;
        typedef T& reference;
        typedef const T& const_reference;
        typedef typename rtt_dsxx::Mat1<T>::pointer pointer;
        typedef typename rtt_dsxx::Mat1<T>::const_pointer const_pointer;
        typedef typename rtt_dsxx::Mat1<T>::iterator iterator;
        typedef typename rtt_dsxx::Mat1<T>::const_iterator const_iterator;
        typedef typename rtt_dsxx::Mat1<T>::difference_type difference_type;
        typedef typename rtt_dsxx::Mat1<T>::size_type size_type;

	nctf( const FieldConstructor& spm_ ) 
          : XYZ_Mapper( spm_->get_Mesh_DB() ),
            data( calcSize(*spm_) ),
            mesh( spm_.bp() ), spm( spm_ ) {}

	nctf( const nctf& n ) 
          : XYZ_Mapper( n.get_Mesh_DB() ),
            data( n.data ), mesh( n.mesh ), spm( n.spm ) {}

        ~nctf() {}

        const Mesh_XYZ& get_Mesh() const { return *mesh; }

        const FieldConstructor& get_FieldConstructor() const { return spm; }

        nctf& operator=( T x )
        {
            data = x;
            return *this;
        }

        nctf& operator=( const nctf& n )
	{
            data = n.data;
            return *this;
	}

        template<class X>
        nctf& operator=( const xm::Xpr< T, X, nctf >& x )
        {
            return assign_from( x );
        }

        const T& operator()( int i ) const { return data(i); }
        T& operator()( int i ) { return data(i); }

	T& operator()( int i, int j, int k )
	{
	    Assert(inRange(i,j,k-zoff));
	    return data( i+(ncx+1)*(j+(ncy+1)*(k-zoff)) );
	}
	const T& operator()( int i, int j, int k ) const 
	{
	    Assert(inRange(i,j,k-zoff));
	    return data( i+(ncx+1)*(j+(ncy+1)*(k-zoff)) );
	}

        const T& operator[]( int i ) const { return data(i); }
        T& operator[]( int i ) { return data(i); }

        iterator begin() { return data.begin(); }
        iterator end() { return data.end(); }

        const_iterator begin() const { return data.begin(); }
        const_iterator end() const { return data.end(); }

        size_type size() const { return data.size(); }
        size_type max_size() const { return data.size(); }
        bool empty() const { return data.empty(); }

        void swap ( nctf& x ) { data.swap(x.data); }

        bool operator==( const nctf& x ) const;
        bool operator!=( const nctf& x ) const;
        bool operator<( const nctf& x ) const;
        bool operator>( const nctf& x ) const;
        bool operator<=( const nctf& x ) const;
        bool operator>=( const nctf& x ) const;

	friend class Mesh_XYZ;
	friend class gnctf<T>;
    };

// Guarded node centered field
// Has a value at each node.

    template<class T>
    class gnctf : private XYZ_Mapper,
		  public xm::Indexable< T, gnctf<T> >
    {
        rtt_dsxx::Mat1<T> data;
        const Mesh_XYZ* mesh;
        FieldConstructor spm;

	static int calcSize(const Mesh_XYZ &m)
	{
	    return (m.get_ncx()+1)*(m.get_ncy()+1)*(m.get_nczp()+1);
	}

	bool inRange(int i, int j, int k) const
	{
	    bool inRng = i >= 0 && i < (ncx + 1);
	    inRng = inRng && j >= 0 && j < (ncy + 1);
	    inRng = inRng && k >= 0 && k < (nczp + 1);
	    return inRng;
	}

      public:
	typedef T value_type;
        typedef T& reference;
        typedef const T& const_reference;
        typedef typename rtt_dsxx::Mat1<T>::pointer pointer;
        typedef typename rtt_dsxx::Mat1<T>::const_pointer const_pointer;
        typedef typename rtt_dsxx::Mat1<T>::iterator iterator;
        typedef typename rtt_dsxx::Mat1<T>::const_iterator const_iterator;
        typedef typename rtt_dsxx::Mat1<T>::difference_type difference_type;
        typedef typename rtt_dsxx::Mat1<T>::size_type size_type;

	gnctf( const FieldConstructor& spm_ ) 
	    : XYZ_Mapper( spm_->get_Mesh_DB() ),
	      data( calcSize(*spm_) ),
	      mesh( spm_.bp() ), spm( spm_ ) {}

        gnctf<T>& operator=( const nctf<T>& n );
        void update_gnctf();

	gnctf( const nctf<T>& n )
	    : XYZ_Mapper( n.get_Mesh_DB() ),
              data( calcSize(*n.mesh) ),
	      mesh( n.mesh ), spm( n.spm )
	{
	    *this = n;
	}

        ~gnctf() {}

        const Mesh_XYZ& get_Mesh() const { return *mesh; }

        const FieldConstructor& get_FieldConstructor() const { return spm; }

        const T& operator()( int i ) const { return data(i); }
        T& operator()( int i ) { return data(i); }

	T& operator()( int i, int j, int k )
	{
	    Assert(inRange(i,j,k-zoff));
	    return data( i+(ncx+1)*(j+(ncy+1)*(k-zoff)) );
	}
	const T& operator()( int i, int j, int k ) const 
	{
	    Assert(inRange(i,j,k-zoff));
	    return data( i+(ncx+1)*(j+(ncy+1)*(k-zoff)) );
	}

        const T& operator[]( int i ) const { return data(i); }
        T& operator[]( int i ) { return data(i); }

        iterator begin() { return data.begin(); }
        iterator end() { return data.end(); }

        const_iterator begin() const { return data.begin(); }
        const_iterator end() const { return data.end(); }

        size_type size() const { return data.size(); }
        size_type max_size() const { return data.size(); }
        bool empty() const { return data.empty(); }

        void swap ( gnctf& x ) { data.swap(x.data); }

	friend class Mesh_XYZ;
    };

// Vertex centered field
// Has a value at each vertex in each cell.

    template<class T>
    class vctf : private XYZ_Mapper,
                 public xm::Indexable< T, vctf<T> > {
	friend class Mesh_XYZ;
	friend class gvctf<T>;

	rtt_dsxx::Mat2<T> data;
        const Mesh_XYZ* mesh;
        FieldConstructor spm;

	vctf( const Mesh_XYZ *m )
	    : XYZ_Mapper( m->get_Mesh_DB() ),
	      data( m->get_ncells(), 8 ), mesh( m ), spm( 0 )
	{}

      public:
	typedef T value_type;
        typedef T& reference;
        typedef const T& const_reference;
        typedef typename rtt_dsxx::Mat2<T>::pointer pointer;
        typedef typename rtt_dsxx::Mat2<T>::const_pointer const_pointer;
        typedef typename rtt_dsxx::Mat2<T>::iterator iterator;
        typedef typename rtt_dsxx::Mat2<T>::const_iterator const_iterator;
        typedef typename rtt_dsxx::Mat2<T>::difference_type difference_type;
        typedef typename rtt_dsxx::Mat2<T>::size_type size_type;

	vctf( const FieldConstructor& spm_ )
	    : XYZ_Mapper( spm_->get_Mesh_DB() ),
	      data( spm_->get_ncells(), 8 ), mesh( spm_.bp() ), spm( spm_ )
	{}

	vctf( const vctf& v ) 
          : XYZ_Mapper( v.get_Mesh_DB() ),
            data( v.data ), mesh( v.mesh ), spm( v.spm ) {}

        ~vctf() {}

        const Mesh_XYZ& get_Mesh() const { return *mesh; }

        const FieldConstructor& get_FieldConstructor() const { return spm; }

	vctf& operator=( T x ) { data = x; return *this; }

        vctf& operator=( const vctf& v )
	{
            data = v.data;
            return *this;
	}

        template<class X>
        vctf& operator=( const xm::Xpr< T, X, vctf >& x )
        {
            return assign_from( x );
        }

    // i, j, k == global xyz cell indicies
    // v == vertex index
    // c == local cell index

	T& operator()( int c, int v ) { return data(c,v); }
	const T& operator()( int c, int v ) const { return data(c,v); }

	T& operator()( int i, int j, int k, int v )
	{
	    return data( local_cell_index(i,j,k), v );
	}
	const T& operator()( int i, int j, int k, int v ) const 
	{
	    return data( local_cell_index(i,j,k), v );
	}

        const T& operator[]( int i ) const { return data[i]; }
        T& operator[]( int i ) { return data[i]; }

        iterator begin() { return data.begin(); }
        iterator end() { return data.end(); }

        const_iterator begin() const { return data.begin(); }
        const_iterator end() const { return data.end(); }

        size_type size() const { return data.size(); }
        size_type max_size() const { return data.size(); }
        bool empty() const { return data.empty(); }

        void swap ( vctf& x ) { data.swap(x.data); }

        bool operator==( const vctf& x ) const;
        bool operator!=( const vctf& x ) const;
        bool operator<( const vctf& x ) const;
        bool operator>( const vctf& x ) const;
        bool operator<=( const vctf& x ) const;
        bool operator>=( const vctf& x ) const;
    };

// Guarded vertex centered field
// Has a value at each vertex in each cell.

    template<class T>
    class gvctf : private XYZ_Mapper,
                  public xm::Indexable< T, gvctf<T> > {
	rtt_dsxx::Mat4<T> data;
        const Mesh_XYZ* mesh;
        FieldConstructor spm;

      public:
        typedef typename rtt_dsxx::Mat4<T>::iterator iterator;
        typedef typename rtt_dsxx::Mat4<T>::const_iterator const_iterator;

	gvctf( const FieldConstructor& spm_ )
	    : XYZ_Mapper( spm_->get_Mesh_DB() ),
	      data( rtt_dsxx::Bounds( 0, 7 ),
                    rtt_dsxx::Bounds( 0, ncx - 1 ),
                    rtt_dsxx::Bounds( 0, ncy - 1 ),
                    rtt_dsxx::Bounds( zoff - 1, zoff + nczp ) ),
              mesh( spm_.bp() ), spm( spm_ )
	{}

        const Mesh_XYZ& get_Mesh() const { return *mesh; }

        gvctf<T>& operator=( const vctf<T>& c );
        void update_gvctf();

        gvctf<T>( const vctf<T>& v )
            : XYZ_Mapper( v.get_Mesh_DB() ),
              data( rtt_dsxx::Bounds( 0, 7 ),
                    rtt_dsxx::Bounds( 0, ncx - 1 ),
                    rtt_dsxx::Bounds( 0, ncy - 1 ),
                    rtt_dsxx::Bounds( zoff - 1, zoff + nczp ) ),
              mesh( v.mesh ), spm( v.spm )
        { *this = v; }

    // i, j, k == global xyz cell indicies
    // v == vertex index

    // Note that the order of the indexing is different than the
    // order of the data layout.

        const T& operator()( int i, int j, int k, int v ) const
        { return data(v,i,j,k); }
        T& operator()( int i, int j, int k, int v )
        { return data(v,i,j,k); }

        const T& operator[]( int i ) const { return data[i]; }
        T& operator[]( int i ) { return data[i]; }

        iterator begin() { return data.begin(); }
        iterator end() { return data.end(); }

        const_iterator begin() const { return data.begin(); }
        const_iterator end() const { return data.end(); }

        int size() const { return data.size(); }
    };

// Boundary specified field
// Has a value on each boundary cell face.

    template <class T>
    class bstf : private XYZ_Mapper,
                 public xm::Indexable<T, bstf<T> > {
        friend class Mesh_XYZ;

        fcdtf<T> data;
        const Mesh_XYZ* mesh;
        FieldConstructor spm;

        bstf ( const Mesh_XYZ *m )
            : XYZ_Mapper( m->get_Mesh_DB() ), data( m ), mesh( m ), spm( 0 )
	{}

        void next_element( int& i, int& j, int& k, int& f) const;
        void get_indexes
             ( int& i, int& j, int& k, int& f, const int index) const;

      public:
	typedef T value_type;
        typedef T& reference;
        typedef const T& const_reference;
        typedef typename fcdtf<T>::pointer pointer;
        typedef typename fcdtf<T>::const_pointer const_pointer;
        typedef typename fcdtf<T>::difference_type difference_type;
        typedef typename fcdtf<T>::size_type size_type;

        class iterator;
        class const_iterator :
          public std::iterator<std::forward_iterator_tag, value_type,
                               difference_type, const_pointer,
                               const_reference>
	{
          private:
            const T* p;
            int i, j, k, f;
            const bstf<T>* bfield;

          public:
            const_iterator( const T* const _p, const int _i, const int _j,
                            const int _k, const int _f,
                            const bstf<T>* const _bfield )
                : p( _p ), i( _i ), j( _j ), k( _k ), f( _f ),
                  bfield( _bfield )
            {}

            const_iterator( const bstf<T>::iterator& iter );

            const_iterator() : p(0), i(0), j(0), k(0), f(0), bfield(0) {}

            bool operator==( const const_iterator& iter ) const
            { return p == iter.p; }
            bool operator!=( const const_iterator& iter ) const
            { return p != iter.p; }
            const_iterator& operator++();
            const_iterator operator++( int dummy );
            const T& operator*() const { return *p; }
            const T* operator->() const { return p; }
            const_iterator& operator=( const const_iterator& iter );
	};

        class iterator :
          public std::iterator<std::forward_iterator_tag, value_type,
                               difference_type, pointer, reference>
	{
          friend class bstf<T>::const_iterator;

          private:
            T* p;
            int i, j, k, f;
            bstf<T>* bfield;

          public:
            iterator( T* const _p, const int _i, const int _j, const int _k,
                      const int _f, bstf<T>* const _bfield )
                : p( _p ), i( _i ), j( _j ), k( _k ), f( _f ),
                  bfield( _bfield )
            {}

            iterator() : p(0), i(0), j(0), k(0), f(0), bfield(0) {}

            bool operator==( const iterator& iter ) const
            { return p == iter.p; }
            bool operator!=( const iterator& iter ) const
            { return p != iter.p; }
            iterator& operator++();
            iterator operator++( int dummy );
            T& operator*() const { return *p; }
            T* operator->() const { return p; }
            iterator& operator=( const iterator& iter );
	};

        bstf( const FieldConstructor& spm_ )
            : XYZ_Mapper( spm_->get_Mesh_DB() ), data( spm_ ),
              mesh( spm_.bp() ), spm( spm_ )
	{}

	bstf( const bstf& b ) 
          : XYZ_Mapper( b.get_Mesh_DB() ),
            data( b.data ), mesh( b.mesh ), spm( b.spm ) {}

        ~bstf() {}

        const Mesh_XYZ& get_Mesh() const { return *mesh; }

        const FieldConstructor& get_FieldConstructor() const { return spm; }

        bstf& operator=( T x );

        bstf& operator=( const bstf& b )
	{
            data = b.data;
            return *this;
	}

        template<class X>
        bstf& operator=( const xm::Xpr<T, X, bstf>& x )
	{
            return assign_from( x );
	}

    // i, j, k == global xyz cell indicies
    // f == face index

        T& operator()( int i, int j, int k, int f );
        const T& operator()( int i, int j, int k, int f ) const;

        const T& operator[]( int index ) const;
        T& operator[]( int index );

        iterator begin();
        iterator end();

        const_iterator begin() const;
        const_iterator end() const;

        size_type size() const;
        size_type max_size() const { return size(); }
        bool empty() const { return data.empty(); }

        void swap ( bstf& x ) { data.swap(x.data); }

        bool operator==( const bstf& x ) const;
        bool operator!=( const bstf& x ) const;
        bool operator<( const bstf& x ) const;
        bool operator>( const bstf& x ) const;
        bool operator<=( const bstf& x ) const;
        bool operator>=( const bstf& x ) const;

        friend class bstf<T>::iterator;
        friend class bstf<T>::const_iterator;
    };

// Small vector class

    template<class T, int N>
    class tiny_vec : public xm::Indexable< T, tiny_vec<T,N> > {

        rtt_dsxx::Mat1<T> data;

      public:
	typedef T value_type;
        typedef T& reference;
        typedef const T& const_reference;
        typedef typename rtt_dsxx::Mat1<T>::pointer pointer;
        typedef typename rtt_dsxx::Mat1<T>::const_pointer const_pointer;
        typedef typename rtt_dsxx::Mat1<T>::iterator iterator;
        typedef typename rtt_dsxx::Mat1<T>::const_iterator const_iterator;
        typedef typename rtt_dsxx::Mat1<T>::difference_type
                                        difference_type;
        typedef typename rtt_dsxx::Mat1<T>::size_type size_type;
        typedef typename rtt_dsxx::Mat1<T>::reverse_iterator
                                        reverse_iterator;
        typedef typename rtt_dsxx::Mat1<T>::const_reverse_iterator
                                        const_reverse_iterator;

        tiny_vec() : data( N ) {}

        tiny_vec( const tiny_vec& v ) : data( v.data ) {}

        tiny_vec& operator=( T x ) { data = x; return *this; }

        static T dot( const tiny_vec<T,N>& x, const tiny_vec<T,N>& y)
        {
            T sum = 0;
            for ( int i = 0; i < N; ++i) {
                sum += x(i)*y(i);
            }
            return sum;
        }

        friend tiny_vec<T,N> operator+( const tiny_vec<T,N>& x,
                                        const tiny_vec<T,N>& y)
        {
            tiny_vec<T,N> results;
            for (int i = 0; i<N; ++i)
                results[i] = x[i] + y[i];
            return results;
        }

        friend tiny_vec<T,N> operator*( const tiny_vec<T,N>& x,
                                        const tiny_vec<T,N>& y)
        {
            tiny_vec<T,N> results;
            for (int i = 0; i<N; ++i)
                results[i] = x[i] * y[i];
            return results;
        }

        friend tiny_vec<T,N> operator-( const tiny_vec<T,N>& x,
                                        const tiny_vec<T,N>& y)
        {
            tiny_vec<T,N> results;
            for (int i = 0; i<N; ++i)
                results[i] = x[i] - y[i];
            return results;
        }

         friend tiny_vec<T,N> operator/( const tiny_vec<T,N>& x,
                                        const tiny_vec<T,N>& y)
        {
            tiny_vec<T,N> results;
            for (int i = 0; i<N; ++i)
                results[i] = x[i] / y[i];
            return results;
        }

        template<class X>
        tiny_vec& operator=( const xm::Xpr< T, X, tiny_vec >& x )
        {
            return assign_from( x );
        }

        tiny_vec& operator=( const tiny_vec& v ) { data = v.data; return *this; }

        T& operator() ( int i ) { return data(i); }
        const T& operator() ( int i ) const { return data(i); }

        const T& operator[] ( int i ) const { return data[i]; }
        T& operator[] ( int i ) { return data[i]; }

        iterator begin() { return data.begin(); }
        iterator end() { return data.end(); }

        const_iterator begin() const { return data.begin(); }
        const_iterator end() const { return data.end(); }

        reverse_iterator rbegin() { return data.rbegin(); }
        reverse_iterator rend() { return data.rend(); }

        const_reverse_iterator rbegin() const { return data.rbegin(); }
        const_reverse_iterator rend() const { return data.rend(); }

        size_type size() const { return data.size(); }
        size_type max_size() const { return size(); }
        bool empty() const { return data.empty(); }

        void swap ( tiny_vec& x ) { data.swap(x.data); }

        bool operator==( const tiny_vec& x ) const;
        bool operator!=( const tiny_vec& x ) const;
        bool operator<( const tiny_vec& x ) const;
        bool operator>( const tiny_vec& x ) const;
        bool operator<=( const tiny_vec& x ) const;
        bool operator>=( const tiny_vec& x ) const;
    };

  private:
    rtt_dsxx::Mat1<double> xc, yc, zc;
    rtt_dsxx::Mat1<double> xf, yf, zf;
    ccsf vc;
    rtt_dsxx::Mat1<double> xA, yA, zA;

    ccsf dX, dY, dZ;
    ccsf xC, yC, zC;
    fcdsf xF, yF, zF;
    fcdvsf face_norms;
    vec3 xhat, yhat, zhat;

  public:
    typedef XYZ_Mapper Coord_Mapper;

    Mesh_XYZ( const Mesh_DB& mdb_in );

    bool operator==( const Mesh_XYZ& m ) const { return this == &m; }
    bool operator!=( const Mesh_XYZ& m ) const { return !(this == &m); }

    const Mesh_DB &get_Mesh_DB() const { return mdb; }

    size_type get_ncx() const { return ncx; }
    size_type get_ncy() const { return ncy; }
    size_type get_ncz() const { return ncz; }
    size_type get_ncells() const { return ncp; }
    size_type get_total_ncells() const { return nct; }
    int get_ncxp() const { return ncx; }
    int get_xoff() const { return 0; }
    int get_ncyp() const { return ncy; }
    int get_yoff() const { return 0; }
    int get_nczp() const { return nczp; }
    int get_zoff() const { return zoff; }
    int get_goff() const { return goff; }

    void get_dx(ccsf& dx) const { dx = dX; }
    void get_dy(ccsf& dy) const { dy = dY; }
    void get_dz(ccsf& dz) const { dz = dZ; }

    const rtt_dsxx::Mat1<double>& get_xc() const { return xc; }
    const rtt_dsxx::Mat1<double>& get_yc() const { return yc; }
    const rtt_dsxx::Mat1<double>& get_zc() const { return zc; }

    void get_xloc(ccsf& xloc) const { xloc = xC; }
    void get_yloc(ccsf& yloc) const { yloc = yC; }
    void get_zloc(ccsf& zloc) const { zloc = zC; }

    const rtt_dsxx::Mat1<double>& get_xf() const { return xf; }
    const rtt_dsxx::Mat1<double>& get_yf() const { return yf; }
    const rtt_dsxx::Mat1<double>& get_zf() const { return zf; }

    void get_xloc(fcdsf& xloc) const { xloc = xF; }
    void get_yloc(fcdsf& yloc) const { yloc = yF; }
    void get_zloc(fcdsf& zloc) const { zloc = zF; }

    void get_face_normals(fcdvsf& fn) const { fn = face_norms; }
    void get_face_areas(fcdsf& fa) const;
    void get_face_lengths(fcdsf& fl) const;

    //    const rtt_dsxx::Mat1<double>& get_xA() const { return xA; }
    //    const rtt_dsxx::Mat1<double>& get_yA() const { return yA; }
    //    const rtt_dsxx::Mat1<double>& get_zA() const { return zA; }

    const ccsf& get_vc() const { return vc; }
    void get_cell_volumes(ccsf &vols) const { vols = vc; }
    void get_vertex_volumes(vcsf &vols) const;
    void get_node_volumes(ncsf &vols) const;

    template <class T1, class T2, class Op>
    static void scatter( fcdtf<T1>& to, const cctf<T2>& from, const Op& op );

    template <class T1, class T2, class Op>
    static void scatter( cctf<T1>& to, const fcdtf<T2>& from, const Op& op );

    template <class T1, class T2, class Op>
    static void scatter( fcdtf<T1>& to, const vctf<T2>& from, const Op& op );

    template <class T1, class T2, class Op>
    static void scatter( vctf<T1>& to, const fcdtf<T2>& from, const Op& op );

    template <class T1, class T2, class Op>
    static void scatter( nctf<T1>& to, const vctf<T2>& from, const Op& op );

    template <class T1, class T2, class Op>
    static void scatter( cctf<T1>& to, const vctf<T2>& from, const Op& op );

    template <class T1, class T2, class Op>
    static void gather( fcdtf<T1>& to, const cctf<T2>& from, const Op& op );

    template <class T1, class T2, class Op>
    static void gather( bstf<T1>& to, const fcdtf<T2>& from, const Op& op );

    template <class T1, class T2, class Op>
    static void gather( fcdtf<T1>& to, const bstf<T2>& from, const Op& op );

    template <class T1, class T2, class Op>
    static void gather( vctf<T1>& to, const nctf<T2>& from, const Op& op );

    template <class T1, class T2, class Op>
    static void gather( vctf<T1>& to, const cctf<T2>& from, const Op& op );

    template <class T>
    static void swap_faces( fcdtf<T>& to, const fcdtf<T>& from,
			    T bndryValue = T());

    template <class T>
    static T sum( const cctf<T>& from );

    template <class T>
    static T sum( const fcdtf<T>& from );

    template <class T>
    static T sum( const nctf<T>& from );

    template <class T>
    static T sum( const vctf<T>& from );

    template <class T>
    static T sum( const bstf<T>& from );

    template <class T>
    static T min( const cctf<T>& from );

    template <class T>
    static T min( const fcdtf<T>& from );

    template <class T>
    static T min( const nctf<T>& from );

    template <class T>
    static T min( const vctf<T>& from );

    template <class T>
    static T min( const bstf<T>& from );

    template <class T>
    static T max( const cctf<T>& from );

    template <class T>
    static T max( const fcdtf<T>& from );

    template <class T>
    static T max( const nctf<T>& from );

    template <class T>
    static T max( const vctf<T>& from );

    template <class T>
    static T max( const bstf<T>& from );

    class OpAssign {
      public:
        template <class T1, class T2>
        void operator() (T1& x, const T2& y) const { x = y; }
    };

    class OpAddAssign {
      public:
        template <class T1, class T2>
        void operator() (T1& x, const T2& y) const { x += y; }
    };

    class OpSubAssign {
      public:
        template <class T1, class T2>
        void operator() (T1& x, const T2& y) const { x -= y; }
    };

    class OpMultAssign {
      public:
        template <class T1, class T2>
        void operator() (T1& x, const T2& y) const { x *= y; }
    };

    class OpMinAssign {
      public:
        template <class T1, class T2>
        void operator() (T1& x, const T2& y) const { x = (x < y) ? x : y; }
    };
    
    class OpMaxAssign {
      public:
        template <class T1, class T2>
        void operator() (T1& x, const T2& y) const { x = (x > y) ? x : y; }
    };
};

template<class T>
void dump( const Mesh_XYZ::cctf<T>& data, char *name );

template<class T>
void dump( const Mesh_XYZ::fcdtf<T>& data, char *name );

template<class T>
void dump( const Mesh_XYZ::nctf<T>& data, char *name );

template<class T>
void dump( const Mesh_XYZ::vctf<T>& data, char *name );

#include "Mesh_XYZ.t.hh"

// Get the ContainerTraits for this class.
#include "Mesh_XYZTraits.hh"

#endif                          // __mesh_Mesh_XYZ_hh__

//---------------------------------------------------------------------------//
//                              end of mesh/Mesh_XYZ.hh
//---------------------------------------------------------------------------//
