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

#include "ds++/Mat.hh"
#include "ds++/SP.hh"

// using dsxx::Mat1;
// using dsxx::Mat2;
// using dsxx::Mat3;

#include "c4/NodeInfo.hh"

#include "xm/xm.hh"

struct XYZ_Mapper : public Mesh_DB, public C4::NodeInfo
{
    int nct;			// # of total cells in problem.
    int nxy;                    // # of cells in an x-y plane.
    int nczp;                   // # of cells in z this processor.
    int ncp;			// # of cells on this processor.

    int zoff;                   // z index offset of this processor's cells.
    int goff;			// global offset of this processor's cells.

// Methods

    XYZ_Mapper( const Mesh_DB& mdb );

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

    const Mesh_DB& get_Mesh_DB() const{ return *this; }

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

    template<class T> class gfcdtf;
    template<class T> class cctf;
    template<class T> class gcctf;
    template<class T> class vctf;
    template<class T, int N> class tiny_vec;

    typedef cctf<double> ccsf;
    typedef cctf<int> ccif;
    typedef gcctf<double> gccsf;
    typedef vctf<double> vcsf;
    typedef vctf<tiny_vec<double, 3> > vcvf;
    typedef tiny_vec<double, 3> vec3;

// Face centered discontinuous field
// Has a value on each face in each cell.

    template<class T>
    class fcdtf : private XYZ_Mapper,
                  public xm::Indexable< T, fcdtf<T> > {
	friend class Mesh_XYZ;
	friend class gfcdtf<T>;

	dsxx::Mat2<T> data;
        const Mesh_XYZ* mesh;

	fcdtf( const Mesh_XYZ *m )
	    : XYZ_Mapper( m->get_Mesh_DB() ),
	      data( m->get_ncp(), 6 ), mesh( m )
	{}

      public:
	typedef T value_type;
        typedef typename dsxx::Mat2<T>::iterator iterator;
        typedef typename dsxx::Mat2<T>::const_iterator const_iterator;

	fcdtf( const dsxx::SP<Mesh_XYZ>& m )
	    : XYZ_Mapper( m->get_Mesh_DB() ),
	      data( m->get_ncp(), 6 ), mesh( m.bp() )
	{}

        const Mesh_XYZ& get_Mesh() const { return *mesh; }

	fcdtf& operator=( T x ) { data = x; return *this; }

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

        int size() const { return data.size(); }
    };

    typedef fcdtf<double> fcdsf;
    typedef fcdtf<int> fcdif;
    typedef fcdtf<tiny_vec<double, 3> > fcdvf;

// Guarded face centered discontinuous field
// Has a value on each face in each cell.

    template<class T>
    class gfcdtf : private XYZ_Mapper,
                   public xm::Indexable< T, gfcdtf<T> > {
	dsxx::Mat4<T> data;
        const Mesh_XYZ* mesh;

      public:
        typedef typename dsxx::Mat4<T>::iterator iterator;
        typedef typename dsxx::Mat4<T>::const_iterator const_iterator;

	gfcdtf( const dsxx::SP<Mesh_XYZ>& m )
	    : XYZ_Mapper( m->get_Mesh_DB() ),
	      data( dsxx::Bounds( 0, 5 ),
                    dsxx::Bounds( 0, ncx - 1 ),
                    dsxx::Bounds( 0, ncy - 1 ),
                    dsxx::Bounds( zoff - 1, zoff + nczp ) ),
              mesh( m.bp() )
	{}

        const Mesh_XYZ& get_Mesh() const { return *mesh; }

        gfcdtf<T>& operator=( const fcdtf<T>& c );
        void update_gfcdtf();

        gfcdtf<T>( const fcdtf<T>& c )
            : XYZ_Mapper( c.get_Mesh_DB() ),
              data( dsxx::Bounds( 0, 5 ),
                    dsxx::Bounds( 0, ncx - 1 ),
                    dsxx::Bounds( 0, ncy - 1 ),
                    dsxx::Bounds( zoff - 1, zoff + nczp ) ),
              mesh( c.mesh )
        { *this = c; }

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
        dsxx::Mat1<T> data;
        const Mesh_XYZ* mesh;

	cctf( const Mesh_XYZ *m ) 
          : XYZ_Mapper( m->get_Mesh_DB() ),
            data( m->get_ncp() ), mesh( m )
        {}

      public:
	typedef T value_type;
        typedef typename dsxx::Mat1<T>::iterator iterator;
        typedef typename dsxx::Mat1<T>::const_iterator const_iterator;

	cctf( const dsxx::SP<Mesh_XYZ>& m ) 
          : XYZ_Mapper( m->get_Mesh_DB() ),
            data( m->get_ncp() ), mesh( m.bp() ) {}

        const Mesh_XYZ& get_Mesh() const { return *mesh; }

        cctf& operator=( T x )
        {
            data = x;
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

        int size() const { return data.size(); }

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
        dsxx::Mat3<T> data;
        const Mesh_XYZ* mesh;

      public:
        typedef typename dsxx::Mat3<T>::iterator iterator;
        typedef typename dsxx::Mat3<T>::const_iterator const_iterator;

        gcctf( const dsxx::SP<Mesh_XYZ>& m )
            : XYZ_Mapper( m->get_Mesh_DB() ),
              data( dsxx::Bounds( 0, ncx - 1 ),
                    dsxx::Bounds( 0, ncy - 1 ),
                    dsxx::Bounds( zoff - 1, zoff + nczp ) ),
              mesh( m.bp() )
        {}

        const Mesh_XYZ& get_Mesh() const { return *mesh; }

        gcctf<T>& operator=( const cctf<T>& c );
        void update_guard_cells();

        gcctf<T>( const cctf<T>& c )
            : XYZ_Mapper( c.get_Mesh_DB() ),
              data( dsxx::Bounds( 0, ncx - 1 ),
                    dsxx::Bounds( 0, ncy - 1 ),
                    dsxx::Bounds( zoff - 1, zoff + nczp ) ),
              mesh( c.mesh )
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

// Vertex centered field
// Has a value at each vertex in each cell.

    template<class T>
    class vctf : private XYZ_Mapper,
                 public xm::Indexable< T, vctf<T> > {
	friend class Mesh_XYZ;

	dsxx::Mat2<T> data;
        const Mesh_XYZ* mesh;

	vctf( const Mesh_XYZ *m )
	    : XYZ_Mapper( m->get_Mesh_DB() ),
	      data( m->get_ncp(), 8 ), mesh( m )
	{}

      public:
	typedef T value_type;
        typedef typename dsxx::Mat2<T>::iterator iterator;
        typedef typename dsxx::Mat2<T>::const_iterator const_iterator;

	vctf( const dsxx::SP<Mesh_XYZ>& m )
	    : XYZ_Mapper( m->get_Mesh_DB() ),
	      data( m->get_ncp(), 8 ), mesh( m.bp() )
	{}

        const Mesh_XYZ& get_Mesh() const { return *mesh; }

	vctf& operator=( T x ) { data = x; return *this; }

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

        bstf ( const Mesh_XYZ *m )
            : XYZ_Mapper( m->get_Mesh_DB() ), data( m ), mesh( m )
	{}

        void next_element( int& i, int& j, int& k, int& f) const;
        void get_indexes
             ( int& i, int& j, int& k, int& f, const int index) const;

      public:
        typedef T value_type;

        class iterator;
        class const_iterator
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

            const_iterator( const iterator& iter );

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

        class iterator
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

        bstf( const dsxx::SP<Mesh_XYZ>& m )
            : XYZ_Mapper( m->get_Mesh_DB() ), data( m ), mesh( m.bp() )
	{}

        const Mesh_XYZ& get_Mesh() const { return *mesh; }

        bstf& operator=( T x );

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

        int size() const;

        friend class bstf<T>::iterator;
        friend class bstf<T>::const_iterator;
    };

    typedef bstf<double> bssf;
    typedef bstf<int> bsif;

// Small vector class

    template<class T, int N>
    class tiny_vec : public xm::Indexable< T, tiny_vec<T,N> > {

        dsxx::Mat1<T> data;

      public:
        typedef T value_type;
        typedef typename dsxx::Mat1<T>::iterator iterator;
        typedef typename dsxx::Mat1<T>::const_iterator const_iterator;

        tiny_vec() : data( N ) {}

        tiny_vec& operator=( T x ) {data = x; return *this;}

        static T dot( const tiny_vec<T,N>& x, const tiny_vec<T,N>& y)
        {
            T sum = 0;
            for ( int i = 0; i < N; ++i) {
                sum += x(i)*y(i);
            }
            return sum;
        }

        template<class X>
        tiny_vec& operator=( const xm::Xpr< T, X, tiny_vec >& x )
        {
            return assign_from( x );
        }

        T& operator() ( int i ) { return data(i); }
        const T& operator() ( int i ) const { return data(i); }

        const T& operator[] ( int i ) const { return data[i]; }
        T& operator[] ( int i ) { return data[i]; }

        iterator begin() { return data.begin(); }
        iterator end() { return data.end(); }

        const_iterator begin() const { return data.begin(); }
        const_iterator end() const { return data.end(); }

        int size() const { return data.size(); }
    };

  private:
    dsxx::Mat1<double> xc, yc, zc;
    dsxx::Mat1<double> xf, yf, zf;
    double       dx, dy, dz;
    ccsf vc;
    dsxx::Mat1<double> xA, yA, zA;

    ccsf dX, dY, dZ;
    ccsf xC, yC, zC;
    fcdsf xF, yF, zF;
    fcdvf face_norms;
    vec3 xhat, yhat, zhat;

    int diags[7];

  public:
    typedef XYZ_Mapper Coord_Mapper;

    Mesh_XYZ( const Mesh_DB& mdb );

    bool operator==( const Mesh_XYZ& m ) const { return this == &m; }

    const Mesh_DB& get_Mesh_DB() const { return XYZ_Mapper::get_Mesh_DB(); }

    int get_ncx() const { return ncx; }
    int get_ncy() const { return ncy; }
    int get_ncz() const { return ncz; }
    int get_ncp() const { return ncp; }
    int get_nct() const { return nct; }
    int get_nczp() const { return nczp; }
    int get_zoff() const { return zoff; }
    int get_goff() const { return goff; }

    double get_dx() const { return dx; }
    double get_dy() const { return dy; }
    double get_dz() const { return dz; }

    void get_dx(ccsf& dx) const { dx = dX; }
    void get_dy(ccsf& dy) const { dy = dY; }
    void get_dz(ccsf& dz) const { dz = dZ; }

    const dsxx::Mat1<double>& get_xc() const { return xc; }
    const dsxx::Mat1<double>& get_yc() const { return yc; }
    const dsxx::Mat1<double>& get_zc() const { return zc; }

    void get_xloc(ccsf& xloc) const { xloc = xC; }
    void get_yloc(ccsf& yloc) const { yloc = yC; }
    void get_zloc(ccsf& zloc) const { zloc = zC; }

    const dsxx::Mat1<double>& get_xf() const { return xf; }
    const dsxx::Mat1<double>& get_yf() const { return yf; }
    const dsxx::Mat1<double>& get_zf() const { return zf; }

    void get_xloc(fcdsf& xloc) const { xloc = xF; }
    void get_yloc(fcdsf& yloc) const { yloc = yF; }
    void get_zloc(fcdsf& zloc) const { zloc = zF; }

    const fcdvf& get_fn() const { return face_norms; }
    void get_face_areas(fcdsf& fa) const;
    void get_face_lengths(fcdsf& fl) const;

    const dsxx::Mat1<double>& get_xA() const { return xA; }
    const dsxx::Mat1<double>& get_yA() const { return yA; }
    const dsxx::Mat1<double>& get_zA() const { return zA; }

    const ccsf& get_vc() const { return vc; }
    void get_cell_volumes(ccsf &vols) const { vols = vc; }

    const int *get_diag_offsets() const { return diags; }

    template <class T1, class T2, class Op>
    static void scatter( fcdtf<T1>& to, const cctf<T2>& from, const Op& op );

    template <class T1, class T2, class Op>
    static void scatter( cctf<T1>& to, const fcdtf<T2>& from, const Op& op );

    template <class T1, class T2, class Op>
    static void scatter( fcdtf<T1>& to, const vctf<T2>& from, const Op& op );

    template <class T1, class T2, class Op>
    static void gather( fcdtf<T1>& to, const cctf<T2>& from, const Op& op );

    template <class T1, class T2, class Op>
    static void gather( bstf<T1>& to, const fcdtf<T2>& from, const Op& op );

    template <class T1, class T2, class Op>
    static void gather( fcdtf<T1>& to, const bstf<T2>& from, const Op& op );

    template <class T>
    static void swap( fcdtf<T>& to, const fcdtf<T>& from );

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
};

template<class T>
void dump( const Mesh_XYZ::cctf<T>& data, char *name );

template<class T>
void dump( const Mesh_XYZ::fcdtf<T>& data, char *name );

#include "Mesh_XYZ.t.cc"

#endif                          // __mesh_Mesh_XYZ_hh__

//---------------------------------------------------------------------------//
//                              end of mesh/Mesh_XYZ.hh
//---------------------------------------------------------------------------//
