//----------------------------------*-C++-*----------------------------------//
// PoomaMesh_XYZ.hh
// Julian C. Cummings
// Mon Sep 21 22:08:34 1998
//---------------------------------------------------------------------------//
// @> A 3-d cartesian structured mesh facility based on POOMA r1.
//---------------------------------------------------------------------------//

#ifndef __mesh_PoomaMesh_XYZ_hh__
#define __mesh_PoomaMesh_XYZ_hh__

// include files

#include "Index/Index.h"
#include "Index/NDIndex.h"
#include "AppTypes/Vektor.h"
#include "FieldLayout/CenteredFieldLayout.h"
#include "Meshes/UniformCartesian.h"
#include "Meshes/Cartesian.h"
#include "Meshes/Centering.h"
#include "Field/Field.h"
#include "Field/LField.h"
#include "Field/Assign.h"
#include "Field/GuardCellSizes.h"
#include "Field/BCond.h"

#include "ds++/config.hh"
#include "ds++/SP.hh"

#include <iterator.h>

//===========================================================================//
// class PoomaMesh_XYZ - A 3-d cartesian mesh class based on POOMA r1
// This is a 3-d cartesisan structured mesh.  The main purpose of having this
// class is in order for it to be instantiated by various test articles
// throughout the Draco system.  It could also be useful to Draco clients as
// a point of reference during the construction of their own mesh classes.
//===========================================================================//

template <class Mesh>
class PoomaMesh_XYZ
{
  public:

// forward declarations
    template <class T> class cctf;
    template <class T> class fcdtf;
    template <class T> class fcdtf_iterator;
    template <class T> class fcdtf_const_iterator;
    template <class T> class bstf;
    template <class T> class bstf_iterator;
    template <class T> class bstf_const_iterator;

// typedefs
    typedef Mesh Mesh_t;
    typedef CenteredFieldLayout<3,Mesh_t,Cell> Layout_t;
    typedef cctf<double> ccsf;
    typedef cctf<int> ccif;
    typedef cctf< Vektor<double,3> > ccvf;  // need for storing mesh info
    typedef fcdtf<double> fcdsf;
    typedef fcdtf<int> fcdif;
    typedef bstf<double> bssf;
    typedef bstf<int> bsif;

  // These refer to PETE operators, which are in global namespace.
    typedef ::OpAssign         OpAssign;
    typedef ::OpAddAssign      OpAddAssign;
    typedef ::OpSubtractAssign OpSubAssign;
    typedef ::OpMultiplyAssign OpMultAssign;

// Face centered discontinuous field
// Has a value on each face in each cell.

    template <class T>
    class fcdtf : public Field<Vektor<T,6>,3,Mesh_t,Cell>
    {
      public:

      // typedefs
	typedef T value_type;
        typedef fcdtf_iterator<T> iterator;
        typedef fcdtf_const_iterator<T> const_iterator;
        typedef Field<Vektor<T,6>,3,Mesh_t,Cell> BaseField_t;
        typedef GuardCellSizes<3> GC_t;
        typedef BConds<Vektor<T,6>,3,Mesh_t,Cell> BC_Container_t;
        typedef ZeroFace<Vektor<T,6>,3,Mesh_t,Cell> BC_t;
        typedef PoomaMesh_XYZ<Mesh_t> MT_t;
        typedef dsxx::SP<MT_t> SPM_t;

      // constructors
	fcdtf(const SPM_t& spm)
          : BaseField_t( const_cast<Mesh_t&>(spm->get_Mesh()),
                         const_cast<Layout_t&>(spm->get_Layout()),
                         GC_t(1), BC_Container_t() ), spm_m(spm)
        {
          // set boundary conditions
          BC_Container_t bconds = getBConds();
          for (int f=0; f<6; ++f)
            bconds[f] = new BC_t(f);
          // compute number of local faces
          computeSize();
        }
	fcdtf(const MT_t& m)
	  : BaseField_t( const_cast<Mesh_t&>(m.get_Mesh()),
                         const_cast<Layout_t&>(m.get_Layout()),
                         GC_t(1), BC_Container_t() ),
            spm_m(const_cast<MT_t*>(&m))
        {
          // set boundary conditions
          BC_Container_t bconds = getBConds();
          for (int f=0; f<6; ++f)
            bconds[f] = new BC_t(f);
          // compute number of local faces
          computeSize();
        }
        fcdtf(const fcdtf<T>& f)
          : BaseField_t( const_cast<Mesh_t&>(f.get_Mesh().get_Mesh()),
                         const_cast<Layout_t&>(f.get_Mesh().get_Layout()),
                         GC_t(1), f.getBConds() ),
            spm_m( f.get_SP_Mesh() ), size_m(f.size_m) {}

      // Assignment from scalar
        const fcdtf<T>& operator=(T x)
        {
          Vektor<T,6> v = x;
          assign(*this,v);
          return *this;
        }
      // Assignment from another field
        const fcdtf<T>& operator=(const fcdtf<T>& x)
        {
          assign(*this,x);
          return *this;
        }
      // Assignment from base field
        const fcdtf<T>& operator=(const BaseField_t& x)
        {
          assign(*this,x);
          return *this;
        }
      // Assignment from an expression
        template <class B>
        const fcdtf<T>& operator=(const PETE_Expr<B>& x)
        {
          assign(*this,x);
          return *this;
        }

      // face accessors
        // i, j, k == global xyz cell indices
        // f == face index
      /*
        cctf<T> operator()(int f) const
        {
          cctf<T> cf(get_Mesh());
          cctf<T>::iterator cfit, cfend = cf.end();
          const BaseField_t& bf = dynamic_cast<const BaseField_t&>(*this);
	  BaseField_t::iterator bfit = bf.begin();
          for (cfit = cf.begin(); cfit != cfend; ++cfit, ++bfit)
            *cfit = (*bfit)(f);
          return cf;
        }
	*/
	value_type& operator()(int i, int j, int k, int f)
	{
          NDIndex<3> loc(Index(i,i),Index(j,j),Index(k,k));
	  return localElement(loc)(f);
	}
	const value_type& operator()(int i, int j, int k, int f) const
	{
          NDIndex<3> loc(Index(i,i),Index(j,j),Index(k,k));
	  return localElement(loc)(f);
	}

      // iterators
        iterator begin()
        {
          return iterator( dynamic_cast<BaseField_t&>(*this) );
        }
        iterator end()
        {
          return iterator( BaseField_t::end() );
        }
        const_iterator begin() const
        {
          return const_iterator( dynamic_cast<const BaseField_t&>(*this) );
        }
        const_iterator end() const
        {
          return const_iterator( BaseField_t::end() );
        }

      // accessors
        const SPM_t& get_SP_Mesh(void) const
        {
          return spm_m;
        }
        const MT_t& get_Mesh(void) const
        {
          return *spm_m;
        }
        int size(void) const { return size_m; }

      private:

      // methods
        void computeSize(void)
        {
          // compute local size of field
          size_m = 0;
          // loop over the lfields and add up elements in each one
          const_iterator_if lfi, lfend = end_if();
          for (lfi = begin_if(); lfi != lfend; ++lfi) {
            LField<Vektor<T,6>,3>& lf(*((*lfi).second));
            size_m += lf.size(0) * lf.size(1) * lf.size(2);
          }
          size_m *= 6;  // multiply by number of faces
        }

      // data
        SPM_t spm_m;
        int size_m;
    };

// Face centered discontinuous field iterators

    template <class T>
    class fcdtf_iterator : public forward_iterator<T,ptrdiff_t>
    {
      friend class fcdtf_const_iterator<T>;

      public:

      // typedefs
        typedef typename fcdtf<T>::BaseField_t BaseField_t;
        typedef typename BaseField_t::iterator BaseField_iterator;

      // constructors
        fcdtf_iterator() : face_m(0) {}
        fcdtf_iterator(BaseField_t& bf)
          : bfi_m(bf.begin()), face_m(0) {}
        fcdtf_iterator(const BaseField_iterator& bfiter)
          : bfi_m(bfiter), face_m(0) {}
        fcdtf_iterator(const fcdtf_iterator<T>& iter)
          : bfi_m(iter.bfi_m), face_m(iter.face_m) {}

      // assignment
        const fcdtf_iterator<T>& operator=(const fcdtf_iterator<T>& x)
        {
          bfi_m = x.bfi_m;
          face_m = x.face_m;
          return *this;
        }
        const fcdtf_iterator<T>& operator=(const bstf_iterator<T>& x)
        {
          // set BaseField iterator and face equal
          FieldLoc<3> loc;
          x.bfi_m.GetCurrentLocation(loc);
          bfi_m.SetCurrentLocation(loc);
          face_m = x.face_m;
          return *this;
        }
        const fcdtf_iterator<T>& operator=(const bstf_const_iterator<T>& x)
        {
          // set BaseField iterator and face equal
          FieldLoc<3> loc;
          x.bfi_m.GetCurrentLocation(loc);
          bfi_m.SetCurrentLocation(loc);
          face_m = x.face_m;
          return *this;
        }

      // pre-increment
        fcdtf_iterator<T>& operator++()
        {
          ++face_m;
          if (face_m>5) nextCell();
          return *this;
        }
      // post-increment
        fcdtf_iterator<T> operator++(int)
        {
          fcdtf_iterator<T> iter(*this);
          ++face_m;
          if (face_m>5) nextCell();
          return iter;
        }
      // dereference
        T& operator*() const { return (*bfi_m)(face_m); }

      // comparisons
        bool operator==(const fcdtf_iterator<T>& rhs) const
        {
          return (bfi_m == rhs.bfi_m && face_m == rhs.face_m);
        }
        bool operator!=(const fcdtf_iterator<T>& rhs) const
        {
          return !(*this == rhs);
        }

      private:

      // methods
        void nextCell()
        {
          ++bfi_m;
          face_m = 0;
        }

      // data
        BaseField_iterator bfi_m;
        int face_m;
    };

    template <class T>
    class fcdtf_const_iterator : public input_iterator<T,ptrdiff_t>
    {
      public:

      // typedefs
        typedef typename fcdtf<T>::BaseField_t BaseField_t;
        typedef typename BaseField_t::iterator BaseField_iterator;

      // constructors
        fcdtf_const_iterator() : face_m(0) {}
        fcdtf_const_iterator(const BaseField_t& bf)
          : bfi_m(bf.begin()), face_m(0) {}
        fcdtf_const_iterator(const BaseField_iterator& bfiter)
          : bfi_m(bfiter), face_m(0) {}
        fcdtf_const_iterator(const fcdtf_const_iterator<T>& iter)
          : bfi_m(iter.bfi_m), face_m(iter.face_m) {}
        fcdtf_const_iterator(const fcdtf_iterator<T>& iter)
          : bfi_m(iter.bfi_m), face_m(iter.face_m) {}

      // assignment
        const fcdtf_const_iterator<T>&
        operator=(const fcdtf_const_iterator<T>& x)
        {
          bfi_m = x.bfi_m;
          face_m = x.face_m;
          return *this;
        }
        const fcdtf_const_iterator<T>&
        operator=(const fcdtf_iterator<T>& x)
        {
          bfi_m = x.bfi_m;
          face_m = x.face_m;
          return *this;
        }
        const fcdtf_const_iterator<T>&
        operator=(const bstf_iterator<T>& x)
        {
          // set BaseField iterator and face equal
          FieldLoc<3> loc;
          x.bfi_m.GetCurrentLocation(loc);
          bfi_m.SetCurrentLocation(loc);
          face_m = x.face_m;
          return *this;
        }
        const fcdtf_const_iterator<T>&
        operator=(const bstf_const_iterator<T>& x)
        {
          // set BaseField iterator and face equal
          FieldLoc<3> loc;
          x.bfi_m.GetCurrentLocation(loc);
          bfi_m.SetCurrentLocation(loc);
          face_m = x.face_m;
          return *this;
        }

      // pre-increment
        fcdtf_const_iterator<T>& operator++()
        {
          ++face_m;
          if (face_m>5) nextCell();
          return *this;
        }
      // post-increment
        fcdtf_const_iterator<T> operator++(int)
        {
          fcdtf_const_iterator<T> iter(*this);
          ++face_m;
          if (face_m>5) nextCell();
          return iter;
        }
      // dereference
        const T& operator*() const { return (*bfi_m)(face_m); }

      // comparisons
        bool operator==(const fcdtf_const_iterator<T>& rhs) const
        {
          return (bfi_m == rhs.bfi_m && face_m == rhs.face_m);
        }
        bool operator!=(const fcdtf_const_iterator<T>& rhs) const
        {
          return !(*this == rhs);
        }

      private:

      // methods
        void nextCell()
        {
          ++bfi_m;
          face_m = 0;
        }

      // data
        BaseField_iterator bfi_m;
        int face_m;
    };

// Cell centered field
// Has a value in each cell.

    template <class T>
    class cctf : public Field<T,3,Mesh_t,Cell>
    {
      public:

      // typedefs
	typedef T value_type;
        typedef Field<T,3,Mesh_t,Cell> BaseField_t;
        // use inherited BaseField iterators (don't have a const BFI yet)
        typedef typename BaseField_t::iterator iterator;
        typedef typename BaseField_t::iterator const_iterator;
        typedef GuardCellSizes<3> GC_t;
        typedef BConds<T,3,Mesh_t,Cell> BC_Container_t;
        typedef ZeroFace<T,3,Mesh_t,Cell> BC_t;
        typedef PoomaMesh_XYZ<Mesh_t> MT_t;
        typedef dsxx::SP<MT_t> SPM_t;

      // constructors
	cctf(const SPM_t& spm)
          : BaseField_t( const_cast<Mesh_t&>(spm->get_Mesh()),
                         const_cast<Layout_t&>(spm->get_Layout()),
                         GC_t(1), BC_Container_t() ), spm_m(spm)
        {
          // set boundary conditions
          BC_Container_t bconds = getBConds();
          for (int f=0; f<6; ++f)
            bconds[f] = new BC_t(f);
          // compute number of local cells
          computeSize();
        }
	cctf(const MT_t& m)
          : BaseField_t( const_cast<Mesh_t&>(m.get_Mesh()),
                         const_cast<Layout_t&>(m.get_Layout()),
                         GC_t(1), BC_Container_t() ),
            spm_m(const_cast<MT_t*>(&m))
        {
          // set boundary conditions
          BC_Container_t bconds = getBConds();
          for (int f=0; f<6; ++f)
            bconds[f] = new BC_t(f);
          // compute number of local cells
          computeSize();
        }
        cctf(const cctf<T>& f)
          : BaseField_t( const_cast<Mesh_t&>(f.get_Mesh().get_Mesh()),
                         const_cast<Layout_t&>(f.get_Mesh().get_Layout()),
                         GC_t(1), f.getBConds() ),
            spm_m( f.get_SP_Mesh() ), size_m(f.size_m) {}

      // Assignment from scalar
        const cctf<T>& operator=(T x)
        {
          assign(*this,x);
          return *this;
        }
      // Assignment from another field
        const cctf<T>& operator=(const cctf<T>& x)
        {
          assign(*this,x);
          return *this;
        }
      // Assignment from base field
        const cctf<T>& operator=(const BaseField_t& x)
        {
          assign(*this,x);
          return *this;
        }
      // Assignment from an expression
        template <class B>
        const cctf<T>& operator=(const PETE_Expr<B>& x)
        {
          assign(*this,x);
          return *this;
        }

      // element accessors
	value_type& operator()(int i, int j, int k)
	{
            NDIndex<3> loc(Index(i,i),Index(j,j),Index(k,k));
            return localElement(loc);
	}
	const value_type& operator()(int i, int j, int k) const
	{
            NDIndex<3> loc(Index(i,i),Index(j,j),Index(k,k));
            return localElement(loc);
	}

      // accessors
        const SPM_t& get_SP_Mesh(void) const
        {
          return spm_m;
        }
        const MT_t& get_Mesh(void) const
        {
          return *spm_m;
        }
        int size(void) const { return size_m; }

      private:

      // methods
        void computeSize(void)
        {
          // compute local size of field
          size_m = 0;
          // loop over the lfields and add up elements in each one
          const_iterator_if lfi, lfend = end_if();
          for (lfi = begin_if(); lfi != lfend; ++lfi) {
            LField<T,3>& lf(*((*lfi).second));
            size_m += lf.size(0) * lf.size(1) * lf.size(2);
          }
        }

      // data
        SPM_t spm_m;
        int size_m;
    };

// Boundary specified field
// Has a value on each boundary face.

    template <class T>
    class bstf : public Field<Vektor<T,6>,3,Mesh_t,Cell>
    {
      friend class bstf_iterator<T>;
      friend class bstf_const_iterator<T>;

      public:

      // typedefs
	typedef T value_type;
        typedef bstf_iterator<T> iterator;
        typedef bstf_const_iterator<T> const_iterator;
        typedef Field<Vektor<T,6>,3,Mesh_t,Cell> BaseField_t;
        typedef PoomaMesh_XYZ<Mesh_t> MT_t;
        typedef dsxx::SP<MT_t> SPM_t;

      // constructors
	bstf(const SPM_t& spm)
	  : BaseField_t( const_cast<Mesh_t&>(spm->get_Mesh()),
                         const_cast<Layout_t&>(spm->get_Layout()) ), spm_m(spm)
        {
          // store global mesh sizes
          ncx_m = spm->get_Mesh().gridSizes[0] - 1;
          ncy_m = spm->get_Mesh().gridSizes[1] - 1;
          ncz_m = spm->get_Mesh().gridSizes[2] - 1;
          // compute local number of boundary faces
          computeSize();
        }
	bstf(const MT_t& m)
	  : BaseField_t( const_cast<Mesh_t&>(m.get_Mesh()),
                         const_cast<Layout_t&>(m.get_Layout()) ),
            spm_m(const_cast<MT_t*>(&m))
        {
          // store global mesh sizes
          ncx_m = m.get_Mesh().gridSizes[0] - 1;
          ncy_m = m.get_Mesh().gridSizes[1] - 1;
          ncz_m = m.get_Mesh().gridSizes[2] - 1;
          // compute local number of boundary faces
          computeSize();
        }
        bstf(const bstf<T>& f)
          : BaseField_t( const_cast<Mesh_t&>(f.get_Mesh().get_Mesh()),
                         const_cast<Layout_t&>(f.get_Mesh().get_Layout()) ),
            spm_m( f.get_SP_Mesh() ), size_m(f.size_m),
            ncx_m(f.ncx_m), ncy_m(f.ncy_m), ncz_m(f.ncz_m) {}

      // Assignment from scalar
        const bstf<T>& operator=(T x) {
          Uncompress();
          iterator it, iend = end();
          for (it = begin(); it != iend; ++it)
            *it = x;
          return *this;
        }
      // Assignment from another field
        const bstf<T>& operator=(const bstf<T>& x) {
          assign(*this,x);
          return *this;
        }
      // Assignment from an expression
        template <class B>
        const bstf<T>& operator=(const PETE_Expr<B>& x) {
          assign(*this,x);
          return *this;
        }

      // face accessors
        // i, j, k == global xyz cell indices
        // f == face index
        // must check that indices refer to a boundary face
	value_type& operator()(int i, int j, int k, int f)
	{
          bool OK = checkIndices(i,j,k,f);
          PInsist(OK, "bstf<T>::operator() detects bad indices!!");
          NDIndex<3> loc(Index(i,i),Index(j,j),Index(k,k));
	  return localElement(loc)(f);
	}
	const value_type& operator()(int i, int j, int k, int f) const
	{
          bool OK = checkIndices(i,j,k,f);
          PInsist(OK, "bstf<T>::operator() detects bad indices!!");
          NDIndex<3> loc(Index(i,i),Index(j,j),Index(k,k));
	  return localElement(loc)(f);
	}

      // iterators
        iterator begin()
        {
          return iterator( dynamic_cast<BaseField_t&>(*this) );
        }
        iterator end()
        {
          return iterator( BaseField_t::end() );
        }
        const_iterator begin() const
        {
          return const_iterator( dynamic_cast<const BaseField_t&>(*this) );
        }
        const_iterator end() const
        {
          return const_iterator( BaseField_t::end() );
        }

      // accessors
        const SPM_t& get_SP_Mesh(void) const
        {
          return spm_m;
        }
        const MT_t& get_Mesh(void) const
        {
          return *spm_m;
        }
        int size(void) const { return size_m; }

      private:

      // methods
        void computeSize(void)
        {
          // compute local size of field
          size_m = 0;
          // iterate over local elements and faces,
          // and count how many are on the boundary
          int loc[3];
          int f;
          bool OK;
          typename BaseField_t::iterator bfi, bfend = BaseField_t::end();
          for (bfi = BaseField_t::begin(); bfi != bfend; ++bfi) {
            bfi.GetCurrentLocation(loc);
            for (f=0; f<6; ++f) {
              OK = checkIndices(loc[0], loc[1], loc[2], f);
              if (OK) ++size_m;
            }
          }
        }
        bool checkIndices(int i, int j, int k, int f) const
        {
          bool result;
          if (i == 0 && f == 0) {
            // left boundary face
            result = true;
          }
          else if (i == ncx_m-1 && f == 1) {
            // right boundary face
            result = true;
          }
          else if (j == 0 && f == 2) {
            // bottom boundary face
            result = true;
          }
          else if (j == ncy_m-1 && f == 3) {
            // top boundary face
            result = true;
          }
          else if (k == 0 && f == 4) {
            // rear boundary face
            result = true;
          }
          else if (k == ncz_m-1 && f == 5) {
            // front boundary face
            result = true;
          }
          else {
            // not a boundary face
            result = false;
          }
          return result;
        }

      // data
        SPM_t spm_m;
        int size_m;
        int ncx_m, ncy_m, ncz_m;
    };

// Boundary specified field iterators

    template <class T>
    class bstf_iterator : public forward_iterator<T,ptrdiff_t>
    {
      friend class bstf_const_iterator<T>;
      friend class fcdtf_iterator<T>;
      friend class fcdtf_const_iterator<T>;

      public:

      // typedefs
        typedef typename bstf<T>::BaseField_t BaseField_t;
        typedef typename BaseField_t::iterator BaseField_iterator;

      // constructors
        bstf_iterator() {}
        bstf_iterator(BaseField_t& bf)
          : bfi_m(bf.begin()), bfend_m(bf.end()), face_m(0)
        {
          // make sure we are at the first valid boundary face
          int loc[3];
          bool OK;
          bfi_m.GetCurrentLocation(loc);
          bstf<T>& bsf = dynamic_cast<bstf<T>&>(bf);
          OK = bsf.checkIndices(loc[0], loc[1], loc[2], face_m);
          if (!OK) ++(*this);  // if not, advance to next valid boundary face
        }
        bstf_iterator(const BaseField_iterator& bfiter)
          : bfi_m(bfiter), bfend_m(bfiter), face_m(0) {}
        bstf_iterator(const bstf_iterator<T>& iter)
          : bfi_m(iter.bfi_m), bfend_m(iter.bfend_m),
            face_m(iter.face_m) {}

      // assignment
        const bstf_iterator<T>& operator=(const bstf_iterator<T>& x)
        {
          bfi_m = x.bfi_m;
          bfend_m = x.bfend_m;
          face_m = x.face_m;
          return *this;
        }

      // pre-increment
        bstf_iterator<T>& operator++()
        {
          bstf<T>& bsf = dynamic_cast<bstf<T>&>(bfi_m.GetBareField());
          // advance to next boundary face
          int loc[3];
          bool OK;
          do {
            ++face_m;
            if (face_m>5) nextCell();
            if (bfi_m != bfend_m) {
              bfi_m.GetCurrentLocation(loc);
              OK = bsf.checkIndices(loc[0], loc[1], loc[2], face_m);
            }
            else {
              // at end, so just mark it OK
              OK = true;
            }
          } while (!OK);
          return *this;
        }
      // post-increment
        bstf_iterator<T> operator++(int)
        {
          bstf<T>& bsf = dynamic_cast<bstf<T>&>(bfi_m.GetBareField());
          // advance to next boundary face
          bstf_iterator<T> iter(*this);
          int loc[3];
          bool OK;
          do {
            ++face_m;
            if (face_m>5) nextCell();
            if (bfi_m != bfend_m) {
              bfi_m.GetCurrentLocation(loc);
              OK = bsf_m->checkIndices(loc[0], loc[1], loc[2], face_m);
            }
            else {
              // at end, so just mark it OK
              OK = true;
            }
          } while (!OK);
          return iter;
        }
      // dereference
        T& operator*() const { return (*bfi_m)(face_m); }

      // comparisons
        bool operator==(const bstf_iterator<T>& rhs) const
        {
          return (bfi_m == rhs.bfi_m && face_m == rhs.face_m);
        }
        bool operator!=(const bstf_iterator<T>& rhs) const
        {
          return !(*this == rhs);
        }

      private:

      // methods
        void nextCell()
        {
          ++bfi_m;
          face_m = 0;
        }

      // data
        BaseField_iterator bfi_m;
        BaseField_iterator bfend_m;
        int face_m;
    };

    template <class T>
    class bstf_const_iterator : public input_iterator<T,ptrdiff_t>
    {
      friend class fcdtf_iterator<T>;
      friend class fcdtf_const_iterator<T>;

      public:

      // typedefs
        typedef typename bstf<T>::BaseField_t BaseField_t;
        typedef typename BaseField_t::iterator BaseField_iterator;

      // constructors
        bstf_const_iterator() {}
        bstf_const_iterator(const BaseField_t& bf)
          : bfi_m(bf.begin()), bfend_m(bf.end()), face_m(0)
        {
          // make sure we are at the first valid boundary face
          int loc[3];
          bool OK;
          bfi_m.GetCurrentLocation(loc);
          const bstf<T>& bsf = dynamic_cast<const bstf<T>&>(bf);
          OK = bsf.checkIndices(loc[0], loc[1], loc[2], face_m);
          if (!OK) ++(*this);  // if not, advance to next valid boundary face
        }
        bstf_const_iterator(const BaseField_iterator& bfiter)
          : bfi_m(bfiter), bfend_m(bfiter), face_m(0) {}
        bstf_const_iterator(const bstf_const_iterator<T>& iter)
          : bfi_m(iter.bfi_m), bfend_m(iter.bfend_m),
            face_m(iter.face_m) {}
        bstf_const_iterator(const bstf_iterator<T>& iter)
          : bfi_m(iter.bfi_m), bfend_m(iter.bfend_m),
            face_m(iter.face_m) {}

      // assignment
        const bstf_const_iterator<T>&
        operator=(const bstf_const_iterator<T>& x)
        {
          bfi_m = x.bfi_m;
          bfend_m = x.bfend_m;
          face_m = x.face_m;
          return *this;
        }
        const bstf_const_iterator<T>&
        operator=(const bstf_iterator<T>& x)
        {
          bfi_m = x.bfi_m;
          bfend_m = x.bfend_m;
          face_m = x.face_m;
          return *this;
        }

      // pre-increment
        bstf_const_iterator<T>& operator++()
        {
          bstf<T>& bsf = dynamic_cast<bstf<T>&>(bfi_m.GetBareField());
          // advance to next boundary face
          int loc[3];
          bool OK;
          do {
            ++face_m;
            if (face_m>5) nextCell();
            if (bfi_m != bfend_m) {
              bfi_m.GetCurrentLocation(loc);
              OK = bsf.checkIndices(loc[0], loc[1], loc[2], face_m);
            }
            else {
              // at end, so just mark it OK
              OK = true;
            }
          } while (!OK);
          return *this;
        }
      // post-increment
        bstf_const_iterator<T> operator++(int)
        {
          bstf<T>& bsf = dynamic_cast<bstf<T>&>(bfi_m.GetBareField());
          // advance to next boundary face
          bstf_const_iterator<T> iter(*this);
          int loc[3];
          bool OK;
          do {
            ++face_m;
            if (face_m>5) nextCell();
            if (bfi_m != bfend_m) {
              bfi_m.GetCurrentLocation(loc);
              OK = bsf.checkIndices(loc[0], loc[1], loc[2], face_m);
            }
            else {
              // at end, so just mark it OK
              OK = true;
            }
          } while (!OK);
          return iter;
        }
      // dereference
        const T& operator*() const { return (*bfi_m)(face_m); }

      // comparisons
        bool operator==(const bstf_const_iterator<T>& rhs) const
        {
          return (bfi_m == rhs.bfi_m && face_m == rhs.face_m);
        }
        bool operator!=(const bstf_const_iterator<T>& rhs) const
        {
          return !(*this == rhs);
        }

      private:

      // methods
        void nextCell()
        {
          ++bfi_m;
          face_m = 0;
        }

      // data
        BaseField_iterator bfi_m;
        BaseField_iterator bfend_m;
        int face_m;
    };

  private:

  // data
    Mesh_t*   Mesh_m;
    Layout_t* FLayout_m;
    ccvf::BaseField_t*     CellSizes_m;
    ccvf::BaseField_t*     CellPositions_m;
    ccsf::BaseField_t*     CellVolumes_m;

  public:

  // constructors
    PoomaMesh_XYZ(int* const ncells, double* const cwidth,
                  e_dim_tag* decomp)
    {
      // create domain for mesh vertices (one larger than no. of cells)
      NDIndex<3> meshDom;
      for (int d=0; d<3; ++d)
        meshDom[d] = Index(ncells[d]+1);
      // create the mesh
      Mesh_m = new Mesh_t(meshDom,cwidth);
      // tell mesh to store mesh spacing fields
      Mesh_m->storeSpacingFields(decomp);
      // create the field layout
      FLayout_m = new Layout_t(*Mesh_m, decomp);
      // initialize internal fields
      initializeCellSizes();
      initializeCellPositions();
      initializeCellVolumes();
    }
    PoomaMesh_XYZ(int* const ncells, double** const cwidth,
                  e_dim_tag* decomp)
    {
      // create domain for mesh vertices (one larger than no. of cells)
      NDIndex<3> meshDom;
      for (int d=0; d<3; ++d)
        meshDom[d] = Index(ncells[d]+1);
      // create the mesh
      Mesh_m = new Mesh_t(meshDom,cwidth);
      // tell mesh to store mesh spacing fields
      Mesh_m->storeSpacingFields(decomp);
      // create the field layout
      FLayout_m = new Layout_t(*Mesh_m, decomp);
      // initialize internal fields
      initializeCellSizes();
      initializeCellPositions();
      initializeCellVolumes();
    }

  // destructor
    ~PoomaMesh_XYZ()
    {
      delete CellSizes_m;
      delete CellPositions_m;
      delete CellVolumes_m;
      delete FLayout_m;
      delete Mesh_m;
    }

  // accessors
    const Mesh_t& get_Mesh() const { return *Mesh_m; }
    Mesh_t&       get_Mesh()       { return *Mesh_m; }
    const Layout_t& get_Layout() const { return *FLayout_m; }
    Layout_t&       get_Layout()       { return *FLayout_m; }

    int get_ncx() const { return (Mesh_m->gridSizes[0] - 1); }
    int get_ncy() const { return (Mesh_m->gridSizes[1] - 1); }
    int get_ncz() const { return (Mesh_m->gridSizes[2] - 1); }

    int get_total_cells() const
    {
      return get_ncx() * get_ncy() * get_ncz();
    }
    int get_cells() const
    {
      // count cells of one of our internal fields
      ccsf::BaseField_t::iterator cvi, cvend = CellVolumes_m->end();
      int count = 0;
      for (cvi = CellVolumes_m->begin(); cvi != cvend; ++cvi)
        ++count;  // increment cell counter
      return count;
    }

    // rmr We no longer need get_ncp()
    int get_ncp() const { return get_cells(); }

    void get_dx(ccsf& f) const
    {
      f.Uncompress();
      ccvf::BaseField_t::iterator csi, csend = CellSizes_m->end();
      ccsf::iterator fi = f.begin();
      for (csi = CellSizes_m->begin(); csi != csend; ++csi, ++fi)
        *fi = (*csi)(0);
    }
    void get_dy(ccsf& f) const
    {
      f.Uncompress();
      ccvf::BaseField_t::iterator csi, csend = CellSizes_m->end();
      ccsf::iterator fi = f.begin();
      for (csi = CellSizes_m->begin(); csi != csend; ++csi, ++fi)
        *fi = (*csi)(1);
    }
    void get_dz(ccsf& f) const
    {
      f.Uncompress();
      ccvf::BaseField_t::iterator csi, csend = CellSizes_m->end();
      ccsf::iterator fi = f.begin();
      for (csi = CellSizes_m->begin(); csi != csend; ++csi, ++fi)
        *fi = (*csi)(2);
    }

    void get_xloc(ccsf& f) const
    {
      f.Uncompress();
      ccvf::BaseField_t::iterator cpi, cpend = CellPositions_m->end();
      ccsf::iterator fi = f.begin();
      for (cpi = CellPositions_m->begin(); cpi != cpend; ++cpi, ++fi)
        *fi = (*cpi)(0);
    }
    void get_yloc(ccsf& f) const
    {
      f.Uncompress();
      ccvf::BaseField_t::iterator cpi, cpend = CellPositions_m->end();
      ccsf::iterator fi = f.begin();
      for (cpi = CellPositions_m->begin(); cpi != cpend; ++cpi, ++fi)
        *fi = (*cpi)(1);
    }
    void get_zloc(ccsf& f) const
    {
      f.Uncompress();
      ccvf::BaseField_t::iterator cpi, cpend = CellPositions_m->end();
      ccsf::iterator fi = f.begin();
      for (cpi = CellPositions_m->begin(); cpi != cpend; ++cpi, ++fi)
        *fi = (*cpi)(2);
    }
    void get_xloc(fcdsf& f) const
    {
      f.Uncompress();
      ccvf::BaseField_t::iterator cpi, cpend = CellPositions_m->end(),
                                  csi = CellSizes_m->begin();
      fcdsf::iterator fi = f.begin();
      double xpos, xsize;
      for (cpi = CellPositions_m->begin(); cpi != cpend; ++cpi, ++csi) {
        xpos = (*cpi)(0);
        xsize = (*csi)(0);
        *fi = xpos - 0.5 * xsize;  ++fi;
        *fi = xpos + 0.5 * xsize;  ++fi;
        *fi = xpos;  ++fi;
        *fi = xpos;  ++fi;
        *fi = xpos;  ++fi;
        *fi = xpos;  ++fi;
      }
    }
    void get_yloc(fcdsf& f) const
    {
      f.Uncompress();
      ccvf::BaseField_t::iterator cpi, cpend = CellPositions_m->end(),
                                  csi = CellSizes_m->begin();
      fcdsf::iterator fi = f.begin();
      double ypos, ysize;
      for (cpi = CellPositions_m->begin(); cpi != cpend; ++cpi, ++csi) {
        ypos = (*cpi)(1);
        ysize = (*csi)(1);
        *fi = ypos;  ++fi;
        *fi = ypos;  ++fi;
        *fi = ypos - 0.5 * ysize;  ++fi;
        *fi = ypos + 0.5 * ysize;  ++fi;
        *fi = ypos;  ++fi;
        *fi = ypos;  ++fi;
      }
    }
    void get_zloc(fcdsf& f) const
    {
      f.Uncompress();
      ccvf::BaseField_t::iterator cpi, cpend = CellPositions_m->end(),
                                  csi = CellSizes_m->begin();
      fcdsf::iterator fi = f.begin();
      double zpos, zsize;
      for (cpi = CellPositions_m->begin(); cpi != cpend; ++cpi, ++csi) {
        zpos = (*cpi)(2);
        zsize = (*csi)(2);
        *fi = zpos;  ++fi;
        *fi = zpos;  ++fi;
        *fi = zpos;  ++fi;
        *fi = zpos;  ++fi;
        *fi = zpos - 0.5 * zsize;  ++fi;
        *fi = zpos + 0.5 * zsize;  ++fi;
      }
    }

    void get_face_lengths(fcdsf& f) const
    {
      f.Uncompress();
      ccvf::BaseField_t::iterator csi, csend = CellSizes_m->end();
      fcdsf::iterator fi = f.begin();
      for (csi = CellSizes_m->begin(); csi != csend; ++csi) {
        *fi = (*csi)(0);  ++fi;
        *fi = (*csi)(0);  ++fi;
        *fi = (*csi)(1);  ++fi;
        *fi = (*csi)(1);  ++fi;
        *fi = (*csi)(2);  ++fi;
        *fi = (*csi)(2);  ++fi;
      }
    }
    void get_face_areas(fcdsf& f) const
    {
      f.Uncompress();
      ccvf::BaseField_t::iterator csi, csend = CellSizes_m->end();
      fcdsf::iterator fi = f.begin();
      for (csi = CellSizes_m->begin(); csi != csend; ++csi) {
        *fi = (*csi)(1) * (*csi)(2);  ++fi;
        *fi = (*csi)(1) * (*csi)(2);  ++fi;
        *fi = (*csi)(2) * (*csi)(0);  ++fi;
        *fi = (*csi)(2) * (*csi)(0);  ++fi;
        *fi = (*csi)(0) * (*csi)(1);  ++fi;
        *fi = (*csi)(0) * (*csi)(1);  ++fi;
      }
    }
    void get_cell_volumes(ccsf& f) const
    {
      f = (*CellVolumes_m);
    }

  // gather/scatter operations
    template <class T1, class T2, class Op>
    static void scatter(fcdtf<T1>& to, const cctf<T2>& from, const Op& op);
    template <class T1, class T2, class Op>
    static void scatter(cctf<T1>& to, const fcdtf<T2>& from, const Op& op);
    template <class T1, class T2, class Op>
    static void gather(fcdtf<T1>& to, const cctf<T2>& from, const Op& op);
    template <class T1, class T2, class Op>
    static void gather(bstf<T1>& to, const fcdtf<T2>& from, const Op& op);
    template <class T1, class T2, class Op>
    static void gather(fcdtf<T1>& to, const bstf<T2>& from, const Op& op);
    template <class T>
    static void swap(fcdtf<T>& f1, fcdtf<T>& f2);

  // reduction operations
    template <class T>
    static T sum(const fcdtf<T>& f);
    template <class T>
    static T min(const fcdtf<T>& f);
    template <class T>
    static T max(const fcdtf<T>& f);
    template <class T>
    static T sum(const cctf<T>& f);
    template <class T>
    static T min(const cctf<T>& f);
    template <class T>
    static T max(const cctf<T>& f);
    template <class T>
    static T sum(const bstf<T>& f);
    template <class T>
    static T min(const bstf<T>& f);
    template <class T>
    static T max(const bstf<T>& f);

  private:

  // methods
    void initializeCellSizes()
    {
      CellSizes_m = new ccvf::BaseField_t(*Mesh_m, *FLayout_m);
      Mesh_m->getDeltaVertexField(*CellSizes_m);
    }
    void initializeCellPositions()
    {
      CellPositions_m = new ccvf::BaseField_t(*Mesh_m, *FLayout_m);
      Mesh_m->getCellPositionField(*CellPositions_m);
    }
    void initializeCellVolumes()
    {
      CellVolumes_m = new ccsf::BaseField_t(*Mesh_m, *FLayout_m);
      Mesh_m->getCellVolumeField(*CellVolumes_m);
    }
};

#include "PoomaMesh_XYZ.t.cc"

#endif                          // __mesh_PoomaMesh_XYZ_hh__

//---------------------------------------------------------------------------//
//                              end of mesh/PoomaMesh_XYZ.hh
//---------------------------------------------------------------------------//
