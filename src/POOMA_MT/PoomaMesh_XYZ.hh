//----------------------------------*-C++-*----------------------------------//
// PoomaMesh_XYZ.hh
// Julian C. Cummings
// Wed Jan 27 1999
//---------------------------------------------------------------------------//
// @> A 3-d cartesian structured mesh facility based on POOMA r1.
//---------------------------------------------------------------------------//


#ifndef __mesh_PoomaMesh_XYZ_hh__
#define __mesh_PoomaMesh_XYZ_hh__

// include files

// Configuration includes

#include <POOMA_MT/config.h>

// POOMA headers
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

// XTM headers
#include "VektorHelper.hh"
#include "ds++/config.hh"
#include "ds++/SP.hh"

#include "traits/MT_traits.hh"

// Standard C++ headers
#include <iterator.h>

// forward declaration
template <class Mesh> class PoomaMesh_XYZ;

// #ifndef __traits_MT_traits_hh__
// #define __traits_MT_traits_hh__

// Vector traits class

namespace rtt_traits
{

// template<class VECTOR>
// class vector_traits
// {
// public:
//	typedef typename VECTOR::value_type value_type;

//	inline static value_type dot(const VECTOR &v1, const VECTOR &v2)
//	{
//	  return VECTOR::dot(v1, v2);
//	}
// };

// specialization for POOMA Vektor class
template <class T, unsigned Dim>
class vector_traits< Vektor<T,Dim> >
{
  public:
    typedef typename Vektor<T,Dim>::Element_t value_type;
    
    inline static value_type dot( const Vektor<T,Dim>& v1,
				  const Vektor<T,Dim>& v2 )
    {
	return ::dot(v1, v2);
    }
};

} // end namespace rtt_traits

// #endif 	// __traits_MT_traits_hh__

//===========================================================================//
// class PoomaMesh_XYZ - A 3-d cartesian mesh class based on POOMA r1
// This is a 3-d cartesian structured mesh.  The main purpose of having this
// class is in order for it to be instantiated by various transport codes.
//===========================================================================//

template <class Mesh>
class PoomaMesh_XYZ
{
  public:

    // forward declarations
    template <class T> class cctf;
    template <class T> class cctf_iterator;
    template <class T> class cctf_const_iterator;

    template <class T> class fcdtf;
    template <class T> class fcdtf_iterator;
    template <class T> class fcdtf_const_iterator;

    template <class T> class bstf;
    template <class T> class bstf_iterator;
    template <class T> class bstf_const_iterator;

    template <class T> class nctf;
    template <class T> class nctf_iterator;
    template <class T> class nctf_const_iterator;

    template <class T> class vctf;
    template <class T> class vctf_iterator;
    template <class T> class vctf_const_iterator;

    template <class FT> class ConnFacesAroundVertices;
    template <class FT> class VertexProxy;
    template <class FT> class CFAV_iterator;
    template <class FT> class CFAV_const_iterator;

    class vec;

    // typedefs
    typedef unsigned long size_type;
    typedef dsxx::SP< PoomaMesh_XYZ<Mesh> > FieldConstructor;

    typedef cctf<double>   ccsf;
    typedef fcdtf<double> fcdsf;
    typedef nctf<double>   ncsf;
    typedef vctf<double>   vcsf;
    typedef bstf<double>   bssf;

    typedef cctf<int>   ccif;
    typedef fcdtf<int> fcdif;
    typedef nctf<int>   ncif;
    typedef vctf<int>   vcif;
    typedef bstf<int>   bsif;
    // typedef cctf<std::vector<int> > ccvif;
    // typedef typename MT:: template cctf<std::vector<int> > ccvif;
    typedef cctf<vec>   ccvsf;
    typedef fcdtf<vec> fcdvsf;
    typedef nctf<vec>   ncvsf;
    typedef vctf<vec>   vcvsf;
    typedef bstf<vec>   bsvsf;

    // These refer to PETE operators, which are in global namespace.
    typedef ::OpAssign OpAssign;
    typedef ::OpAddAssign OpAddAssign;
    typedef ::OpSubtractAssign OpSubAssign;
    typedef ::OpMultiplyAssign OpMultAssign;
    typedef ::OpMinAssign OpMinAssign;
    typedef ::OpMaxAssign OpMaxAssign;

    // Class vec, a model of a 3D tiny_vec of doubles

    class vec : public Vektor<double,3>
    {
      public:

	// typedefs for container
	typedef double value_type;
	typedef double& reference;
	typedef const double& const_reference;
	typedef double* pointer;
	typedef const double* const_pointer;
	typedef pointer iterator;
	typedef const_pointer const_iterator;
	typedef int difference_type;
	typedef unsigned int size_type;

	// typedefs for reverse_iterators
	typedef std::reverse_iterator<iterator> reverse_iterator;
	typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
	
	// constructors
	vec() {}
	vec(const double& x) : Vektor<double,3>(x) {}
	vec(const double& x00, const double& x01, const double& x02)
	    : Vektor<double,3>(x00,x01,x02) {}
	vec(const vec& rhs) : Vektor<double,3>(rhs) {}
	vec(const Vektor<double,3>& rhs) : Vektor<double,3>(rhs) {}

	// destructor
	~vec() {}

        // iterators
        iterator begin(void)
        {
	    vec &self = *this;
            iterator p = &(self[0]);
            return p;
        }
        iterator end(void)
        {
            iterator p = begin() + 3;
            return p;
        }
        const_iterator begin(void) const
        {
	    vec &self = const_cast<vec&>(*this);
            const_iterator p = &(self[0]);
            return p;
        }
        const_iterator end(void) const
        {
            const_iterator p = begin() + 3;
            return p;
        }

        // reverse iterators
        reverse_iterator rbegin(void)
        {
            return reverse_iterator(end());
        }
        reverse_iterator rend(void)
        {
            return reverse_iterator(begin());
        }
        const_reverse_iterator rbegin(void) const
        {
            return const_reverse_iterator(end());
        }
        const_reverse_iterator rend(void) const
        {
            return const_reverse_iterator(begin());
        }

        // accessors
        size_type size(void) const { return 3; }
        size_type max_size(void) const { return 3; }
        bool empty(void) const { return false; }
        void swap(vec& rhs)
        {
            vec temp = *this;
            *this = rhs;
            rhs = temp;
        }

        // LessThan Comparable
        bool operator<(const vec& rhs) const
        {
            const_iterator i, iend = end(), rhsi = rhs.end();
            for (i = begin(); i != iend; ++i, ++rhsi) {
                if (*i < *rhsi) return true;
                if (*rhsi < *i) return false;
            }
            return false;
        }
        bool operator>(const vec& rhs) const { return (rhs < *this); }
        bool operator<=(const vec& rhs) const { return !(*this > rhs); }
        bool operator>=(const vec& rhs) const { return !(*this < rhs); }

        // Operator =
        vec& operator=(const Vektor<double,3>& rhs)
        {
            for (int d=0; d<3; ++d)
                (*this)(d) = rhs(d);
            return *this;
        }

        // Dot product
        static double dot(const vec& v1, const vec& v2)
        {
            return ::dot(v1,v2);
        }
    };

    // Cell centered field
    // Has a value in each cell.

    template <class T>
    class cctf : public Field<T,3,Mesh,Cell>
    {
        friend class PoomaMesh_XYZ<Mesh>;
        friend class cctf_iterator<T>;
        friend class cctf_const_iterator<T>;

      public:

        // typedefs for container
        typedef T value_type;
        typedef T& reference;
        typedef const T& const_reference;
        typedef T* pointer;
        typedef const T* const_pointer;
        typedef cctf_iterator<T> iterator;
        typedef cctf_const_iterator<T> const_iterator;
        typedef long difference_type;
        typedef unsigned long size_type;
        // typedefs for MT field
        typedef PoomaMesh_XYZ<Mesh> MT_t;

      private:

    // Some additional typedefs for internal use
        typedef CenteredFieldLayout<3,Mesh,Cell> CLayout_t;
        typedef Field<T,3,Mesh,Cell> BaseField_t;
        typedef GuardCellSizes<3> GC_t;

      public:

    // constructors
        cctf(const MT_t::FieldConstructor& spm)
            : BaseField_t( const_cast<Mesh&>(spm->get_Mesh()),
                           const_cast<CLayout_t&>(spm->get_CLayout()),
                           GC_t(1) ),
              spm_m(spm) {}
        cctf(const MT_t& m)
            : BaseField_t( const_cast<Mesh&>(m.get_Mesh()),
                           const_cast<CLayout_t&>(m.get_CLayout()),
                           GC_t(1) ),
              spm_m(const_cast<MT_t*>(&m)) {}
        cctf(const cctf<T>& f)
            : BaseField_t( const_cast<Mesh&>(f.get_Mesh().get_Mesh()),
                           const_cast<CLayout_t&>(f.get_Mesh().get_CLayout()),
                           GC_t(1) ),
              spm_m( f.get_FieldConstructor() ) {}

        // Assignment from scalar
        cctf<T>& operator=(T x)
        {
            assign(*this,x);
            return *this;
        }
        // Assignment from another field
        cctf<T>& operator=(const cctf<T>& x)
        {
            PAssert( size() == x.size() );
            assign(*this,x);
            return *this;
        }
        // Assignment from base field
        cctf<T>& operator=(const BaseField_t& x)
        {
            assign(*this,x);
            return *this;
        }
        // Assignment from an expression
        template <class B>
        cctf<T>& operator=(const PETE_Expr<B>& x)
        {
            assign(*this,x);
            return *this;
        }

        // element accessors
        value_type& operator()(unsigned int i, unsigned int j, unsigned int k)
        {
            NDIndex<3> loc(Index(i,i),Index(j,j),Index(k,k));
            return localElement(loc);
        }
        const value_type&
        operator()(unsigned int i, unsigned int j, unsigned int k) const
        {
            NDIndex<3> loc(Index(i,i),Index(j,j),Index(k,k));
            return localElement(loc);
        }
        value_type& operator()(size_type c)
        {
            unsigned int i,j,k;
            spm_m->get_cell_indices(c,i,j,k);
            NDIndex<3> loc(Index(i,i),Index(j,j),Index(k,k));
            return localElement(loc);
        }
        const value_type& operator()(size_type c) const
        {
            unsigned int i,j,k;
            spm_m->get_cell_indices(c,i,j,k);
            NDIndex<3> loc(Index(i,i),Index(j,j),Index(k,k));
            return localElement(loc);
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
        const MT_t::FieldConstructor& get_FieldConstructor(void) const
        {
            return spm_m;
        }
        const MT_t& get_Mesh(void) const
        {
            return *spm_m;
        }
        size_type size(void) const { return spm_m->get_ncells(); }
        size_type max_size(void) const { return size(); }
        bool empty(void) const { return (size() == 0); }
        
        void swap(cctf<T>& f)
        {
            // just exchange vmap of pointers to local fields
            Locals_ac.swap(f.Locals_ac);
        }

        // comparisons
        bool operator==(const cctf<T>& rhs) const
        {
            if (size() != rhs.size()) return false;
            const_iterator i, iend = end(), rhsi = rhs.begin();
            bool result = true;
            for (i = begin(); i != iend; ++i, ++rhsi)
                if (*i != *rhsi) result = false;
            return result;
        }
        bool operator!=(const cctf<T>& rhs) const { return !(*this==rhs); }
        
        bool operator<(const cctf<T>& rhs) const
        {
            const_iterator i, iend, rhsi;
            if (size() <= rhs.size()) {
                // use this iterator as loop bounds
                iend = end();
                rhsi = rhs.begin();
                for (i = begin(); i != iend; ++i, ++rhsi) {
                    if (*i < *rhsi) return true;
                    if (*rhsi < *i) return false;
                }
                // all elements equal up to this point
                if (size() < rhs.size())
                    return true;
                else
                    return false;
            }
            else {
                // use rhs iterator as loop bounds
                iend = rhs.end();
                i = begin();
                for (rhsi = rhs.begin(); rhsi != iend; ++i, ++rhsi) {
                    if (*i < *rhsi) return true;
                    if (*rhsi < *i) return false;
                }
                // all elements equal up to this point
                return false;
            }
        }
        bool operator>(const cctf<T>& rhs) const { return (rhs < *this); }
        bool operator<=(const cctf<T>& rhs) const { return !(*this > rhs); }
        bool operator>=(const cctf<T>& rhs) const { return !(*this < rhs); }

      private:
        // data
        MT_t::FieldConstructor spm_m;

    };
        
    // Cell centered field iterators

    template <class T>
    class cctf_iterator : public forward_iterator<T,ptrdiff_t>
    {
        friend class cctf_const_iterator<T>;

        // typedefs
        typedef typename cctf<T>::BaseField_t BaseField_t;
        typedef typename BaseField_t::iterator BaseField_iterator;

      public:

        // constructors
        cctf_iterator() {}
        cctf_iterator(BaseField_t& bf)
            : bfi_m(bf.begin()) {}
        cctf_iterator(const BaseField_iterator& bfiter)
            : bfi_m(bfiter) {}
        cctf_iterator(const cctf_iterator<T>& iter)
            : bfi_m(iter.bfi_m) {}

        // assignment
        cctf_iterator<T>& operator=(const cctf_iterator<T>& x)
        {
            bfi_m = x.bfi_m;
            return *this;
        }

        // pre-increment
        cctf_iterator<T>& operator++()
        {
            ++bfi_m;
            return *this;
        }
        // post-increment
        cctf_iterator<T> operator++(int)
        {
            cctf_iterator<T> iter(*this);
            ++bfi_m;
            return iter;
        }
        // dereference
        T& operator*() const
        {
            return *bfi_m;
        }
	// member
	T* operator->() const
	{
	    return &(*bfi_m);
	}

        // comparisons
        bool operator==(const cctf_iterator<T>& rhs) const
        {
            return (bfi_m == rhs.bfi_m);
        }
        bool operator!=(const cctf_iterator<T>& rhs) const
        {
            return !(*this == rhs);
        }

      private:
        // data
        BaseField_iterator bfi_m;
    };
    
    template <class T>
    class cctf_const_iterator : public input_iterator<T,ptrdiff_t>
    {
        // typedefs
        typedef typename cctf<T>::BaseField_t BaseField_t;
        typedef typename BaseField_t::iterator BaseField_iterator;

      public:

        // constructors
        cctf_const_iterator() {}
        cctf_const_iterator(const BaseField_t& bf)
            : bfi_m(bf.begin()) {}
        cctf_const_iterator(const BaseField_iterator& bfiter)
            : bfi_m(bfiter) {}
        cctf_const_iterator(const cctf_const_iterator<T>& iter)
            : bfi_m(iter.bfi_m) {}
        cctf_const_iterator(const cctf_iterator<T>& iter)
            : bfi_m(iter.bfi_m) {}

        // assignment
        cctf_const_iterator<T>&
        operator=(const cctf_const_iterator<T>& x)
        {
            bfi_m = x.bfi_m;
            return *this;
        }
        cctf_const_iterator<T>&
        operator=(const cctf_iterator<T>& x)
        {
            bfi_m = x.bfi_m;
            return *this;
        }
        
        // pre-increment
        cctf_const_iterator<T>& operator++()
        {
            ++bfi_m;
            return *this;
        }
        // post-increment
        cctf_const_iterator<T> operator++(int)
        {
            cctf_const_iterator<T> iter(*this);
            ++bfi_m;
            return iter;
        }
        // dereference
        const T& operator*() const { return (*bfi_m); }
	// member
	const T* operator->() const { return &(*bfi_m); }

        // comparisons
        bool operator==(const cctf_const_iterator<T>& rhs) const
        {
            return (bfi_m == rhs.bfi_m);
        }
        bool operator!=(const cctf_const_iterator<T>& rhs) const
        {
            return !(*this == rhs);
        }

      private:

        // data
        BaseField_iterator bfi_m;
    };
    
    // Face centered discontinuous field
    // Has a value on each face in each cell.

    template <class T>
    class fcdtf : public Field<Vektor<T,6>,3,Mesh,Cell>
    {
        friend class PoomaMesh_XYZ<Mesh>;
        friend class fcdtf_iterator<T>;
        friend class fcdtf_const_iterator<T>;

      public:

        // typedefs for container
        typedef T value_type;
        typedef T& reference;
        typedef const T& const_reference;
        typedef T* pointer;
        typedef const T* const_pointer;
        typedef fcdtf_iterator<T> iterator;
        typedef fcdtf_const_iterator<T> const_iterator;
        typedef long difference_type;
        typedef unsigned long size_type;
        // typedefs for MT field
        typedef PoomaMesh_XYZ<Mesh> MT_t;

      private:

        // Some additional typedefs for internal use
        typedef CenteredFieldLayout<3,Mesh,Cell> CLayout_t;
        typedef Field<Vektor<T,6>,3,Mesh,Cell> BaseField_t;
        typedef GuardCellSizes<3> GC_t;

      public:

        // constructors
        fcdtf(const MT_t::FieldConstructor& spm)
            : BaseField_t( const_cast<Mesh&>(spm->get_Mesh()),
                           const_cast<CLayout_t&>(spm->get_CLayout()),
                           GC_t(1) ),
              spm_m(spm) {}
        fcdtf(const MT_t& m)
            : BaseField_t( const_cast<Mesh&>(m.get_Mesh()),
                           const_cast<CLayout_t&>(m.get_CLayout()),
                           GC_t(1) ),
              spm_m(const_cast<MT_t*>(&m)) {}
        fcdtf(const fcdtf<T>& f)
            : BaseField_t( const_cast<Mesh&>(f.get_Mesh().get_Mesh()),
                           const_cast<CLayout_t&>(f.get_Mesh().get_CLayout()),
                           GC_t(1) ),
              spm_m( f.get_FieldConstructor() ) {}

        // Assignment from scalar
        fcdtf<T>& operator=(T x)
        {
            Vektor<T,6> v = x;
            assign(*this,v);
            return *this;
        }
        // Assignment from another field
        fcdtf<T>& operator=(const fcdtf<T>& x)
        {
            PAssert( size() == x.size() );
            assign(*this,x);
            return *this;
        }
        // Assignment from base field
        fcdtf<T>& operator=(const BaseField_t& x)
        {
            assign(*this,x);
            return *this;
        }
        // Assignment from an expression
        template <class B>
        fcdtf<T>& operator=(const PETE_Expr<B>& x)
        {
            assign(*this,x);
            return *this;
        }
	// AddAssign
	template <class Expr>
	fcdtf<T>& operator+=(const Expr &expr)
	{
	    *this = *this + expr;
	    return *this;
	}
	// SubAssign
	template <class Expr>
	fcdtf<T>& operator-=(const Expr &expr)
	{
	    *this = *this - expr;
	    return *this;
	}
	// MulAssign
	template <class Expr>
	fcdtf<T>& operator*=(const Expr &expr)
	{
	    *this = *this * expr;
	    return *this;
	}
	// DivAssign
	template <class Expr>
	fcdtf<T>& operator/=(const Expr &expr)
	{
	    *this = *this / expr;
	    return *this;
	}
        
        // face accessors
        // i, j, k == global xyz cell indices
        // c == cell number
        // f == face index
        cctf<T> operator()(unsigned int f) const
        {
            cctf<T> cf(get_FieldConstructor());
            typename cctf<T>::iterator cfit, cfend = cf.end();
            const BaseField_t& bf = dynamic_cast<const BaseField_t&>(*this);
            typename BaseField_t::iterator bfit = bf.begin();
            for (cfit = cf.begin(); cfit != cfend; ++cfit, ++bfit)
                *cfit = (*bfit)(f);
            return cf;
        }
        value_type&
        operator()(unsigned int i, unsigned int j, unsigned int k,
                   unsigned int f)
        {
            NDIndex<3> loc(Index(i,i),Index(j,j),Index(k,k));
            return localElement(loc)(f);
        }
        const value_type&
        operator()(unsigned int i, unsigned int j, unsigned int k,
                   unsigned int f) const
        {
            NDIndex<3> loc(Index(i,i),Index(j,j),Index(k,k));
            return localElement(loc)(f);
        }
        value_type& operator()(size_type c, unsigned int f)
        {
            unsigned int i,j,k;
            spm_m->get_cell_indices(c,i,j,k);
            NDIndex<3> loc(Index(i,i),Index(j,j),Index(k,k));
            return localElement(loc)(f);
        }
        const value_type& operator()(size_type c, unsigned int f) const
        {
            unsigned int i,j,k;
            spm_m->get_cell_indices(c,i,j,k);
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
        const MT_t::FieldConstructor& get_FieldConstructor(void) const
        {
            return spm_m;
        }
        const MT_t& get_Mesh(void) const
        {
            return *spm_m;
        }
        size_type size(void) const { return (spm_m->get_ncells() * 6); }
        size_type max_size(void) const { return size(); }
        bool empty(void) const { return (size() == 0); }

        void swap(fcdtf<T>& f)
        {
            // just exchange vmap of pointers to local fields
            Locals_ac.swap(f.Locals_ac);
        }

        // comparisons
        bool operator==(const fcdtf<T>& rhs) const
        {
            if (size() != rhs.size()) return false;
            const_iterator i, iend = end(), rhsi = rhs.begin();
            bool result = true;
            for (i = begin(); i != iend; ++i, ++rhsi)
                if (*i != *rhsi) result = false;
            return result;
        }
        bool operator!=(const fcdtf<T>& rhs) const { return !(*this==rhs); }

        bool operator<(const fcdtf<T>& rhs) const
        {
            const_iterator i, iend, rhsi;
            if (size() <= rhs.size()) {
                // use this iterator as loop bounds
                iend = end();
                rhsi = rhs.begin();
                for (i = begin(); i != iend; ++i, ++rhsi) {
                    if (*i < *rhsi) return true;
                    if (*rhsi < *i) return false;
                }
                // all elements equal up to this point
                if (size() < rhs.size())
                    return true;
                else
                    return false;
            }
            else {
                // use rhs iterator as loop bounds
                iend = rhs.end();
                i = begin();
                for (rhsi = rhs.begin(); rhsi != iend; ++i, ++rhsi) {
                    if (*i < *rhsi) return true;
                    if (*rhsi < *i) return false;
                }
                // all elements equal up to this point
                return false;
            }
        }
        bool operator>(const fcdtf<T>& rhs) const { return (rhs < *this); }
        bool operator<=(const fcdtf<T>& rhs) const { return !(*this > rhs); }
        bool operator>=(const fcdtf<T>& rhs) const { return !(*this < rhs); }

      private:
        // data
        MT_t::FieldConstructor spm_m;

    };
    
    // Face centered discontinuous field iterators

    template <class T>
    class fcdtf_iterator : public forward_iterator<T,ptrdiff_t>
    {
        friend class fcdtf_const_iterator<T>;

        // typedefs
        typedef typename fcdtf<T>::BaseField_t BaseField_t;
        typedef typename BaseField_t::iterator BaseField_iterator;

      public:

        // constructors
        fcdtf_iterator() : face_m(0) {}
        fcdtf_iterator(BaseField_t& bf)
            : bfi_m(bf.begin()), face_m(0) {}
        fcdtf_iterator(const BaseField_iterator& bfiter)
            : bfi_m(bfiter), face_m(0) {}
        fcdtf_iterator(const fcdtf_iterator<T>& iter)
            : bfi_m(iter.bfi_m), face_m(iter.face_m) {}

        // assignment
        fcdtf_iterator<T>& operator=(const fcdtf_iterator<T>& x)
        {
            bfi_m = x.bfi_m;
            face_m = x.face_m;
            return *this;
        }
        fcdtf_iterator<T>& operator=(const bstf_iterator<T>& x)
        {
            // set BaseField iterator and face equal
            FieldLoc<3> loc;
            x.bfi_m.GetCurrentLocation(loc);
            bfi_m.SetCurrentLocation(loc);
            face_m = x.face_m;
            return *this;
        }
        fcdtf_iterator<T>& operator=(const bstf_const_iterator<T>& x)
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
            if (face_m==6) nextCell();
            return *this;
        }
        // post-increment
        fcdtf_iterator<T> operator++(int)
        {
            fcdtf_iterator<T> iter(*this);
            ++face_m;
            if (face_m==6) nextCell();
            return iter;
        }
        // dereference
        T& operator*() const { return (*bfi_m)(face_m); }
        // member
        T* operator->() const { return &((*bfi_m)(face_m)); }

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
        unsigned int face_m;
    };
    
    template <class T>
    class fcdtf_const_iterator : public input_iterator<T,ptrdiff_t>
    {
        // typedefs
        typedef typename fcdtf<T>::BaseField_t BaseField_t;
        typedef typename BaseField_t::iterator BaseField_iterator;

      public:

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
        fcdtf_const_iterator<T>&
        operator=(const fcdtf_const_iterator<T>& x)
        {
            bfi_m = x.bfi_m;
            face_m = x.face_m;
            return *this;
        }
        fcdtf_const_iterator<T>&
        operator=(const fcdtf_iterator<T>& x)
        {
            bfi_m = x.bfi_m;
            face_m = x.face_m;
            return *this;
        }
        fcdtf_const_iterator<T>&
        operator=(const bstf_iterator<T>& x)
        {
            // set BaseField iterator and face equal
            FieldLoc<3> loc;
            x.bfi_m.GetCurrentLocation(loc);
            bfi_m.SetCurrentLocation(loc);
            face_m = x.face_m;
            return *this;
        }
        fcdtf_const_iterator<T>&
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
            if (face_m==6) nextCell();
            return *this;
        }
        // post-increment
        fcdtf_const_iterator<T> operator++(int)
        {
            fcdtf_const_iterator<T> iter(*this);
            ++face_m;
            if (face_m==6) nextCell();
            return iter;
        }
        // dereference
        const T& operator*() const { return (*bfi_m)(face_m); }
	// member
        const T* operator->() const { return &((*bfi_m)(face_m)); }

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
        unsigned int face_m;
    };

    // Boundary specified field
    // Has a value on each boundary face.

    template <class T>
    class bstf : public Field<Vektor<T,6>,3,Mesh,Cell>
    {
        friend class PoomaMesh_XYZ<Mesh>;
        friend class bstf_iterator<T>;
        friend class bstf_const_iterator<T>;

      public:

        // typedefs for container
        typedef T value_type;
        typedef T& reference;
        typedef const T& const_reference;
        typedef T* pointer;
        typedef const T* const_pointer;
        typedef bstf_iterator<T> iterator;
        typedef bstf_const_iterator<T> const_iterator;
        typedef long difference_type;
        typedef unsigned long size_type;
        // typedefs for MT field
        typedef PoomaMesh_XYZ<Mesh> MT_t;

      private:

	// Some additional typedefs for internal use
	typedef CenteredFieldLayout<3,Mesh,Cell> CLayout_t;
	typedef Field<Vektor<T,6>,3,Mesh,Cell> BaseField_t;

      public:

	// constructors
	bstf(const MT_t::FieldConstructor& spm)
	    : BaseField_t( const_cast<Mesh&>(spm->get_Mesh()),
			   const_cast<CLayout_t&>(spm->get_CLayout()) ),
	      spm_m(spm)
	{
	    // compute local number of boundary faces
	    computeSize();
	}
	bstf(const MT_t& m)
	    : BaseField_t( const_cast<Mesh&>(m.get_Mesh()),
			   const_cast<CLayout_t&>(m.get_CLayout()) ),
	      spm_m(const_cast<MT_t*>(&m))
	{
	    // compute local number of boundary faces
	    computeSize();
	}
	bstf(const bstf<T>& f)
	    : BaseField_t( const_cast<Mesh&>(f.get_Mesh().get_Mesh()),
			   const_cast<CLayout_t&>(f.get_Mesh().get_CLayout()) ),
	      spm_m( f.get_FieldConstructor() ), size_m(f.size_m)
	{
	}

	// Assignment from scalar
	bstf<T>& operator=(T x) {
	    iterator it, iend = end();
	    for (it = begin(); it != iend; ++it)
		*it = x;
	    return *this;
	}
	// Assignment from another field
	bstf<T>& operator=(const bstf<T>& x) {
	    PAssert( size() == x.size() );
	    assign(*this,x);
	    return *this;
	}
	// Assignment from an expression
	template <class B>
	bstf<T>& operator=(const PETE_Expr<B>& x) {
	    assign(*this,x);
	    return *this;
	}
	// AddAssign
	template <class Expr>
	bstf<T>& operator+=(const Expr &expr)
	{
	    *this = *this + expr;
	    return *this;
	}
	// SubAssign
	template <class Expr>
	bstf<T>& operator-=(const Expr &expr)
	{
	    *this = *this - expr;
	    return *this;
	}
	// MulAssign
	template <class Expr>
	bstf<T>& operator*=(const Expr &expr)
	{
	    *this = *this * expr;
	    return *this;
	}
	// DivAssign
	template <class Expr>
	bstf<T>& operator/=(const Expr &expr)
	{
	    *this = *this / expr;
	    return *this;
	}

	// face accessors
	// i, j, k == global xyz cell indices
	// c == cell number
	// f == face index
	// must check that indices refer to a boundary face
	value_type&
	operator()(unsigned int i, unsigned int j, unsigned int k,
		   unsigned int f)
	{
	    bool OK = checkIndices(i,j,k,f);
	    PInsist(OK, "bstf<T>::operator() detects bad indices!!");
            NDIndex<3> loc(Index(i,i),Index(j,j),Index(k,k));
            return localElement(loc)(f);
        }
        const value_type&
        operator()(unsigned int i, unsigned int j, unsigned int k,
                   unsigned int f) const
        {
            bool OK = checkIndices(i,j,k,f);
            PInsist(OK, "bstf<T>::operator() detects bad indices!!");
            NDIndex<3> loc(Index(i,i),Index(j,j),Index(k,k));
            return localElement(loc)(f);
        }
        value_type& operator()(size_type c, unsigned int f)
        {
            unsigned int i,j,k;
            spm_m->get_cell_indices(c,i,j,k);
            bool OK = checkIndices(i,j,k,f);
            PInsist(OK, "bstf<T>::operator() detects bad indices!!");
            NDIndex<3> loc(Index(i,i),Index(j,j),Index(k,k));
            return localElement(loc)(f);
        }
        const value_type& operator()(size_type c, unsigned int f) const
        {
            unsigned int i,j,k;
            spm_m->get_cell_indices(c,i,j,k);
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
        const MT_t::FieldConstructor& get_FieldConstructor(void) const
        {
            return spm_m;
        }
        const MT_t& get_Mesh(void) const
        {
            return *spm_m;
        }
        size_type size(void) const { return size_m; }
        size_type max_size(void) const { return size(); }
        bool empty(void) const { return (size() == 0); }

        void swap(bstf<T>& f)
        {
            // just exchange vmap of pointers to local fields
            Locals_ac.swap(f.Locals_ac);
        }
        
        // comparisons
        bool operator==(const bstf<T>& rhs) const
        {
            if (size() != rhs.size()) return false;
            const_iterator i, iend = end(), rhsi = rhs.begin();
            bool result = true;
            for (i = begin(); i != iend; ++i, ++rhsi)
                if (*i != *rhsi) result = false;
            return result;
        }
        bool operator!=(const bstf<T>& rhs) const { return !(*this==rhs); }

        bool operator<(const bstf<T>& rhs) const
        {
            const_iterator i, iend, rhsi;
            if (size() <= rhs.size()) {
                // use this iterator as loop bounds
                iend = end();
                rhsi = rhs.begin();
                for (i = begin(); i != iend; ++i, ++rhsi) {
                    if (*i < *rhsi) return true;
                    if (*rhsi < *i) return false;
                }
                // all elements equal up to this point
                if (size() < rhs.size())
                    return true;
                else
                    return false;
            }
            else {
                // use rhs iterator as loop bounds
                iend = rhs.end();
                i = begin();
                for (rhsi = rhs.begin(); rhsi != iend; ++i, ++rhsi) {
                    if (*i < *rhsi) return true;
                    if (*rhsi < *i) return false;
                }
                // all elements equal up to this point
                return false;
            }
        }
        bool operator>(const bstf<T>& rhs) const { return (rhs < *this); }
        bool operator<=(const bstf<T>& rhs) const { return !(*this > rhs); }
        bool operator>=(const bstf<T>& rhs) const { return !(*this < rhs); }

      private:

        // methods
        void computeSize(void)
        {
            // compute local size of field
            size_m = 0;
            // iterate over local elements and faces,
            // and count how many are on the boundary
            int loc[3];
            unsigned int f;
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
        bool
        checkIndices(unsigned int i, unsigned int j, unsigned int k,
                     unsigned int f) const
        {
            bool result;
            if (i == 0 && f == 0) {
                // left boundary face
                result = true;
            }
            else if (i == spm_m->get_ncx()-1 && f == 1) {
                // right boundary face
                result = true;
            }
            else if (j == 0 && f == 2) {
                // bottom boundary face
                result = true;
            }
            else if (j == spm_m->get_ncy()-1 && f == 3) {
                // top boundary face
                result = true;
            }
            else if (k == 0 && f == 4) {
                // rear boundary face
                result = true;
            }
            else if (k == spm_m->get_ncz()-1 && f == 5) {
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
        MT_t::FieldConstructor spm_m;
        size_type size_m;

    };
    
    // Boundary specified field iterators

    template <class T>
    class bstf_iterator : public forward_iterator<T,ptrdiff_t>
    {
        friend class bstf_const_iterator<T>;
        friend class fcdtf_iterator<T>;
        friend class fcdtf_const_iterator<T>;

        // typedefs
        typedef typename bstf<T>::BaseField_t BaseField_t;
        typedef typename BaseField_t::iterator BaseField_iterator;

      public:

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
            if (!OK) ++(*this); // if not, advance to next valid boundary face
        }
        bstf_iterator(const BaseField_iterator& bfiter)
            : bfi_m(bfiter), bfend_m(bfiter), face_m(0) {}
        bstf_iterator(const bstf_iterator<T>& iter)
            : bfi_m(iter.bfi_m), bfend_m(iter.bfend_m),
              face_m(iter.face_m) {}

        // assignment
        bstf_iterator<T>& operator=(const bstf_iterator<T>& x)
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
                if (face_m==6) nextCell();
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
                if (face_m==6) nextCell();
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
        T& operator*() const { return (*bfi_m)(face_m); }
        // member
        T* operator->() const { return &((*bfi_m)(face_m)); }

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
        unsigned int face_m;
    };

    template <class T>
    class bstf_const_iterator : public input_iterator<T,ptrdiff_t>
    {
        friend class fcdtf_iterator<T>;
        friend class fcdtf_const_iterator<T>;

        // typedefs
        typedef typename bstf<T>::BaseField_t BaseField_t;
        typedef typename BaseField_t::iterator BaseField_iterator;

      public:

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
            if (!OK) ++(*this); // if not, advance to next valid boundary face
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
        bstf_const_iterator<T>&
        operator=(const bstf_const_iterator<T>& x)
        {
            bfi_m = x.bfi_m;
            bfend_m = x.bfend_m;
            face_m = x.face_m;
            return *this;
        }
        bstf_const_iterator<T>&
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
                if (face_m==6) nextCell();
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
                if (face_m==6) nextCell();
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
	// member
	const T* operator->() const { return &((*bfi_m)(face_m)); }

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
	unsigned int face_m;
    };

    // Node centered field
    // Has a value at each node.

    template <class T>
    class nctf : public Field<T,3,Mesh,Vert>
    {
	friend class PoomaMesh_XYZ<Mesh>;
	friend class nctf_iterator<T>;
	friend class nctf_const_iterator<T>;

      public:

	// typedefs for container
	typedef T value_type;
	typedef T& reference;
	typedef const T& const_reference;
	typedef T* pointer;
	typedef const T* const_pointer;
	typedef nctf_iterator<T> iterator;
	typedef nctf_const_iterator<T> const_iterator;
	typedef long difference_type;
	typedef unsigned long size_type;
	// typedefs for MT field
	typedef PoomaMesh_XYZ<Mesh> MT_t;

      private:

	// Some additional typedefs for internal use
	typedef CenteredFieldLayout<3,Mesh,Vert> VLayout_t;
	typedef Field<T,3,Mesh,Vert> BaseField_t;
	typedef GuardCellSizes<3> GC_t;

      public:

	// constructors
	nctf(const MT_t::FieldConstructor& spm)
	    : BaseField_t( const_cast<Mesh&>(spm->get_Mesh()),
			   const_cast<VLayout_t&>(spm->get_VLayout()),
			   GC_t(1) ), spm_m(spm)
	{
	    // compute number of local nodes
	    computeSize();
	}
	nctf(const MT_t& m)
	    : BaseField_t( const_cast<Mesh&>(m.get_Mesh()),
			   const_cast<VLayout_t&>(m.get_VLayout()),
			   GC_t(1) ), spm_m(const_cast<MT_t*>(&m))
	{
	    // compute number of local nodes
	    computeSize();
	}
	nctf(const nctf<T>& f)
	    : BaseField_t( const_cast<Mesh&>(f.get_Mesh().get_Mesh()),
			   const_cast<VLayout_t&>(f.get_Mesh().get_VLayout()),
			   GC_t(1) ),
	      spm_m( f.get_FieldConstructor() ), size_m(f.size_m) {}

	// Assignment from scalar
	nctf<T>& operator=(T x)
	{
	    assign(*this,x);
	    return *this;
	}
	// Assignment from another field
	nctf<T>& operator=(const nctf<T>& x)
	{
	    PAssert( size() == x.size() );
	    assign(*this,x);
	    return *this;
	}
	// Assignment from base field
	nctf<T>& operator=(const BaseField_t& x)
	{
	    assign(*this,x);
	    return *this;
	}
	// Assignment from an expression
	template <class B>
	nctf<T>& operator=(const PETE_Expr<B>& x)
	{
	    assign(*this,x);
	    return *this;
	}
	
	// element accessors
	value_type&
	operator()(unsigned int i, unsigned int j, unsigned int k)
	{
	    NDIndex<3> loc(Index(i,i),Index(j,j),Index(k,k));
	    return localElement(loc);
	}
	const value_type&
	operator()(unsigned int i, unsigned int j, unsigned int k) const
	{
	    NDIndex<3> loc(Index(i,i),Index(j,j),Index(k,k));
	    return localElement(loc);
	}
	value_type& operator()(size_type n)
	{
	    unsigned int i,j,k;
	    spm_m->get_node_indices(n,i,j,k);
	    NDIndex<3> loc(Index(i,i),Index(j,j),Index(k,k));
	    return localElement(loc);
	}
	const value_type& operator()(size_type n) const
	{
	    unsigned int i,j,k;
	    spm_m->get_node_indices(n,i,j,k);
	    NDIndex<3> loc(Index(i,i),Index(j,j),Index(k,k));
	    return localElement(loc);
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
	const MT_t::FieldConstructor& get_FieldConstructor(void) const
	{
	    return spm_m;
	}
	const MT_t& get_Mesh(void) const
	{
	    return *spm_m;
	}
	size_type size(void) const { return size_m; }
	size_type max_size(void) const { return size(); }
	bool empty(void) const { return (size() == 0); }

	void swap(nctf<T>& f)
	{
	    // just exchange vmap of pointers to local fields
	    Locals_ac.swap(f.Locals_ac);
	}
	
	// comparisons
	bool operator==(const nctf<T>& rhs) const
	{
	    if (size() != rhs.size()) return false;
	    const_iterator i, iend = end(), rhsi = rhs.begin();
	    bool result = true;
	    for (i = begin(); i != iend; ++i, ++rhsi)
		if (*i != *rhsi) result = false;
	    return result;
	}
	bool operator!=(const nctf<T>& rhs) const { return !(*this==rhs); }

	bool operator<(const nctf<T>& rhs) const
	{
	    const_iterator i, iend, rhsi;
	    if (size() <= rhs.size()) {
		// use this iterator as loop bounds
		iend = end();
		rhsi = rhs.begin();
		for (i = begin(); i != iend; ++i, ++rhsi) {
		    if (*i < *rhsi) return true;
		    if (*rhsi < *i) return false;
		}
		// all elements equal up to this point
		if (size() < rhs.size())
		    return true;
		else
		    return false;
	    }
	    else {
		// use rhs iterator as loop bounds
		iend = rhs.end();
		i = begin();
		for (rhsi = rhs.begin(); rhsi != iend; ++i, ++rhsi) {
		    if (*i < *rhsi) return true;
		    if (*rhsi < *i) return false;
		}
		// all elements equal up to this point
		return false;
	    }
	}
	bool operator>(const nctf<T>& rhs) const { return (rhs < *this); }
	bool operator<=(const nctf<T>& rhs) const { return !(*this > rhs); }
	bool operator>=(const nctf<T>& rhs) const { return !(*this < rhs); }

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
	MT_t::FieldConstructor spm_m;
	size_type size_m;

    };
    
    // Node centered field iterators

    template <class T>
    class nctf_iterator : public forward_iterator<T,ptrdiff_t>
    {
	friend class nctf_const_iterator<T>;

	// typedefs
	typedef typename nctf<T>::BaseField_t BaseField_t;
	typedef typename BaseField_t::iterator BaseField_iterator;

      public:

	// constructors
	nctf_iterator() {}
	nctf_iterator(BaseField_t& bf)
	    : bfi_m(bf.begin()) {}
	nctf_iterator(const BaseField_iterator& bfiter)
	    : bfi_m(bfiter) {}
	nctf_iterator(const nctf_iterator<T>& iter)
	    : bfi_m(iter.bfi_m) {}

	// assignment
	nctf_iterator<T>& operator=(const nctf_iterator<T>& x)
	{
	    bfi_m = x.bfi_m;
	    return *this;
	}

	// pre-increment
	nctf_iterator<T>& operator++()
	{
	    ++bfi_m;
	    return *this;
	}
	// post-increment
	nctf_iterator<T> operator++(int)
	{
	    nctf_iterator<T> iter(*this);
	    ++bfi_m;
	    return iter;
	}
	// dereference
	T& operator*() const
	{
	    return *bfi_m;
	}
	// member
	T* operator->() const
	{
	    return &(*bfi_m);
	}
	
	// comparisons
	bool operator==(const nctf_iterator<T>& rhs) const
	{
	    return (bfi_m == rhs.bfi_m);
	}
	bool operator!=(const nctf_iterator<T>& rhs) const
	{
	    return !(*this == rhs);
	}
	
      private:
	// data
	BaseField_iterator bfi_m;
    };

    template <class T>
    class nctf_const_iterator : public input_iterator<T,ptrdiff_t>
    {
	// typedefs
	typedef typename nctf<T>::BaseField_t BaseField_t;
	typedef typename BaseField_t::iterator BaseField_iterator;

      public:

	// constructors
	nctf_const_iterator() {}
	nctf_const_iterator(const BaseField_t& bf)
	    : bfi_m(bf.begin()) {}
	nctf_const_iterator(const BaseField_iterator& bfiter)
	    : bfi_m(bfiter) {}
	nctf_const_iterator(const nctf_const_iterator<T>& iter)
	    : bfi_m(iter.bfi_m) {}
	nctf_const_iterator(const nctf_iterator<T>& iter)
	    : bfi_m(iter.bfi_m) {}

	// assignment
	nctf_const_iterator<T>&
	operator=(const nctf_const_iterator<T>& x)
	{
	    bfi_m = x.bfi_m;
	    return *this;
	}
	nctf_const_iterator<T>&
	operator=(const nctf_iterator<T>& x)
	{
	    bfi_m = x.bfi_m;
	    return *this;
	}
	
	// pre-increment
	nctf_const_iterator<T>& operator++()
	{
	    ++bfi_m;
	    return *this;
	}
	// post-increment
	nctf_const_iterator<T> operator++(int)
	{
	    nctf_const_iterator<T> iter(*this);
	    ++bfi_m;
	    return iter;
	}
	// dereference
	const T& operator*() const { return (*bfi_m); }
	// member
	const T* operator->() const { return &(*bfi_m); }

	// comparisons
	bool operator==(const nctf_const_iterator<T>& rhs) const
	{
	    return (bfi_m == rhs.bfi_m);
	}
	bool operator!=(const nctf_const_iterator<T>& rhs) const
	{
	    return !(*this == rhs);
	}
	
      private:

	// data
	BaseField_iterator bfi_m;
    };
    
    // Vertex centered field
    // Has a value at each vertex.

    template <class T>
    class vctf : public Field<Vektor<T,8>,3,Mesh,Cell>
    {
	friend class PoomaMesh_XYZ<Mesh>;
	friend class vctf_iterator<T>;
	friend class vctf_const_iterator<T>;

      public:

	// typedefs for container
	typedef T value_type;
	typedef T& reference;
	typedef const T& const_reference;
	typedef T* pointer;
	typedef const T* const_pointer;
	typedef vctf_iterator<T> iterator;
	typedef vctf_const_iterator<T> const_iterator;
	typedef long difference_type;
	typedef unsigned long size_type;
	// typedefs for MT field
	typedef PoomaMesh_XYZ<Mesh> MT_t;

      private:

	// Some additional typedefs for internal use
	typedef CenteredFieldLayout<3,Mesh,Cell> CLayout_t;
	typedef Field<Vektor<T,8>,3,Mesh,Cell> BaseField_t;
	typedef GuardCellSizes<3> GC_t;

      public:

	// constructors
	vctf(const MT_t::FieldConstructor& spm)
	    : BaseField_t( const_cast<Mesh&>(spm->get_Mesh()),
			   const_cast<CLayout_t&>(spm->get_CLayout()),
			   GC_t(1) ),
	      spm_m(spm) {}
	vctf(const MT_t& m)
	    : BaseField_t( const_cast<Mesh&>(m.get_Mesh()),
			   const_cast<CLayout_t&>(m.get_CLayout()),
			   GC_t(1) ),
	      spm_m(const_cast<MT_t*>(&m)) {}
	vctf(const vctf<T>& f)
	    : BaseField_t( const_cast<Mesh&>(f.get_Mesh().get_Mesh()),
			   const_cast<CLayout_t&>(f.get_Mesh().get_CLayout()),
			   GC_t(1) ),
	      spm_m( f.get_FieldConstructor() ) {}

	// Assignment from scalar
	vctf<T>& operator=(T x)
	{
	    Vektor<T,8> v = x;
	    assign(*this,v);
	    return *this;
	}
	// Assignment from another field
	vctf<T>& operator=(const vctf<T>& x)
	{
	    PAssert( size() == x.size() );
	    assign(*this,x);
	    return *this;
	}
	// Assignment from base field
	vctf<T>& operator=(const BaseField_t& x)
	{
	    assign(*this,x);
	    return *this;
	}
	// Assignment from an expression
	template <class B>
	vctf<T>& operator=(const PETE_Expr<B>& x)
	{
	    assign(*this,x);
	    return *this;
	}
	// AddAssign
	template <class Expr>
	vctf<T>& operator+=(const Expr &expr)
	{
	    *this = *this + expr;
	    return *this;
	}
	// SubAssign
	template <class Expr>
	vctf<T>& operator-=(const Expr &expr)
	{
	    *this = *this - expr;
	    return *this;
	}
	// MulAssign
	template <class Expr>
	vctf<T>& operator*=(const Expr &expr)
	{
	    *this = *this * expr;
	    return *this;
	}
	// DivAssign
	template <class Expr>
	vctf<T>& operator/=(const Expr &expr)
	{
	    *this = *this / expr;
	    return *this;
	}
	
	// vertex accessors
	// i, j, k == global xyz cell indices
	// c == cell number
	// v == vertex index
	cctf<T> operator()(unsigned int v) const
	{
	    cctf<T> cf(get_FieldConstructor());
	    typename cctf<T>::iterator cfit, cfend = cf.end();
	    const BaseField_t& bf = dynamic_cast<const BaseField_t&>(*this);
	    typename BaseField_t::iterator bfit = bf.begin();
	    for (cfit = cf.begin(); cfit != cfend; ++cfit, ++bfit)
		*cfit = (*bfit)(v);
	    return cf;
	}
	value_type&
	operator()(unsigned int i, unsigned int j, unsigned int k,
		   unsigned int v)
	{
	    NDIndex<3> loc(Index(i,i),Index(j,j),Index(k,k));
	    return localElement(loc)(v);
	}
	const value_type&
	operator()(unsigned int i, unsigned int j, unsigned int k,
		   unsigned int v) const
	{
	    NDIndex<3> loc(Index(i,i),Index(j,j),Index(k,k));
	    return localElement(loc)(v);
	}
	value_type& operator()(size_type c, unsigned int v)
	{
	    unsigned int i,j,k;
	    spm_m->get_cell_indices(c,i,j,k);
	    NDIndex<3> loc(Index(i,i),Index(j,j),Index(k,k));
	    return localElement(loc)(v);
	}
	const value_type& operator()(size_type c, unsigned int v) const
	{
	    unsigned int i,j,k;
	    spm_m->get_cell_indices(c,i,j,k);
	    NDIndex<3> loc(Index(i,i),Index(j,j),Index(k,k));
	    return localElement(loc)(v);
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
	const MT_t::FieldConstructor& get_FieldConstructor(void) const
	{
	    return spm_m;
	}
	const MT_t& get_Mesh(void) const
	{
	    return *spm_m;
	}
	size_type size(void) const { return (spm_m->get_ncells() * 8); }
	size_type max_size(void) const { return size(); }
	bool empty(void) const { return (size() == 0); }

	void swap(vctf<T>& f)
	{
	    // just exchange vmap of pointers to local fields
	    Locals_ac.swap(f.Locals_ac);
	}
	
	// comparisons
	bool operator==(const vctf<T>& rhs) const
	{
	    if (size() != rhs.size()) return false;
	    const_iterator i, iend = end(), rhsi = rhs.begin();
	    bool result = true;
	    for (i = begin(); i != iend; ++i, ++rhsi)
		if (*i != *rhsi) result = false;
	    return result;
	}
	bool operator!=(const vctf<T>& rhs) const { return !(*this==rhs); }

	bool operator<(const vctf<T>& rhs) const
	{
	    const_iterator i, iend, rhsi;
	    if (size() <= rhs.size()) {
		// use this iterator as loop bounds
		iend = end();
		rhsi = rhs.begin();
		for (i = begin(); i != iend; ++i, ++rhsi) {
		    if (*i < *rhsi) return true;
		    if (*rhsi < *i) return false;
		}
		// all elements equal up to this point
		if (size() < rhs.size())
		    return true;
		else
		    return false;
	    }
	    else {
		// use rhs iterator as loop bounds
		iend = rhs.end();
		i = begin();
		for (rhsi = rhs.begin(); rhsi != iend; ++i, ++rhsi) {
		    if (*i < *rhsi) return true;
		    if (*rhsi < *i) return false;
		}
		// all elements equal up to this point
		return false;
	    }
	}
	bool operator>(const vctf<T>& rhs) const { return (rhs < *this); }
	bool operator<=(const vctf<T>& rhs) const { return !(*this > rhs); }
	bool operator>=(const vctf<T>& rhs) const { return !(*this < rhs); }

      private:
	// data
	MT_t::FieldConstructor spm_m;

    };
    
    // Vertex centered field iterators

    template <class T>
    class vctf_iterator : public forward_iterator<T,ptrdiff_t>
    {
	friend class vctf_const_iterator<T>;

	// typedefs
	typedef typename vctf<T>::BaseField_t BaseField_t;
	typedef typename BaseField_t::iterator BaseField_iterator;

      public:

	// constructors
	vctf_iterator() : vertex_m(0) {}
	vctf_iterator(BaseField_t& bf)
	    : bfi_m(bf.begin()), vertex_m(0) {}
	vctf_iterator(const BaseField_iterator& bfiter)
	    : bfi_m(bfiter), vertex_m(0) {}
	vctf_iterator(const vctf_iterator<T>& iter)
	    : bfi_m(iter.bfi_m), vertex_m(iter.vertex_m) {}

	// assignment
	vctf_iterator<T>& operator=(const vctf_iterator<T>& x)
	{
	    bfi_m = x.bfi_m;
	    vertex_m = x.vertex_m;
	    return *this;
	}
	
	// pre-increment
	vctf_iterator<T>& operator++()
	{
	    ++vertex_m;
	    if (vertex_m==8) nextCell();
	    return *this;
	}
	// post-increment
	vctf_iterator<T> operator++(int)
	{
	    vctf_iterator<T> iter(*this);
	    ++vertex_m;
	    if (vertex_m==8) nextCell();
	    return iter;
	}
	// dereference
	T& operator*() const { return (*bfi_m)(vertex_m); }
	// member
	T* operator->() const { return &((*bfi_m)(vertex_m)); }

	// comparisons
	bool operator==(const vctf_iterator<T>& rhs) const
	{
	    return (bfi_m == rhs.bfi_m && vertex_m == rhs.vertex_m);
	}
	bool operator!=(const vctf_iterator<T>& rhs) const
	{
	    return !(*this == rhs);
	}

      private:

	// methods
	void nextCell()
	{
	    ++bfi_m;
	    vertex_m = 0;
	}
	
	// data
	BaseField_iterator bfi_m;
	unsigned int vertex_m;
    };

    template <class T>
    class vctf_const_iterator : public input_iterator<T,ptrdiff_t>
    {
	// typedefs
	typedef typename vctf<T>::BaseField_t BaseField_t;
	typedef typename BaseField_t::iterator BaseField_iterator;

      public:

	// constructors
	vctf_const_iterator() : vertex_m(0) {}
	vctf_const_iterator(const BaseField_t& bf)
	    : bfi_m(bf.begin()), vertex_m(0) {}
	vctf_const_iterator(const BaseField_iterator& bfiter)
	    : bfi_m(bfiter), vertex_m(0) {}
	vctf_const_iterator(const vctf_const_iterator<T>& iter)
	    : bfi_m(iter.bfi_m), vertex_m(iter.vertex_m) {}
	vctf_const_iterator(const vctf_iterator<T>& iter)
	    : bfi_m(iter.bfi_m), vertex_m(iter.vertex_m) {}

	// assignment
	vctf_const_iterator<T>&
	operator=(const vctf_const_iterator<T>& x)
	{
	    bfi_m = x.bfi_m;
	    vertex_m = x.vertex_m;
	    return *this;
	}
	vctf_const_iterator<T>&
	operator=(const vctf_iterator<T>& x)
	{
	    bfi_m = x.bfi_m;
	    vertex_m = x.vertex_m;
	    return *this;
	}

	// pre-increment
	vctf_const_iterator<T>& operator++()
	{
	    ++vertex_m;
	    if (vertex_m==8) nextCell();
	    return *this;
	}
	// post-increment
	vctf_const_iterator<T> operator++(int)
	{
	    vctf_const_iterator<T> iter(*this);
	    ++vertex_m;
	    if (vertex_m==8) nextCell();
	    return iter;
	}
	// dereference
	const T& operator*() const { return (*bfi_m)(vertex_m); }
	// member
	const T* operator->() const { return &((*bfi_m)(vertex_m)); }
	
	// comparisons
	bool operator==(const vctf_const_iterator<T>& rhs) const
	{
	    return (bfi_m == rhs.bfi_m && vertex_m == rhs.vertex_m);
	}
	bool operator!=(const vctf_const_iterator<T>& rhs) const
	{
	    return !(*this == rhs);
	}
	
      private:

	// methods
	void nextCell()
	{
	    ++bfi_m;
	    vertex_m = 0;
	}

	// data
	BaseField_iterator bfi_m;
	unsigned int vertex_m;
    };

    // Special accessor classes for face-centered discontinuous field

    template <class FT>
    class ConnFacesAroundVertices
    {
      public:

	// typedefs
	typedef VertexProxy<FT> value_type;
	typedef CFAV_iterator<FT> iterator;
	typedef CFAV_const_iterator<FT> const_iterator;

	// constructors
	ConnFacesAroundVertices(const FT& field)
	    : field_m(const_cast<FT&>(field))
	{
	    numCells_m = field_m.size() / 6;
	}
	ConnFacesAroundVertices(const ConnFacesAroundVertices<FT>& model)
	    : field_m(model.field_m)
	{
	    numCells_m = field_m.size() / 6;
	}

	// destructor
	~ConnFacesAroundVertices() { }

	// assignment
	ConnFacesAroundVertices<FT>&
	operator=(const ConnFacesAroundVertices<FT>& rhs)
	{
	    field_m = rhs.field_m;
	    numCells_m = rhs.numCells_m;
	    return *this;
	}
	
	// accessors
	iterator begin()
	{
	    return iterator(field_m, 0, 0);
	}
	iterator end()
	{
	    return iterator(field_m, numCells_m, 0);
	}
	const_iterator begin() const
	{
	    return const_iterator(field_m, 0, 0);
	}
	const_iterator end() const
	{
	    return const_iterator(field_m, numCells_m, 0);
	}
	
      private:
	FT& field_m;
	typename FT::size_type numCells_m;
    };
    
    // CFAV iterators

    template <class FT>
    class CFAV_iterator
	: public forward_iterator<VertexProxy<FT>,typename FT::difference_type>
    {
	friend class CFAV_const_iterator<FT>;

      public:

	// constructors
	CFAV_iterator(FT& field,
		      typename FT::size_type cell,
		      unsigned int vertex)
	    : field_m(field),
	      cell_m(cell),
	      vertex_m(vertex),
	      proxy_m(field, cell, vertex) { }
	CFAV_iterator(const CFAV_iterator<FT>& model)
	    : field_m(model.field_m),
	      cell_m(model.cell_m),
	      vertex_m(model.vertex_m),
	      proxy_m(model.proxy_m) { }

	// assignment
	CFAV_iterator<FT>& operator=(const CFAV_iterator<FT>& rhs)
	{
	    field_m = rhs.field_m;
	    cell_m = rhs.cell_m;
	    vertex_m = rhs.vertex_m;
	    proxy_m = rhs.proxy_m;
	    return *this;
	}

	// pre-increment
	CFAV_iterator<FT>& operator++()
	{
	    ++vertex_m;
	    if (vertex_m==8) nextCell();
	    proxy_m.assign(cell_m, vertex_m);
	    return *this;
	}
	
	// post-increment
	CFAV_iterator<FT> operator++(int)
	{
	    CFAV_iterator<FT> iter(*this);
	    ++vertex_m;
	    if (vertex_m==8) nextCell();
	    proxy_m.assign(cell_m, vertex_m);
	    return iter;
	}

	// dereference
	VertexProxy<FT>& operator*() const { return proxy_m; }
	// member
	VertexProxy<FT>* operator->() const { return &proxy_m; }

	// comparisons
	bool operator==(const CFAV_iterator<FT>& rhs) const
	{
	    return (cell_m == rhs.cell_m && vertex_m == rhs.vertex_m);
	}
	bool operator!=(const CFAV_iterator<FT>& rhs) const
	{
	    return !(*this == rhs);
	}

      private:

	// methods
	void nextCell()
	{
	    ++cell_m;
	    vertex_m = 0;
	}
	
	// data
	FT& field_m;
	typename FT::size_type cell_m;
	unsigned int vertex_m;
	mutable VertexProxy<FT> proxy_m;
    };

    template <class FT>
    class CFAV_const_iterator
	: public input_iterator<VertexProxy<FT>,typename FT::difference_type>
    {
      public:

	// constructors
	CFAV_const_iterator(FT& field,
			    typename FT::size_type cell,
			    unsigned int vertex)
	    : field_m(field),
	      cell_m(cell),
	      vertex_m(vertex),
	      proxy_m(field, cell, vertex) { }
	CFAV_const_iterator(const CFAV_iterator<FT>& model)
	    : field_m(model.field_m),
	      cell_m(model.cell_m),
	      vertex_m(model.vertex_m),
	      proxy_m(model.proxy_m) { }
	CFAV_const_iterator(const CFAV_const_iterator<FT>& model)
	    : field_m(model.field_m),
	      cell_m(model.cell_m),
	      vertex_m(model.vertex_m),
	      proxy_m(model.proxy_m) { }

	// assignment
	CFAV_const_iterator<FT>& operator=(const CFAV_const_iterator<FT>& rhs)
	{
	    field_m = rhs.field_m;
	    cell_m = rhs.cell_m;
	    vertex_m = rhs.vertex_m;
	    proxy_m = rhs.proxy_m;
	    return *this;
	}
	CFAV_const_iterator<FT>& operator=(const CFAV_iterator<FT>& rhs)
	{
	    field_m = rhs.field_m;
	    cell_m = rhs.cell_m;
	    vertex_m = rhs.vertex_m;
	    proxy_m = rhs.proxy_m;
	    return *this;
	}

	// pre-increment
	CFAV_const_iterator<FT>& operator++()
	{
	    ++vertex_m;
	    if (vertex_m==8) nextCell();
	    proxy_m.assign(cell_m, vertex_m);
	    return *this;
	}
	
	// post-increment
	CFAV_const_iterator<FT> operator++(int)
	{
	    CFAV_const_iterator<FT> iter(*this);
	    ++vertex_m;
	    if (vertex_m==8) nextCell();
	    proxy_m.assign(cell_m, vertex_m);
	    return iter;
	}

	// dereference
	const VertexProxy<FT>& operator*() const { return proxy_m; }
	// member
	const VertexProxy<FT>* operator->() const { return &proxy_m; }

	// comparisons
	bool operator==(const CFAV_const_iterator<FT>& rhs) const
	{
	    return (cell_m == rhs.cell_m && vertex_m == rhs.vertex_m);
	}
	bool operator!=(const CFAV_const_iterator<FT>& rhs) const
	{
	    return !(*this == rhs);
	}
	
      private:

	// methods
	void nextCell()
	{
	    ++cell_m;
	    vertex_m = 0;
	}
	
	// data
	FT& field_m;
	typename FT::size_type cell_m;
	unsigned int vertex_m;
	mutable VertexProxy<FT> proxy_m;
    };
    
    template <class FT>
    class VertexProxy
    {
	friend class CFAV_iterator<FT>;
	friend class CFAV_const_iterator<FT>;
	struct FacesAtVertex;

      public:

	// typedefs
	typedef typename FT::value_type value_type;

	// nested classes

	class iterator
	    : public forward_iterator<value_type,typename FT::difference_type>
	{
	  public:

	    // constructors
	    iterator(VertexProxy<FT>& proxy, unsigned int face)
		: proxy_m(proxy), face_m(face) { }
	    iterator(const iterator& model)
		: proxy_m(model.proxy_m), face_m(model.face_m) { }

	    // comparisons
	    bool operator==(const iterator& rhs) const
	    {
		return &proxy_m == &rhs.proxy_m && face_m == rhs.face_m;
	    }
	    bool operator!=(const iterator& rhs) const
	    {
		return !(*this == rhs);
	    }

	    // dereference
	    typename iterator::reference operator*() const
	    {
		return proxy_m.field_m(proxy_m.cell_m,
				       proxy_m.faces_m.dat[face_m]);
	    }
	    // member
	    typename iterator::pointer operator->() const
	    {
		return &(proxy_m.field_m(proxy_m.cell_m,
					 proxy_m.faces_m.dat[face_m]));
	    }

	    // increment
	    iterator& operator++()
	    {
		face_m++;
		return *this;
	    }
	    iterator operator++(int)
	    {
		iterator tmp(*this);
		face_m++;
		return tmp;
	    }
	    
	  private:
	    VertexProxy<FT>& proxy_m;
	    unsigned int face_m;
	};

	class const_iterator
	    : public input_iterator<value_type,typename FT::difference_type>
	{
	  public:

	    // constructors
	    const_iterator(VertexProxy<FT>& proxy, unsigned int face)
		: proxy_m(proxy), face_m(face) { }
	    const_iterator(const const_iterator& model)
		: proxy_m(model.proxy_m), face_m(model.face_m) { }

	    // comparisons
	    bool operator==(const const_iterator& rhs) const
	    {
		return &proxy_m == &rhs.proxy_m && face_m == rhs.face_m;
	    }
	    bool operator!=(const const_iterator& rhs) const
	    {
		return !(*this == rhs);
	    }
	    
	    // dereference
	    const value_type& operator*() const
	    {
		return proxy_m.field_m(proxy_m.cell_m,
				       proxy_m.faces_m.dat[face_m]);
	    }
	    // member
	    const value_type* operator->() const
	    {
		return &(proxy_m.field_m(proxy_m.cell_m,
					 proxy_m.faces_m.dat[face_m]));
	    }
	
	    // increment
	    const_iterator& operator++()
	    {
		face_m++;
		return *this;
	    }
	    const_iterator operator++(int)
	    {
		const_iterator tmp(*this);
		face_m++;
		return tmp;
	    }

	  private:
	    VertexProxy<FT>& proxy_m;
	    unsigned int face_m;
	};
	
	// constructors
	VertexProxy(FT& field,
		    typename FT::size_type cell,
		    unsigned int vertex)
	    : field_m(field), cell_m(cell), vertex_m(vertex), faces_m(vertex)
	{
	}
	
	VertexProxy(const VertexProxy<FT>& model)
	    : field_m(model.field_m),
	      cell_m(model.cell_m),
	      vertex_m(model.vertex_m),
	      faces_m(model.faces_m)
	{
	}
	
	// assignment
	VertexProxy<FT>& operator=(const VertexProxy<FT>& rhs)
	{
	    field_m = rhs.field_m;
	    cell_m = rhs.cell_m;
	    vertex_m = rhs.vertex_m;
	    faces_m = rhs.faces_m;
	    return *this;
	}

	// accessors
	iterator begin()
	{
	    return iterator(*this, 0);
	}
	iterator end()
	{
	    return iterator(*this, 3);
	}

	const_iterator begin() const
	{
	    return const_iterator(const_cast<VertexProxy<FT>&>(*this), 0);
	    
	}
	const_iterator end() const
	{
	    return const_iterator(const_cast<VertexProxy<FT>&>(*this), 3);
	    
	}
	
      private:
	friend class iterator;
	friend class const_iterator;

	// nested class
	struct FacesAtVertex
	{
	    static const int vf[8][3];
	    int dat[3];

	    FacesAtVertex(int iv)
	    {
		std::copy(vf[iv], vf[iv+1], dat);
	    }
	};

	//  const int FacesAtVertex::vf = {
	//	{0,2,4}, {1,2,4}, {0,3,4}, {1,3,4},
	//	{0,2,5}, {1,2,5}, {0,3,5}, {1,3,5}
	//  };

	// methods
	void assign(typename FT::size_type cell, unsigned int vertex)
	{
	    cell_m = cell;
	    vertex_m = vertex;
	    faces_m = FacesAtVertex(vertex_m);
	}
	
	// data
	FT& field_m;
	typename FT::size_type cell_m;
	unsigned int vertex_m;
	FacesAtVertex faces_m;
    };
    
  private:

    // Some additional typedefs for internal use
    typedef CenteredFieldLayout<3,Mesh,Cell> CLayout_t;
    typedef CenteredFieldLayout<3,Mesh,Vert> VLayout_t;

    // data
    Mesh* Mesh_m;
    CLayout_t* CLayout_m;
    VLayout_t* VLayout_m;
    Field<Vektor<double,3>,3,Cartesian<3,double>,Vert>*	VertSizes_m;
    Field<Vektor<double,3>,3,Cartesian<3,double>,Cell>*	CellSizes_m;
    Field<Vektor<double,3>,3,Cartesian<3,double>,Cell>*	CellPositions_m;
    //	typename ncvsf::BaseField_t* VertSizes_m;
    //	typename ccvsf::BaseField_t* CellSizes_m;
    //	typename ccvsf::BaseField_t* CellPositions_m;
    typename ccsf::BaseField_t* CellVolumes_m;
    Field<Vektor<double,3>,3,Cartesian<3,double>,Cell>**	CellSurfaceNormals_m;
    //	typename ccvsf::BaseField_t** CellSurfaceNormals_m;
    unsigned int ncx_m, ncy_m, ncz_m;
    size_type ncp_m;

  public:

    // constructors
#if 0
    PoomaMesh_XYZ(int* const ncells,
		  typename Mesh::MeshValue_t* const cwidth,
		  e_dim_tag* decomp)
    {
	// store global mesh sizes
	ncx_m = ncells[0];
	ncy_m = ncells[1];
	ncz_m = ncells[2];
	// create domain for mesh vertices (one larger than no. of cells)
	NDIndex<3> meshDom;
	for (unsigned int d=0; d<3; ++d)
	    meshDom[d] = Index(ncells[d]+1);
	// create the mesh
	Mesh_m = new Mesh(meshDom,cwidth);
	// tell mesh to store mesh spacing fields
	Mesh_m->storeSpacingFields(decomp);
	// create the field layouts
	CLayout_m = new CLayout_t(*Mesh_m, decomp);
	VLayout_m = new VLayout_t(*Mesh_m, decomp);
	// compute number of cells on this processor
	ncp_m = 0;
	typename CLayout_t::iterator_iv cli, clend = CLayout_m->end_iv();
	for (cli = CLayout_m->begin_iv(); cli != clend; ++cli)
	    ncp_m += (*cli).second->getDomain().size();
	// initialize internal fields
	initializeVertSizes();
	initializeCellSizes();
	initializeCellPositions();
	initializeCellVolumes();
	initializeCellSurfaceNormals();
    }
#endif
    PoomaMesh_XYZ(int* const ncells,
		  typename Mesh::MeshValue_t** const cwidth,
		  e_dim_tag* decomp)
    {
	// store global mesh sizes
	ncx_m = ncells[0];
	ncy_m = ncells[1];
	ncz_m = ncells[2];
	// create domain for mesh vertices (one larger than no. of cells)
	NDIndex<3> meshDom;
	for (unsigned int d=0; d<3; ++d)
	    meshDom[d] = Index(ncells[d]+1);
	// create the mesh
	Mesh_m = new Mesh(meshDom,cwidth);
	// tell mesh to store mesh spacing fields
	Mesh_m->storeSpacingFields(decomp);
	// create the field layouts
	CLayout_m = new CLayout_t(*Mesh_m, decomp);
	VLayout_m = new VLayout_t(*Mesh_m, decomp);
	// compute number of cells on this processor
	ncp_m = 0;
	typename CLayout_t::iterator_iv cli, clend = CLayout_m->end_iv();
	for (cli = CLayout_m->begin_iv(); cli != clend; ++cli)
	    ncp_m += (*cli).second->getDomain().size();
	// initialize internal fields
	initializeVertSizes();
	initializeCellSizes();
	initializeCellPositions();
	initializeCellVolumes();
	initializeCellSurfaceNormals();
    }
    
    // destructor
    ~PoomaMesh_XYZ()
    {
	delete VertSizes_m;
	delete CellSizes_m;
	delete CellPositions_m;
	delete CellVolumes_m;
	for (unsigned int f=0; f<6; ++f) delete CellSurfaceNormals_m[f];
	delete [] CellSurfaceNormals_m;
	delete CLayout_m;
	delete VLayout_m;
	delete Mesh_m;
    }

    // Equality Comparable
    bool operator==(const PoomaMesh_XYZ<Mesh>& rhs) const
    {
	return (this == &rhs);
    }
    bool operator!=(const PoomaMesh_XYZ<Mesh>& rhs) const
    {
	return !(*this == rhs);
    }

    // accessors
    size_type get_ncells() const
    {
	return ncp_m;
    }
    size_type get_total_ncells() const
    {
	return (ncx_m * ncy_m * ncz_m);
    }

    void get_dx(ccsf& f) const
    {
	PAssert(f.get_Mesh() == *this);
	Field<Vektor<double,3>,3,Cartesian<3,double>,Cell>::iterator csi, csend =
	    CellSizes_m->end();
	typename ccsf::iterator fi = f.begin();
	for (csi = CellSizes_m->begin(); csi != csend; ++csi, ++fi)
	    *fi = (*csi)(0);
    }
    void get_dy(ccsf& f) const
    {
	PAssert(f.get_Mesh() == *this);
	Field<Vektor<double,3>,3,Cartesian<3,double>,Cell>::iterator csi, csend =
	    CellSizes_m->end();
	typename ccsf::iterator fi = f.begin();
	for (csi = CellSizes_m->begin(); csi != csend; ++csi, ++fi)
	    *fi = (*csi)(1);
    }
    void get_dz(ccsf& f) const
    {
	PAssert(f.get_Mesh() == *this);
	Field<Vektor<double,3>,3,Cartesian<3,double>,Cell>::iterator csi, csend =
	    CellSizes_m->end();
	typename ccsf::iterator fi = f.begin();
	for (csi = CellSizes_m->begin(); csi != csend; ++csi, ++fi)
	    *fi = (*csi)(2);
    }
    
    void get_xloc(ccsf& f) const
    {
	PAssert(f.get_Mesh() == *this);
	Field<Vektor<double,3>,3,Cartesian<3,double>,Cell>::iterator cpi,
	    cpend = CellPositions_m->end();
	typename ccsf::iterator fi = f.begin();
	for (cpi = CellPositions_m->begin(); cpi != cpend; ++cpi, ++fi)
	    *fi = (*cpi)(0);
    }
    void get_yloc(ccsf& f) const
    {
	PAssert(f.get_Mesh() == *this);
	Field<Vektor<double,3>,3,Cartesian<3,double>,Cell>::iterator cpi,
	    cpend = CellPositions_m->end();
	typename ccsf::iterator fi = f.begin();
	for (cpi = CellPositions_m->begin(); cpi != cpend; ++cpi, ++fi)
	    *fi = (*cpi)(1);
    }
    void get_zloc(ccsf& f) const
    {
	PAssert(f.get_Mesh() == *this);
	Field<Vektor<double,3>,3,Cartesian<3,double>,Cell>::iterator cpi,
	    cpend = CellPositions_m->end();
	typename ccsf::iterator fi = f.begin();
	for (cpi = CellPositions_m->begin(); cpi != cpend; ++cpi, ++fi)
	    *fi = (*cpi)(2);
    }
    void get_xloc(fcdsf& f) const
    {
	PAssert(f.get_Mesh() == *this);
	Field<Vektor<double,3>,3,Cartesian<3,double>,Cell>::iterator cpi,
	    cpend = CellPositions_m->end(), csi = CellSizes_m->begin();
	typename fcdsf::iterator fi = f.begin();
	double xpos, xsize;
	for (cpi = CellPositions_m->begin(); cpi != cpend; ++cpi, ++csi) {
	    xpos = (*cpi)(0);
	    xsize = (*csi)(0);
	    *fi = xpos - 0.5 * xsize; ++fi;
	    *fi = xpos + 0.5 * xsize; ++fi;
	    *fi = xpos; ++fi;
	    *fi = xpos; ++fi;
	    *fi = xpos; ++fi;
	    *fi = xpos; ++fi;
	}
    }
    void get_yloc(fcdsf& f) const
    {
	PAssert(f.get_Mesh() == *this);
	Field<Vektor<double,3>,3,Cartesian<3,double>,Cell>::iterator cpi,
	    cpend = CellPositions_m->end(), csi = CellSizes_m->begin();
	typename fcdsf::iterator fi = f.begin();
	double ypos, ysize;
	for (cpi = CellPositions_m->begin(); cpi != cpend; ++cpi, ++csi) {
	    ypos = (*cpi)(1);
	    ysize = (*csi)(1);
	    *fi = ypos; ++fi;
	    *fi = ypos; ++fi;
	    *fi = ypos - 0.5 * ysize; ++fi;
	    *fi = ypos + 0.5 * ysize; ++fi;
	    *fi = ypos; ++fi;
	    *fi = ypos; ++fi;
	}
    }
    void get_zloc(fcdsf& f) const
    {
	PAssert(f.get_Mesh() == *this);
	Field<Vektor<double,3>,3,Cartesian<3,double>,Cell>::iterator cpi,
	    cpend = CellPositions_m->end(), csi = CellSizes_m->begin();
	typename fcdsf::iterator fi = f.begin();
	double zpos, zsize;
	for (cpi = CellPositions_m->begin(); cpi != cpend; ++cpi, ++csi) {
	    zpos = (*cpi)(2);
	    zsize = (*csi)(2);
	    *fi = zpos; ++fi;
	    *fi = zpos; ++fi;
	    *fi = zpos; ++fi;
	    *fi = zpos; ++fi;
	    *fi = zpos - 0.5 * zsize; ++fi;
	    *fi = zpos + 0.5 * zsize; ++fi;
	}
    }
    
    void get_face_normals(fcdvsf& f) const
    {
	PAssert(f.get_Mesh() == *this);
	Field<Vektor<double,3>,3,Cartesian<3,double>,Cell>::iterator csn0,
	    csn1 = CellSurfaceNormals_m[1]->begin(),
	    csn2 = CellSurfaceNormals_m[2]->begin(),
	    csn3 = CellSurfaceNormals_m[3]->begin(),
	    csn4 = CellSurfaceNormals_m[4]->begin(),
	    csn5 = CellSurfaceNormals_m[5]->begin(),
	    csnend = CellSurfaceNormals_m[0]->end();
	typename fcdvsf::iterator fi = f.begin();
	for (csn0 = CellSurfaceNormals_m[0]->begin(); csn0 != csnend;
	     ++csn0, ++csn1, ++csn2, ++csn3, ++csn4, ++csn5) {
	    *fi = *csn0; ++fi;
	    *fi = *csn1; ++fi;
	    *fi = *csn2; ++fi;
	    *fi = *csn3; ++fi;
	    *fi = *csn4; ++fi;
	    *fi = *csn5; ++fi;
	}
    }
    void get_face_areas(fcdsf& f) const
    {
	PAssert(f.get_Mesh() == *this);
	Field<Vektor<double,3>,3,Cartesian<3,double>,Cell>::iterator csi, csend =
	    CellSizes_m->end();
	typename fcdsf::iterator fi = f.begin();
	for (csi = CellSizes_m->begin(); csi != csend; ++csi) {
	    *fi = (*csi)(1) * (*csi)(2); ++fi;
	    *fi = (*csi)(1) * (*csi)(2); ++fi;
	    *fi = (*csi)(2) * (*csi)(0); ++fi;
	    *fi = (*csi)(2) * (*csi)(0); ++fi;
	    *fi = (*csi)(0) * (*csi)(1); ++fi;
	    *fi = (*csi)(0) * (*csi)(1); ++fi;
	}
    }
    void get_face_lengths(fcdsf& f) const
    {
	PAssert(f.get_Mesh() == *this);
	Field<Vektor<double,3>,3,Cartesian<3,double>,Cell>::iterator csi, csend =
	    CellSizes_m->end();
	typename fcdsf::iterator fi = f.begin();
	for (csi = CellSizes_m->begin(); csi != csend; ++csi) {
	    *fi = (*csi)(0); ++fi;
	    *fi = (*csi)(0); ++fi;
	    *fi = (*csi)(1); ++fi;
	    *fi = (*csi)(1); ++fi;
	    *fi = (*csi)(2); ++fi;
	    *fi = (*csi)(2); ++fi;
	}
    }

    void get_cell_volumes(ccsf& f) const
    {
	PAssert(f.get_Mesh() == *this);
	f = (*CellVolumes_m);
    }
    void get_vertex_volumes(vcsf& f) const
    {
	PAssert(f.get_Mesh() == *this);
	typename ccsf::BaseField_t::iterator cvi, cvend = CellVolumes_m->end();
	typename vcsf::iterator fi = f.begin();
	for (cvi = CellVolumes_m->begin(); cvi != cvend; ++cvi) {
	    *fi = *cvi / 8.0; ++fi;
	    *fi = *cvi / 8.0; ++fi;
	    *fi = *cvi / 8.0; ++fi;
	    *fi = *cvi / 8.0; ++fi;
	    *fi = *cvi / 8.0; ++fi;
	    *fi = *cvi / 8.0; ++fi;
	    *fi = *cvi / 8.0; ++fi;
	    *fi = *cvi / 8.0; ++fi;
	}
    }
    void get_node_volumes(ncsf& f) const
    {
	PAssert(f.get_Mesh() == *this);
	Field<Vektor<double,3>,3,Cartesian<3,double>,Vert>::iterator vsi, vsend =
	    VertSizes_m->end();
	typename ncsf::iterator fi = f.begin();
	for (vsi = VertSizes_m->begin(); vsi != vsend; ++vsi, ++fi) {
	    *fi = (*vsi)(0) * (*vsi)(1) * (*vsi)(2);
	}
    }
    
    // additional accessor functions
    const Mesh& get_Mesh() const { return *Mesh_m; }
    Mesh& get_Mesh() { return *Mesh_m; }
    const CLayout_t& get_CLayout() const { return *CLayout_m; }
    CLayout_t& get_CLayout() { return *CLayout_m; }
    const VLayout_t& get_VLayout() const { return *VLayout_m; }
    VLayout_t& get_VLayout() { return *VLayout_m; }

    unsigned int get_ncx() const { return ncx_m; }
    unsigned int get_ncy() const { return ncy_m; }
    unsigned int get_ncz() const { return ncz_m; }

    // conversion function for cell numbering
    void
    get_cell_indices(size_type ncell,
		     unsigned int& i, unsigned int& j, unsigned int& k) const
    {
	size_type res = ncell;
	k = res % (ncx_m * ncy_m);
	res -= k * (ncx_m * ncy_m);
	j = res % ncx_m;
	res -= j * ncx_m;
	i = res;
	return;
    }
    
    // conversion function for node numbering
    void
    get_node_indices(size_type nnode,
		     unsigned int& i, unsigned int& j, unsigned int& k) const
    {
	size_type res = nnode;
	k = res % ( (ncx_m+1) * (ncy_m+1) );
	res -= k * ( (ncx_m+1) * (ncy_m+1) );
	j = res % (ncx_m+1);
	res -= j * (ncx_m+1);
	i = res;
	return;
    }

    // gather/scatter operations

    template <class T1, class T2, class Op>
    static void scatter(fcdtf<T1>& to, const cctf<T2>& from, const Op& op);

    template <class T1, class T2, class Op>
    static void scatter(cctf<T1>& to, const fcdtf<T2>& from, const Op& op);

    template <class T1, class T2, class Op>
    static void scatter(fcdtf<T1>& to, const vctf<T2>& from, const Op& op);

    template <class T1, class T2, class Op>
    static void scatter(vctf<T1>& to, const fcdtf<T2>& from, const Op& op);

    template <class T1, class T2, class Op>
    static void scatter(nctf<T1>& to, const vctf<T2>& from, const Op& op);

    template <class T1, class T2, class Op>
    static void scatter(cctf<T1>& to, const vctf<T2>& from, const Op& op);

    template <class T1, class T2, class Op>
    static void gather(fcdtf<T1>& to, const cctf<T2>& from, const Op& op);

    template <class T1, class T2, class Op>
    static void gather(bstf<T1>& to, const fcdtf<T2>& from, const Op& op);

    template <class T1, class T2, class Op>
    static void gather(fcdtf<T1>& to, const bstf<T2>& from, const Op& op);

    template <class T1, class T2, class Op>
    static void gather(vctf<T1>& to, const nctf<T2>& from, const Op& op);

    template <class T1, class T2, class Op>
    static void gather(vctf<T1>& to, const cctf<T2>& from, const Op& op);

    template <class T>
    static void swap_faces(fcdtf<T>& to, const fcdtf<T>& from);


    // reduction operations
    template <class T>
    static T sum(const cctf<T>& f);
    template <class T>
    static T sum(const fcdtf<T>& f);
    template <class T>
    static T sum(const bstf<T>& f);
    template <class T>
    static T sum(const nctf<T>& f);
    template <class T>
    static T sum(const vctf<T>& f);

    template <class T>
    static T min(const cctf<T>& f);
    template <class T>
    static T min(const fcdtf<T>& f);
    template <class T>
    static T min(const bstf<T>& f);
    template <class T>
    static T min(const nctf<T>& f);
    template <class T>
    static T min(const vctf<T>& f);

    template <class T>
    static T max(const cctf<T>& f);
    template <class T>
    static T max(const fcdtf<T>& f);
    template <class T>
    static T max(const bstf<T>& f);
    template <class T>
    static T max(const nctf<T>& f);
    template <class T>
    static T max(const vctf<T>& f);

  private:

    // methods
    void initializeVertSizes()
    {
	VertSizes_m =
	    new Field<Vektor<double,3>,3,Cartesian<3,double>,Vert>(*Mesh_m, *VLayout_m);
	Mesh_m->getDeltaCellField(*VertSizes_m);
    }
    void initializeCellSizes()
    {
	CellSizes_m =
	    new Field<Vektor<double,3>,3,Cartesian<3,double>,Cell>(*Mesh_m, *CLayout_m);
	Mesh_m->getDeltaVertexField(*CellSizes_m);
    }
    void initializeCellPositions()
    {
	CellPositions_m =
	    new Field<Vektor<double,3>,3,Cartesian<3,double>,Cell>(*Mesh_m, *CLayout_m);
	Mesh_m->getCellPositionField(*CellPositions_m);
    }
    void initializeCellVolumes()
    {
	CellVolumes_m = new typename ccsf::BaseField_t(*Mesh_m, *CLayout_m);
	Mesh_m->getCellVolumeField(*CellVolumes_m);
    }
    void initializeCellSurfaceNormals()
    {
	CellSurfaceNormals_m =
	    new Field<Vektor<double,3>,3,Cartesian<3,double>,Cell>*[6];
	for (unsigned int f=0; f<6; ++f) {
	    CellSurfaceNormals_m[f] =
		new Field<Vektor<double,3>,3,Cartesian<3,double>,Cell>(*Mesh_m,
								       *CLayout_m);
	}
	Mesh_m->getSurfaceNormalFields(CellSurfaceNormals_m);
    }
};

template<class Mesh>
template<class FT>
const int PoomaMesh_XYZ<Mesh>::VertexProxy<FT>::FacesAtVertex::vf[8][3] = {
    {0,2,4}, {1,2,4}, {0,3,4}, {1,3,4},
    {0,2,5}, {1,2,5}, {0,3,5}, {1,3,5}
};

#include "PoomaMesh_XYZ.t.hh"

#endif		// __mesh_PoomaMesh_XYZ_hh__

//---------------------------------------------------------------------------//
//		end of mesh/PoomaMesh_XYZ.hh
//---------------------------------------------------------------------------//
