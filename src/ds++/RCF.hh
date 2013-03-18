//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/RCF.hh
 * \author Thomas M. Evans, Rob Lowrie
 * \date   Mon Jan 26 15:12:22 2004
 * \brief  Reference Counted Field class definition file.
 * \note   Copyright (C) 2003-2013 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_ds_RCF_hh
#define rtt_ds_RCF_hh

#include "Assert.hh"
#include "SP.hh"

namespace rtt_dsxx
{

//===========================================================================//
/*!
 * \class RCF
 * 
 * \brief Adds reference counting to a field.
 *
 * The RCF provides reference counting to a field type.  In a sense it is an
 * analog of rtt_dsxx::SP, except it is tailored to field types.  The field
 * types are expected to be random access containers, and functions that are
 * part of the Container Concept (with the exception of swap) are also
 * included (Austern, Sec. 9.1.1).  This means that the following operations
 * are provided:
 * - operator[]()
 * - size()
 * - empty()
 * - begin() and end()
 * .
 *
 * The RCF is templated on Field_t (Field Type).  The field must
 * provide the following functions
 * - Field_t::Field_t(int n, value_type v);
 * - Field_t::Field_t(const_iterator b, const_iterator e);
 * - const T& Field_t::operator[] const
 * - T& Field_t::operator[]
 * - size_type size() const
 * - bool empty() const
 * - const_iterator begin() const
 * - iterator begin()
 * - const_iterator end() const
 * - iterator end()
 * .
 * 
 * Finally, the Field_t must supply the following type definitions:
 * - Field_t::size_type
 * - Field_t::value_type
 * - Field_t::iterator
 * - Field_t::const_iterator
 * .
 *
 * An example of usage is:
 * \code
 *     RCF<vector<double> > x(new vector<double>(10, 0.0));
 *     for (int i = 0; i < x.size(); i++)
 *     {
 *          x[i] += 10.0;
 *     }
 * \endcode
 * 
 * or
 *
 * \code
 *     RCF<vector<double> > x(10, 0.0); // shorthand
 *     for (vector<double>::iterator i = x.begin(); 
 *          i != x.end(); i++)
 *     {
 *          *i += 10.0;
 *     }
 * \endcode
 *
 * When a field type pointer is given to the RCF, the RCF takes ownership.
 * The RCF is reference counted; thus, when the last RCF goes out of scope,
 * the field is deleted.  It is dangerous to attempt to access the underlying
 * field owned by a RCF through an address, ie.
 * \code
 *     RCF<vector<double> > x(new vector<double>());
 *     vector<double> *y = &(x.get_field()); // VERY DANGEROUS
 *     delete y;                             // YIKES!!!!!!!!!
 * \endcode
 * 
 * \sa rtt_dsxx::SP for more details on reference counting and Austern,
 * "Generic Programming and the STL" Sec. 9.1.1 for container concepts.
 */
/*!
 * \example ds++/test/tstRCF.cc
 *
 * rtt_dsxx::RCF usage example.
 */
//===========================================================================//

template<typename Field_t>
class RCF
{
  public:
    // Useful typedefs.
    typedef typename Field_t::value_type     value_type;
    typedef typename Field_t::size_type      size_type;
    typedef typename Field_t::iterator       iterator;
    typedef typename Field_t::const_iterator const_iterator;
    typedef Field_t                          Field;

  private: 
    // >>> DATA

    //! Smart pointer to field.
    SP<Field_t> sp_field;
    
    //! Dumb pointer to the field
    /*! This member is used to implement all of the member functions.  There
      tends to be overhead associated with accessing data through an SP: this
      approach avoids that overhead */
    Field_t* ptr_field;

    // Friendship for const specialization.
    template<typename X> friend class RCF;

  public:

    //! Default constructor.
    RCF(void) : sp_field(NULL), ptr_field(NULL) { /* */ }

    //! Constructor.
    inline explicit RCF(int        const n,
			value_type const v = value_type());

    // Explicit constructor for type Field_t *.
    inline explicit RCF(Field_t *p_in);

    // Range constructor.
    inline RCF(const_iterator b, const_iterator e);

    // copy construction and assignament (class has pointer members). 
    RCF( RCF<Field_t> const & );
    RCF<Field_t> & operator=( RCF<Field_t> const & rhs );
    
    // Assignment operator for type Field_t *.
    inline RCF<Field_t>& operator=(Field_t *p_in);

    //! Get the field (const).
    const Field_t& get_field() const 
    { 
	Require(assigned()); 
	return *ptr_field; 
    }

    //! Get the field (l-value).
    Field_t& get_field() 
    { 
	Require(assigned()); 
	return *ptr_field; 
    }

    //! Determine if field is assigned.
    bool assigned() const { return bool(sp_field); }

    // Expose operator[] on underlying Field_t (const).
    inline const value_type& operator[](const size_type) const;

    // Expose operator[] on underlying Field_t.
    inline value_type& operator[](const size_type);

    // Expose begin() (const).
    inline const_iterator begin() const;

    // Expose begin().
    inline iterator begin();

    // Expose end() (const).
    inline const_iterator end() const;

    // Expose end().
    inline iterator end();

    // Expose size().
    size_type size() const { Require(assigned()); return ptr_field->size(); }

    // Expose empty().
    bool empty() const { Require(assigned()); return ptr_field->empty(); }

};

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Construct an RCF from its size and initial value.
 *
 * This constructor has the following usage:
 * \code
 *     // make an RCF field from a vector<double>(10, 0.0)
 *     RCF<vector<double> > x(10, 0.0);
 * \endcode
 */
template<typename Field_t>
RCF<Field_t>::RCF(const int        n,
		  const value_type v)
    : sp_field( new Field_t(n,v) ),
      ptr_field( sp_field.bp() )
{
    // empty
}

//---------------------------------------------------------------------------//
/*!
 * \brief Construct a RCF from a pointer to the field.
 *
 * Once a pointer is given to the RCF it \b owns the field.  It would be very
 * dangerous (although legal) to try to access the pointer to the native
 * field.  The field is deleted when the last copy of RCF goes out of scope.
 *
 * This constructor has the following usage:
 * \code
 *     // make a RCF field from a vector<double>
 *     RCF<vector<double> > x(new vector<double>(10, 0.0));
 * \endcode
 */
template<typename Field_t>
RCF<Field_t>::RCF(Field_t *p_in)
    : sp_field(p_in),
      ptr_field(sp_field.bp())
{
    // nothing to check because this could be a NULL field pointer
}

//---------------------------------------------------------------------------//
/*!
 * \brief Range constructor.
 *
 * This constructor creates a new field, of size \a e - \a b,
 * and initializes its values using the iterator range [\a b, \a e).
 *
 * \param b Starting iterator.
 * \param e Ending iterator.
 */
template<typename Field_t>
RCF<Field_t>::RCF(const_iterator b,
		  const_iterator e)
    : sp_field( new Field_t(b,e) ),
      ptr_field( sp_field.bp() )
{
    /* empty */
}

//---------------------------------------------------------------------------//
/*!
 * \brief Assignment operator for type Field_t *.
 *
 * As in copy construction, the RCF owns the pointer after assignment. Here
 * is an example of usage:
 * \code
 *     RCF<Mat> f;         // has reference to NULL field
 *     f      = new Mat;   // now has 1 count of Mat
 *     Mat *g = new Mat; 
 *     f      = g;        // f's original pointer to Mat is deleted
 *                        // because count goes to zero; f now has
 *                        // 1 reference to g
 * \endcode
 */
template<typename Field_t>
RCF<Field_t> & RCF<Field_t>::operator=(Field_t *p_in)
{
    // check if we already own this field
    if (sp_field.bp() == p_in)
        return *this;

    // reassign the existing smart pointer
    sp_field = p_in;
    ptr_field = sp_field.bp();
    return *this;
}

//! Assignment operator
template<typename Field_t>
RCF<Field_t> & RCF<Field_t>::operator=( RCF<Field_t> const & rhs )
{
    // check if we already own this field
    if (this == &rhs)
        return *this;
    
    // reassign the existing smart pointer
    sp_field  = rhs.sp_field;
    ptr_field = sp_field.bp();
    return *this;
}
//! Copy constructor
template<typename Field_t>
RCF<Field_t>::RCF( RCF<Field_t> const & rhs)
    : sp_field( rhs.sp_field ),
      ptr_field( sp_field.bp() )
{
    /* empty */
}

//---------------------------------------------------------------------------//
/*!
 * \brief Expose operator[] on underlying Field_t (const).
 */
template<typename Field_t>
const typename RCF<Field_t>::value_type& 
RCF<Field_t>::operator[](const size_type i) const
{
    Require (assigned());
    return ptr_field->operator[](i);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Expose operator[] on underlying Field_t.
 */
template<typename Field_t>
typename RCF<Field_t>::value_type& RCF<Field_t>::operator[](const size_type i)
{
    Require (assigned());
    return ptr_field->operator[](i);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Expose begin() (const).
 */
template<typename Field_t>
typename RCF<Field_t>::const_iterator RCF<Field_t>::begin() const
{
    Require (assigned());
    return ptr_field->begin();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Expose begin().
 */
template<typename Field_t>
typename RCF<Field_t>::iterator RCF<Field_t>::begin()
{
    Require (assigned());
    return ptr_field->begin();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Expose end() (const).
 */
template<typename Field_t>
typename RCF<Field_t>::const_iterator RCF<Field_t>::end() const
{
    Require (assigned());
    return ptr_field->end();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Expose end().
 */
template<typename Field_t>
typename RCF<Field_t>::iterator RCF<Field_t>::end()
{
    Require (assigned());
    return ptr_field->end();
}

//===========================================================================//
/*
 * \class RCF<const Field_t>
 *
 * \brief Specialization of RCF on const Field_t.
 *
 * This class specializes RCF (reference counted field) for const Field_t.
 * You can copy a RCF<T> to a RCF<const T> but not the other way.
 *
 * \sa rtt_dsxx::RCF for details.
 */
//===========================================================================//

template<typename Field_t>
class RCF<const Field_t>
{
  public:
    // Useful typedefs.
    typedef typename Field_t::value_type     value_type;
    typedef typename Field_t::size_type      size_type;
    typedef typename Field_t::const_iterator const_iterator;
    typedef Field_t                          Field;

  private: 
    // >>> DATA

    //! Smart pointer to field.
    SP<const Field_t> sp_field;

    //! Dumb pointer to the field
    /*! This member is used to implement all of the member functions.  There
      tends to be overhead associated with accessing data through an SP: this
      approach avoids that overhead */
    Field_t const * ptr_field;

  public:
    //! Default constructor.
    RCF(void) : sp_field(NULL), ptr_field(NULL) { /* empty */ }
    
    //! Constructor.
    inline explicit RCF(const int n,
			const value_type v = value_type());

    //! Constructor from non-const Field_t.
    RCF(const RCF<Field_t> &x) 
	: sp_field(x.sp_field),
          ptr_field(sp_field.bp())
    {/*...*/}

    // Explicit constructor for type Field_t *.
    inline explicit RCF(Field_t *p_in);

    // Range constructor.
    inline RCF(const_iterator b, const_iterator e);

    // Assignment operator for type Field_t *.
    inline RCF<const Field_t>& operator=(Field_t *p_in);

    // Provide copy construction and assignemnt operations
    RCF( RCF<const Field_t> const & rhs);
    RCF<const Field_t> & operator=( RCF<const Field_t> const & rhs);
    
    // Const member functions apparently must be defined within this class
    // definition for the IBM compiler to work.

    //! Get the field (const).
    const Field_t& get_field() const
    {
        Require(assigned()); 
	return *ptr_field;
    }

    //! Determine if field is assigned.
    bool assigned() const { return bool(sp_field); }

    // Expose operator[] on underlying Field_t (const).
    const value_type& operator[](const size_type i) const
    {
       Require (assigned());
       return ptr_field->operator[](i);
    }

    // Expose begin() (const).
    const_iterator begin() const
    {
        Require (assigned());
        return ptr_field->begin();
    }

    // Expose end() (const).
    const_iterator end() const
    {
      Require (assigned());
      return ptr_field->end();
    }

    // Expose size().
    size_type size() const { Require(assigned()); return ptr_field->size(); }

    // Expose empty().
    bool empty() const { Require(assigned()); return ptr_field->empty(); }
};

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Construct an RCF from its size and initial value.
 *
 * This constructor has the following usage:
 * \code
 *     // make an RCF field from a vector<double>(10, 0.0)
 *     RCF<const vector<double> > x(10, 0.0);
 * \endcode
 */
template<typename Field_t>
RCF<const Field_t>::RCF(const int        n,
			const value_type v)
    : sp_field( new Field_t(n,v) ),
      ptr_field( sp_field.bp() )
{
    /* empty */
}

//---------------------------------------------------------------------------//
/*!
 * \brief Construct a RCF to a const Field_t from a pointer to the field.
 *
 * Once a pointer is given to the RCF it \b owns the field.  It would be very
 * dangerous (although legal) to try to access the pointer to the native
 * field.  The field is deleted when the last copy of RCF goes out of scope.
 *
 * This constructor has the following usage:
 * \code
 *     // make a RCF field from a vector<double>
 *     RCF<const vector<double> > x(new vector<double>(10, 0.0));
 * \endcode
 */
template<typename Field_t>
RCF<const Field_t>::RCF(Field_t *p_in)
    : sp_field(p_in)
    , ptr_field(sp_field.bp())
{
    // nothing to check because this could be a NULL field pointer
}

//---------------------------------------------------------------------------//
/*!
 * \brief Range constructor.
 *
 * This constructor creates a new field, of size \a e - \a b,
 * and initializes its values using the iterator range [\a b, \a e).
 *
 * \param b Starting iterator.
 * \param e Ending iterator.
 */
template<typename Field_t>
RCF<const Field_t>::RCF(const_iterator b,
			const_iterator e)
    : sp_field( new Field_t(b,e) ),
      ptr_field( sp_field.bp() )
{
   /* empty */
}

//---------------------------------------------------------------------------//
/*!
 * \brief Assignment operator for type Field_t * to RCF<const Field_t>.
 *
 * As in copy construction, the RCF owns the pointer after assignment. Here
 * is an example of usage:
 * \code
 *     RCF<const Mat> f;         // has reference to NULL field
 *     f      = new Mat;         // now has 1 count of Mat
 *     Mat *g = new Mat; 
 *     f      = g;               // f's original pointer to Mat is deleted
 *                               // because count goes to zero; f now has
 *                               // 1 reference to g
 * \endcode
 */
template<typename Field_t>
RCF<const Field_t>& RCF<const Field_t>::operator=(Field_t *p_in)
{
    // check if we already own this field
    if (sp_field.bp() == p_in)
        return *this;

    // reassign the existing smart pointer
    sp_field = p_in;
    ptr_field = sp_field.bp();
    return *this;
}

//! Assignment operator
template<typename Field_t>
RCF<const Field_t> & RCF<const Field_t>::operator=( RCF<const Field_t> const & rhs )
{
    // check if we already own this field
    if (this == &rhs)
        return *this;
    
    // reassign the existing smart pointer
    sp_field  = rhs.sp_field;
    ptr_field = sp_field.bp();
    return *this;
}
//! Copy constructor
template<typename Field_t>
RCF<const Field_t>::RCF( RCF<const Field_t> const & rhs)
    : sp_field( rhs.sp_field ),
      ptr_field( sp_field.bp() )
{
    /* empty */
}

//---------------------------------------------------------------------------//
/*!
 * \brief Expose operator[] on underlying Field_t (const).
 */
//
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// NOTE: Moving the following function definition outside the class
// definition, as below, does not compile with the IBM compiler on Purple.
// Apparently it is the const qualifier, because the other definitions above
// compile just fine.  
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// template<typename Field_t>
// const typename RCF<const Field_t>::value_type& 
// RCF<const Field_t>::operator[](const size_type i) const
// {
//     Require (assigned());
//     return ptr_field->operator[](i);
// }

} // end namespace rtt_dsxx

#endif // rtt_ds_RCF_hh

//---------------------------------------------------------------------------//
// end of ds++/RCF.hh
//---------------------------------------------------------------------------//
