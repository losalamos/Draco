//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/RCF.hh
 * \author Thomas M. Evans, Rob Lowrie
 * \date   Mon Jan 26 15:12:22 2004
 * \brief  Reference Counted Field class definition file.
 * \note   Copyright © 2003 The Regents of the University of California.
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
 * part of the Basic Container Concept (with the exception of swap) are also
 * included.  This means that the following operations are provided:
 * - operator[]()
 * - size()
 * - empty()
 * - begin() and end()
 * .
 *
 * The RCF is templated on Field_t (Field Type).  The field must
 * provide the following functions
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
 *     RCF<vector<double> > x(new vector<double>(10, 0.0));
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
 * \sa rtt_dsxx::SP for more details on reference counting.
 */
/*!
 * \example ds++/test/tstRCF.cc
 *
 * rtt_dsxx::RCF usage example.
 */
// revision history:
// -----------------
// 0) (Mon Jan 26 15:12:22 2004) Thomas M. Evans: original
// 
//===========================================================================//

template<class Field_t>
class RCF
{
  public:
    // Useful typedefs.
    typedef typename Field_t::value_type     value_type;
    typedef typename Field_t::size_type      size_type;
    typedef typename Field_t::iterator       iterator;
    typedef typename Field_t::const_iterator const_iterator;

  private: 
    // >>> DATA

    //! Smart pointer to field.
    SP<Field_t> sp_field;

  public:
    //! Default constructor.
    RCF() {/*...*/}

    // Explicit constructor for type Field_t *.
    inline explicit RCF(Field_t *p_in);

    // Assignment operator for type Field_t *.
    inline RCF<Field_t>& operator=(Field_t *p_in);

    //! Get the field (const).
    const Field_t& get_field() const { Require(assigned()); return *sp_field; }

    //! Get the field (l-value).
    Field_t& get_field() { Require(assigned()); return *sp_field; }

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
    size_type size() const { Require(assigned()); return sp_field->size(); }

    // Expose empty().
    bool empty() const { Require(assigned()); return sp_field->empty(); }
};

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
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
template<class Field_t>
RCF<Field_t>::RCF(Field_t *p_in)
    : sp_field(p_in)
{
    // nothing to check because this could be a NULL field pointer
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
template<class Field_t>
RCF<Field_t>& RCF<Field_t>::operator=(Field_t *p_in)
{
    // check if we already own this field
    if (sp_field.bp() == p_in)
        return *this;

    // reassign the existing smart pointer
    sp_field = p_in;
    return *this;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Expose operator[] on underlying Field_t (const).
 */
template<class Field_t>
const typename RCF<Field_t>::value_type& 
RCF<Field_t>::operator[](const size_type i) const
{
    Require (assigned());
    return sp_field->operator[](i);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Expose operator[] on underlying Field_t.
 */
template<class Field_t>
typename RCF<Field_t>::value_type& RCF<Field_t>::operator[](const size_type i)
{
    Require (assigned());
    return sp_field->operator[](i);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Expose begin() (const).
 */
template<class Field_t>
typename RCF<Field_t>::const_iterator RCF<Field_t>::begin() const
{
    Require (assigned());
    return sp_field->begin();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Expose begin().
 */
template<class Field_t>
typename RCF<Field_t>::iterator RCF<Field_t>::begin()
{
    Require (assigned());
    return sp_field->begin();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Expose end() (const).
 */
template<class Field_t>
typename RCF<Field_t>::const_iterator RCF<Field_t>::end() const
{
    Require (assigned());
    return sp_field->end();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Expose end().
 */
template<class Field_t>
typename RCF<Field_t>::iterator RCF<Field_t>::end()
{
    Require (assigned());
    return sp_field->end();
}

} // end namespace rtt_dsxx

#endif // rtt_ds_RCF_hh

//---------------------------------------------------------------------------//
//              end of ds++/RCF.hh
//---------------------------------------------------------------------------//
