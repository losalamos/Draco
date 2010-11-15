//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/value.hh
 * \author Kent Budge
 * \brief  Definition of function template value
 * \note   Copyright 2007 LANSLLC All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef dsxx_value_hh
#define dsxx_value_hh

#include "Field_Traits.hh"

namespace rtt_dsxx
{

//---------------------------------------------------------------------------//
/*! Strip a field type of any labeling.
 *
 * Implicit conversion of a labeled field type (such as a
 * value-plus-derivatives class) to an underlying unlabeled field type can
 * be dangerous. However, unlike conversion constructors, conversion operators
 * cannot be flagged explicit. Our alternative is to templatize a conversion
 * function from labeled field type to unlabeled field type.
 *
 * The default implementation assumes there is no labeling to strip.
 *
 * \arg \a Field A field type
 */
template<class Field>
inline
typename Field_Traits<Field>::unlabeled_type &value(Field &x)
{
    return x;
}

//---------------------------------------------------------------------------//
//! A version of the value function template for const arguments.
template<class Field>
inline
typename Field_Traits<Field>::unlabeled_type const &value(Field const &x)
{
    return x;
}

} // end namespace rtt_dsxx

#endif // dsxx_value_hh

//---------------------------------------------------------------------------//
//              end of ds++/value.hh
//---------------------------------------------------------------------------//
