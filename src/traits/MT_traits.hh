//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   traits/MT_traits.hh
 * \author Randy M. Roberts
 * \date   Wed May  6 14:50:14 1998
 * \brief  
 * \note   Copyright (C) 1998-2014 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __traits_MT_traits_hh__
#define __traits_MT_traits_hh__

#include <iosfwd>

namespace rtt_traits
{
 //===========================================================================//
 // class MT_traits - 
 //
 // Date created :
 // Purpose      :
 //
 // revision history:
 // -----------------
 // 0) original
 // 
 //===========================================================================//

 template<class VECTOR>
 class vector_traits
 {

     // NESTED CLASSES AND TYPEDEFS

   public:

     typedef typename VECTOR::value_type value_type;
    

     // DATA
    
     // CREATORS
    
     // MANIPULATORS
    
     // ACCESSORS

     // STATIC UTILITIES

   public:
     
     inline static value_type dot(const VECTOR &v1, const VECTOR &v2)
     {
	 return VECTOR::dot(v1, v2);
     }
     
   private:
    
     // IMPLEMENTATION
 };

template <class MT>
struct MT_Traits
{
    static std::ostream &print(std::ostream &os, const MT &mesh)
    {
        os << "No Mesh Output Defined." << std::endl;
        return os;
    }
};

} // end namespace rtt_traits

#endif // __traits_MT_traits_hh__

//---------------------------------------------------------------------------//
// end of traits/MT_traits.hh
//---------------------------------------------------------------------------//
