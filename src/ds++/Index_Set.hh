//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/Index_Set.hh
 * \author Mike Buksas
 * \date   Thu Feb  2 10:01:46 2006
 * \brief  
 * \note   © Copyright 2006 LANSLLC All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef dsxx_Index_Set_hh
#define dsxx_Index_Set_hh

#include <functional>
#include <numeric>
#include <algorithm>

#include "Assert.hh"

namespace rtt_dsxx
{

//===========================================================================//
/*!
 * \class Index_Set
 * \brief Represents a D-dimensional set if indices.
 *
 * \sa Index_Set.cc for detailed descriptions.
 *
 */
/*! 
 * \example ds++/test/tstIndex_Set.cc 
 * 
 * description of example
 */
//===========================================================================//
template <unsigned D, int OFFSET>
class Index_Set 
{
  public:

    // NESTED CLASSES AND TYPEDEFS

    // CREATORS
    
    //! Default constructors.
    Index_Set() { /* ... */ }

    //! Construct with pointer to sizes
    Index_Set(unsigned const * const dimensions) { set_size(dimensions); }

    //! Construct with all dimensions equal
    Index_Set(const unsigned dimension) { set_size(dimension); }

    //! Destructor
    virtual ~Index_Set() { /* ... */ }
    
    //! Comparison operator
    bool operator==(const Index_Set &rhs) const;

    //! Negative comparison operator
    bool operator!=(const Index_Set &rhs) const { return !(*this==rhs); }

    //! Re-assignment operator
    void set_size(unsigned const* const dimensions);

    //! Uniform size re-assignment operator
    void set_size(const unsigned size);

    bool index_in_range(int index) const {
        return (index >= OFFSET) && (index < array_size + OFFSET);
    }
    bool index_in_range(int index, unsigned dimension) const;

    template <typename IT> bool indices_in_range(IT indices) const;

    int get_size()     const { return array_size; }
    int min_of_index() const { return OFFSET; }
    int max_of_index() const { return OFFSET + array_size - 1; }
    int limit_of_index(const bool positive) const {
        return positive ? max_of_index() : min_of_index();
    }

    int get_size(const int d) const { Check(dimension_okay(d)); return dimensions[d]; }
    int min_of_index(const unsigned d) const { Check(dimension_okay(d)); return OFFSET; }
    int max_of_index(const unsigned d) const {
        Check(dimension_okay(d)); return OFFSET+dimensions[d]-1;
    }
     int limit_of_index(const unsigned d, const bool positive) const {
         return positive ? max_of_index(d) : min_of_index(d);
    }
            
    static bool direction_okay(const int d) { return (d >  0) && (d <= 2*D); }
    static bool dimension_okay(const int d) { return (d >= 0) && (d <  D);   }

  private:

    void compute_size();
    
    unsigned array_size;    //!< Sizes of the whole index range
    unsigned dimensions[D]; //!< Sizes of each dimension

  protected:

    // Make sure the index sizes are all positive when creating or resizing:
    bool sizes_okay() const {
        return (std::find(dimensions, dimensions+D, 0) == dimensions+D);
    }

    // Allow derived classes const access to the dimensions.
    unsigned const * get_dimensions() const { return dimensions; }

};


//---------------------------------------------------------------------------//
/**
 * \brief Set the size of the Index_Set. Discards old size information
 *
 * \arg sizes Pointer to unsigned integers for the index set sizes.
 * 
 */
template <unsigned D, int OFFSET>
void Index_Set<D,OFFSET>::set_size(unsigned const* const dimensions_)
{

    std::copy(dimensions_, dimensions_+D, dimensions);

    Require(sizes_okay());

    compute_size();

}

//---------------------------------------------------------------------------//
/**
 * \brief Set the size of the Index_Set to a uniform dimension. Discards old
 * size information
 *
 * \arg dimension The uniform dimension of the index set.
 * 
 */
template <unsigned D, int OFFSET>
void Index_Set<D,OFFSET>::set_size(const unsigned dimension)
{

    for (unsigned* dim = dimensions; dim < dimensions + D; ++dim)
        *dim = dimension;

    compute_size();
}

//---------------------------------------------------------------------------//
/**
 * \brief Comparison routine
 *
 * \arg The Index_Set object to compare to.
 * 
 */
template <unsigned D, int OFFSET>
inline bool Index_Set<D,OFFSET>::operator==(const Index_Set &rhs) const
{

    if (array_size != rhs.array_size) return false;

    return std::equal(dimensions, dimensions+D, rhs.dimensions);

}


//---------------------------------------------------------------------------//
/**
 * \brief Make sure the indices are with the range for each dimension
 *
 * \arg iterator An itertator to a range of indices.
 * 
 */
template <unsigned D, int OFFSET>
template <typename IT>
bool Index_Set<D,OFFSET>::indices_in_range(IT indices) const
{

    int dimension = 0;
    for (IT index = indices; index != indices + D; ++index, ++dimension)
        if (!index_in_range(*index, dimension)) return false;

    return true;
}


//---------------------------------------------------------------------------//
/**
 * \brief Return true iff the given index is within the range for the given
 * dimension
 *
 * \arg index The index value
 * \arg dimension The dimension of the index
 * 
 */
template <unsigned D, int OFFSET>
inline bool Index_Set<D,OFFSET>::index_in_range(int index, unsigned dimension) const
{
    Check(dimension_okay(dimension));

    return ((index >= OFFSET) && (index < dimensions[dimension] + OFFSET));
}



//---------------------------------------------------------------------------//
// IMPLEMENTAION
//---------------------------------------------------------------------------//
template <unsigned D, int OFFSET>
inline void Index_Set<D,OFFSET>::compute_size()
{

    array_size =
        std::accumulate<unsigned*>(dimensions, dimensions+D, 1, std::multiplies<unsigned>());

    Ensure(array_size > 0);

}


} // end namespace rtt_dsxx

#endif // dsxx_Index_Set_hh

//---------------------------------------------------------------------------//
//              end of ds++/Index_Set.hh
//---------------------------------------------------------------------------//
