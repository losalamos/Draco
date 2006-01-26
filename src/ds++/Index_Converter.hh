//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/Index_Converter.hh
 * \author Mike Buksas
 * \date   Fri Jan 20 14:51:51 2006
 * \brief  
 * \note   Copyright 2006 The Regents of the University of California.
 *
 * Long description.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef dsxx_Index_Converter_hh
#define dsxx_Index_Converter_hh

#include <vector>
#include <functional>
#include <algorithm>
#include <numeric>

namespace rtt_dsxx
{

//===========================================================================//
/*!
 * \class Index_Converter
 * \brief Utiltity class for converting one dimension indicies to and from
 * N-dimensional ones.
 *
 * 
 * 
 * \sa Index_Converter.cc for detailed descriptions.
 *
 */
/*! 
 * \example ds++/test/tstIndex_Converter.cc 
 * 
 * description of example
 */
//===========================================================================//
template<unsigned D, int OFFSET>
class Index_Converter 
{
  public:

    //! Default constructor
    Index_Converter() {/* ... */}

    //! Construct with just a pointer to the sizes
    Index_Converter(const unsigned* dimensions);

    //! Construct a with all dimensions equal
    Index_Converter(const unsigned dimension);

    //! Destructor.
    ~Index_Converter() {/* ... */}

    //! Assignment operator for Index_Converter.
    Index_Converter& operator=(const Index_Converter &rhs);

    //! Comparison operator
    bool operator==(const Index_Converter &rhs) const;

    bool operator!=(const Index_Converter &rhs) const { return !(*this==rhs);}

    //! Re-assignment operator
    void resize(const unsigned* dimensions);

    // ACCESSORS

    //! Convert N-index to 1-index
    template <typename IT> int operator()(IT indices) const;

    //! Convert 1-index to N-index
    std::vector<int> operator()(int index) const;

    //! Convert 1-index to N-index and store in provided iterator.
    template <typename IT> void operator()(int index, IT begin) const;

    int get_size() const { return array_size; }
    int get_size(int dim) const;

    int get_next_index(int index, int direction) const;


  private:

    // DATA

    unsigned array_size;

    unsigned dimensions[D]; //!< Sizes of each dimension
    unsigned sub_sizes [D]; //!< Sizes of sub-grids of increasing dimension.



    // IMPLEMENTATION
    
    bool sizes_okay() const;

    template <typename IT>
    bool indices_in_range(IT indices) const;
    bool index_in_range(int index) const;

    
};

//---------------------------------------------------------------------------//
// Function Definitions
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/**
 * \brief Construct from a pointer to unsigned values for the dimensions.
 * 
 */
template <unsigned D, int OFFSET>
Index_Converter<D,OFFSET>::Index_Converter(const unsigned* dimensions_)
{
    resize(dimensions_);
}

//---------------------------------------------------------------------------//
/**
 * \brief Resize the index converter object with new dimensions.
 * 
 */
template <unsigned D, int OFFSET>
void Index_Converter<D,OFFSET>::resize(const unsigned* dimensions_)
{
    std::copy(dimensions_, dimensions_+D, dimensions);

    array_size = std::accumulate<unsigned*>(
        dimensions, dimensions+D, 1, std::multiplies<unsigned>());

    sub_sizes[0] = 1;
    
    unsigned* end = std::partial_sum<unsigned*>(
        dimensions, dimensions+D-1, sub_sizes+1, std::multiplies<unsigned>());

    Ensure(end == sub_sizes +D);
    Ensure(sizes_okay());
}

//---------------------------------------------------------------------------//
/**
 * \brief Comparison routine
 *
 * \arg The Index_Converter object to compare to.
 * 
 */
template <unsigned D, int OFFSET>
inline bool Index_Converter<D,OFFSET>::operator==(const Index_Converter &rhs) const
{

    if (array_size != rhs.array_size) return false;

    return std::equal(dimensions, dimensions+D, rhs.dimensions);

}

//---------------------------------------------------------------------------//
/**
 * \brief Convert an N-index to a 1-index
 * 
 */
template <unsigned D, int OFFSET>
template <typename IT>
int Index_Converter<D,OFFSET>::operator()(IT indices) const
{

    Check(indices_in_range(indices));

    int one_index_value = 0;
    int dimension = 0;
    for (IT index = indices; index != indices + D; ++index, ++dimension)
    {
        one_index_value += (*index - OFFSET) * sub_sizes[dimension];
    }
                           
    one_index_value += OFFSET;

    Ensure(index_in_range(one_index_value));

    return one_index_value;
    
}


//---------------------------------------------------------------------------//
/**
 * \brief Convert a 1-index to an N-index
 *
 * \arg index The 1-index value
 *
 * This function dispatches to the write-in-place version of the function and
 * stores the result in a local int[] array. It then constructs the return
 * vector in the return statement in order to allow the compiler to perform
 * return value optimization (RVO). This can potentially eliminate the
 * creation of a temporary return object.
 * 
 */
template <unsigned D, int OFFSET>
std::vector<int> Index_Converter<D,OFFSET>::operator()(int index) const
{

    Check(index_in_range(index));

    static int indices[D];

    operator()(index, indices);

    // Construct in return statement for RVO.
    return std::vector<int>(indices, indices+D);
    
}


//---------------------------------------------------------------------------//
/**
 * \brief Convert a 1-index to an N-index. Store in provided pointer
 *
 * \arg index The index
 * \arg iterator The iterator pointing to the place to store the results.
 * 
 */
template <unsigned D, int OFFSET>
template <typename IT>
void Index_Converter<D,OFFSET>::operator()(int index, IT iter) const
{

    Check(index_in_range(index));

    IT point(iter+D);

    index -= OFFSET;

    for (int dimension = D-1; dimension >= 0; --dimension)
    {
        *(--point) = index / sub_sizes[dimension] + OFFSET;
        index %= sub_sizes[dimension];
    }

    Ensure (point == iter);

}

//---------------------------------------------------------------------------//
/**
 * \brief Return the next index in a given direction. Return -1 if this
 * direction is outside the range of indices
 *
 * \arg index   The index in question
 * \arg direction The direction, 1-based numbered (negative,positive) by
 * dimension.
 *
 */
template <unsigned D, int OFFSET>
int Index_Converter<D,OFFSET>::get_next_index(int index, int direction) const
{

    int indices[D];

    return 0;
}

//---------------------------------------------------------------------------//
// IMPLEMENTATION ROUTINES
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/**
 * \brief Internal routine to make sure all of the dimensions are positive.
 *
 * Since the dimensions are stored as unsigned integers, we need only check
 * for zeros.
 * 
 */
template <unsigned D, int OFFSET>
inline bool Index_Converter<D,OFFSET>::sizes_okay() const
{
    
    return (std::find(dimensions, dimensions+D, 0) == dimensions+D);

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
bool Index_Converter<D,OFFSET>::indices_in_range(IT indices) const
{

    int dimension = 0;
    for (IT index = indices; index != indices + D; ++index, ++dimension)
        if (*index < OFFSET || *index >= dimensions[dimension]+OFFSET) return false;

    return true;
}


//---------------------------------------------------------------------------//
/**
 * \brief Make sure the 1-index is in range
 *
 * \arg index The value of the 1-index
 * 
 */
template <unsigned D, int OFFSET>
inline bool Index_Converter<D,OFFSET>::index_in_range(int index) const
{

    return (index >= OFFSET) && (index < array_size + OFFSET);

}

} // end namespace rtt_dsxx

#endif // dsxx_Index_Converter_hh

//---------------------------------------------------------------------------//
//              end of ds++/Index_Converter.hh
//---------------------------------------------------------------------------//
