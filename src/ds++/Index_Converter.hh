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
 * Code Sample:
 * \code
 *     cout << "Hello, world." << endl;
 * \endcode
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

    //! Construct with just a pointer to the sizes
    Index_Converter(unsigned* dimensions);

    //! Construct a with all dimensions equal
    Index_Converter(unsigned dimension);

    //! Destructor.
    ~Index_Converter() {/* ... */}

    //! Assignment operator for Index_Converter.
    Index_Converter& operator=(const Index_Converter &rhs);

    //! Comparison operator
    bool operator==(const Index_Converter &rhs) const;

    bool operator!=(const Index_Converter &rhs) const { return !(*this==rhs);}


    // ACCESSORS

    //! Convert N-index to 1-index
    template <typename IT>
    int operator()(IT indices) const;

    //! Convert 1-index to N-index
    std::vector<int> operator()(int index) const;

    int get_size() const { return array_size; }
    int get_size(int dim) const;


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
Index_Converter<D,OFFSET>::Index_Converter(unsigned* dimensions_)
{

    std::copy(dimensions_, dimensions_+D, dimensions);

    array_size = std::accumulate<unsigned*>(
        dimensions, dimensions+D, 1, std::multiplies<unsigned>());

    sub_sizes[0] = 1;
    
    std::partial_sum<unsigned*>(
        dimensions, dimensions+D-1, sub_sizes+1, std::multiplies<unsigned>());


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

    return std::equal(dimensions, dimensions+D, rhs.dimension);

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
 */
template <unsigned D, int OFFSET>
std::vector<int> Index_Converter<D,OFFSET>::operator()(int index) const
{

    Check(index_in_range(index));

    std::vector<int> indices(D);

    index -= OFFSET;

    for (int dimension = D; dimension > 0; --dimension)
    {

        indices[dimension-1] = index / sub_sizes[dimension-1] + OFFSET;

        index = index % sub_sizes[dimension-1];

    }

    return indices;
    
}


//---------------------------------------------------------------------------//
// IMPLEMENTATION ROUTINES
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/**
 * \brief Internal routine to make sure all of the dimensions are positive.
 * 
 */
template <unsigned D, int OFFSET>
inline bool Index_Converter<D,OFFSET>::sizes_okay() const
{
    
    return (std::find_if(dimensions, dimensions+D,
                         std::bind2nd(std::less_equal<int>(), 0)) == dimensions+D);

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
bool Index_Converter<D,OFFSET>::index_in_range(int index) const
{

    return (index >= OFFSET) && (index < array_size + OFFSET);

}

} // end namespace rtt_dsxx

#endif // dsxx_Index_Converter_hh

//---------------------------------------------------------------------------//
//              end of ds++/Index_Converter.hh
//---------------------------------------------------------------------------//
