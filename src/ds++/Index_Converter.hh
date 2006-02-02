//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/Index_Converter.hh
 * \author Mike Buksas
 * \date   Fri Jan 20 14:51:51 2006
 * \brief  Decleration and Definition of Index_Converter
 * \note   Copyright 2006 The Regents of the University of California.
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
 * \sa Index_Converter.cc for detailed descriptions.
 *
 */
/*! 
 * \example ds++/test/tstIndex_Converter.cc 
 * 
 */
//===========================================================================//
template<unsigned D, int OFFSET>
class Index_Converter 
{
  public:

    //! Default constructor
    Index_Converter() { /* ... */ }

    //! Construct with just a pointer to the sizes
    Index_Converter(const unsigned* dimensions);

    //! Construct a with all dimensions equal
    Index_Converter(const unsigned dimension);

    //! Destructor.
    virtual ~Index_Converter() {/* ... */}

    //! Assignment operator for Index_Converter.
    Index_Converter& operator=(const Index_Converter &rhs);

    //! Comparison operator
    bool operator==(const Index_Converter &rhs) const;

    bool operator!=(const Index_Converter &rhs) const { return !(*this==rhs);}

    //! Re-assignment operator
    void resize(const unsigned* dimensions);

    //! Uniform size re-assignment operator
    void resize(unsigned size);

    // ACCESSORS

    //! Convert N-index to 1-index
    template <typename IT> int get_index(IT indices) const;

    //! Convert 1-index to N-index
    std::vector<int> get_indices(int index) const;

    //! Convert 1-index to N-index and store in provided iterator.
    template <typename IT> void get_indices(int index, IT begin) const;

    int get_size() const { return array_size; }
    int get_size(int d) const { Check(dimension_okay(d)); return dimensions[d]; }

    int get_next_index(int index, int direction) const;

    bool index_in_range(int index) const;
    bool index_in_range(int index, unsigned dimension) const;
    int min_of_index() const { return OFFSET; }
    int max_of_index() const { return OFFSET + array_size - 1; }
    int min_of_index(unsigned d) const { Check(dimension_okay(d)); return OFFSET; }
    int max_of_index(unsigned d) const {
        Check(dimension_okay(d)); return OFFSET+dimensions[d]-1;
    }
            
    
    bool direction_okay(int face) const;

  private:

    // DATA

    unsigned array_size;

    unsigned dimensions[D]; //!< Sizes of each dimension
    unsigned sub_sizes [D]; //!< Sizes of sub-grids of increasing dimension.



    // IMPLEMENTATION
    
    template <typename IT> bool indices_in_range(IT indices) const;
    bool sizes_okay() const;
    bool dimension_okay(int d) const { return (d >= 0) && (d < D); }
    void compute_sizes();

    
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
Index_Converter<D,OFFSET>::Index_Converter(const unsigned* dimensions)
{
    resize(dimensions);
}


//---------------------------------------------------------------------------//
/**
 * \brief Construct from a single size
 * \arg dimension  The size of all dimensions of the Index_Converter
 * 
 */
template <unsigned D, int OFFSET>
Index_Converter<D,OFFSET>::Index_Converter(unsigned dimension)
{
    resize(dimension);
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

    compute_sizes();

}

//---------------------------------------------------------------------------//
/**
 * \brief Resize the index converter with a uniform size
 *
 * \arg dimension The new size
 */
template <unsigned D, int OFFSET>
void Index_Converter<D,OFFSET>::resize(unsigned dimension)
{
    for (unsigned* it = dimensions; it != dimensions+D; ++it) *it = dimension;

    compute_sizes();
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
int Index_Converter<D,OFFSET>::get_index(IT indices) const
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
std::vector<int> Index_Converter<D,OFFSET>::get_indices(int index) const
{

    Check(index_in_range(index));

    static int indices[D];

    get_indices(index, indices);

    // Construct in return statement for RVO.
    return std::vector<int>(indices, indices+D);
    
}


//---------------------------------------------------------------------------//
/**
 * \brief Convert a 1-index to an N-index. Store in provided pointer
 *
 * \arg index The index
 * \arg iterator The iterator pointing to the place to store the results. Must
 * be a reversible iterator, e.g. implement '--'
 * 
 */
template <unsigned D, int OFFSET>
template <typename IT>
void Index_Converter<D,OFFSET>::get_indices(int index, IT iter) const
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

    Check(index_in_range(index));
    Check(direction_okay(direction));

    --direction;

    unsigned direction_axis = direction / 2;
    int      direction_sign = 2*(direction % 2) - 1;

    int indices[D];

    get_indices(index, indices);

    indices[direction_axis] += direction_sign;

    if (indices[direction_axis] < OFFSET ||
        indices[direction_axis] >= dimensions[direction_axis]+OFFSET)
        return -1;

    return get_index(indices);

}

//---------------------------------------------------------------------------//
// IMPLEMENTATION ROUTINES
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/**
 * \brief Assign the internal data members.
 *
 * Used once the dimensions array has been filled.
 * 
 */
template <unsigned D, int OFFSET>
void Index_Converter<D,OFFSET>::compute_sizes()
{

    Require(sizes_okay());

    array_size = std::accumulate<unsigned*>(
        dimensions, dimensions+D, 1, std::multiplies<unsigned>()
        );

    sub_sizes[0] = 1;
    
    unsigned* end = std::partial_sum<unsigned*>(
        dimensions, dimensions+D-1, sub_sizes+1, std::multiplies<unsigned>()
        );

    Ensure(end == sub_sizes+D);

}


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


//---------------------------------------------------------------------------//
/**
 * \brief Make sure the direction index is valid
 *
 * \arg direction The direcition index.
 * 
 */
template <unsigned D, int OFFSET>
inline bool Index_Converter<D,OFFSET>::direction_okay(int direction) const
{
    return (direction >= 1) && (direction <= 2*D);
}

} // end namespace rtt_dsxx

#endif // dsxx_Index_Converter_hh

//---------------------------------------------------------------------------//
//              end of ds++/Index_Converter.hh
//---------------------------------------------------------------------------//
