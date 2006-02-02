//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/Index_Iterator.hh
 * \author Mike Buksas
 * \date   Tue Jan 31 16:45:39 2006
 * \brief  
 * \note   Copyright 2006 The Regents of the University of California.
 *
 * Long description.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef dsxx_Index_Iterator_hh
#define dsxx_Index_Iterator_hh

#include "Assert.hh"
#include "Index_Converter.hh"

namespace rtt_dsxx
{

//===========================================================================//
/*!
 * \class Index_Iterator
 * \brief Facilitates iterating over a multi-dimensional range of indices.
 *
 *
 * \sa Index_Iterator.cc for detailed descriptions.
 *
 */
/*! 
 * \example ds++/test/tstIndex_Iterator.cc 
 * 
 */
//===========================================================================//
template <unsigned D, int OFFSET>
class Index_Iterator 
{
  public:

    // NESTED CLASSES AND TYPEDEFS

    // CREATORS
    
    //! Default constructors.
    Index_Iterator(const Index_Converter<D,OFFSET>& index_converter);

    //! Destructor.
    ~Index_Iterator() { /* ... */ }

    // MANIPULATORS
    
    //! Assignment operator for Index_Iterator.
    Index_Iterator& operator=(const Index_Iterator &rhs);

    // ACCESSORS

    Index_Iterator& operator++() { increment(); return *this; }
    Index_Iterator& operator--() { decrement(); return *this; }

    // Accessors for the 1-index
    int get_index() const { return index; }
    operator int()  const { return index; }

    // Accessors for the N-indices
    int get_index(unsigned d) const { Check(dimension_okay(d)); return indices[d]; }

    std::vector<int> get_indices() const { return std::vector<int>(indices,indices+D);}

    template <typename IT>
    void get_indices(IT out) const { std::copy(indices, indices+D, out); }

    bool is_in_range() const { return in_range; }

  private:

    // DATA

    const Index_Converter<D,OFFSET>& index_converter;

    int indices[D];
    int index;
    bool in_range;
    

    // IMPLEMENTATION

    void increment();
    void decrement();

    bool dimension_okay(int d) const { return (d >= 0) && (d < D); }


};

//---------------------------------------------------------------------------//
/**
 * \brief Construct from an Index_Converter object
 * 
 */
template <unsigned D, int OFFSET>
Index_Iterator<D,OFFSET>::Index_Iterator(const Index_Converter<D,OFFSET>& converter)
    : index_converter(converter),
      in_range(true),
      index(OFFSET)
{
    for (int d=0; d < D; ++d) indices[d] = OFFSET;
}

//---------------------------------------------------------------------------//
// IMPLEMENTATION
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/**
 * \brief Increment the iterator
 * 
 */
template <unsigned D, int OFFSET>
void Index_Iterator<D,OFFSET>::increment()
{

    Require(in_range);

    ++indices[0];
    ++index;

    for (int d = 0; d < D-1; ++d)
    {
        if (indices[d] > index_converter.max_of_index(d))
        {
            ++indices[d+1];
            indices[d] = index_converter.min_of_index(d);
        }
        else
            break;
    }

    if (indices[D-1] > index_converter.max_of_index(D-1))
    {
        indices[D-1] = index_converter.min_of_index(D-1);
        in_range = false;
    }

    Ensure (!in_range || (index_converter.get_index(indices) == index));
    
}

//---------------------------------------------------------------------------//
/**
 * \brief Decrement the iterator
 * 
 */
template <unsigned D, int OFFSET>
void Index_Iterator<D,OFFSET>::decrement()
{

    Require(in_range);

    --indices[0];
    --index;

    for (int d = 0; d < D-1; ++d)
    {
        if (indices[d] < index_converter.min_of_index(d))
        {
            indices[d] = index_converter.max_of_index(d);
            --indices[d+1];
        }
        else
            break;
    
    }

    if (indices[D-1] < index_converter.min_of_index(D-1))
    {
        indices[D-1] = index_converter.max_of_index(D-1);
        in_range = false;
    }

    Ensure ( !in_range || (index_converter.get_index(indices) == index));

}


} // end namespace rtt_dsxx

#endif // dsxx_Index_Iterator_hh

//---------------------------------------------------------------------------//
//              end of ds++/Index_Iterator.hh
//---------------------------------------------------------------------------//
