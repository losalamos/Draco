//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Parallel_Data_Operator.hh
 * \author Thomas M. Evans
 * \date   Fri Dec 10 09:58:21 1999
 * \brief  Parallel_Data_Operator header file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __mc_Parallel_Data_Operator_hh__
#define __mc_Parallel_Data_Operator_hh__

#include "Topology.hh"
#include "c4/global.hh"
#include "ds++/Assert.hh"
#include "ds++/SP.hh"

namespace rtt_mc
{
 
//===========================================================================//
/*!
 * \class Parallel_Data_Operator
 *
 * \brief Class to operate on field-type data across processors.
 *
 * Parallel_Data_Operator takes field-type arguments and performs operations
 * on the data across processor space.  The fields must have random access
 * iterators.  Thus, vectors, c-style arrays, deques, etc. will work.
 * Additional functionality is required for containers and fields that have
 * iterators with less functionality than random access.
 *
 * The Parallel_Data_Operator uses the rtt_mc::Topology class to determine
 * mapping between processors.
 *
 * While most functions in this class are not inlined, they are contained in
 * the header file (template functions have internal linkage).  Thus,
 * explicit instantiations on different field types are not necessary.
 *
 * This class also contains two equivalence-type functions to check data
 * across processor space.  These are intended to be used for DBC checking.
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class Parallel_Data_Operator 
{
  public:
    // some typedefs used in this class
    typedef dsxx::SP<Topology> SP_Topology;

  private:
    //! Smart pointer to a topology class.
    SP_Topology topology;

  public:
    // Constructor.
    Parallel_Data_Operator(SP_Topology);

    // Communication functions for DBC checks.
    bool check_global_equiv(int) const;
    bool check_global_equiv(double, double = 1.0e-8) const;

    // Do a global sum of a global mesh-sized field.
    template<class IT> void global_sum(IT begin, IT end);
    template<class T>  void global_sum(T *, T *);
};

//---------------------------------------------------------------------------//
// PARALLEL_DATA_OPERATOR GLOBAL SUMMING FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Sum up global mesh-sized fields on all processors.
 *
 * This function takes field type iterators and sums up the data on all
 * processors.  The summed data is redistributed to all processors.  The size
 * of the field must be equal to the number of global mesh cells.  The
 * function is templated on Field Type (FT).
 *
 * This function works on any topology with the following provisions: \arg
 * each processor gives a globally sized field (equal to the number of global
 * mesh cells), \arg only straight summations are performed (all data in like
 * cells is added).
 *
 * This function assumes that the iterator is not a pointer; thus, a copy to
 * and from a c-style array is performed.  This operation may incur some
 * performance penalty; however, communication requires c-style arrays.
 *
 * \param begin random access iterator to the first element in the field
 * \param end random access iterator to the last element in the field.  
 */
template<class IT>
void Parallel_Data_Operator::global_sum(IT begin, IT end)
{
    // get the value_type from the iterator class
    typedef typename std::iterator_traits<IT>::value_type T;

    // define iterator
    IT itr;

    // field size
    int size = end - begin;
    Require (size == topology->num_cells());
    
    // copy data to a c-style array
    T *data = new T[size];
    int count;

    count = 0;
    for (itr = begin; itr != end; itr++)
	data[count++] = *itr;

    // do a gsum
    C4::gsum(data, size);

    // write data back to field
    count = 0;
    for (itr = begin; itr != end; itr++)
	*itr = data[count++];

    Check (count == size);

    // reclaim memory
    delete [] data;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Sum up global mesh-sized fields on all processors.
 *
 * This function is an optimized version of global_sum.  If the field
 * iterators are pointers, then copying data to c-style arrays is not
 * necessary and is bypassed.  For example, some implementations of
 * vector::iterator are simply T * pointers.  Thus, containers that have
 * iterators that are pointers can run much more efficiently.  Otherwise,
 * this function has the same functionality as template<class IT>
 * Parallel_Data_Operator::global_sum(IT, IT).
 *
 * T* pointers are models of the random access iterator concept.
 *
 * \param begin T* pointer to beginning of field
 * \param end T* pointer to end of field
 */
template<class T>
void Parallel_Data_Operator::global_sum(T *begin, T *end)
{
    // field size
    int size = end - begin;
    Require (size == topology->num_cells());

    // do a gsum
    C4::gsum(begin, size);
}

} // end namespace rtt_mc

#endif                          // __mc_Parallel_Data_Operator_hh__

//---------------------------------------------------------------------------//
//                              end of mc/Parallel_Data_Operator.hh
//---------------------------------------------------------------------------//
