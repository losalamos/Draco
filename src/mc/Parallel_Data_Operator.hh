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

#include <iterator>

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
/*!
 * \example mc/test/tstParallel_Data_Op.cc
 *
 * Example usages of the Parallel_Data_Operator class.  In particular, note
 * how in functions test_mapping_DD() and test_mapping_replication() the
 * nested operation classes are used to determine the mapping for various
 * data fields in parallel topologies.  
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
    typedef rtt_dsxx::SP<Topology> SP_Topology;

    // nested operations classes
    struct Data_Replicated;
    struct Data_Decomposed;
    struct Data_Distributed;

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

    // Map local data fields to global data fields.
    template<class FT, class Op>
    void local_to_global(FT &local, FT &global, Op);
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

//---------------------------------------------------------------------------//
// LOCAL-TO-GLOBAL MAPPING OPERATIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Map local data fields to global data fields.

 * This function takes a local field and maps the data to a global field
 * across all processors.  The function is templated on field type (FT) and
 * Op type.  The Op type provides a functor class (Data_Decomposed,
 * Data_Distributed, Data_Replicated) that describes how data should be
 * mapped across processors in certain topologies.  The following list
 * describes the data combinations: \arg \b "replication"/Data_Replicated
 * local_field is copied directly into the global_field, it is assumed that
 * the local_fields are the same on each processor; \arg \b
 * "replication"/Data_Decomposed the local_field is copied into the global
 * field on each processor, the global_fields are subsequently summed across
 * all processors; \arg \b "DD"/Data_Distributed local_fields are mapped into
 * the appropriate cells on the global_field on each processor, the
 * global_fields are then summed across all processors; \arg \b
 * "DD/replication"/Data_Decomposed first the local_field data is mapped into
 * local copies of the global_field on each processor, the global_fields are
 * then summed-up across all processors; \arg \b
 * "DD/replication"/Data_Replicated first the local_field data is mapped into
 * local copies of the global_field on each processor, the global_fields are
 * then mapped across all processors with replicated cells being assigned the
 * proper value and non-replicated cells communicating its results to all
 * other processors.
 *
 * The Parallel_Data_Operator knows the problem topology.  The user must
 * specify an appropriate nested operations functor class when calling this
 * function.  The choices are Parallel_Data_Operator::Data_Replicated,
 * Parallel_Data_Operator::Data_Decomposed, and
 * Parallel_Data_Operator::Data_Distributed.  See the examples for info.
 *
 * \param local_field local field data on a processor 
 * \param global_field mutable global field data sized to the number of
 * global cells
 * \param mapping nested operations functor class (see above) 
 */
template<class FT, class Op>
void Parallel_Data_Operator::local_to_global(FT &local_field,
					     FT &global_field,
					     Op mapping)
{
    // the global field should be equal to the number of global cells and
    // exists on every processor
    Require(global_field.size() == topology->num_cells());

    // if topology is replication we do a simple copy
    if (topology->get_parallel_scheme() == "replication")
    {
	// the local field and global field should be the same size for
	// replication topologies
	Check(local_field.size() == global_field.size());

	// perform the local-global mapping depending upon the data-type
	mapping(local_field, global_field, *this);
    }
    else if (topology->get_parallel_scheme() == "DD")
    {
	Check(local_field.size() == topology->num_cells(C4::node()));

	// iterators for local field
	typename FT::iterator begin  = local_field.begin();
	typename FT::iterator end    = local_field.end();
	
	// iterators for global field
	typename FT::iterator global_begin = global_field.begin();
	typename FT::iterator global_end   = global_field.end();
	
	// iterators for reading/writing
	typename FT::iterator itr;
	typename FT::iterator global;

	// local and global cell indices
	int local_cell;
	int global_cell;

	// sweep through local cells and add them to the appropriate point in 
	// the local processor version of the global field
	for (itr = begin; itr != end; itr++)
	{
	    // calculate the local cell index
	    local_cell = (itr - begin) + 1;
	    Check(local_cell > 0 && local_cell <= local_field.size());

	    // calculate the global cell index
	    global_cell = topology->global_cell(local_cell);
	    
	    // calculate the global cell iterator position corresponding to
	    // the global_cell -> we require a random access iterator type in 
	    // order to perform iterator::difference_type + iterator addition
	    global = (global_cell - 1) + global_begin;
	    Check(global >= global_begin && global < global_end);

	    // assign the local cell value to the global cell value
	    *global = *itr;
	}

	// now do a summation over all processors of the local versions of
	// the global fields
	global_sum(global_begin, global_end);
    }
    else 
    {
	Insist(0, "Haven't done the rest yet!");
    }
}

//---------------------------------------------------------------------------//
// NESTED OPERATIONS CLASSES
//---------------------------------------------------------------------------//
/*!
 * \brief Nested operations functor class for data-replicated data.
 */
struct Parallel_Data_Operator::Data_Replicated
{
    // operator for performing data replicated operations in full replication 
    // topoogies
    template<class FT>
    void operator()(FT &local_field, FT &global_field, 
		    Parallel_Data_Operator data_op)
    {	
	Require (local_field.size() == global_field.size());

	typename FT::iterator begin  = local_field.begin();
	typename FT::iterator end    = local_field.end();
	typename FT::iterator global = global_field.begin();

	// on each processor assign local data into global data
	for (typename FT::iterator itr = begin; itr != end; itr++)
	{
	    *global = *itr;
	    global++;
	}
    }
    
    template<class FT>
    void operator()(FT &local_field, FT &global_field, SP_Topology top, 
		    Parallel_Data_Operator data_op)
    {
	// this will be a combination of data_distributed and data_replicated 
	// operations

	// first do a data distributed (full DD-like) operation locally

	// second cannot do a blind global sum: for cells that only live on
	// one processor, do a gsum(), for cells that are replicated, assign
	// if they are on the processor, communicate to processors that don't 
	// have the cell --> then check on the global_field
    }
};

//---------------------------------------------------------------------------//
/*!
 * \brief Nested operations functor class for data-decomposed data.
 */
struct Parallel_Data_Operator::Data_Decomposed
{
    // operator for performing data decomposed operations in full replication
    // topologies
    template<class FT>
    void operator()(FT &local_field, FT &global_field, 
		    Parallel_Data_Operator data_op)
    {	
	Require (local_field.size() == global_field.size());

	typename FT::iterator local_begin  = local_field.begin();
	typename FT::iterator local_end    = local_field.end();
	
	typename FT::iterator global       = global_field.begin();
	typename FT::iterator global_begin = global_field.begin();
	typename FT::iterator global_end   = global_field.end();

	typename FT::iterator itr;

	// on each processor assign local data into global data
	for (itr = local_begin; itr != local_end; itr++) 
	{
	    *global = *itr;
	    global++;
	}

	// sum up global data
	data_op.global_sum(global_begin, global_end);
    }

    template<class FT>
    void operator()(FT &local_field, FT &global_field, SP_Topology top, 
		    Parallel_Data_Operator data_op)
    {
	// this will be a combination of data_distributed and data_decomposed 
	// operations

	// first do a data distributed (full DD-like) operation locally

	// second do a blind global_sum
    }
};

//---------------------------------------------------------------------------//
/*!
 * \brief Nested operations functor class for data-distributed data.
 */
struct Parallel_Data_Operator::Data_Distributed
{
    template<class FT>
    void operator()(FT &local, FT &global, Parallel_Data_Operator data_op)
    {
	// Do nothing
    }
};

} // end namespace rtt_mc

#endif                          // __mc_Parallel_Data_Operator_hh__

//---------------------------------------------------------------------------//
//                              end of mc/Parallel_Data_Operator.hh
//---------------------------------------------------------------------------//
