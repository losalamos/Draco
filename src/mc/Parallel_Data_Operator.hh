//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Parallel_Data_Operator.hh
 * \author Thomas M. Evans
 * \date   Fri Dec 10 09:58:21 1999
 * \brief  Parallel_Data_Operator header file.
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_mc_Parallel_Data_Operator_hh
#define rtt_mc_Parallel_Data_Operator_hh

#include "Topology.hh"
#include "Comm_Patterns.hh"
#include "c4/global.hh"
#include "ds++/Assert.hh"
#include "ds++/SP.hh"

#include <iterator>
#include <vector>

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
// 1) 05-03-00 : modified local_to_global mapping to allow different local
//               and global field types (e.g. vector and ccsf)
// 
//===========================================================================//

class Parallel_Data_Operator 
{
  public:
    // some typedefs used in this class
    typedef rtt_dsxx::SP<Topology>      SP_Topology;
    typedef rtt_dsxx::SP<Comm_Patterns> SP_Comm_Patterns;

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
    template<class FT, class GT, class Op>
    void local_to_global(FT &local, GT &global, Op);

    // Create boundary cell data fields in a spatially decomposed topology.
    template<class T> 
    void gather_bnd_cell_data(SP_Comm_Patterns, const std::vector<T> &, 
			      std::vector<T> &);
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
 * describes the data combinations:

 * \arg \b "replication"/Data_Replicated local_field is copied directly into
 * the global_field, it is assumed that the local_fields are the same on each
 * processor;

 * \arg \b "replication"/Data_Decomposed the local_field is copied into the
 * global field on each processor, the global_fields are subsequently summed
 * across all processors;

 * \arg \b "DD"/Data_Distributed local_fields are mapped into the appropriate
 * cells on the global_field on each processor, the global_fields are then
 * summed across all processors;

 * \arg \b "DD/replication"/Data_Decomposed first the local_field data is
 * mapped into local copies of the global_field on each processor, the
 * global_fields are then summed-up across all processors;

 * \arg \b "DD/replication"/Data_Replicated first the local_field data is
 * mapped into local copies of the global_field on each processor, the
 * global_fields are then mapped across all processors with replicated cells
 * being assigned the proper value and non-replicated cells communicating its
 * results to all other processors.
 
 * The Parallel_Data_Operator knows the problem topology.  The user must
 * specify an appropriate nested operations functor class when calling this
 * function.  The choices are Parallel_Data_Operator::Data_Replicated,
 * Parallel_Data_Operator::Data_Decomposed, and
 * Parallel_Data_Operator::Data_Distributed.  See the examples for info.
 
 * \param local_field local field data on a processor

 * \param global_field mutable global field data sized to the number of
 * global cells

 * \param mapping nested operations functor class (see above)

 */
template<class FT, class GT, class Op>
void Parallel_Data_Operator::local_to_global(FT &local_field,
					     GT &global_field,
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
	Check(global_field.size() == topology->num_cells());

	// iterators for local field
	const typename FT::iterator begin  = local_field.begin();
	const typename FT::iterator end    = local_field.end();
	
	// iterators for global field
	const typename GT::iterator global_begin = global_field.begin();
	const typename GT::iterator global_end   = global_field.end();
	
	// iterators for reading/writing
	typename FT::iterator itr;
	typename GT::iterator global;

	// local and global cell indices
	int local_cell;
	int global_cell;

	// sweep through local cells and add them to the appropriate point in 
	// the local processor version of the global field
	for (itr = begin; itr != end; itr++)
	{
	    // calculate the local cell index
	    local_cell = std::distance(begin, itr) + 1;
	    Check(local_cell > 0 && local_cell <= local_field.size());

	    // calculate the global cell index
	    global_cell = topology->global_cell(local_cell);
	    
	    // calculate the global cell iterator position corresponding to
	    // the global_cell index -> the advance function is much faster
	    // when global is a random access iterator
	    global = global_begin;
	    std::advance(global, (global_cell-1));
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

    template<class FT, class GT>
    void operator()(FT &local_field, GT &global_field, SP_Topology top, 
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
    template<class FT, class GT>
    void operator()(FT &local, GT &global, Parallel_Data_Operator data_op)
    {
	// Do nothing
    }
};

//---------------------------------------------------------------------------//
// Boundary Cell Parallel Operations (gather/scatter type)
//---------------------------------------------------------------------------//
/*!

 * \brief Fill boundary cell data fields in a spatially decomposed topology.

 * This function utilizes the rtt_mc::Comm_Patterns class to "gather" data
 * into boundary cell fields.  What this entails is the following:

 * \arg using the Comm_Patterns, each processor posts requests to processors
 * it needs data from;

 * \arg using the Comm_Patterns and local data fields, each processor sends
 * data out to processors that need data;

 * \arg each processor receives the data for its boundary cells and writes it
 * into a local boundary cell field

 * The Topology must be spatially decomposed (not full replication).

 * \param pattern active Comm_Patterns object, the patterns must be set

 * \param local_data vector of local data on processor that will be sent to
 * other processors according to the Comm_Patterns

 * \param bc_data empty boundary cell field that will be filled

 */
template<class T> void
Parallel_Data_Operator::gather_bnd_cell_data(SP_Comm_Patterns pattern,
					     const std::vector<T> &local_data, 
					     std::vector<T> &bc_data)
{
    using C4::C4_Req;
    using std::vector;

    Require (bc_data.empty());
    Require (local_data.size() == topology->num_cells(C4::node()));
    Require (*pattern);

    // first make a vector of data that will be received from each processor
    vector<C4_Req>   data_requests(pattern->get_num_recv_procs());
    vector<double *> recv_data(pattern->get_num_recv_procs());

    // iterators to the recv and send data maps
    Comm_Patterns::const_iterator itor;
    Comm_Patterns::const_iterator send_begin = pattern->get_send_begin();
    Comm_Patterns::const_iterator send_end   = pattern->get_send_end();
    Comm_Patterns::const_iterator recv_begin = pattern->get_recv_begin();
    Comm_Patterns::const_iterator recv_end   = pattern->get_recv_end();

    // post the receives to the data
    int processor;
    int size;
    int index = 0;
    for (itor = recv_begin; itor != recv_end; itor++)	
    {
	// determine the processor and the size of the data field
	processor = itor->first;
	size      = itor->second.size();

	// allocate the dataerature arrays
	recv_data[index] = new double[size];

	// post the receive to the data
	C4::RecvAsync(data_requests[index], recv_data[index], size,
		      processor, 605);
	index++;
    }
    Check (itor == recv_end);
    Check (index == pattern->get_num_recv_procs());

    // send out the data
    int local_cell;
    for (itor = send_begin; itor != send_end; itor++)
    {
	// determine the processor we are sending to
	processor = itor->first;
	size      = itor->second.size();

	// make a data field to send
	double *send_data = new double[size];

	// fill it up with dataerature data
	for (int i = 0; i < size; i++)
	{
	    // get the local cell index for a global cell that is needed by
	    // processor and lives on this processor
	    local_cell = topology->local_cell(itor->second[i]);
	    Check (local_cell);
	    
	    // write the data into the temporary sending field
	    send_data[i] = local_data[local_cell-1];
	}

	// send out data and reclaim memory
	C4::Send<double>(send_data, size, processor, 605);
	delete [] send_data;
    }
    Check (itor == send_end);

    // receive the data and write out bc_data field

    // size bc_data data
    int num_bcells = topology->get_boundary_cells(C4::node());
    bc_data.resize(num_bcells);

    // receive data from the processors we communicate with and write their
    // data into the bc_data field
    int global_cell;
    int bound_cell;
    int finished = 0;
    while (finished < pattern->get_num_recv_procs())
    {
	// define an iterator to the recv processor map
	itor = recv_begin;

	// loop through the C4 requests for each processor that we
	// communicate with and see if the data is here
	for (int req = 0; req < data_requests.size(); req++)
	{
	    // determine the processor and size of the message we expect to
	    // receive
	    processor = itor->first;
	    size      = itor->second.size();

	    // see if the message is compete
	    if (data_requests[req].complete())
	    {
		// indicate that this data has been received
		finished++;
		Check (!data_requests[req].inuse());

		// write the data to the bc_data field
		for (int c = 0; c < size; c++)
		{
		    // determine the global_cell
		    global_cell = itor->second[c];

		    // determine the boundary cell
		    bound_cell = topology->global_to_boundary(global_cell,
							      C4::node());
		    Check (bound_cell);

		    // add the data to the boundary data field
		    bc_data[bound_cell-1] = recv_data[req][c];
		}
		
		// reclaim memory
		delete [] recv_data[req];
	    }
	
	    // advance the iterator
	    itor++;
	}
	Check (itor == recv_end);
    }
    Check (finished == pattern->get_num_recv_procs());

    // sync everything
    C4::gsync();
}

} // end namespace rtt_mc

#endif                          // rtt_mc_Parallel_Data_Operator_hh

//---------------------------------------------------------------------------//
//                              end of mc/Parallel_Data_Operator.hh
//---------------------------------------------------------------------------//
