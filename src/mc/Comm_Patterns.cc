//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Comm_Patterns.cc
 * \author Thomas M. Evans
 * \date   Mon May  1 19:50:44 2000
 * \brief  Comm_Patterns implementation file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Comm_Patterns.hh"
#include "c4/global.hh"
#include "ds++/Assert.hh"
#include <iostream>

namespace rtt_mc
{

// Draco namespace includes
using C4::C4_Req;
using C4::RecvAsync;
using C4::Send;

// STD namespace includes
using std::pair;
using std::vector;

//---------------------------------------------------------------------------//
// Build The Communication Patterns
//---------------------------------------------------------------------------//
/*

 * \brief Build communication patterns between processor's boundary cells.

 * This function uses the Topology class to determine what processors a
 * processor interfacts with across its processor boundary. It calculates two
 * information maps:

 * \arg map of processor/boundary cell data that the processor sends out to
 * other processors.

 * \arg map of processor/boundary cell data that the processor receives from
 * other processors.

 * The former data map is a list of boundary cells that live on other
 * processors.  For the current processor to perform work, it needs to know
 * where it can get data on the boundary cells.  This field tells it what
 * data will come from what processors.  The later field tells the current
 * processor what data it owns that it needs to send out to other processors
 * so that they may do their boundary work.  The former field requires
 * communication to calculate.  The later field can be calculated without
 * communication.

 * We support three basic topology types: replication, DD, DD/replication.
 * If the topology is replication then no data is set because there are no
 * boundary cells.  If the topology is DD/replication then an assertion is
 * thrown because we do not support that topology fully yet.  Data maps for
 * full DD topologies are calculated.

 */
void Comm_Patterns::calc_patterns(SP_Topology topology)
{
    Require (topology);

    // empty the query maps
    empty_query_maps();

    // check the topology and launch the correct calculator
    if (topology->get_parallel_scheme() == "replication")
    {
	// there are no boundary cells so we don't do anything here
	Ensure (send_query_map.empty());
	Ensure (send_query_map.empty());

	return;
    }
    else if (topology->get_parallel_scheme() == "DD")
    {
	// we can use a simple function to calculate the boundary cell maps
	// because we have exact one-one communication
	calc_recv_map_DD(topology);

	// do the requisite communication to calculate the send maps after
	// the recv_query_map has been calculated
	calc_send_map_DD(topology);

	Ensure (recv_query_map.size() == send_query_map.size());;
    }
    else if (topology->get_parallel_scheme() == "DD/replication")
    {
	// we don't support this yet!
	Insist(0, "We don't support general DD schemes yet!");
    }
    else
    {
	Insist(0, "Don't know your topology!");
    }
}

//---------------------------------------------------------------------------//
// Comm_Patterns Implementation
//---------------------------------------------------------------------------//
/*!
 * \brief Calculate recv_query_map patterns for DD topologies.
 */
void Comm_Patterns::calc_recv_map_DD(SP_Topology topology)
{
    Require (topology->get_parallel_scheme() == "DD");
    Require (recv_query_map.empty());

    // get number of boundary cells
    int num_bcells = topology->get_boundary_cells(C4::node());
    
    // insertion notification for recv_query_map
    pair<iterator, bool> insert_indicator;

    // entries into recv_query_map
    proc_pair entry;

    // global cell and processor indices
    int global_cell;
    int processor;

    // loop through boundary cells for this processor and determine the
    // processor they live on and their global cell index.  Insert them into
    // the recv_query_map.  These are the cells that this processor will
    // require information about to do calculations
    for (int bc = 1; bc <= num_bcells; bc++)
    {
	// determine global cell index of the boundary cell
	global_cell = topology->boundary_to_global(bc, C4::node());

	// find out what processor its on; for full DD each cell is only on
	// one processor so the vector list of cells that is returned from
	// topology will only have one element in it
	processor = topology->get_procs(global_cell).front();

	// insert the processor/cell pair into the map
	entry.first      = processor;
	entry.second     = sf_int(1, global_cell);
	insert_indicator = recv_query_map.insert(entry);

	// if the data wasn't inserted because this processor key already
	// exists in the map then add teh cell to the processor cell data
	// listing
	if (!insert_indicator.second)
	    recv_query_map[processor].push_back(global_cell);
    }

    // the recv_query_map should not be empty because there must be at least
    // one boundary cell
    Ensure (!recv_query_map.empty());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Calcualte send_query_map patterns for DD topologies.
 */
void Comm_Patterns::calc_send_map_DD(SP_Topology topology)
{
    Require (topology->get_parallel_scheme() == "DD");
    Require (send_query_map.empty());
    Require (!recv_query_map.empty());

    // number of boundary processors we interact with, it is the number of
    // processors we receive messages from because full DD is one-to-one
    // communication
    int num_bproc = recv_query_map.size();

    // C4 Requests for the number of cells that other processors will query
    // us about
    vector<C4_Req> size_recv(num_bproc);
    vector<C4_Req> cells_recv(num_bproc);

    // data that we receive from the processors we communicate with telling
    // us the number of cells that are on this processor that they need
    // information about and what (global)cells those are
    vector<int>   number_of_cells(num_bproc);
    vector<int *> list_of_cells(num_bproc);

    // iterators to beginning and end of recv_query_map
    const_iterator recv_begin = recv_query_map.begin();
    const_iterator recv_end   = recv_query_map.end();

    // post receives from other processors for the number of cells that they
    // require information about and send out to those processors the number
    // of cells that we need information from them about
    int index = 0;
    int processor;
    int query_size;
    for (const_iterator i = recv_begin; i != recv_end; i++)
    {
	// determine the processor that we are communicating with and the
	// number of cells that reside on that processor
	processor  = i->first;
	query_size = i->second.size();

	// post receives from that processor
	RecvAsync(size_recv[index], &number_of_cells[index], 1, processor,
		  600);

	index++;

	// send the number of cells that the processor needs to send us
	// information about
	Send(query_size, processor, 600);
    }
    Check (index == recv_query_map.size());

    // receive the size data from sent from other processors->we only receive 
    // one message from each of those processors
    int finished = 0;
    while (finished < size_recv.size())
    {
	// loop through size C4 requests and see if the message has completed
	for (int i = 0; i < size_recv.size(); i++)
	{
	    // check to see if the message has been completed
	    if (size_recv[i].complete())
	    {
		finished++;
		Check (!size_recv[i].inuse());
		Check (!size_recv[i].complete());
	    }
	}
    }
    Check (finished == size_recv.size());

    // post receives for the list of cells that we will receive from each
    // processor; send out the list of cell data that we require from each
    // processor
    index = 0;
    for (const_iterator i = recv_begin; i != recv_end; i++)
    {
	// determine the processor that we are communicating with and the
	// number of cells that reside on that processor
	processor  = i->first;
	query_size = i->second.size();
	

	// allocate list of cell data arrays
	list_of_cells[index] = new int[number_of_cells[index]];

	// post receives for the cell lists that we will receive from other
	// processors
	RecvAsync(cells_recv[index], list_of_cells[index],
		  number_of_cells[index], processor, 601);

	index++;

	// send the list of cells that this processor needs information about 
	// to the other processor
	int *query_data = new int[query_size];
	for (int cell = 0; cell < query_size; cell++)
	    query_data[cell] = i->second[cell];
	Send<int>(query_data, query_size, processor, 601);
	delete [] query_data;
    }
    Check (index == recv_query_map.size());

    // data we will need to make the send_query_map
    const_iterator       iter;
    proc_pair            entry;   
    pair<iterator, bool> insert_indicator;

    // receive the cell lists from the processors that we communicate with
    // and build the send_query_map
    finished = 0;
    while (finished < cells_recv.size())
    {
	// assign the iterator to the recv_query_map to the beginning
	// location
	iter = recv_begin;

	// loop through cell list C4 requests and see if the message is
	// complete
	for (int i = 0; i < cells_recv.size(); i++)
	{
	    // determine the processor that we are communicating with and
	    // make it the key entry for the send_query_map; because this is
	    // a DD topology the send_query_map and recv_query_map have the
	    // same prcocessor keys
	    entry.first = iter->first;
	    
	    // see if the message is complete
	    if (cells_recv[i].complete())
	    {
		// indicate that the message is complete
		finished++;
		Check (!cells_recv[i].inuse());
		Check (!cells_recv[i].complete());

		// write the list of cells to the entry
		entry.second.resize(number_of_cells[i]);
		for (int cell = 0; cell < entry.second.size(); cell++)
		    entry.second[cell] = list_of_cells[i][cell];
		delete [] list_of_cells[i];

		// add the entry to the send_query_map
		insert_indicator = send_query_map.insert(entry);
		Check (insert_indicator.second);
	    }

	    // advance processor map
	    iter++;
	}
	Check (iter == recv_end);
    }
    Check (finished == cells_recv.size());

    // the send_query_map should not be empty because there must be at least
    // one boundary cell
    Ensure (!send_query_map.empty());

    // get all processors to the same place before proceeding
    C4::gsync();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Empty the query maps.
 */
void Comm_Patterns::empty_query_maps()
{
    // empty send_query_map
    if (!send_query_map.empty())
	send_query_map.erase(send_query_map.begin(), send_query_map.end());
    
    // empty recv_query_map
    if (!recv_query_map.empty())
	recv_query_map.erase(recv_query_map.begin(), recv_query_map.end());

    Ensure (send_query_map.empty());
    Ensure (recv_query_map.empty());
}

//---------------------------------------------------------------------------//
// Comm_Patterns Interface
//---------------------------------------------------------------------------//
/*!
  
 * \brief Determine if the Comm_Patterns are set.

 * The overloaded boolean operator determines if the Comm_Patterns are set or
 * not.  It has no way of determining if the data is correct.

 */
Comm_Patterns::operator bool() const
{
    return !send_query_map.empty() && !recv_query_map.empty();
}

} // end namespace rtt_mc

//---------------------------------------------------------------------------//
//                              end of Comm_Patterns.cc
//---------------------------------------------------------------------------//
