//----------------------------------*-C++-*----------------------------------//
// CAR_CU_Parallel_Builder.t.hh
// Thomas M. Evans
// Tue Apr 14 14:50:22 1998
//---------------------------------------------------------------------------//
// @> Parallel_Builder implementation file for CAR_CU_Mesh
//---------------------------------------------------------------------------//

#include "Parallel_Builder.hh"
#include "XYCoord_sys.hh"
#include "XYZCoord_sys.hh"
#include "imc/Global.hh"
#include <string>
#include <iomanip>
#include <cstdio>
#include <fstream>

namespace rtt_imc 
{

// draco necessities
using C4::Send;
using C4::Recv;
using C4::node;
using C4::nodes;
using rtt_mc::global::max;
using rtt_mc::global::min;
using rtt_mc::Coord_sys;
using rtt_mc::Layout;
using rtt_mc::XYCoord_sys;
using rtt_mc::XYZCoord_sys;

// std necessities
using std::string;
using std::fill;
using std::remove;
using std::endl;
using std::cout;
using std::setw;
using std::ios;

//---------------------------------------------------------------------------//
// constructors
//---------------------------------------------------------------------------//
// host node constructor that determines the toplogy on the IMC-nodes

template<class MT>
Parallel_Builder<MT>::Parallel_Builder(const MT &mesh, 
				       const Source_Init<MT> &sinit)
{
    Require (!node());

    // calculate the parameters for splitting the problem amongst many
    // processors 
    parallel_topology(mesh, sinit);
    Check (cells_per_proc.size() != 0);
    Check (procs_per_cell.size() != 0);

    // if we have more than 1 processor, send out the global_cells index
    for (int np = 1; np < nodes(); np++)
    {
	// assign global cell lists for each processor
	int num_cells = cells_per_proc[np].size();
	int *send_global = new int[num_cells];
	for (int lcell = 0; lcell < num_cells; lcell++)
	    send_global[lcell] = cells_per_proc[np][lcell];

	// send out the global cell lists
	Send (num_cells, np, 48);
	Send <int>(send_global, num_cells, np, 49);

	// reclaim storage
	delete [] send_global;
    }

    // assign the global cell list on the host processor
    global_cells = cells_per_proc[0];
}

//---------------------------------------------------------------------------//
// default constructor for IMC-nodes

template<class MT>
Parallel_Builder<MT>::Parallel_Builder()
{
    Require (node());

    // receive the global_cell lists from the host node
    int num_cells;
    Recv (num_cells, 0, 48);
    int *recv_global = new int[num_cells];
    Recv (recv_global, num_cells, 0, 49);

    // assign the global cell list to global_cells
    global_cells.resize(num_cells);
    for (int i = 0; i < num_cells; i++)
	global_cells[i] = recv_global[i];

    // reclaim storage
    delete [] recv_global;
}

//---------------------------------------------------------------------------//
// IMC Topology calculation
//---------------------------------------------------------------------------//
// calculate topology parameters to determine what cells go where

template<class MT>
void Parallel_Builder<MT>::parallel_topology(const MT &mesh,
					     const Source_Init<MT> &sinit)
{
    // make sure we have a Source_Init and Mesh
    Require (&sinit);
    Require (&mesh);

    // number of cells in the global (master) mesh
    int num_cells = mesh.num_cells();

    // size the partitioning vectors
    cells_per_proc.resize(nodes());
    procs_per_cell.resize(num_cells);

    // calculate the total capacity of all processors
    int total_capacity = sinit.get_capacity() * nodes();
    Check (total_capacity >= num_cells);

    // do full replication / full DD / rep-DD
    if (num_cells <= sinit.get_capacity())
    {
	// do full replication as we can fit the whole mesh on each processor
	parallel_scheme = "replication";

	for (int proc = 0; proc < nodes(); proc++)
	{
	    for (int cell = 1; cell <= num_cells; cell++)
	    {
		procs_per_cell[cell-1].push_back(proc);
		cells_per_proc[proc].push_back(cell);
		Ensure (procs_per_cell[cell-1].size() <= nodes());
	    }
	    Ensure (cells_per_proc[proc].size() == num_cells);
	}
    }
    else if (num_cells == total_capacity)
    {
	// do full DD
	parallel_scheme = "DD";

	int cell = 1;
	for (int proc = 0; proc < nodes(); proc++)
	{
	    int last_cell = 1;
	    while (last_cell <= sinit.get_capacity())
	    {
		cells_per_proc[proc].push_back(cell);
		procs_per_cell[cell-1].push_back(proc);
		Ensure (procs_per_cell[cell-1].size() == 1);
		last_cell++;
		cell++;
	    }
	    Ensure (procs_per_cell[proc].size() <= sinit.get_capacity());
	}
    }
    else
    {
	// do DD/replication
	parallel_scheme = "DD/replication";

	// set up variables for cell replication
	vector<int> cellrep(num_cells);
	fill(cellrep.begin(), cellrep.end(), 0);
	double sfrac = 0;
	int total_rep = 0;

	// loop for first estimate of cell replication
	for (int cell = 1; cell <= num_cells; cell++)
	{
	    // calculate the source fraction
	    sfrac = static_cast<double>(sinit.get_ncen(cell) +
					sinit.get_nvol(cell) +
					sinit.get_nss(cell)) / 
		(sinit.get_ncentot() + sinit.get_nsstot() + 
		 sinit.get_nvoltot());

	    // estimate the number of cell replicates where the number must be
	    // >=1 and <=number of processors
	    cellrep[cell-1] = min(nodes(), max(1, static_cast<int>
					       (sfrac * total_capacity)));
	    total_rep += cellrep[cell-1];
	}

	// calculate amount of extra capacity still available
	int xtra_reps = total_capacity - total_rep;

	// loop until xtra_capacity is gone
	while (xtra_reps != 0)
	{
	    // reset counters
	    total_rep = 0;

	    if (xtra_reps > 0)
	    {
		// loop over cells and add one replicate at a time
		int cell = 1;
		while (cell <= num_cells && total_rep < xtra_reps)
		{
		    if (cellrep[cell-1] < nodes())
			cellrep[cell-1]++;
		    total_rep++;
		    cell++;
		}
	    }
	    else if (xtra_reps < 0)
	    {
		// loop over cells and subtract one replicate at a time
		int cell = 1;
		while (cell <= num_cells && total_rep > xtra_reps)
		{
		    if (cellrep[cell-1] > 1)
			cellrep[cell-1]--;
		    total_rep--;
		    cell++;
		}
	    }
	    xtra_reps -= total_rep;
	}

	// now we can do our DD/replication partitioning
	for (int cell = 1; cell <= num_cells; cell++)
	{	
	    // marker to determine if we have put a cell on this node
	    vector<bool> flag(nodes(), false);

	    // fill up the modes with cell
	    int replicates = 1;
	    int old_replicates = 0;
	    do
	    {
		for (int proc = 0; proc < nodes(); proc++)
		{
		    if (replicates <= cellrep[cell-1] && 
			cells_per_proc[proc].size() < sinit.get_capacity()
			&& !flag[proc])
		    {
			procs_per_cell[cell-1].push_back(proc);
			cells_per_proc[proc].push_back(cell);
			flag[proc] = true;
			replicates++;
		    }
		    Check (cells_per_proc[proc].size() <=
			   sinit.get_capacity());
		}

		// test on old replicates to prevent infinite loops
		if (old_replicates == replicates)
		    break;
		old_replicates = replicates;
	    } while (replicates <= cellrep[cell-1]);
 	    Ensure (procs_per_cell[cell-1].size() > 0);
 	    Ensure (procs_per_cell[cell-1].size() <= nodes());
	}
    }
}

//---------------------------------------------------------------------------//
// Source Distribution interface
//---------------------------------------------------------------------------//
// send the source to the IMC-topology processors

template<class MT>
template<class PT> SP<Source<MT> >
Parallel_Builder<MT>::send_Source(SP<MT> mesh, SP<Mat_State<MT> > mat,
				  SP<Rnd_Control> rcon,
				  const Source_Init<MT> &sinit, 
				  const Particle_Buffer<PT> &buffer)
{
    // check that we are on the host node only
    Require (!node());
    Require (mesh);
    Require (mat);
    Require (mesh->num_cells() == mat->num_cells());

    // data necessary to build Source on host
    SP<Source<MT> > host_source;
    typename MT::CCSF_int volrn(mesh);
    typename MT::CCSF_int nvol(mesh);
    typename MT::CCSF_double ew_vol(mesh);
    typename MT::CCVF_double t4_slope(mesh);
    typename MT::CCSF_int ssrn(mesh);
    typename MT::CCSF_int nss(mesh);
    typename MT::CCSF_int fss(mesh);
    typename MT::CCSF_double ew_ss(mesh);
    Particle_Buffer<PT>::Census census;

    // send the numbers of each type of source to the other processors
    for (int i = 1; i < nodes(); i++)
    {
	Send (sinit.get_ncentot(), i, 34);
	Send (sinit.get_nvoltot(), i, 35);
	Send (sinit.get_nsstot(), i, 36);
    }

    // first distribute the census
    if (sinit.get_ncentot() > 0) 
	dist_census(sinit, buffer, census);

    // next do the volume source
    if (sinit.get_nvoltot() > 0)
	dist_vol(sinit, volrn, nvol, ew_vol, t4_slope);

    // finally do the surface source
    if (sinit.get_nsstot() > 0)
	dist_ss(sinit, ssrn, nss, fss, ew_ss);

    // get the number of volume, census, and surface sources for this processor
    int nvoltot = 0;
    int nsstot = 0;
    for (int i = 1; i <= mesh->num_cells(); i++)
    {
	nvoltot += nvol(i);
	nsstot  += nss(i);
    }

    // do the source distributions

    // surface source distribution
    string ss_dist = sinit.get_ss_dist();
    for (int i = 1; i < nodes(); i++)
    {
	if (ss_dist == "none")
	    Send (0, i, 41);
	else if (ss_dist == "normal")
	    Send (-1, i, 41);
	else if (ss_dist == "cosine")
	    Send (1, i, 41);
	else
	    Check (0);
    }

    // make the source
    host_source = new Source<MT>(volrn, nvol, ew_vol, t4_slope, ssrn, nss, 
				 fss, ew_ss, census, ss_dist, nvoltot,
				 nsstot, rcon, buffer, mat); 

    // return source
    return host_source;
}

//---------------------------------------------------------------------------//
// receive the source from the host processor

template<class MT>
template<class PT> SP<Source<MT> > 
Parallel_Builder<MT>::recv_Source(SP<MT> mesh, SP<Mat_State<MT> > mat, 
				  SP<Rnd_Control> rcon,
				  const Particle_Buffer<PT> &buffer)
{
    // require that we are on an IMC node only
    Require (node());
    Require (mesh);
    Require (mat);

    // declare the source
    SP<Source<MT> > imc_source;

    // receive the total numbers of particles for each type of source
    int global_ncentot;
    int global_nvoltot;
    int global_nsstot;
    Recv (global_ncentot, 0, 34);
    Recv (global_nvoltot, 0, 35);
    Recv (global_nsstot, 0, 36);

    // define the variables for Source construction
    typename MT::CCSF_int volrn(mesh);
    typename MT::CCSF_int nvol(mesh);
    typename MT::CCSF_double ew_vol(mesh);
    typename MT::CCVF_double t4_slope(mesh);
    typename MT::CCSF_int ssrn(mesh);
    typename MT::CCSF_int nss(mesh);
    typename MT::CCSF_int fss(mesh);
    typename MT::CCSF_double ew_ss(mesh);
    Particle_Buffer<PT>::Census census;

    // receive the census
    if (global_ncentot)
	recv_census(buffer, census);

    // receive the volume source
    if (global_nvoltot)
	recv_vol(volrn, nvol, ew_vol, t4_slope);

    // receive the surface source
    if (global_nsstot)
	recv_ss(ssrn, nss, fss, ew_ss);

    // get the number of volume, census, and surface sources for this processor
    int nvoltot = 0;
    int nsstot = 0;
    for (int i = 1; i <= mesh->num_cells(); i++)
    {
	nvoltot += nvol(i);
	nsstot  += nss(i);
    }

    // get and calculate source distributions
    
    // surface source distributions
    string ss_dist;
    int ss_int;
    Recv (ss_int, 0, 41);
    if (ss_int == 0)
	ss_dist = "none";
    else if (ss_int == -1)
	ss_dist = "normal";
    else if (ss_int == 1)
	ss_dist = "cosine";
    else
	Check (0);

    // make the source
    imc_source = new Source<MT>(volrn, nvol, ew_vol, t4_slope, ssrn, nss, 
				fss, ew_ss, census, ss_dist, nvoltot, nsstot, 
				rcon, buffer, mat);

    // return source
    return imc_source;
}

//---------------------------------------------------------------------------//
// source passing implementation
//---------------------------------------------------------------------------//
// distribute the census particles

template<class MT>
template<class PT>
void Parallel_Builder<MT>::dist_census(const Source_Init<MT> &sinit, 
				       const Particle_Buffer<PT> &buffer,
				       typename Particle_Buffer<PT>::Census
				       &census_bank) 
{
    // get number of cells on the global mesh
    int num_cells = procs_per_cell.size();
    Require (num_cells > 0);
    Require (census_bank.size() == 0);

    // calculate census particle-to-processor data
    vector<int> ncen2proc(num_cells);
    vector<int> ncenleft(num_cells);
    for (int cell = 1; cell <= num_cells; cell++)
    {
	ncen2proc[cell-1] = sinit.get_ncen(cell) /
	    procs_per_cell[cell-1].size();
	ncenleft[cell-1] = sinit.get_ncen(cell) %
	    procs_per_cell[cell-1].size();
	Check (ncenleft[cell-1] < procs_per_cell[cell-1].size());
    }

    // calculate and send census to processors

    // get the census buffer
    SP<Particle_Buffer<PT>::Census> old_census = sinit.get_census();
    Check (old_census);
    Check (old_census->size() == sinit.get_ncentot());

    // initialize counters for census read/write
    int total_placed = 0;
    int total_read = 0;
    vector<int> num_placed_per_cell(num_cells);
    vector<int> num_placed_per_proc(nodes());
    vector<int> num_to_send(nodes());

    // read census file and distribute to IMC_Procs

    // necessary variables
    int cell;
    int offset;
    int proc_index;
    int proc_goto;

    // Census Particle Buffers
    vector<Particle_Buffer<PT>::Comm_Buffer> cen_buffer(nodes());

    // loop to read census particles and send them
    while (old_census->size())
    {
	// get a particle from the old census
	SP<PT> particle = old_census->top();
	old_census->pop();
	total_read++;

	// determine the cell for this particle
	cell = particle->get_cell();

	// in a cell, total number of census particles after which the
	// leftover particles are placed
	offset = (ncen2proc[cell-1] + 1) * ncenleft[cell-1];

	// determine the processor to go to
	if (num_placed_per_cell[cell-1] < offset)
	    proc_index = num_placed_per_cell[cell-1] / 
		(ncen2proc[cell-1] + 1); 
	else
	    proc_index = (num_placed_per_cell[cell-1] - offset) / 
		ncen2proc[cell-1] + ncenleft[cell-1];
	proc_goto = procs_per_cell[cell-1][proc_index];

	// update the particle cell with the local index of the cell on the
	// IMC_processor
	particle->set_cell(imc_cell(cell, proc_goto));
	buffer.buffer_particle(cen_buffer[proc_goto], *particle);

	// increment some counters
	num_placed_per_cell[cell-1]++;
	num_placed_per_proc[proc_goto]++;
	total_placed++;
	num_to_send[proc_goto]++;

	// send these guys out
	if (num_to_send[proc_goto] ==
	    Particle_Buffer<PT>::get_buffer_s())
	{
	    if (!proc_goto)
	    {
		// if we are on the host fill up the census bank
		buffer.add_to_bank(cen_buffer[proc_goto], census_bank);
		Check (cen_buffer[proc_goto].n_part == 0);
		num_to_send[proc_goto] = 0;
	    }
	    else
	    {
		// if we are on an IMC processor send these guys
		buffer.send_buffer(cen_buffer[proc_goto], proc_goto);
		cen_buffer[proc_goto].n_part = 0;
		num_to_send[proc_goto] = 0;
	    }
	}
    }
	    
    // some Checks
    Check (total_placed == sinit.get_ncentot());
    Check (total_read == sinit.get_ncentot());

    // send the guys off that haven't been sent yet

    // on the host
    if (num_to_send[0] > 0)
    {
	buffer.add_to_bank(cen_buffer[0], census_bank);
	Check (cen_buffer[0].n_part == 0);
    }

    // on the IMC processors
    for (int proc = 1; proc < nodes(); proc++)
    {
	if (num_to_send[proc] > 0)
	{
	    buffer.send_buffer(cen_buffer[proc], proc);
	    cen_buffer[proc].n_part = 0;
	    num_to_send[proc] = 0;  
	}
    }

    // send a final message to indicate completion
    for (int proc = 1; proc < nodes(); proc++)
    {
	Check (cen_buffer[proc].n_part == 0);
	cen_buffer[proc].n_part = -1;
	buffer.send_buffer(cen_buffer[proc], proc);
    }

    // the old census should now be zero
    Ensure (old_census->size() == 0);
}

//---------------------------------------------------------------------------//
// receive the census buffers and write out the census files

template<class MT>
template<class PT>
void Parallel_Builder<MT>::recv_census(const Particle_Buffer<PT> &buffer,
				       typename Particle_Buffer<PT>::Census
				       &census_bank) 
{
    Require (census_bank.size() == 0);

    // condition to keep receiving
    bool receive = true;
    
    // receive the Comm_Buffers from the host
    while (receive)
    {
	// receive the buffer
	SP<Particle_Buffer<PT>::Comm_Buffer> cen_buffer =
	    buffer.recv_buffer(0); 

	// if the buffer has particles bankit, else stop this nonsense
	if (cen_buffer->n_part >= 0)
	    buffer.add_to_bank(*cen_buffer, census_bank);
	else if (cen_buffer->n_part < 0)
	    receive = false;
	else
	    Insist(0, "Buffers are < 0 during census remap!");
    }	
}

//---------------------------------------------------------------------------//
// distribute the volume source stuff
	    
template<class MT> 
void Parallel_Builder<MT>::dist_vol(const Source_Init<MT> &sinit,
				    typename MT::CCSF_int &volrn,
				    typename MT::CCSF_int &numvol,
				    typename MT::CCSF_double &ew,
				    typename MT::CCVF_double &t4_slope)
{
    // find the number of cells on the global mesh
    int num_cells = procs_per_cell.size();
    int dimension = numvol.get_Mesh().get_Coord().get_dim();
    Require (num_cells > 0);

    // make the volume source counters
    vector<int> nvol(num_cells);
    vector<int> nvol_xtra(num_cells);
    vector<int> streamnum(num_cells);
    int counter = rtt_rng::rn_stream;
    for (int cell = 1; cell <= num_cells; cell++)
    {
	nvol[cell-1] = sinit.get_nvol(cell) / procs_per_cell[cell-1].size(); 
	nvol_xtra[cell-1] = sinit.get_nvol(cell) % 
	    procs_per_cell[cell-1].size();
	streamnum[cell-1] = counter;
	counter += sinit.get_nvol(cell);
	Check (nvol_xtra[cell-1] < procs_per_cell[cell-1].size());
    }

    // update the Global rn_stream counter
    rtt_rng::rn_stream += sinit.get_nvoltot();
    Check (counter == rtt_rng::rn_stream);

    // loop over processors and make the send info
    int voltot = 0;
    for (int i = 0; i < nodes(); i++)
    {
	// calculate the number of cells on proc i mesh
	int num_cells = cells_per_proc[i].size();

	// make a c-style array for passing to other nodes
	int *nvol_send   = new int[num_cells];
	int *stream_send = new int[num_cells];
	double *ew_send  = new double[num_cells];
	double *t4_send  = new double[num_cells * dimension];
	
	// loop through cells on this processor and send stuff out
	for (int cell = 1; cell <= num_cells; cell++)
	{
	    // determine the global cell index
	    int global_cell = cells_per_proc[i][cell-1];

	    // calculate the nvol source, ew, and slope
	    nvol_send[cell-1] = nvol[global_cell-1];
	    ew_send[cell-1]   = sinit.get_ew_vol(global_cell);
	    for (int d = 0; d < dimension; d++)
	    {
		int ptr               = d * num_cells;
		t4_send[ptr + cell-1] = sinit.get_t4_slope(d+1, global_cell); 
	    } 

	    // check to see if we need to add extra sources to this cell on
	    // this processor
	    if (nvol_xtra[global_cell-1] >= i+1)
		nvol_send[cell-1]++;

	    // calculate the stream number for proc i
	    stream_send[cell-1] = streamnum[global_cell-1];
	    streamnum[global_cell-1] += nvol_send[cell-1];

	    // do some check counters
	    voltot += nvol_send[cell-1];
	}

	// send out the data 
	if (i)
	{
	    // send to proc i if we are not on the host
	    Send (num_cells, i, 24);
	    Send <int>(nvol_send, num_cells, i, 25);
	    Send <int>(stream_send, num_cells, i, 26);
	    Send <double>(ew_send, num_cells, i, 27);
	    Send <double>(t4_send, num_cells * dimension, i, 28);
	}
	else
	{
	    // return the needed source stuff for the host
	  
	    // some assertions	    
	    Check (num_cells == numvol.get_Mesh().num_cells());

	    // assign the data
	    for (int n = 1; n <= num_cells; n++)
	    {
		volrn(n)  = stream_send[n-1];
		numvol(n) = nvol_send[n-1];
		ew(n)     = ew_send[n-1];
		for (int d = 0; d < dimension; d++)
		{
		    int ptr          = d * num_cells;
		    t4_slope(d+1, n) = t4_send[ptr + n-1];
		}
	    }
	}

	// reclaim dynamic memory
	delete [] nvol_send;
	delete [] stream_send;
	delete [] ew_send;
	delete [] t4_send;
    }

    // final assertion
    Ensure (voltot == sinit.get_nvoltot());
}

//---------------------------------------------------------------------------//
// receive the volume source stuff

template<class MT>
void Parallel_Builder<MT>::recv_vol(typename MT::CCSF_int &volrn,
				    typename MT::CCSF_int &nvol, 
				    typename MT::CCSF_double &ew,
				    typename MT::CCVF_double &t4)
{
    Require (volrn.get_Mesh() == nvol.get_Mesh());
    Require (ew.get_Mesh()    == nvol.get_Mesh());
    Require (t4.get_Mesh()    == t4.get_Mesh());

    // lets get the size of the thing from the host processor
    int num_cells;
    int dimension = volrn.get_Mesh().get_Coord().get_dim();
    Recv (num_cells, 0, 24);
    Check (num_cells == volrn.get_Mesh().num_cells());

    // receive the volume source data from the host processor
    int *nvol_recv   = new int[num_cells];
    int *stream_recv = new int[num_cells];
    double *ew_recv  = new double[num_cells];
    double *t4_recv  = new double[num_cells * dimension];
    Recv (nvol_recv, num_cells, 0, 25);
    Recv (stream_recv, num_cells, 0, 26);
    Recv (ew_recv, num_cells, 0, 27);
    Recv (t4_recv, num_cells * dimension, 0, 28);

    // assign data to the CCSFs
    for (int cell = 1; cell <= num_cells; cell++)
    {
	nvol(cell)  = nvol_recv[cell-1];
	volrn(cell) = stream_recv[cell-1];
	ew(cell)    = ew_recv[cell-1];
	for (int d = 0; d < dimension; d++)
	{
	    int ptr       = d * num_cells;
	    t4(d+1, cell) = t4_recv[ptr + cell-1];
	}
    }

    // release the dynamic storage
    delete [] nvol_recv;
    delete [] stream_recv;
    delete [] ew_recv;
    delete [] t4_recv;
}

//---------------------------------------------------------------------------//
// distribute the surface source stuff
	    
template<class MT>
void Parallel_Builder<MT>::dist_ss(const Source_Init<MT> &sinit,
				   typename MT::CCSF_int &ssrn,
				   typename MT::CCSF_int &numss,
				   typename MT::CCSF_int &fss,
				   typename MT::CCSF_double &ew)
{
    // find the number of cells on the global mesh
    int num_cells = procs_per_cell.size();
    Require (num_cells > 0);

    // make the surface source counters
    vector<int> nss(num_cells);
    vector<int> nss_xtra(num_cells);
    vector<int> streamnum(num_cells);
    int counter = rtt_rng::rn_stream;
    for (int cell = 1; cell <= num_cells; cell++)
    {
	nss[cell-1] = sinit.get_nss(cell) / procs_per_cell[cell-1].size();
	nss_xtra[cell-1] = sinit.get_nss(cell) % 
	    procs_per_cell[cell-1].size();
	streamnum[cell-1] = counter;
	counter += sinit.get_nss(cell);
	Check (nss_xtra[cell-1] < procs_per_cell[cell-1].size());
    }
    
    // update the Global rn_stream counter
    rtt_rng::rn_stream += sinit.get_nsstot();
    Check (counter == rtt_rng::rn_stream);

    // loop over processors and make the send info
    int sstot = 0;
    for (int i = 0; i < nodes(); i++)
    {
	// calculate the number of cells on proc i mesh
	int num_cells = cells_per_proc[i].size();

	// make a c-style array for passing to other nodes
	int *nss_send    = new int[num_cells];
	int *stream_send = new int[num_cells];
	int *fss_send    = new int[num_cells];
	double *ew_send  = new double[num_cells];
	
	// loop through cells on this processor and send stuff out
	for (int cell = 1; cell <= num_cells; cell++)
	{
	    // determine the global cell index
	    int global_cell = cells_per_proc[i][cell-1];

	    // calculate the nss source, fss, and ew
	    nss_send[cell-1] = nss[global_cell-1];
	    fss_send[cell-1] = sinit.get_fss(global_cell);
	    ew_send[cell-1]  = sinit.get_ew_ss(global_cell);

	    // check to see if we need to add extra sources to this cell on
	    // this processor
	    if (nss_xtra[global_cell-1] >= i+1)
		nss_send[cell-1]++;
	    
	    // calculate the stream number for proc i
	    stream_send[cell-1] = streamnum[global_cell-1];
	    streamnum[global_cell-1] += nss_send[cell-1];

	    // do some check counters
	    sstot += nss_send[cell-1];
	}

	// send out the data 
	if (i)
	{
	    // send to proc i if we are not on the host
	    Send (num_cells, i, 29);
	    Send <int>(nss_send, num_cells, i, 30);
	    Send <int>(stream_send, num_cells, i, 31);
	    Send <int>(fss_send, num_cells, i, 32);
	    Send <double>(ew_send, num_cells, i, 33);
	}
	else
	{
	    // return the needed source stuff for the host
	  
	    // some assertions
	    Check (num_cells == numss.get_Mesh().num_cells());

	    // assign the vectors to our point array
	    for (int n = 1; n <= num_cells; n++)
	    {
		ssrn(n)  = stream_send[n-1];
		numss(n) = nss_send[n-1];
		fss(n)   = fss_send[n-1];
		ew(n)    = ew_send[n-1];
	    }
	}

	// reclaim dynamic memory
	delete [] nss_send;
	delete [] stream_send;
	delete [] fss_send;
	delete [] ew_send;
    }

    // final assertion
    Ensure (sstot == sinit.get_nsstot());
}

//---------------------------------------------------------------------------//
// receive the surface source stuff

template<class MT>
void Parallel_Builder<MT>::recv_ss(typename MT::CCSF_int &ssrn,
				   typename MT::CCSF_int &nss,
				   typename MT::CCSF_int &fss,
				   typename MT::CCSF_double &ew)
{
    Require (ssrn.get_Mesh() == nss.get_Mesh());
    Require (ew.get_Mesh()   == nss.get_Mesh());
    Require (fss.get_Mesh()  == nss.get_Mesh());

    // lets get the size of the thing from the host processor
    int num_cells;
    Recv (num_cells, 0, 29);
    Check (num_cells == ssrn.get_Mesh().num_cells());

    // receive the surface source data from the host processor
    int *nss_recv    = new int[num_cells];
    int *stream_recv = new int[num_cells];
    int *fss_recv    = new int[num_cells];
    double *ew_recv  = new double[num_cells];
    Recv (nss_recv, num_cells, 0, 30);
    Recv (stream_recv, num_cells, 0, 31);
    Recv (fss_recv, num_cells, 0, 32);
    Recv (ew_recv, num_cells, 0, 33);

    // assign data to the CCSFs
    for (int cell = 1; cell <= num_cells; cell++)
    {
	nss(cell)  = nss_recv[cell-1];
	ssrn(cell) = stream_recv[cell-1];
	fss(cell)  = fss_recv[cell-1];
	ew(cell)   = ew_recv[cell-1];
    }

    // release the dynamic storage
    delete [] nss_recv;
    delete [] stream_recv;
    delete [] fss_recv;
    delete [] ew_recv;
}

//---------------------------------------------------------------------------//
// Mesh passing interface
//---------------------------------------------------------------------------//
// send out the Mesh components and return a new mesh to the host process

template<class MT>
SP<MT> Parallel_Builder<MT>::send_Mesh(const MT &mesh)
{
    // send the mesh components out to the other processors

    // assure that we are on the host node only
    Check (!node());

    // new host mesh objects
    SP<Coord_sys> host_coord =  mesh.get_SPCoord();
    Layout host_layout = mesh.get_Layout();
    typename MT::NCVF_d host_vertex(mesh.get_Coord().get_dim());
    typename MT::CCVF_i host_cell_pair(cells_per_proc[0].size());
    typename MT::CCSF_i host_generation(cells_per_proc[0].size());
    typename MT::CCSF_b host_has_kids(cells_per_proc[0].size());
    SP<MT> host_mesh;

    // let us pass a coordinate system, shall we
    send_Coord(host_coord);

    // let us pass the Layout
    send_Layout(host_layout);

    // let us pass the vertex, cell_pairings, and generations, ie. the cells
    send_cells(mesh, host_vertex, host_cell_pair, host_generation);

  // build the new mesh on the host, we assume it is a submesh
    host_mesh = new MT(host_coord, host_layout, host_vertex,
		       host_cell_pair, host_generation, true);

    return host_mesh;
}

//---------------------------------------------------------------------------//
// receive the Mesh components

template<class MT>
SP<MT> Parallel_Builder<MT>::recv_Mesh()
{
    // receive the Mesh components from the host and build the new mesh on 
    // this topology

    // assure that we are not on the host node
    Check (node());

    // rebuilt Mesh SP
    SP<MT> return_mesh;

    // get coordinate system
    SP<Coord_sys> coord = recv_Coord();

    // get Layout
    Layout layout = recv_Layout();

    // get vertices and cell_pair
    typename MT::NCVF_d vertex    = recv_vertex();
    typename MT::CCVF_i cell_pair = recv_cellpair();
    typename MT::CCSF_i generation = recv_generation();
    typename MT::CCSF_b has_kids;

  // build mesh, we assume that this is a submesh!!!
  // SOLARIS MPICH fails at this point
    return_mesh = new MT(coord, layout, vertex, cell_pair, generation, true);

    // return mesh
    return return_mesh;
}

//---------------------------------------------------------------------------//
// Mesh passing implementation
//---------------------------------------------------------------------------//
// pass the Coord_sys object

template<class MT>
void Parallel_Builder<MT>::send_Coord(SP<Coord_sys> coord)
{
    // send out the coordinate system designator

    // send variables
    string cs    = coord->get_Coord();
    const char *sendcs = cs.c_str();
    int cs_size  = cs.size() + 1;

    // send the Coord_sys string
    for (int np = 1; np < nodes(); np++)
    {
	Send (cs_size, np, 1);
	Send (sendcs, cs_size, np, 2);
    }
}

//---------------------------------------------------------------------------//
// receive the Coord_sys

template<class MT>
SP<Coord_sys> Parallel_Builder<MT>::recv_Coord()
{
    // receive and build the coordinate system

    // define Coord_sys SP
    SP<Coord_sys> coord;

    // receive the designation
    int cs_size;
    Recv (cs_size, 0, 1);
    char *reccs = new char[cs_size];
    Recv (reccs, cs_size, 0, 2);
    string cs = reccs;
    delete [] reccs;

    // build the coordinate system
    if (cs == "xy")
    {
	SP<XYCoord_sys> xycoord(new XYCoord_sys);
	coord = xycoord;
    }
    else if (cs == "xyz")
    {
	SP<XYZCoord_sys> xyzcoord(new XYZCoord_sys);
	coord = xyzcoord;
    }
    
    // return base class SP to a derived Coord_sys
    return coord;
}

//---------------------------------------------------------------------------//
// pass the Layout

template<class MT>
void Parallel_Builder<MT>::send_Layout(Layout &host_layout)
{
    // make sure we have the global layout at this junction
    Require (host_layout.num_cells() == procs_per_cell.size());

    // size the boundary cells object
    bound_cells.resize(nodes());
    for (int i = 0; i < bound_cells.size(); i++)
	bound_cells[i].resize(0);

    // loop through processors and send the layouts
    for (int np = 1; np < nodes(); np++)
    {
	// get the Layout for processor np
	Layout layout = build_Layout(host_layout, np);

	// set the Layout size
	int num_cells = layout.num_cells();
	Check (num_cells == cells_per_proc[np].size());

	// calculate the number of faces on each cell, the number of adjacent 
	// cells on each face, and the total size of the Layout (layout_size =
	//  SUM_(i)SUM_(f)^(num_cells) (num_faces[i]) (num_adj[i][f])
	int * num_faces = new int[num_cells];
	int layout_size = 0;
	int faces_size = 0;
	for (int c = 1; c <= num_cells; c++)
	{
	    num_faces[c-1] = layout.num_faces(c);
	    for (int f = 1; f <= num_faces[c-1]; f++)
	    {
	        layout_size += num_faces[c-1]*layout.num_adj(c,f);
	    }
	    faces_size += layout.num_faces(c);
	}

	// write the faces and num_adj arrays for passing
	int * num_adj = new int[faces_size];
	int * faces = new int[layout_size];
	int face_index  = 0;
	int adj_index  = 0;
	for (int c = 1; c <= num_cells; c++)
	{
	    for (int f = 1; f <= layout.num_faces(c); f++)
	    {
	        num_adj[face_index] = layout.num_adj(c,f);
	        for (int a = 1; a <= layout.num_adj(c,f); a++)
		{
		    faces[adj_index] = layout(c,f,a);
		    adj_index++;
		}
		face_index++;
	    }
	}

	// pass the Layout cell size
	Send (num_cells, np, 3);

	// pass the Layout number of adjacent cells size
	Send (faces_size, np, 4);

	// pass the Layout total size
	Send (layout_size, np, 5);

	// pass the num_faces array
	Send <int>(&num_faces[0], num_cells, np, 6);

	// pass the num_adj array
	Send <int>(&num_adj[0], faces_size, np, 7);

	// pass the face-values array
	Send <int>(faces, layout_size, np, 8);

	// delete dynamically allocated faces array
	delete [] faces;
	delete [] num_adj;
	delete [] num_faces;
    }

    // now rebuild the host layout
    host_layout = build_Layout(host_layout, 0);
    Ensure (host_layout.num_cells() == cells_per_proc[0].size());
}

//---------------------------------------------------------------------------//
// receive the Layout

template<class MT>
Layout Parallel_Builder<MT>::recv_Layout()
{
    // receive and build the Layout

    // get the number of cells and size of the Layout
    int num_cells;
    int faces_size;
    int layout_size;
    Recv (num_cells, 0, 3);
    Recv (faces_size, 0, 4);
    Recv (layout_size, 0, 5);

    // make the Layout
    Layout layout = num_cells;

    // receive the Layout data arrays
    int * num_faces = new int[num_cells];
    int * num_adj = new int[faces_size];
    int * faces = new int[layout_size];
    Recv (num_faces, num_cells, 0, 6);
    Recv (num_adj, faces_size, 0, 7);
    Recv (faces, layout_size, 0, 8);

    // rebuild the Layout
    int face_index = 0;
    int adj_index = 0;
    for (int cell = 1; cell <= num_cells; cell++)
    {
	layout.set_size(cell, num_faces[cell-1]);
	for (int face = 1; face <= num_faces[cell-1]; face++)
	{
	    layout.set_adj_size(cell, face, num_adj[face_index]);
	    for (int adj = 1; adj <= num_adj[face_index]; adj++)
	    {
	        layout(cell, face, adj) = faces[adj_index];
		adj_index++;
	    }
	    face_index++;
	}
    }

    // get rid of dynamic arrays
    delete [] num_faces;
    delete [] num_adj;
    delete [] faces;
    
    // return the new Layout on this node
    return layout;
}

//---------------------------------------------------------------------------//
// build a new layout on each processor

template<class MT>
Layout Parallel_Builder<MT>::build_Layout(const Layout &layout, int proc)
{
    // assure that we have the global Layout
    Require (layout.num_cells() == procs_per_cell.size());

    // determine the number of cells on this processor
    int num_cells = cells_per_proc[proc].size();
    Layout imc_layout(num_cells);

    // loop over cells on IMC processor and calculate new Layout
    int global_cell;
    int next_global_cell;
    int next_imc_cell;
    for (int cell = 1; cell <= num_cells; cell++)
    {
	// determine the global cell index
	global_cell = cells_per_proc[proc][cell-1];
	
	// size the new layout to the proper number of faces
	imc_layout.set_size(cell, layout.num_faces(global_cell));
	
	// loop through faces and build the Layout for this IMC cell
	for (int face = 1; face <= imc_layout.num_faces(cell); face++)
	{
	    // size the new layout to the proper number of adjacent cells
	    imc_layout.set_adj_size(cell, face,
				    layout.num_adj(global_cell,face));
	    for (int adj = 1; adj <= imc_layout.num_adj(cell, face); adj++)
	    {
	        next_global_cell = layout(global_cell, face, adj);
		next_imc_cell    = imc_cell(next_global_cell, proc);
		if (!next_imc_cell && next_global_cell)
		{
		    bound_cells[proc].push_back(next_global_cell);
		    next_imc_cell = -bound_cells[proc].size();
		}
		imc_layout(cell, face, adj) = next_imc_cell;
	    }
	}
    }
    
    // return the new Layout
    return imc_layout;
}

//---------------------------------------------------------------------------//
// send out the cells in the form of cell vertices and cell_pair arrays

template<class MT>
void Parallel_Builder<MT>::send_cells(const MT &host_mesh,
				      typename MT::NCVF_d &host_vertex,
				      typename MT::CCVF_i &host_cellpair,
				      typename MT::CCSF_i &host_generation)
{
    // loop over processors and build new vertices and cell_pairs and send them 
    // out
    for (int np = 1; np < nodes(); np++)
    {
	// build new vertices, cell_pairs, and generations on this node
	typename MT::NCVF_d imc_vertex(host_mesh.get_Coord().get_dim());
	typename MT::CCVF_i imc_cellpair(cells_per_proc[np].size());
	typename MT::CCSF_i imc_generation(cells_per_proc[np].size());
	build_cells(host_mesh, imc_vertex, imc_cellpair,imc_generation, np);

	// send them to the other processors
	send_vertex(imc_vertex, np);
	send_cellpair(imc_cellpair, np);
 	send_generation(imc_generation, np);
   }

    // now rebuild the host vertices and cell_pair
    build_cells(host_mesh, host_vertex, host_cellpair, host_generation, 0);
    Ensure (host_cellpair.size() == cells_per_proc[0].size());
}

//---------------------------------------------------------------------------//
// build a new vertex and cell_pair for each processor

template<class MT>
void Parallel_Builder<MT>::build_cells(const MT &mesh,
				       typename MT::NCVF_d &vertex, 
				       typename MT::CCVF_i &cellpair,
				       typename MT::CCSF_i &generation,
				       int proc)
{
    Require (cells_per_proc[proc].size() == cellpair.size());
    Require (mesh.get_Coord().get_dim() == vertex.size());
    Require (cells_per_proc[proc].size() == generation.size());

    // determine dimension and number of cells
    int num_cells = cells_per_proc[proc].size();
    int dim       = vertex.size();

    // loop over cells on IMC processor and calculate new vertex
    int global_cell;
    for (int cell = 1; cell <= num_cells; cell++)
    {
	// determine the global_cell index
	global_cell = cells_per_proc[proc][cell-1];
	
	// get vertices from the global mesh
	typename MT::NCVF_d cell_vert = mesh.get_vertices(global_cell);
	Check (cell_vert.size() == dim);
	int num_vert = cell_vert[0].size();

	// assign these vertices to the new IMC vertices and cell_pair
	for (int i = 0; i < num_vert; i++)
	{
	    for (int d = 0; d < dim; d++)
	    {
		vertex[d].push_back(cell_vert[d][i]);
		Check (num_vert == cell_vert[d].size());
	    }
	    cellpair[cell-1].push_back(vertex[0].size());
	}
	Check (cellpair[cell-1].size() == num_vert);
	// assign the generation to this cell
	generation[cell-1] = mesh.get_generation(global_cell);
    }
}

//---------------------------------------------------------------------------//
// pass the mesh vertices

template<class MT>
void Parallel_Builder<MT>::send_vertex(const typename MT::NCVF_d &vertex, 
				       int proc)
{
    // send the vertices to another processor
    
    // get all necessary vertex dimensions
    int vertex_dim  = vertex.size();
    int vertex_size = vertex[0].size();
    int total_size  = vertex_dim * vertex_size;

    // push all vertices onto one array for communication
    double *vert = new double[total_size];
    int index = 0;
    for (int d = 0; d < vertex_dim; d++)
    {
	// some assertions to check that the vertex size is constant over all
	// dimensions 
	Check (vertex_size == vertex[d].size());

	// now load up the communication vertex array
	for (int i = 0; i < vertex_size; i++)
	{
	    vert[index] = vertex[d][i];
	    index++;
	}
    }

    // send the dimensionality of the vertex-array
    Send (vertex_dim, proc, 9);

    // send the total size of the vertex-array
    Send (total_size, proc, 10);

    // send the vertex array
    Send <double>(vert, total_size, proc, 11);

    // delete dynamically allocated arrays
    delete [] vert;
}
	   
//---------------------------------------------------------------------------//
// receive the vertex

template<class MT>
typename MT::NCVF_d Parallel_Builder<MT>::recv_vertex()
{
    // receive and rebuild the vertices

    // first get the sizes
    int vertex_dim;
    int total_size;
    Recv (vertex_dim, 0, 9);
    Recv (total_size, 0, 10);

    // find the size of each dimensional vertex array
    Check (!(total_size % vertex_dim));
    int vertex_size = total_size / vertex_dim;

    // make a new NCVF_d vertex array
    typename MT::NCVF_d vertex(vertex_dim);

    // get the vertices from the host
    double *vert = new double[total_size];
    Recv (vert, total_size, 0, 11);

    // assign values to the new vertex object
    int index = 0;
    for (int d = 0; d < vertex_dim; d++)
    {
	vertex[d].resize(vertex_size);
	for (int i = 0; i < vertex_size; i++)
	{
	    vertex[d][i] = vert[index];
	    index++;
	}
    }
    Check (index == total_size);

    // delete dynamic allocated arrays
    delete [] vert;

    // return the new vertex guy
    return vertex;
}

//---------------------------------------------------------------------------//
// pass the cell_pair

template<class MT>
void Parallel_Builder<MT>::send_cellpair(const typename MT::CCVF_i &cellpair, 
					 int proc)  
{
    // send the cell_pair array to another processor
    
    Require (cellpair.size() == cells_per_proc[proc].size())

	// set the cell_pair size
	int num_cells = cellpair.size();

    // calculate the number of vertices per cell and the total size of the
    // cell_pair object
    int *num_vert = new int[num_cells];
    int size = 0;
    for (int i = 0; i < num_cells; i++)
    {
	num_vert[i] = cellpair[i].size();
	size += num_vert[i];
    }

    // write the cell_pair array for passing
    int *vertices = new int[size];
    int index = 0;
    for (int i = 0; i < num_cells; i++)
	for (int j = 0; j < cellpair[i].size(); j++)
	{
	    vertices[index] = cellpair[i][j];
	    index++;
	}
    Check (index == size);

    // pass the size
    Send (num_cells, proc, 12);

    // pass the total size of the cell_pair object
    Send (size, proc, 13);

    // pass the num_vert array
    Send <int>(num_vert, num_cells, proc, 14);
	
    // pass the vertices-values array
    Send <int>(vertices, size, proc, 15);

    // delete the dynamically allocated arrays
    delete [] num_vert;
    delete [] vertices;
}

//---------------------------------------------------------------------------//
// receive the cell_pair object

template<class MT>
typename MT::CCVF_i Parallel_Builder<MT>::recv_cellpair()
{
    // receive and rebuild the cell_pair array
    
    // first get the sizes 
    int num_cells;
    int size;
    Recv (num_cells, 0, 12);
    Recv (size, 0, 13);

    // make a new cell_pair object
    typename MT::CCVF_i cell_pair(num_cells);

    // receive the cell_pair data
    int *num_vert = new int[num_cells];
    int *vertices = new int[size];
    Recv (num_vert, num_cells, 0, 14);
    Recv (vertices, size, 0, 15);

    // rebuild the cell_pair object
    int index = 0;
    for (int i = 0; i < num_cells; i++)
    {
	cell_pair[i].resize(num_vert[i]);
	for (int j = 0; j < num_vert[i]; j++)
	{
	    cell_pair[i][j] = vertices[index];
	    index++;
	}
    }
    Check (index == size);

    // delete dynamic arrays
    delete [] num_vert;
    delete [] vertices;

    // return the cell_pair object
    return cell_pair;
}

//---------------------------------------------------------------------------//
// pass the mesh generations

template<class MT>
void Parallel_Builder<MT>::send_generation(const typename MT::CCSF_i & 
					   generation, int proc)
{
    // send the generations to another processor
    
    // get all necessary generation dimensions
    int num_cells = generation.size();
    int *gen = new int[num_cells];
    for (int i = 0; i < num_cells; i++)
	gen[i] = generation[i];

    // send the dimensionality of the generation-array
    Send(num_cells, proc, 16);

    // send the generation array
    Send <int>(gen, num_cells, proc, 17);

    delete [] gen;
}

//---------------------------------------------------------------------------//
// receive the generation

template<class MT>
typename MT::CCSF_i Parallel_Builder<MT>::recv_generation()
{
    // receive and rebuild the generations

    // first get the sizes
    int num_cells;
    Recv (num_cells, 0, 16);
    int *gen = new int[num_cells];

    // make a new CCSF_i generation array
    typename MT::CCSF_i generation(num_cells);

    // get the generations from the host
    Recv (gen, num_cells, 0, 17);

    for (int i = 0; i < num_cells; i++)
	generation[i] = gen[i];

    delete [] gen;

    // return the new generation
    return generation;
}

//---------------------------------------------------------------------------//
// Opacity passing interface
//---------------------------------------------------------------------------//
// send the Opacity object

template<class MT> SP<Opacity<MT> > 
Parallel_Builder<MT>::send_Opacity(SP<MT> mesh, const Opacity<MT> &opacity) 
{
    // send out the Opacities, one component at a time

    // assure that we are on the host node
    Require (!node());

    // we must have a local mesh and a global opacity
    Require (mesh);
    Require (opacity.num_cells() == procs_per_cell.size());

    // data necessary to build Opacity on host
    SP<Opacity<MT> > host_opacity;
    typename MT::CCSF_double sigma(mesh);
    typename MT::CCSF_double sigma_thomson(mesh);
    typename MT::CCSF_double planck(mesh);
    typename MT::CCSF_double fleck(mesh);

    // loop over procs and send out the Opacities
    for (int np = 0; np < nodes(); np++)
    {
	// determine the number of cells
	int num_cells = cells_per_proc[np].size();

	// assign the Opacity data
	double *sigma_send      = new double[num_cells];
	double *sigma_thom_send = new double[num_cells];
	double *planck_send     = new double[num_cells];
	double *fleck_send      = new double[num_cells];
    
	int global_cell;
	for (int cell = 1; cell <= num_cells; cell++)
	{
	    // get the global_cell index
	    global_cell = cells_per_proc[np][cell-1];

	    // assign the data
	    sigma_send[cell-1]      = opacity.get_sigma_abs(global_cell); 
	    sigma_thom_send[cell-1] = opacity.get_sigma_thomson(global_cell); 
	    planck_send[cell-1]     = opacity.get_planck(global_cell); 
	    fleck_send[cell-1]      = opacity.get_fleck(global_cell);
	}

	// send the Opacity data to the IMC nodes
	if (np)
	{
	    Send (num_cells, np, 20);
	    Send <double>(sigma_send, num_cells, np, 21);
	    Send <double>(sigma_thom_send, num_cells, np, 47);
	    Send <double>(planck_send, num_cells, np, 22);
	    Send <double>(fleck_send, num_cells, np, 23);
	}
	else
	{
	    // build the host processor CCSF
	    
	    // check to make sure the local mesh is the right size
	    Check (mesh->num_cells() == num_cells);
	    
	    // assign the data
	    for (int cell = 1; cell <= num_cells; cell++)
	    {
		sigma(cell)         = sigma_send[cell-1];
		sigma_thomson(cell) = sigma_thom_send[cell-1];
		planck(cell)        = planck_send[cell-1];
		fleck(cell)         = fleck_send[cell-1];
	    }
	}

	// delete dynamic allocation
	delete [] sigma_send;
	delete [] sigma_thom_send;
	delete [] planck_send;
	delete [] fleck_send;
    }
    
    // make and return the Opacity to the host
    host_opacity = new Opacity<MT>(sigma, sigma_thomson, planck, fleck);
    Ensure (host_opacity->num_cells() == mesh->num_cells());
    return host_opacity;
}

//---------------------------------------------------------------------------//
// receive the Opacity object

template<class MT>
SP<Opacity<MT> > Parallel_Builder<MT>::recv_Opacity(SP<MT> mesh)
{
    // receive and rebuild the Opacity object

    // assure we are on receive nodes and have a valid mesh
    Require (node());
    Require (mesh);

    // declare return opacity object
    SP< Opacity<MT> > imc_opacity;
    typename MT::CCSF_double sigma(mesh);
    typename MT::CCSF_double sigma_thomson(mesh);
    typename MT::CCSF_double planck(mesh);
    typename MT::CCSF_double fleck(mesh);

    // receive the size of this guy
    int num_cells;
    Recv (num_cells, 0, 20);
    Check (num_cells == mesh->num_cells());

    // receive data from host
    double *rsigma      = new double[num_cells];
    double *rsigma_thom = new double[num_cells];
    double *rplanck     = new double[num_cells];
    double *rfleck      = new double[num_cells];
    Recv (rsigma, num_cells, 0, 21);
    Recv (rsigma_thom, num_cells, 0, 47);
    Recv (rplanck, num_cells, 0, 22);
    Recv (rfleck, num_cells, 0, 23);

    // assign to new Opacity objects
    for (int cell = 1; cell <= num_cells; cell++)
    {
	sigma(cell)         = rsigma[cell-1];
	sigma_thomson(cell) = rsigma_thom[cell-1];
	planck(cell)        = rplanck[cell-1];
	fleck(cell)         = rfleck[cell-1];
    }

    // reclaim dynamic memory
    delete [] rsigma;
    delete [] rsigma_thom;
    delete [] rplanck;
    delete [] rfleck;

    // build and return new opacity object
    imc_opacity = new Opacity<MT>(sigma, sigma_thomson, planck, fleck);
    Ensure (imc_opacity->num_cells() == num_cells);
    return imc_opacity;
}

//---------------------------------------------------------------------------//
// Mat_State passing interface
//---------------------------------------------------------------------------//
// send out the Mat_State and build it on this processor

template<class MT> SP<Mat_State<MT> >
Parallel_Builder<MT>::send_Mat(SP<MT> mesh, const Mat_State<MT> &mat) 
{
    // assure that we are on the host node
    Require (!node());
    
    // we must have a local mesh and a global Mat_State
    Require (mesh);
    Require (mat.num_cells() == procs_per_cell.size());

    // data necessary to build Mat_State on host
    // specific heat is not needed here, since dedt is current and valid
    SP<Mat_State<MT> > host_mat;
    typename MT::CCSF_double density(mesh);
    typename MT::CCSF_double T(mesh);
    typename MT::CCSF_double dedt(mesh);
    typename MT::CCSF_double sp_heat(mesh);
    string analytic_sp_heat = "straight";

    // loop over procs and send out the Mat_States
    for (int np = 0; np < nodes(); np++)
    {
	// determine the number of cells on this processor
	int num_cells = cells_per_proc[np].size();

	// data for sending/receiving Mat_State
	double *density_send = new double[num_cells];
	double *T_send       = new double[num_cells];
	double *dedt_send    = new double[num_cells];

	// loop over on-proc cells and assign the data
	int global_cell;
	for (int cell = 1; cell <= num_cells; cell++)
	{
	    // get the global cell index
	    global_cell = cells_per_proc[np][cell-1];

	    // assign the data
	    density_send[cell-1] = mat.get_rho(global_cell);
	    T_send[cell-1]       = mat.get_T(global_cell);
	    dedt_send[cell-1]    = mat.get_dedt(global_cell);
	}
	
	// send the Mat_State data to the IMC nodes
	if (np)
	{
	    Send (num_cells, np, 37);
	    Send <double>(density_send, num_cells, np, 38);
	    Send <double>(T_send, num_cells, np, 39);
	    Send <double>(dedt_send, num_cells, np, 40);
	}
	else
	{
	    // build the host-processor CCSFs
	 
	    // check to make sure the local mesh is the right size
	    Check (mesh->num_cells() == num_cells);

	    // assign the data
	    for (int cell = 1; cell <= num_cells; cell++)
	    {
		density(cell) = density_send[cell-1];
		T(cell)       = T_send[cell-1];
		dedt(cell)    = dedt_send[cell-1];
		sp_heat(cell) = 0.0;
	    }
	}

	// delete dynamic allocation
	delete [] density_send;
	delete [] T_send;
	delete [] dedt_send;
    }

    // make and return the Mat_State to the host
    host_mat = new Mat_State<MT>(density, T, dedt, sp_heat,
				 analytic_sp_heat); 
    Ensure (host_mat->num_cells() == mesh->num_cells());
    return host_mat;
}

//---------------------------------------------------------------------------//
// receive the Mat_State object

template<class MT> SP<Mat_State<MT> > 
Parallel_Builder<MT>::recv_Mat(SP<MT> mesh)
{
    // insure that we are on an IMC-node and that we have a Mesh
    Require (node());
    Require (mesh);

    // declare the return Mat_State
    // specific heat is needed only to construct Mat_State
    SP<Mat_State<MT> > imc_mat;
    typename MT::CCSF_double density(mesh);
    typename MT::CCSF_double T(mesh);
    typename MT::CCSF_double dedt(mesh);
    typename MT::CCSF_double sp_heat(mesh);
    string analytic_sp_heat = "straight";

    // get the num_cells from the host
    int num_cells;
    Recv (num_cells, 0, 37);
    Check (num_cells == mesh->num_cells());
    
    // now receive the Mat_State data from the host
    double *density_recv = new double[num_cells];
    double *T_recv       = new double[num_cells];
    double *dedt_recv    = new double[num_cells];
    Recv (density_recv, num_cells, 0, 38);
    Recv (T_recv, num_cells, 0, 39);
    Recv (dedt_recv, num_cells, 0, 40);

    // assign the CCSFs
    for (int cell = 1; cell <= num_cells; cell++)
    {
	density(cell) = density_recv[cell-1];
	T(cell)       = T_recv[cell-1];
	dedt(cell)    = dedt_recv[cell-1];
	sp_heat(cell) = 0.0;
    }

    // release the dynamic storage
    delete [] density_recv;
    delete [] T_recv;
    delete [] dedt_recv;

    // build and return the Mat_state
    imc_mat = new Mat_State<MT>(density, T, dedt, sp_heat, analytic_sp_heat); 
    Ensure (imc_mat->num_cells() == num_cells);
    return imc_mat;
}

//---------------------------------------------------------------------------//
// Communicator passing interface
//---------------------------------------------------------------------------//
// build and send out a Communicator on each Processor

template<class MT>
template<class PT>
SP<Communicator<PT> > Parallel_Builder<MT>::send_Communicator()
{
    Require (!node());

    // first let's build a Communicator on each processor
    for (int np = 1; np < nodes(); np++)
    {
	// build a communicator on an IMC processor amd get its components 
	SP<Communicator<PT> > comm  = build_Communicator<PT>(np);
	vector<vector<int> > b_node = comm->get_b_node();
	vector<vector<int> > b_cell = comm->get_b_cell();
	vector<int> recv_nodes      = comm->get_recv_nodes();
	vector<int> send_nodes      = comm->get_send_nodes();

	// define c-style arrays for sending out
	Check (b_node.size() == b_cell.size());
	int bound_size = 0;
	for (int i = 0; i < b_node.size(); i++)
	{
	    bound_size += b_node[i].size();
	    Check (b_node[i].size() == b_cell[i].size());
	}
	int *b_node_send = new int[bound_size];
	int *b_cell_send = new int[bound_size];
	int *b_num       = new int[b_node.size()];
	int *nodes_send  = new int[recv_nodes.size() + send_nodes.size()];
	int *sizes       = new int[4];

	// assign sending arrays

	// assign nodes
	for (int i = 0; i < recv_nodes.size(); i++)
	    nodes_send[i] = recv_nodes[i];
	for (int i = 0; i < send_nodes.size(); i++)
	    nodes_send[i + recv_nodes.size()] = send_nodes[i];

	// assign boundary cell info
	int spacer = 0;
	for (int i = 0; i < b_node.size(); i++)
	{
	    b_num[i] = b_node[i].size();
	    for (int j = 0; j < b_node[i].size(); j++)
	    {
		b_node_send[spacer] = b_node[i][j];
		b_cell_send[spacer] = b_cell[i][j];
		spacer++;
	    }
	}
	Check (spacer == bound_size);

	// assign sizes of things
	sizes[0] = bound_size;
	sizes[1] = recv_nodes.size();
	sizes[2] = send_nodes.size();
	sizes[3] = b_node.size();

	// send the data out to the IMC processors
	Send <int>(sizes, 4, np, 42);
	Send <int>(b_node_send, bound_size, np, 43);
	Send <int>(b_cell_send, bound_size, np, 44);
	Send <int>(b_num, b_node.size(), np, 45);
	Send <int>(nodes_send, recv_nodes.size()+send_nodes.size(), np, 46);

	// reclaim storage
	delete [] b_node_send;
	delete [] b_cell_send;
	delete [] b_num;
	delete [] nodes_send;
	delete [] sizes;
    }

    // return the host communicator
    SP<Communicator<PT> > host_comm = build_Communicator<PT>(0);
    return host_comm;
}

//---------------------------------------------------------------------------//
// receive a Communicator on each processor

template<class MT>
template<class PT>
SP<Communicator<PT> > Parallel_Builder<MT>::recv_Communicator()
{
    Require (node());

    // receive the sizes of object data needed for the Communicator
    int *sizes = new int[4];
    Recv (sizes, 4, 0, 42);
    int bound_size = sizes[0];
    int recv_size  = sizes[1];
    int send_size  = sizes[2];
    int num_bcells = sizes[3];
    delete [] sizes;

    // make objects necessary for a communicator
    vector<vector<int> > b_node(num_bcells);
    vector<vector<int> > b_cell(num_bcells);
    vector<int> recv_nodes(recv_size);
    vector<int> send_nodes(send_size);
    int *b_node_recv = new int[bound_size];
    int *b_cell_recv = new int[bound_size];
    int *b_num       = new int[num_bcells];
    int *nodes_recv  = new int[recv_size + send_size];
    Recv (b_node_recv, bound_size, 0, 43);
    Recv (b_cell_recv, bound_size, 0, 44);
    Recv (b_num, num_bcells, 0, 45);
    Recv (nodes_recv, recv_size+send_size, 0, 46);

    // assign data to these objects

    // first do the receive and send nodes
    for (int i = 0; i < recv_size; i++)
	recv_nodes[i] = nodes_recv[i];
    for (int i = 0; i < send_size; i++)
	send_nodes[i] = nodes_recv[recv_size + i];

    // now do the boundary_nodes and cells
    int index = 0;
    for (int i = 0; i < num_bcells; i++)
    {
	b_node[i].resize(b_num[i]);
	b_cell[i].resize(b_num[i]);
	for (int j = 0; j < b_num[i]; j++)
	{
	    b_node[i][j] = b_node_recv[index];
	    b_cell[i][j] = b_cell_recv[index];
	    index++;
	}
    }
    Check (index == bound_size);

    // reclaim storage
    delete [] b_node_recv;
    delete [] b_cell_recv;
    delete [] b_num;
    delete [] nodes_recv;

    // make new communicator
    SP<Communicator<PT> > comm(new Communicator<PT>(recv_nodes, send_nodes, 
						    b_node, b_cell));
    // comm->print(std::cout);
    return comm;
}


//---------------------------------------------------------------------------//
// Communicator passing implementation
//---------------------------------------------------------------------------//
// build the communicator on processor np, we may have to change this when we 
// goto more general DD/rep modes

template<class MT>
template<class PT>
SP<Communicator<PT> > Parallel_Builder<MT>::build_Communicator(int np)
{
    // return value
    SP<Communicator<PT> > return_com;

    // get number of boundary cells on this processor
    int num_cells = bound_cells[np].size();
    vector<vector<int> > b_node(num_cells);
    vector<vector<int> > b_cell(num_cells);
    vector<int> com_nodes;
    vector<bool> procs(nodes(), false);

    // find the nodes this processor communicates with
    int global_cell;
    for (int i = 0; i < num_cells; i++)
    {
	// global cell index of the boundary cell
	global_cell = bound_cells[np][i];

	// the procs that this cell is on
	b_node[i] = procs_per_cell[global_cell-1];
	b_cell[i].resize(b_node[i].size());

	// loop over the nodes for this boundary cell and assign processors 
	// and local cell stuff
	int local_cell;
	int send_proc;
	for (int n = 0; n < b_node[i].size(); n++)
	{
	    send_proc        = b_node[i][n];
	    local_cell       = imc_cell(global_cell, send_proc);
	    procs[send_proc] = true;
	    b_cell[i][n]     = local_cell;
	}
    }

    // loop over the nodes on the problem and see if we communicate with it
    for (int n = 0; n < nodes(); n++)
	if (procs[n])
	    com_nodes.push_back(n);

    // return the communicator, at this point we assume one-to-one
    // communication
    return_com = new Communicator<PT>(com_nodes, com_nodes, b_node, b_cell);
    return return_com;
}

//---------------------------------------------------------------------------//
// diagnostics
//---------------------------------------------------------------------------//

template<class MT>
void Parallel_Builder<MT>::print(ostream &out) const
{
    out << endl;
    out << ">>> PARALLEL BUILD DATA <<<" << endl;
    out << "===========================" << endl;

    // check to see if we are on an IMC processor
    if (cells_per_proc.size() == 0 && procs_per_cell.size() == 0)
    {
	out << " ** You are on an IMC processor, no topology here!" << endl;
	return;
    }
}
    
} // end namespace rtt_imc

//---------------------------------------------------------------------------//
//                              end of CAR_CU_Parallel_Builder.t.hh
//---------------------------------------------------------------------------//
