//----------------------------------*-C++-*----------------------------------//
// Parallel_Builder.cc
// Thomas M. Evans
// Tue Apr 14 14:50:22 1998
//---------------------------------------------------------------------------//
// @> Parallel_Builder implementation file
//---------------------------------------------------------------------------//

#include "imc/Parallel_Builder.hh"
#include "imc/XYCoord_sys.hh"
#include "imc/XYZCoord_sys.hh"
#include "imc/Global.hh"
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cstdio>

IMCSPACE

// draco necessities
using C4::Send;
using C4::Recv;
using Global::max;
using Global::min;

// std necessities
using std::string;
using std::fill;
using std::ofstream;
using std::ifstream;
using std::remove;
using std::ostringstream;

//---------------------------------------------------------------------------//
// constructors
//---------------------------------------------------------------------------//
// host node constructor that determines the toplogy on the IMC-nodes

template<class MT>
Parallel_Builder<MT>::Parallel_Builder(const MT &mesh, 
				       const Source_Init<MT> &sinit)
{
  // calculate the parameters for splitting the problem amongst many
  // processors 
    parallel_topology(mesh, sinit);
}

//---------------------------------------------------------------------------//
// default constructor for IMC-nodes

template<class MT>
Parallel_Builder<MT>::Parallel_Builder()
{
  // the IMC-nodes don't use this data
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
	std::cout << " ** Doing full replication" << std::endl;
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
	std::cout << " ** Doing full DD" << std::endl;
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

	std::cout << " ** Doing general DD/rep" << std::endl;
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

  // send the numbers of each type of source to the other processors
    for (int i = 1; i < nodes(); i++)
    {
	Send (sinit.get_ncentot(), i, 34);
	Send (sinit.get_nvoltot(), i, 35);
	Send (sinit.get_nsstot(), i, 36);
    }

  // first distribute the census
    if (sinit.get_ncentot() > 0) 
	dist_census(sinit, buffer);

  // next do the volume source
    if (sinit.get_nvoltot() > 0)
	dist_vol(sinit, volrn, nvol, ew_vol, t4_slope);

  // finally do the surface source
    if (sinit.get_nsstot() > 0)
	dist_ss(sinit, ssrn, nss, fss, ew_ss);

  // get the number of volume, census, and surface sources for this processor
    int ncentot = 0;
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
				 fss, ew_ss, ss_dist, "census.0", nvoltot, 
				 nsstot, ncentot, rcon, buffer, mat); 

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

  // receive the census
    if (global_ncentot)
	recv_census(buffer);

  // receive the volume source
    if (global_nvoltot)
	recv_vol(volrn, nvol, ew_vol, t4_slope);

  // receive the surface source
    if (global_nsstot)
	recv_ss(ssrn, nss, fss, ew_ss);

  // get the number of volume, census, and surface sources for this processor
    int ncentot = 0;
    int nvoltot = 0;
    int nsstot = 0;
    for (int i = 1; i <= mesh->num_cells(); i++)
    {
	nvoltot += nvol(i);
	nsstot  += nss(i);
    }

  // get the census name
    ostringstream cenfile;
    cenfile << "census." << node();

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
				fss, ew_ss, ss_dist, cenfile.str(), nvoltot, 
				nsstot, ncentot, rcon, buffer, mat);

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
				       const Particle_Buffer<PT> &buffer)
{
  // get number of cells
    int num_cells = procs_per_cell.size();
    Require (num_cells > 0);

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

  // open host census file
    ifstream host_census("census");
    if (!host_census)
	Insist(0, "The host did not provide an initial census file!");

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
    SP<Particle_Buffer<PT>::Census_Buffer> cenpart;
    vector<Particle_Buffer<PT>::Comm_Buffer> cen_buffer(nodes());

  // loop to read census particles and send them
    do
    {
	cenpart = buffer.read_census(host_census);
	if (cenpart)
	{
	    total_read++;
	    cell = cenpart->cell;

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
	    buffer.buffer_census(cen_buffer[proc_goto], *cenpart);

	  // increment some counters
	    num_placed_per_cell[cell-1]++;
	    num_placed_per_proc[proc_goto]++;
	    total_placed++;
	    num_to_send[proc_goto]++;

	  // send these guys out
	    if (num_to_send[proc_goto] ==
		Particle_Buffer<PT>::get_buffer_s())
	    {
		buffer.send_buffer(cen_buffer[proc_goto], proc_goto);
		cen_buffer[proc_goto].n_part = 0;
		num_to_send[proc_goto] = 0;
	    }
	}
    } while (cenpart);
	    
  // some Checks
    Check (total_placed == sinit.get_ncentot());
    Check (total_read == sinit.get_ncentot());

  // send the guys off that haven't been sent yet
    for (int proc = 0; proc < nodes(); proc++)
    {
	if (num_to_send[proc] > 0)
	{
	    buffer.send_buffer(cen_buffer[proc], proc_goto);
	    cen_buffer[proc].n_part = 0;
	    num_to_send[proc] = 0;  
	}
    }
}

//---------------------------------------------------------------------------//
// receive the census buffers and write out the census files

template<class MT>
template<class PT>
void Parallel_Builder<MT>::recv_census(const Particle_Buffer<PT> &buffer)
{
  // do something here
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
    int counter = RNG::rn_stream;
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
    RNG::rn_stream += sinit.get_nvoltot();
    Check (counter == RNG::rn_stream);

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
	    Send (nvol_send, num_cells, i, 25);
	    Send (stream_send, num_cells, i, 26);
	    Send (ew_send, num_cells, i, 27);
	    Send (t4_send, num_cells * dimension, i, 28);
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
    int counter = RNG::rn_stream;
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
    RNG::rn_stream += sinit.get_nsstot();
    Check (counter == RNG::rn_stream);

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
	    Send (nss_send, num_cells, i, 30);
	    Send (stream_send, num_cells, i, 31);
	    Send (fss_send, num_cells, i, 32);
	    Send (ew_send, num_cells, i, 33);
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
// send out the Mesh components

template<class MT>
void Parallel_Builder<MT>::send_Mesh(const MT &mesh)
{
  // send the mesh components out to the other processors

  // assure that we are on the host node only
    Check (!node());

  // let us pass a coordinate system, shall we
    send_Coord(mesh.get_Coord());

  // let us pass the Layout
    send_Layout(mesh.get_Layout());

  // let us pass the vertex
    send_vertex(mesh.get_vertex());

  // let us pass the cell_pair
    send_cellpair(mesh.get_cell_pair());
}

//---------------------------------------------------------------------------//
// receive the Mesh components

template<class MT>
SP<MT> Parallel_Builder<MT>::recv_Mesh()
{
  // receive the Mesh components from the host and build the new mesh on this 
  // topology

  // assure that we are not on the host node
    Check (node());

  // rebuilt Mesh SP
    SP<MT> return_mesh;

  // get coordinate system
    SP<Coord_sys> coord = recv_Coord();

  // get Layout
    Layout layout = recv_Layout();

  // get vertices and cell_pair
    typename MT::CCVF_a vertex = recv_vertex();
    typename MT::CCVF_i cell_pair = recv_cellpair();

  // build mesh
  // SOLARIS MPICH fails at this point
    return_mesh = new MT(coord, layout, vertex, cell_pair);

  // return mesh
    return return_mesh;
}

//---------------------------------------------------------------------------//
// Mesh passing implementation
//---------------------------------------------------------------------------//
// pass the Coord_sys object

template<class MT>
void Parallel_Builder<MT>::send_Coord(const Coord_sys &coord)
{
  // send out the coordinate system designator

  // send variables
    string cs          = coord.get_Coord();
    const char *sendcs = cs.c_str();
    int cs_size        = cs.size();

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
	SP<XYCoord_sys> xycoord = new XYCoord_sys;
	coord = xycoord;
    }
    else if (cs == "xyz")
    {
	SP<XYZCoord_sys> xyzcoord = new XYZCoord_sys;
	coord = xyzcoord;
    }
    
  // return base class SP to a derived Coord_sys
    return coord;
}

//---------------------------------------------------------------------------//
// pass the Layout

template<class MT>
void Parallel_Builder<MT>::send_Layout(const Layout &layout)
{
  // set the Layout size
    int num_cells = layout.num_cells();

  // calculate the number of faces on each cell and total size of the Layout
  // (layout_size = SUM_(i)^(num_cells) (num_faces[i]))
    int *num_faces = new int[num_cells];
    int layout_size = 0;
    for (int i = 1; i <= num_cells; i++)
    {
	num_faces[i-1] = layout.num_faces(i);
	layout_size += num_faces[i-1];
    }

  // write the faces array for passing
    int *faces = new int[layout_size];
    int index = 0;
    for (int i = 1; i <= num_cells; i++)
	for (int j = 1; j <= layout.num_faces(i); j++)
	{
	    faces[index] = layout(i,j);
	    index++;
	}
    
  // pass all this good stuff to the other processors (for now)
    for (int np = 1; np < nodes(); np++)
    {
      // pass the Layout cell size
	Send (num_cells, np, 3);

      // pass the Layout total size
	Send (layout_size, np, 4);

      // pass the num_faces array
	Send (&num_faces[0], num_cells, np, 5);

      // pass the face-values array
	Send (faces, layout_size, np, 6);
    }

  // delete dynamically allocated faces array
    delete [] faces;
    delete [] num_faces;
}

//---------------------------------------------------------------------------//
// receive the Layout

template<class MT>
Layout Parallel_Builder<MT>::recv_Layout()
{
  // receive and build the Layout

  // get the number of cells and size of the Layout
    int num_cells;
    int layout_size;
    Recv (num_cells, 0, 3);
    Recv (layout_size, 0, 4);

  // make the Layout
    Layout layout = num_cells;

  // receive the Layout data arrays
    int *num_faces = new int[num_cells];
    int *faces = new int[layout_size];
    Recv (num_faces, num_cells, 0, 5);
    Recv (faces, layout_size, 0, 6);

  // rebuild the Layout
    int index = 0;
    for (int cell = 1; cell <= num_cells; cell++)
    {
	layout.set_size(cell, num_faces[cell-1]);
	for (int i = 1; i <= num_faces[cell-1]; i++)
	{
	    layout(cell, i) = faces[index];
	    index++;
	}
    }

  // get rid of dynamic arrays
    delete [] num_faces;
    delete [] faces;
    
  // return the new Layout on this node
    return layout;
}

//---------------------------------------------------------------------------//
// pass the mesh vertices

template<class MT>
void Parallel_Builder<MT>::send_vertex(const typename MT::CCVF_a &vertex)
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

  // send out the goodies
    for (int np = 1; np < nodes(); np++)
    {
      // send the dimensionality of the vertex-array
	Send (vertex_dim, np, 7);

      // send the total size of the vertex-array
	Send (total_size, np, 8);

      // send the vertex array
	Send (vert, total_size, np, 9);
    }

  // delete dynamically allocated arrays
    delete [] vert;
}
	   
//---------------------------------------------------------------------------//
// receive the vertex

template<class MT>
typename MT::CCVF_a Parallel_Builder<MT>::recv_vertex()
{
  // receive and rebuild the vertices

  // first get the sizes
    int vertex_dim;
    int total_size;
    Recv (vertex_dim, 0, 7);
    Recv (total_size, 0, 8);

  // find the size of each dimensional vertex array
    Check (!(total_size % vertex_dim));
    int vertex_size = total_size / vertex_dim;

  // make a new CCVF_a vertex array
    typename MT::CCVF_a vertex(vertex_dim);

  // get the vertices from the host
    double *vert = new double[total_size];
    Recv (vert, total_size, 0, 9);

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

  // delete dynamic allocated arrays
    delete [] vert;

  // return the new vertex guy
    return vertex;
}

//---------------------------------------------------------------------------//
// pass the cell_pair

template<class MT>
void Parallel_Builder<MT>::send_cellpair(const typename MT::CCVF_i
					 &cell_pair)
{
  // send the cell_pair array to another processor
    
  // set the cell_pair size
    int num_cells = cell_pair.size();

  // calculate the number of vertices per cell and the total size of the
  // cell_pair object
    int *num_vert = new int[num_cells];
    int size = 0;
    for (int i = 0; i < num_cells; i++)
    {
	num_vert[i] = cell_pair[i].size();
	size += num_vert[i];
    }

  // write the cell_pair array for passing
    int *vertices = new int[size];
    int index = 0;
    for (int i = 0; i < num_cells; i++)
	for (int j = 0; j < cell_pair[i].size(); j++)
	{
	    vertices[index] = cell_pair[i][j];
	    index++;
	}

  // pass the important stuff to other processors
    for (int np = 1; np < nodes(); np++)
    {
      // pass the size
	Send (num_cells, np, 10);

      // pass the total size of the cell_pair object
	Send (size, np, 11);

      // pass the num_vert array
	Send (num_vert, num_cells, np, 12);
	
      // pass the vertices-values array
	Send (vertices, size, np, 13);
    }

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
    Recv (num_cells, 0, 10);
    Recv (size, 0, 11);

  // make a new cell_pair object
    typename MT::CCVF_i cell_pair(num_cells);

  // receive the cell_pair data
    int *num_vert = new int[num_cells];
    int *vertices = new int[size];
    Recv (num_vert, num_cells, 0, 12);
    Recv (vertices, size, 0, 13);

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
// Opacity passing interface
//---------------------------------------------------------------------------//
// send the Opacity object

template<class MT>
void Parallel_Builder<MT>::send_Opacity(const Opacity<MT> &opacity)
{
  // send out the Opacities, one component at a time

  // assure that we are on the host node
    Check (!node());

  // determine the number of cells
    int num_cells = opacity.num_cells();

  // assign the Opacity data
    double *sigma  = new double[num_cells];
    double *planck = new double[num_cells];
    double *fleck  = new double[num_cells];
    
    for (int cell = 1; cell <= num_cells; cell++)
    {
	sigma[cell-1]  = opacity.get_sigma_abs(cell);
	planck[cell-1] = opacity.get_planck(cell);
	fleck[cell-1]  = opacity.get_fleck(cell);
    }

  // send the Opacity data
    for (int np = 1; np < nodes(); np++)
    {
	Send (num_cells, np, 20);
	Send (sigma, num_cells, np, 21);
	Send (planck, num_cells, np, 22);
	Send (fleck, num_cells, np, 23);
    }

  // delete dynamic allocation
    delete [] sigma;
    delete [] planck;
    delete [] fleck;
}

//---------------------------------------------------------------------------//
// receive the Opacity object

template<class MT>
SP<Opacity<MT> > Parallel_Builder<MT>::recv_Opacity(SP<MT> mesh)
{
  // receive and rebuild the Opacity object

  // assure we are on receive nodes
    Check (node());

  // declare return opacity object
    SP< Opacity<MT> > return_opacity;

  // check to make sure we have a valid mesh pointer
    Check (mesh);

  // receive the size of this guy
    int num_cells;
    Recv (num_cells, 0, 20);

  // check to make sure our meshes are of proper size
    Check (num_cells == mesh->num_cells());

  // make new Opacity objects
    typename MT::CCSF_double sigma(mesh);
    typename MT::CCSF_double planck(mesh);
    typename MT::CCSF_double fleck(mesh);

  // receive data from host
    double *rsigma  = new double[num_cells];
    double *rplanck = new double[num_cells];
    double *rfleck  = new double[num_cells];
    Recv (rsigma, num_cells, 0, 21);
    Recv (rplanck, num_cells, 0, 22);
    Recv (rfleck, num_cells, 0, 23);

  // assign to new Opacity objects
    for (int cell = 1; cell <= num_cells; cell++)
    {
	sigma(cell)  = rsigma[cell-1];
	planck(cell) = rplanck[cell-1];
	fleck(cell)  = rfleck[cell-1];
    }

  // reclaim dynamic memory
    delete [] rsigma;
    delete [] rplanck;
    delete [] rfleck;

  // build and return new opacity object
    return_opacity = new Opacity<MT>(sigma, planck, fleck);
    return return_opacity;
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
    SP<Mat_State<MT> > host_mat;
    typename MT::CCSF_double density(mesh);
    typename MT::CCSF_double T(mesh);
    typename MT::CCSF_double Cv(mesh);

  // loop over procs and send out the Mat_States
    for (int np = 0; np < nodes(); np++)
    {
      // determine the number of cells on this processor
	int num_cells = cells_per_proc[np].size();

      // data for sending/receiving Mat_State
	double *density_send = new double[num_cells];
	double *T_send       = new double[num_cells];
	double *Cv_send      = new double[num_cells];

      // loop over on-proc cells and assign the data	
	for (int cell = 1; cell <= num_cells; cell++)
	{
	  // get the global cell index
	    int global_cell = cells_per_proc[np][cell-1];

	  // assign the data
	    density_send[cell-1] = mat.get_rho(global_cell);
	    T_send[cell-1]       = mat.get_T(global_cell);
	    Cv_send[cell-1]      = mat.get_Cv(global_cell);
	}
	
      // send the Opacity data to the IMC nodes
	if (np)
	{
	    Send (num_cells, np, 37);
	    Send (density_send, num_cells, np, 38);
	    Send (T_send, num_cells, np, 39);
	    Send (Cv_send, num_cells, np, 40);
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
		Cv(cell)      = Cv_send[cell-1];
	    }
	}

      // delete dynamic allocation
	delete [] density_send;
	delete [] T_send;
	delete [] Cv_send;
    }

  // make and return the Mat_State to the host
    host_mat = new Mat_State<MT>(density, T, Cv);
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
    SP<Mat_State<MT> > imc_mat;
    typename MT::CCSF_double density(mesh);
    typename MT::CCSF_double T(mesh);
    typename MT::CCSF_double Cv(mesh);

  // get the num_cells from the host
    int num_cells;
    Recv (num_cells, 0, 37);
    Check (num_cells == mesh->num_cells());
    
  // now receive the Mat_State data from the host
    double *density_recv = new double[num_cells];
    double *T_recv       = new double[num_cells];
    double *Cv_recv      = new double[num_cells];
    Recv (density_recv, num_cells, 0, 38);
    Recv (T_recv, num_cells, 0, 39);
    Recv (Cv_recv, num_cells, 0, 40);

  // assign the CCSFs
    for (int cell = 1; cell <= num_cells; cell++)
    {
	density(cell) = density_recv[cell-1];
	T(cell)       = T_recv[cell-1];
	Cv(cell)      = Cv_recv[cell-1];
    }

  // release the dynamic storage
    delete [] density_recv;
    delete [] T_recv;
    delete [] Cv_recv;

  // build and return the Mat_state
    imc_mat = new Mat_State<MT>(density, T, Cv);
    Ensure (imc_mat->num_cells() == num_cells);
    return imc_mat;
}

CSPACE

//---------------------------------------------------------------------------//
//                              end of Parallel_Builder.cc
//---------------------------------------------------------------------------//
