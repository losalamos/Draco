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

  // number of cells in the mesh
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
	std::cout << "** Doing full replication" << std::endl;
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
	std::cout << "** Doing full DD" << std::endl;
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

	std::cout << "** Doing general DD/rep" << std::endl;
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
Parallel_Builder<MT>::send_Source(SP<MT> mesh, const Source_Init<MT> &sinit, 
				  const Particle_Buffer<PT> &buffer,
				  SP<Rnd_Control> rcon)
{
  // check that we are on the host node only
    Require (!node());
    Require (mesh);

  // data necessary to build Source on host
    SP<Source<MT> > host_source;
    vector<vector<int> > vol(mesh->num_cells());
    vector<vector<int> > ss(mesh->num_cells());

  // first distribute the census
    if (sinit.get_ncentot() > 0) 
	dist_census(sinit, buffer);

  // next do the volume source
    if (sinit.get_nvoltot() > 0)
	vol = dist_vol(sinit);

  // finally do the surface source
    if (sinit.get_nsstot() > 0)
	ss = dist_ss(sinit);

  // now let's build the source on the host
    typename MT::CCSF_int volrn(mesh, vol[1]);
    typename MT::CCSF_int nvol(mesh, vol[0]);
    typename MT::CCSF_int ssrn(mesh, ss[1]);
    typename MT::CCSF_int nss(mesh, ss[0]);

  // get the number of volume, census, and surface sources for this processor
    int ncentot = 0;
    int nvoltot = 0;
    int nsstot = 0;
    for (int i = 1; i <= mesh->num_cells(); i++)
    {
	nvoltot += nvol(i);
	nsstot  += nss(i);
    }

  // make the source
    host_source = new Source<MT>(volrn, nvol, ssrn, nss, "census.0", 
				 nvoltot, nsstot, ncentot, rcon, buffer);

  // return source
    return host_source;
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
	    if (num_to_send[proc_goto] == Global::buffer_s)
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
// distribute the volume source stuff
	    
template<class MT> vector<vector<int> > 
Parallel_Builder<MT>::dist_vol(const Source_Init<MT> &sinit)
{
  // return data for the host source
    vector<vector<int> > host_vol(2);

  // find the number of cells on the global mesh
    int num_cells = procs_per_cell.size();
    Require (num_cells > 0);

  // make the volume source counters
    vector<int> nvol(num_cells);
    vector<int> nvol_xtra(num_cells);
    vector<int> streamnum(num_cells);
    int counter = Global::rn_stream;
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
    Global::rn_stream += sinit.get_nvoltot();
    Check (counter == Global::rn_stream);

  // loop over processors and make the send info
    int voltot = 0;
    for (int i = 0; i < nodes(); i++)
    {
      // calculate the number of cells on proc i mesh
	int num_cells = cells_per_proc[i].size();

      // make a c-style array for passing to other nodes
	int *nvol_send   = new int[num_cells];
	int *stream_send = new int[num_cells];
	
      // loop through cells on this processor and send stuff out
	for (int cell = 1; cell <= num_cells; cell++)
	{
	  // determine the global cell index
	    int global_cell = cells_per_proc[i][cell-1];

	  // check to see if we need to add extra sources to this cell on
	  // this processor
	    if (nvol_xtra[global_cell-1] <= i+1)
		nvol[global_cell-1]++;
	    
	  // calculate the nvol source and stream number for proc i
	    nvol_send[cell-1]   = nvol[global_cell-1];
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
	}
	else
	{
	  // return the needed source stuff for the host
	  
	  // nvol and stream num data
	    host_vol[0].resize(num_cells);
	    host_vol[1].resize(num_cells);

	  // assign the vectors to our point array
	    for (int n = 0; n < num_cells; n++)
	    {
		host_vol[0][n] = nvol_send[n];
		host_vol[1][n] = stream_send[n];
	    }
	}

      // reclaim dynamic memory
	delete [] nvol_send;
	delete [] stream_send;
    }

  // final assertion
    std::cout << voltot << " " << sinit.get_nvoltot() << std::endl;
    Ensure (voltot == sinit.get_nvoltot());

  // return the host volume source to the host processor
    return host_vol;
}

//---------------------------------------------------------------------------//
// distribute the surface source stuff
	    
template<class MT> vector<vector<int> > 
Parallel_Builder<MT>::dist_ss(const Source_Init<MT> &sinit)
{
  // return data for the host source
    vector<vector<int> > host_ss(2);

  // find the number of cells on the global mesh
    int num_cells = procs_per_cell.size();
    Require (num_cells > 0);

  // make the surface source counters
    vector<int> nss(num_cells);
    vector<int> nss_xtra(num_cells);
    vector<int> streamnum(num_cells);
    int counter = Global::rn_stream;
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
    Global::rn_stream += sinit.get_nsstot();
    Check (counter == Global::rn_stream);

  // loop over processors and make the send info
    int sstot = 0;
    for (int i = 0; i < nodes(); i++)
    {
      // calculate the number of cells on proc i mesh
	int num_cells = cells_per_proc[i].size();

      // make a c-style array for passing to other nodes
	int *nss_send   = new int[num_cells];
	int *stream_send = new int[num_cells];
	
      // loop through cells on this processor and send stuff out
	for (int cell = 1; cell <= num_cells; cell++)
	{
	  // determine the global cell index
	    int global_cell = cells_per_proc[i][cell-1];

	  // check to see if we need to add extra sources to this cell on
	  // this processor
	    if (nss_xtra[global_cell-1] <= i)
		nss[global_cell-1]++;
	    
	  // calculate the nss source and stream number for proc i
	    nss_send[cell-1]    = nss[global_cell-1];
	    stream_send[cell-1] = streamnum[global_cell-1];
	    streamnum[global_cell-1] += nss_send[cell-1];

	  // do some check counters
	    sstot += nss_send[cell-1];
	}

      // send out the data 
	if (node())
	{
	  // send to proc i if we are not on the host
	    Send (num_cells, i, 27);
	    Send (nss_send, num_cells, i, 28);
	    Send (stream_send, num_cells, i, 29);
	}
	else
	{
	  // return the needed source stuff for the host
	  
	  // nss and stream num data
	    host_ss[0].resize(num_cells);
	    host_ss[1].resize(num_cells);

	  // assign the vectors to our point array
	    for (int n = 0; n < num_cells; n++)
	    {
		host_ss[0][n] = nss_send[n];
		host_ss[1][n] = stream_send[n];
	    }
	}

      // reclaim dynamic memory
	delete [] nss_send;
	delete [] stream_send;
    }

  // final assertion
    Ensure (sstot == sinit.get_nsstot());

  // return the host surface source to the host processor
    return host_ss;
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
	sigma[cell-1]  = opacity.get_sigma(cell);
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
    return_opacity = new Opacity<MT>(sigma);
    return return_opacity;
}
    
CSPACE

//---------------------------------------------------------------------------//
//                              end of Parallel_Builder.cc
//---------------------------------------------------------------------------//
