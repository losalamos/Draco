//----------------------------------*-C++-*----------------------------------//
// Particle_Buffer.cc
// Thomas M. Evans
// Tue May 12 14:34:34 1998
//---------------------------------------------------------------------------//
// @> Particle_Buffer class implementation file
//---------------------------------------------------------------------------//

#include "imc/Particle_Buffer.hh"
#include <cstdlib>

IMCSPACE

//---------------------------------------------------------------------------//
// constructors
//---------------------------------------------------------------------------//
// templated constructor for Particle_Buffer

template<class PT>
template<class MT>
Particle_Buffer<PT>::Particle_Buffer(const MT &mesh, const Rnd_Control &rcon)
{
    
  // determine size of double info from Particles;
  // 5 = omega(3) + ew + fraction
    dsize = mesh.get_Coord().get_dim() + 5;

  // determine size of integer info from Particles;
  // 2 = cell + streamnum
    isize = 2;

  // determine size of character (RN state) from Particles;
  // in bytes
    csize = rcon.get_size();
    Check (csize <= RNG::max_buffer);
}

//---------------------------------------------------------------------------//
// constructors for Particle_Buffer<PT>::Census Buffer struct

template<class PT>
Particle_Buffer<PT>::Census_Buffer::Census_Buffer(vector<double> &r_,
						  vector<double> &omega_,
						  double ew_,
						  double fraction_,
						  int cell_,
						  Sprng random_)
    : r(r_), omega(omega_), ew(ew_), fraction(fraction_), cell(cell_),
      random(random_)
{
  // constructor for abbreviated particle data that comes back from census
  // files
}

template<class PT>
Particle_Buffer<PT>::Census_Buffer::Census_Buffer()
    : random(0, 0)
{
  // constructor for use with STL, this cannot be used
    Insist (0, "You tried to default construct a Census_Buffer!");
}

//---------------------------------------------------------------------------//
// IO FUNCTIONS
//---------------------------------------------------------------------------//
// write a single particle to an output

template<class PT>
void Particle_Buffer<PT>::write_census(ostream &cenfile,
				       const PT &particle) const
{
  // dynamicaly assign arrays for output, types double, int
    double *ddata = new double[dsize];
    int *idata = new int[isize];

  // assign all data of type double about the particle
    int index = 0;
    ddata[index++] = particle.ew;
    ddata[index++] = particle.fraction;
    for (int i = 0; i < particle.omega.size(); i++)
	ddata[index++] = particle.omega[i];
    for (int i = 0; i < particle.r.size(); i++)
	ddata[index++] = particle.r[i];
    Check (index == dsize);

  // assign integer data
    idata[0] = particle.cell;
    idata[1] = particle.random.get_num();

  // set the size of dynamic storage for the Random number state and pack it
    char *rdata;
    int size = pack_sprng(particle.random.get_id(), &rdata);
    Check (size == csize);

  // now dump particle data to the census file

  // make sure census file exists
    Check (cenfile);

  // write the output
    cenfile.write(reinterpret_cast<const char *>(ddata), dsize *
		  sizeof(double));
    cenfile.write(reinterpret_cast<const char *>(idata), isize *
		  sizeof(int)); 
    cenfile.write(reinterpret_cast<const char *>(rdata), csize);

  // reclaim dynamic memory
    delete [] ddata;
    delete [] idata;
    std::free(rdata);
}

//---------------------------------------------------------------------------//
// write a full Comm_Buffer to an output

template<class PT>
void Particle_Buffer<PT>::write_census(ostream &cenfile,
				       Comm_Buffer &buffer) const
{
  // check for output file
    if (!cenfile)
	Insist(0, 
	       "You tried to write census particles to a non-existent file!");

  // determine number of particles in the buffer
    int num_particles = buffer.n_part;
    Check (num_particles > 0);

  // make dynamically allocatable arrays
    double *ddata = new double[dsize];
    int    *idata = new int[isize];
    char   *cdata = new char[csize];

  // define indices for data
    int id = 0;
    int ii = 0;
    int ic = 0;

  // loop through comm_buffer and write the particle data to a census file
    for (int n = 1; n <= num_particles; n++)
    {
      // get the double array data
	for (int d = 0; d < dsize; d++)
	    ddata[d] = buffer.array_d[d+id];

      // get the int array data
	for (int i = 0; i < isize; i++)
	    idata[i] = buffer.array_i[i+ii];

      // get the char array data
	for (int c = 0; c < csize; c++)
	    cdata[c] = buffer.array_c[c+ic];

      // write the particle data
	cenfile.write(reinterpret_cast<const char *>(ddata), dsize *
		      sizeof(double));
	cenfile.write(reinterpret_cast<const char *>(idata), isize *
		      sizeof(int));
	cenfile.write(reinterpret_cast<const char *>(cdata), csize);

      // update the counters
	id += dsize;
	ii += isize;
	ic += csize;
	buffer.n_part--;

      // asserts to make sure we haven't gone over
	Check (id <= Global::buffer_d);
	Check (ii <= Global::buffer_i);
	Check (ic <= Global::buffer_c);
    }

  // recover storage
    delete [] ddata;
    delete [] idata;
    delete [] cdata;

  // reset the Comm_Buffer to 0
    Ensure (buffer.n_part == 0);
}

//---------------------------------------------------------------------------//
// read a single particle from an output and return a Census_Buffer

template<class PT>
SP<typename Particle_Buffer<PT>::Census_Buffer> 
Particle_Buffer<PT>::read_census(istream &cenfile) const
{
  // make sure file exists
    Check (cenfile);

  // cast smart pointer to Census_Buffer
    SP<Census_Buffer> return_part;
    
  // set pointers for dynamic memory storage
    double *ddata = new double[dsize];
    int    *idata = new int[isize];
    char   *rdata = new char[csize];

  // read in data
    cenfile.read(reinterpret_cast<char *>(ddata), dsize * sizeof(double));
    if (!cenfile.eof())
    {
      // read in integer data
	cenfile.read(reinterpret_cast<char *>(idata), isize * sizeof(int));
	Check (!cenfile.eof());

      // read in random number state
	cenfile.read(reinterpret_cast<char *>(rdata), csize);
	Check (!cenfile.eof());

      // assign data to proper structures for Census Particle
	double ew   = ddata[0];
	double frac = ddata[1];
	vector<double> r;
	vector<double> omega;
	for (int i = 2; i < 5; i++)
	    omega.push_back(ddata[i]);
	for (int i = 5; i < dsize; i++)
	    r.push_back(ddata[i]);
	int cell = idata[0];

      // make new random number
	int *id = unpack_sprng(rdata);
	Sprng random(id, idata[1]);

      // make new Census_Buffer
	return_part = new Census_Buffer(r, omega, ew, frac, cell, random);
    }

  // reclaim dynamic memory
    delete [] idata;
    delete [] ddata;
    delete [] rdata;

  // return Census_Buffer
    return return_part;
}

//---------------------------------------------------------------------------//
// Particle send and receives
//---------------------------------------------------------------------------//
// do a block send of a Comm_Buffer

template<class PT>
void Particle_Buffer<PT>::send_buffer(Comm_Buffer &buffer, int proc) const
{
  // <<CONTINUE HERE>>
}

//---------------------------------------------------------------------------//
// do a block recv of a Comm_Buffer

template<class PT>
SP<typename Particle_Buffer<PT>::Comm_Buffer> 
Particle_Buffer<PT>::recv_buffer(int proc) const
{
  // return Comm_Buffer declaration
    SP<Comm_Buffer> return_buffer;

  // <<CONTINUE HERE>>

  // return SP
    return return_buffer;
}

//---------------------------------------------------------------------------//
// Do an asyncronous send using C4

template<class PT>
void Particle_Buffer<PT>::asend_bank(Comm_Buffer &buffer, int proc, 
				     Comm_Bank &bank) const
{
  // find out the number of Particles
    int num_part = bank.size();
    Check (num_part > 0 && num_part <= Global::buffer_s);
    
  // define indices for data
    int id = 0;
    int ii = 0;
    int ic = 0;

  // define arrays for buffering
    double array_d[Global::buffer_d];
    int    array_i[Global::buffer_i];
    char   array_c[Global::buffer_c];

  // loop through particles and get the goods
    while (bank.size())
    {
      // get the double info from the particle
	array_d[id++] = bank.top().ew;
	array_d[id++] = bank.top().fraction;
	array_d[id++] = bank.top().time_left;
	for (int j = 0; j < bank.top().omega.size(); j++)
	    array_d[id++] = bank.top().omega[j];
	for (int j = 0; j < bank.top().r.size(); j++)
	    array_d[id++] = bank.top().r[j];
	Check (id <= Global::buffer_d);

      // get the int info from the particle
	array_i[ii++] = bank.top().cell;
	array_i[ii++] = bank.top().random.get_num();
	Check (ii <= Global::buffer_i);

      // get the char info from the particle
	char *bytes;
	int size = pack_sprng(bank.top().random.get_id(), &bytes);
	Check (size == csize);
	for (int j = 0; j < csize; j++)
	    array_c[ic++] = bytes[j];
	std::free(bytes);
	Check (ic <= Global::buffer_c);

      // pop the highest level particle
	bank.pop();
    }

  // send the particle buffers
    SendAsync(buffer.comm_n, &num_part, 1, proc, 100);
    SendAsync(buffer.comm_d, &array_d[0], Global::buffer_d, proc, 101);
    SendAsync(buffer.comm_i, &array_i[0], Global::buffer_i, proc, 102);
    SendAsync(buffer.comm_c, &array_c[0], Global::buffer_c, proc, 103);
}
    
//---------------------------------------------------------------------------//
// post async recives

template<class PT>
void Particle_Buffer<PT>::arecv_bank(Comm_Buffer &buf, int proc) const
{
  // post c4 async receives
    RecvAsync(buf.comm_n, &(buf.n_part), 1, proc, 100);
    RecvAsync(buf.comm_d, &(buf.array_d[0]), Global::buffer_d, proc, 101);
    RecvAsync(buf.comm_i, &(buf.array_i[0]), Global::buffer_i, proc, 102);
    RecvAsync(buf.comm_c, &(buf.array_c[0]), Global::buffer_c, proc, 103);
}

//---------------------------------------------------------------------------//
// use the filled buffer to create a new bank of particles

template<class PT> void 
Particle_Buffer<PT>::add_to_bank(Comm_Buffer &buffer, Comm_Bank &bank) const
{
  // wait on recieve buffers to make sure they are full
    buffer.comm_n.wait();
    buffer.comm_d.wait();
    buffer.comm_i.wait();
    buffer.comm_c.wait();

  // make the new particle bank from the buffer
  // ...
}

//---------------------------------------------------------------------------//
// Buffering functions
//---------------------------------------------------------------------------//
// add a Census_Buffer to a Comm_Buffer

template<class PT>
void Particle_Buffer<PT>::buffer_census(Comm_Buffer &comm, 
					const Census_Buffer &census) const
{
  // <<CONTINUE HERE>>
}

//---------------------------------------------------------------------------//
// add a Particle to a Comm_Buffer

template<class PT>
void Particle_Buffer<PT>::buffer_particle(Comm_Buffer &comm,
					  const PT &particle) const
{
  // <<CONTINUE HERE>>
}

CSPACE

//---------------------------------------------------------------------------//
//                              end of Particle_Buffer.cc
//---------------------------------------------------------------------------//
