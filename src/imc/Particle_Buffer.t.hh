//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Particle_Buffer.t.hh
 * \author Thomas M. Evans
 * \date   Tue May 12 14:34:34 1998
 * \brief  Particle_Buffer member function definitions.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Particle_Buffer.hh"
#include <cstdlib>
#include <algorithm>
#include <numeric>

namespace rtt_imc 
{

//---------------------------------------------------------------------------//
// CONSTRUCTORS
//---------------------------------------------------------------------------//
// templated constructor for Particle_Buffer

template<class PT>
template<class MT>
Particle_Buffer<PT>::Particle_Buffer(const MT &mesh, 
				     const Rnd_Type_Control &rcon) 
{
    // Note: dsize is the size of the double data for the census. This
    // differs from the complete set of particle data by one double value
    // (time_left). Hence, we set dsize to _one_less_than_ the value returned
    // by Particle::Pack. Since, the communication buffers include all of the
    // double data in the particle, dsize+1 is passed to set_buffer

    // Call setup for Particle::Pack
    PT_Pack::setup_buffer_sizes(mesh, rcon);

    // Get sizes from Particle::Pack
    dsize = PT_Pack::get_double_size()-1;  // -1 : See note above
    isize = PT_Pack::get_int_size();
    csize = PT_Pack::get_char_size();

    Check(dsize); Check(isize); Check(csize);

    // set the static buffer variables 
    set_buffer(dsize+1, isize, csize);      // +1 : See note above
}

//---------------------------------------------------------------------------//
// Constructor that gets the random number state size from Rnd_Control; but
// allows the user to enter the dimensionality instead of relying on a fully
// built mesh.

template<class PT>
Particle_Buffer<PT>::Particle_Buffer(int dimension, 
				     const Rnd_Type_Control &rcon)
{
    // Note: See comments for Particle_Buffer(const MT&, const Rnd_Control&)
    // above

    // Call setup for Particle::Pack
    PT_Pack::setup_buffer_sizes(dimension, rcon);

    // Get sizes from Particle::Pack
    dsize = PT_Pack::get_double_size()-1;
    isize = PT_Pack::get_int_size();
    csize = PT_Pack::get_char_size();

    // set the static buffer variables
    set_buffer(dsize+1, isize, csize);  
}

//---------------------------------------------------------------------------//
// constructor for Particle_Buffer which allows the user to enter the 
// Particle sizes

// This constructor is dangerous because it does not sync Particle::Pack
// buffer sizes

template<class MT>
Particle_Buffer<MT>::Particle_Buffer(int d, int i, int c)
    : dsize(d), isize(i), csize(c) 
{
    // set the static buffer variables
    set_buffer(dsize+1, isize, csize);
}

//---------------------------------------------------------------------------//
// constructors for Particle_Buffer<PT>::Census Buffer struct

template<class PT>
Particle_Buffer<PT>::Census_Buffer::Census_Buffer(sf_double &r_,
						  sf_double &omega_,
						  double ew_,
						  double fraction_,
						  int cell_,
						  Rnd_Type random_)
    : r(r_), 
      omega(omega_), 
      ew(ew_), 
      fraction(fraction_), 
      cell(cell_),
      random(random_)
{
    // constructor for abbreviated particle data that comes back from census
    // files
}

//---------------------------------------------------------------------------//
// STATIC PRIVATE DATA MEMBERS
//---------------------------------------------------------------------------//

// buffer size
template<class PT> int Particle_Buffer<PT>::buffer_s = 1000;

// size of doubles, ints, and chars in buffer
template<class PT> 
int Particle_Buffer<PT>::buffer_d = buffer_s * PT_Pack::get_double_size();

template<class PT> 
int Particle_Buffer<PT>::buffer_i = buffer_s * PT_Pack::get_int_size();

template<class PT>
int Particle_Buffer<PT>::buffer_c = buffer_s * PT_Pack::get_char_size();

//---------------------------------------------------------------------------//
// BUFFER SIZE FUNCTIONS
//---------------------------------------------------------------------------//
// set the buffer size (public interface)

template<class PT>
void Particle_Buffer<PT>::set_buffer_size(int s)
{
    Require(s);

    // reset the buffer sizes
    buffer_d = s*PT_Pack::get_double_size();
    buffer_i = s*PT_Pack::get_int_size();
    buffer_c = s*PT_Pack::get_char_size();
    buffer_s = s;
}

//---------------------------------------------------------------------------//
// set the buffers (private implementation)

template<class PT>
void Particle_Buffer<PT>::set_buffer(int d, int i, int c)
{
    // This function is dangerous because it resets buffer sizes without
    // consulting Particle::Pack

    // reset the double, int, and char buffer sizes
    buffer_d = buffer_s * d;
    buffer_i = buffer_s * i;
    buffer_c = buffer_s * c;
}

template<class PT>
void Particle_Buffer<PT>::set_buffer(int d, int i, int c, int s)
{    
    // This function is dangerous because it resets buffer sizes without
    // consulting Particle::Pack

    // reset the buffer size, double, int, and char buffer sizes
    buffer_s = s;
    buffer_d = buffer_s * d;
    buffer_i = buffer_s * i;
    buffer_c = buffer_s * c;
}

//---------------------------------------------------------------------------//
// IO FUNCTIONS
//---------------------------------------------------------------------------//
// write a single particle to an output

template<class PT>
void Particle_Buffer<PT>::write_census(std_ostream &cenfile,
				       const PT &particle) const
{
    // create a packed version of the particle
    PT_Pack packed(particle);
    
    // now dump particle data to the census file

    // make sure census file exists
    Check (cenfile);

    // Require that buffer sizes are set correctly
    Check (packed.get_double_size() == dsize+1);
    Check (packed.get_int_size()    == isize  );
    Check (packed.get_char_size()   == csize  );

    // write the output
    // Note that we skip the first element of the double data. This is
    // time_left, which is not included in the census
    cenfile.write(reinterpret_cast<const char *>(packed.double_begin()+1),
		  dsize * sizeof(double));
    cenfile.write(reinterpret_cast<const char *>(packed.int_begin()),
		  isize * sizeof(int)); 
    cenfile.write(reinterpret_cast<const char *>(packed.char_begin()),
		  csize);
}

//---------------------------------------------------------------------------//
// write a full Comm_Buffer to an output

template<class PT>
void Particle_Buffer<PT>::write_census(std_ostream &cenfile,
				       Comm_Buffer &buffer) const
{
    // check for output file
    if (!cenfile)
	Insist(0, 
	       "You tried to write census particles to a non-existent file!");

    // determine number of particles in the buffer
    int num_particles = buffer.n_part;
    if (num_particles == 0)
	return;

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
	// get the double array data, add 1 because Comm_Buffer includes
	// time_left in its info
	for (int d = 0; d < dsize; d++)
	    ddata[d] = buffer.array_d[d+id+1];

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
	id += dsize + 1;
	ii += isize;
	ic += csize;
	buffer.n_part--;

	// asserts to make sure we haven't gone over
	Check (id <= buffer_d);
	Check (ii <= buffer_i);
	Check (ic <= buffer_c);
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
rtt_dsxx::SP<typename Particle_Buffer<PT>::Census_Buffer> 
Particle_Buffer<PT>::read_census(std_istream &cenfile) const
{
    using rtt_dsxx::SP;

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
	sf_double r;
	sf_double omega;
	for (int i = 2; i < 5; i++)
	    omega.push_back(ddata[i]);
	for (int i = 5; i < dsize; i++)
	    r.push_back(ddata[i]);
	int cell = idata[0];

	// make new random number
	int *id = unpack_sprng(rdata);
	Rnd_Type random(id, idata[1]);

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
// PARTICLE COMMUNICATION FUNCTIONS
//---------------------------------------------------------------------------//
// do a block send of a Comm_Buffer

template<class PT>
void Particle_Buffer<PT>::send_buffer(Comm_Buffer &buffer, int proc) const
{
    using C4::Send;

    // send out a Comm_Buffer, add explicit template parameter so that highly 
    // standard compliant compilers will not face ambiguity in determining
    // which Send to use
    Send <int>(buffer.n_part, proc, 200);
    Send <double>(&buffer.array_d[0], buffer_d, proc, 201);
    Send <int>(&buffer.array_i[0], buffer_i, proc, 202);
    Send <char>(&buffer.array_c[0], buffer_c, proc, 203);

    // empty the buffer
    buffer.n_part = 0;
}

//---------------------------------------------------------------------------//
// do a block recv of a Comm_Buffer

template<class PT>
rtt_dsxx::SP<typename Particle_Buffer<PT>::Comm_Buffer> 
Particle_Buffer<PT>::recv_buffer(int proc) const
{
    using C4::Recv;
    using rtt_dsxx::SP;

    // return Comm_Buffer declaration
    SP<Comm_Buffer> buffer(new Comm_Buffer());

    // receive the n_part, as above, explicitly give the template function
    // arguments 
    Recv (buffer->n_part, proc, 200);
    Recv (&buffer->array_d[0], buffer_d, proc, 201);
    Recv (&buffer->array_i[0], buffer_i, proc, 202);
    Recv (&buffer->array_c[0], buffer_c, proc, 203);

    // return SP
    return buffer;
}
    
//---------------------------------------------------------------------------//
// do an async send of a Comm_Buffer

template<class PT>
void Particle_Buffer<PT>::asend_buffer(Comm_Buffer &buffer, int proc) const
{
    using C4::SendAsync;

    // async send this Comm_Buffer
    SendAsync(buffer.comm_n, &buffer.n_part, 1, proc, 100);
    SendAsync(buffer.comm_d, &buffer.array_d[0], buffer_d, proc, 101);
    SendAsync(buffer.comm_i, &buffer.array_i[0], buffer_i, proc, 102);
    SendAsync(buffer.comm_c, &buffer.array_c[0], buffer_c, proc, 103);

    // empty the buffer
    buffer.n_part = 0;
}

//---------------------------------------------------------------------------//
// post async receives

template<class PT>
void Particle_Buffer<PT>::post_arecv(Comm_Buffer &buffer, int proc) const
{
    using C4::RecvAsync;

    // post c4 async receives
    RecvAsync(buffer.comm_n, &buffer.n_part, 1, proc, 100);
    RecvAsync(buffer.comm_d, &buffer.array_d[0], buffer_d, proc, 101);
    RecvAsync(buffer.comm_i, &buffer.array_i[0], buffer_i, proc, 102);
    RecvAsync(buffer.comm_c, &buffer.array_c[0], buffer_c, proc, 103);
}

//---------------------------------------------------------------------------//
// post waits on the receives and sends to fill up a Comm_Buffer

template<class PT> 
void Particle_Buffer<PT>::async_wait(Comm_Buffer &buffer) const
{
    // wait on recieve buffers to make sure they are full
    buffer.comm_n.wait();
    buffer.comm_d.wait();
    buffer.comm_i.wait();
    buffer.comm_c.wait();
}

//---------------------------------------------------------------------------//
// check to see if our ship has come in or gone out

template<class PT>
bool Particle_Buffer<PT>::async_check(Comm_Buffer &buffer) const
{
    using std::vector;
    using std::fill;
    using std::accumulate;

    // tag to check what is in; we want all or nothing
    vector<int> arrived(4);
    fill(arrived.begin(), arrived.end(), 0);
    int total = 0;
    int count = 0;
    
    // check to see if the buffers have been received
    do
    {
	// check comm_n
	if (arrived[0] == 0)
	    if (buffer.comm_n.complete())
		arrived[0] = 1;
	
	// check comm_d
	if (arrived[1] == 0)
	    if (buffer.comm_d.complete())
		arrived[1] = 1;

	// check comm_i
	if (arrived[2] == 0)
	    if (buffer.comm_i.complete())
		arrived[2] = 1;

	// check comm_c
	if (arrived[3] == 0)
	    if (buffer.comm_c.complete())
		arrived[3] = 1;

	// accumulate total
	total = accumulate(arrived.begin(), arrived.end(), 0);
	count++;
	
    } while (total > 0 && total < 4);

    Ensure (total == 0 || total == 4);
    return total;
}

//---------------------------------------------------------------------------//
// free the Comm_Buffer

template<class PT>
void Particle_Buffer<PT>::async_free(Comm_Buffer &buffer) const
{
    // free the C4_Req objects from Async posts, note that these must be
    // reassigned with a post request before they can be received or tested 
    buffer.comm_n.free();
    buffer.comm_d.free();
    buffer.comm_i.free();
    buffer.comm_c.free();
}

//---------------------------------------------------------------------------//
// check to see if the C4_Req objects are in use or not

template<class PT>
bool Particle_Buffer<PT>::comm_status(Comm_Buffer &buffer) const
{
    if (!buffer.comm_n.inuse())
	return false;
    else if (!buffer.comm_d.inuse())
	return false;
    else if (!buffer.comm_i.inuse())
	return false;
    else if (!buffer.comm_c.inuse())
	return false;

    // if we haven't returned then these are still active
    return true;
}

//---------------------------------------------------------------------------//
// BUFFERING FUNCTIONS
//---------------------------------------------------------------------------//
// add a Census_Buffer to a Comm_Buffer

template<class PT>
void Particle_Buffer<PT>::buffer_census(Comm_Buffer &comm, 
					const Census_Buffer &census) const
{
    // check to make sure this Comm_Buffer isn't full
    Require (comm.n_part < buffer_s);
    Require (comm.n_part >= 0);

    // calculate indices for the buffer
    int id = comm.n_part * (dsize + 1);
    int ii = comm.n_part * isize;
    int ic = comm.n_part * csize;

    // add one Census_Buffer particle to the buffer
  
    // start with the double data, time_left is not part of the census, but is
    // part of the Comm_Buffer, any Census_Particle will have a time_left of 1
    comm.array_d[id++] = 1;
    comm.array_d[id++] = census.ew;
    comm.array_d[id++] = census.fraction;
    for (int i = 0; i < census.omega.size(); i++)
	comm.array_d[id++] = census.omega[i];
    for (int i = 0; i < census.r.size(); i++)
	comm.array_d[id++] = census.r[i];
    Check (id - comm.n_part * (dsize+1) == (dsize+1));

    // do the int data
    comm.array_i[ii++] = census.cell;
    comm.array_i[ii++] = census.random.get_num();
    Check (ii - comm.n_part * isize == isize);

    // do the char data
    char *bytes;
    int size = pack_sprng(census.random.get_id(), &bytes);
    Check (size == csize);
    for (int i = 0; i < csize; i++)
	comm.array_c[ic++] = bytes[i];
    std::free(bytes);
    Check (ic - comm.n_part * csize == csize);
  
    // update the n_particles
    comm.n_part++;

    // do some assertions
    Ensure (id <= buffer_d);
    Ensure (ii <= buffer_i);
    Ensure (ic <= buffer_c);
}

//---------------------------------------------------------------------------//
// add a Particle to a Comm_Buffer

template<class PT>
void Particle_Buffer<PT>::buffer_particle(Comm_Buffer &buffer,
					  const PT &particle) const
{
    // check to make sure this Comm_Buffer isn't full
    Require (buffer.n_part < buffer_s);
    Require (buffer.n_part >= 0);

    // calculate indices for the buffer, remember the particle info required
    // during a timestep must also include the time left, hence dsize+1.
    int id = buffer.n_part * (dsize+1);
    int ii = buffer.n_part * isize;
    int ic = buffer.n_part * csize;
    
    // create a packed version of the particle;
    PT_Pack packed(particle);

    // Check that the data sizes are consistent
    Check (packed.get_double_size() == dsize+1);
    Check (packed.get_int_size()    == isize  );
    Check (packed.get_char_size()   == csize  );

    // add the double data
    std::copy(packed.double_begin(), packed.double_end(), buffer.array_d+id);
    std::copy(packed.int_begin(),    packed.int_end(),    buffer.array_i+ii);
    std::copy(packed.char_begin(),   packed.char_end(),   buffer.array_c+ic);

    // update the n_particles
    buffer.n_part++;

    // do some assertions
    Ensure (id+dsize+1 <= buffer_d);
    Ensure (ii+isize   <= buffer_i);
    Ensure (ic+csize   <= buffer_c);
}

//---------------------------------------------------------------------------//
// make a Particle bank out of a Comm_Buffer

template<class PT>
void Particle_Buffer<PT>::add_to_bank(Comm_Buffer &buffer, Bank &bank) const 
{
    using rtt_dsxx::SP;

    // get the number of Particles in the Comm_Buffer
    int num_part = buffer.n_part;
    Require (num_part <= buffer_s);
    Require (num_part >= 0);

    // loop through Comm_Buffer and add the Particles to the Bank.
    // Particles are numbered 0 through num_part-1
    for (int n = 0; n < num_part; n++)
    {

	// find data locations
	double* d_data = buffer.array_d + (dsize+1)*n;
	int*    i_data = buffer.array_i +  isize   *n;
	char*   c_data = buffer.array_c +  csize   *n;

	// Check that enough data remains to be read
	Check(d_data + dsize+1 <= buffer.array_d + buffer_d);
	Check(i_data + isize   <= buffer.array_i + buffer_i);
	Check(c_data + csize   <= buffer.array_c + buffer_c);

	// unpack the particle, make a pointer to it
	SP<PT> particle(PT_Pack::unpack(dsize+1, d_data, isize, i_data, 
					csize, c_data)) ;

	// add the Particle to the Bank and update the num_part counter
	bank.push(particle);
	buffer.n_part--;
    }

    // the buffer should now be empty
    Ensure (buffer.n_part == 0);
}

//===========================================================================//
// CENSUS BUFFER FUNCTIONS
//===========================================================================//
/*!
 * \brief Make a particle out of a Census_Buffer.
 */
template<class PT>
rtt_dsxx::SP<PT> Particle_Buffer<PT>::Census_Buffer::make_Particle() const
{
    // make a particle
    rtt_dsxx::SP<PT> particle;

    // load particle data, set time_left to 1.0-->this will be adjusted in
    // the source
    particle = new PT(r, omega, ew, cell, random, fraction, 1.0, PT::CENSUS);

    return particle;
}

} // end namespace rtt_imc

//---------------------------------------------------------------------------//
//                              end of Particle_Buffer.t.hh
//---------------------------------------------------------------------------//
