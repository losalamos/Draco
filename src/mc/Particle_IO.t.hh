//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Particle_IO.t.hh
 * \author Thomas M. Evans
 * \date   Fri Dec 21 10:10:55 2001
 * \brief  Particle_IO class template members.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __mc_Particle_IO_t_hh__
#define __mc_Particle_IO_t_hh__

#include "Particle_IO.hh"
#include "Particle_Buffer.hh"
#include "ds++/SP.hh"
#include "ds++/Assert.hh"
#include <vector>

namespace rtt_mc
{

//---------------------------------------------------------------------------//
// IO FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Write a particle to an output stream.
 *
 * The particle is packed and written to the output stream in binary format.
 *
 * \param out output stream of type std::ostream, this will usually refer to
 * an ofstream file object
 *
 * \param particle particle that has a pack function
 */
template<class PT>
void Particle_IO<PT>::write_particle(std_ostream &out, const PT &particle)
{
    using std::vector;

    // create a packed version of the particle
    vector<char> packed_PT = particle.pack();
    int          size      = packed_PT.size();
    int          number    = 1;
    
    // now dump particle data to the output stream

    // make sure output object exists
    Check (out);

    // write the output
    out.write(reinterpret_cast<const char *>(&number), sizeof(int));
    out.write(reinterpret_cast<const char *>(&size), sizeof(int));
    out.write(&packed_PT[0], packed_PT.size());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Write a Particle_Buffer to an output stream.
 *
 * The Particle_Buffer is packed and written to the output stream in binary
 * format.
 *
 * \param out output stream of type std::ostream, this will usually refer to
 * an ofstream file object
 *
 * \param buffer a Particle_Buffer object
 */
template<class PT>
void Particle_IO<PT>::write_Particle_Buffer(std_ostream &out,
					    const Particle_Buffer<PT> &buffer)
{
    using std::vector;

    // get buffer size and number of particles data
    int size   = buffer.int_data[1];
    int number = buffer.int_data[0];
    
    // now dump particle data to the output stream

    // make sure output object exists
    Check (out);

    // write the output
    out.write(reinterpret_cast<const char *>(&number), sizeof(int));
    out.write(reinterpret_cast<const char *>(&size), sizeof(int));
    out.write(&buffer.packed_particles[0], size);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Read all the particles from an input stream into a bank.
 *
 * This function reads all particles from the input stream into the bank.  It
 * expects the input stream to consist entirely of particles.
 * 
 * \param in input stream of type std::istream, this will usually refer to an
 * ifstream file object
 *
 * \param bank bank to add particles to
 */
template<class PT>
void Particle_IO<PT>::read_particles(std_istream &in, Bank &bank)
{
    using std::vector;
    using rtt_dsxx::SP;

    // make sure file exists
    Check (in);

    // size of particle(s) and number of particles
    int size                 = 0;
    int number               = 0;
    int size_packed_particle = 0;
    SP<PT> particle;

    // boolean for finishing reading the file
    bool done = false;

    while (!done)
    {
	// read the number of particles
	in.read(reinterpret_cast<char *>(&number), sizeof(int));

	// if we are at the end of the file then quit
	if (in.eof())
	{
	    done = true;
	}

	// else continue reading particles
	else
	{
	    // read the size of the characters
	    in.read(reinterpret_cast<char *>(&size), sizeof(int));
	    Check (!in.eof());

	    // determine the size of a single particle
	    if (number > 0)
	    {
		// determine the size of the packed particle
		size_packed_particle = size / number;
		Check (size % number == 0);

		// make a packed particle
		vector<char> packed_PT(size_packed_particle);
		
		for (int i = 0; i < number; i++)
		{
		    // read in the particle
		    in.read(&packed_PT[0], size_packed_particle);
		    Check (!in.eof());
		    
		    // build the particle
		    particle = new PT(packed_PT);
		    Check (particle);
		    
		    // add it to the bank
		    bank.push(particle);
		}
	    }
	}
    }
}

} // end namespace rtt_mc

#endif                          // __mc_Particle_IO_t_hh__

//---------------------------------------------------------------------------//
//                        end of mc/Particle_IO.t.hh
//---------------------------------------------------------------------------//
