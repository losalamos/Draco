//----------------------------------*-C++-*----------------------------------//
// Particle_Buffer.cc
// Thomas M. Evans
// Tue May 12 14:34:34 1998
//---------------------------------------------------------------------------//
// @> Particle_Buffer class implementation file
//---------------------------------------------------------------------------//

#include "imctest/Particle_Buffer.hh"

IMCSPACE

//---------------------------------------------------------------------------//
// constructors
//---------------------------------------------------------------------------//
// templated constructor

template<class PT>
template<class MT>
Particle_Buffer<PT>::Particle_Buffer(const MT &mesh)
{
  // determine size of double info from Particles;
  // 5 = omega(3) + ew + fraction
    dsize = mesh.get_Coord().get_dim() + 5;

  // determine size of integer info from Particles;
  // 3 = cell + streamnum + random_state_size
    isize = 3;
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
    idata[2] = pack_sprng(particle.random.get_id(), &rdata);
    Check (idata[2] < RNG::max_buffer);

  // now dump particle data to the census file

  // make sure census file exists
    Check (cenfile);

  // write the output
    cenfile.write(reinterpret_cast<const char *>(ddata), dsize *
		  sizeof(double));
    cenfile.write(reinterpret_cast<const char *>(idata), isize * sizeof(int));
    cenfile.write(reinterpret_cast<const char *>(rdata), idata[2]);

  // reclaim dynamic memory
    delete [] ddata;
    delete [] idata;
    std::free(rdata);
}

CSPACE

//---------------------------------------------------------------------------//
//                              end of Particle_Buffer.cc
//---------------------------------------------------------------------------//
