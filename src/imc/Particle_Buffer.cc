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
// templated constructor for Particle_Buffer

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
// constructors for Particle_Buffer<PT>::Census Particle struct

template<class PT>
Particle_Buffer<PT>::Census_Particle::Census_Particle(vector<double> &r_,
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
Particle_Buffer<PT>::Census_Particle::Census_Particle()
    : random(0, 0)
{
  // constructor for use with STL, this cannot be used
    Insist (0, "You tried to default construct a Census_Particle!");
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

//---------------------------------------------------------------------------//
// read a single particle from an output

template<class PT>
SP<Particle_Buffer<PT>::Census_Particle> 
Particle_Buffer<PT>::read_census(istream &cenfile)
{
  // make sure file exists
    Check (cenfile);

  // cast smart pointer to Census_Particle
    SP<Census_Particle> return_part;
    
  // set pointers for dynamic memory storage
    double *ddata = new double[dsize];
    int    *idata = new int[isize];
    char   *rdata;

  // read in data
    cenfile.read(reinterpret_cast<char *>(ddata), dsize * sizeof(double));
    if (!cenfile.eof())
    {
      // read in integer data
	cenfile.read(reinterpret_cast<char *>(idata), isize * sizeof(int));
	Check (!cenfile.eof());

      // get size of random number state
	rdata = new char[idata[2]];
	cenfile.read(reinterpret_cast<char *>(rdata), idata[2]);
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

      // make new Census_Particle
	return_part = new Census_Particle(r, omega, ew, frac, cell, random);

      // reclaim dynamic memory
	delete [] rdata;
    }

  // reclaim dynamic memory
    delete [] idata;
    delete [] ddata;

  // return Census_Particle
    return return_part;
}

CSPACE

//---------------------------------------------------------------------------//
//                              end of Particle_Buffer.cc
//---------------------------------------------------------------------------//
