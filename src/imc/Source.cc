//----------------------------------*-C++-*----------------------------------//
// Source.cc
// Thomas M. Evans
// Thu May 14 08:45:49 1998
//---------------------------------------------------------------------------//
// @> Source class implementation file
//---------------------------------------------------------------------------//

#include "imc/Source.hh"
#include "imc/Global.hh"
#include "ds++/Assert.hh"
#include <iomanip>
#include <cmath>

IMCSPACE

// Draco functions
using Global::pi;

// STL functions
using std::ios;
using std::cos;
using std::sin;
using std::sqrt;
using std::pow;
using std::endl;
using std::setiosflags;
using std::setw;


//---------------------------------------------------------------------------//
// constructor
//---------------------------------------------------------------------------//
// standard constructor for the Source class

template<class MT, class PT>
Source<MT, PT>::Source(typename MT::CCSF_int &vol_rnnum_, 
		       typename MT::CCSF_int &nvol_,
		       typename MT::CCSF_double &ew_vol_,
		       typename MT::CCVF_double &t4_,
		       typename MT::CCSF_int &ss_rnnum_, 
		       typename MT::CCSF_int &nss_, 
		       typename MT::CCSF_int &fss_,
		       typename MT::CCSF_double &ew_ss_,
		       string title,
		       int nvoltot_, int nsstot_, int ncentot_,
		       SP<Rnd_Control> rcon_, 
		       const Particle_Buffer<PT> &buffer_,
		       SP<Mat_State<MT> > mat_state) 
    : vol_rnnum(vol_rnnum_), nvol(nvol_), ew_vol(ew_vol_), t4_slope(t4_),
      ss_rnnum(ss_rnnum_), nss(nss_), fss(fss_), ew_ss(ew_ss_),
      census(title.c_str(), ios::in), nvoltot(nvoltot_), nsstot(nsstot_),
      ncentot(ncentot_), rcon(rcon_), buffer(buffer_), material(mat_state) 
{
  // some assertions
    Check (vol_rnnum.get_Mesh() == nvol.get_Mesh());
    Check (nvol.get_Mesh()      == ss_rnnum.get_Mesh());
    Check (nvol.get_Mesh()      == nss.get_Mesh());
    Check (nvol.get_Mesh()      == ew_vol.get_Mesh());
    Check (nvol.get_Mesh()      == t4_slope.get_Mesh());
    Check (nvol.get_Mesh()      == fss.get_Mesh());
    Check (nvol.get_Mesh()      == ew_ss.get_Mesh());

  // nsrcdone_cell is the running number of source particles completed for a
  // particular source type in a particular cell.  
    nsrcdone_cell = 0;

  // Begin with first cell
    current_cell = 1;

  // running totals of completed source particles, by type
    nssdone  = 0;
    nvoldone = 0;
    ncendone = 0;
}

//---------------------------------------------------------------------------//
// Source interface functions
//---------------------------------------------------------------------------//
// get a source particle

template<class MT, class PT>
SP<PT> Source<MT, PT>::get_Source_Particle(double delta_t)
{
    bool sampled = false;

  // instantiate particle to return
    SP<PT> source_particle;

  // do all surface source particles, one-by-one
    while (!sampled && nssdone < nsstot)
    {
	if (nsrcdone_cell < nss(current_cell))
	{
	    rcon->set_num(ss_rnnum(current_cell) + nsrcdone_cell);
	    source_particle = get_ss(delta_t);
	    sampled = true;
	    nsrcdone_cell++;
	    nssdone++;
	}
	else
	{
	    current_cell++;
	    nsrcdone_cell = 0;
	    if (current_cell > nss.get_Mesh().num_cells())
	    {
		Check (nssdone == nsstot);
		current_cell = 1;
	    }
	}
    }

  // do all volume emission particles, one-by-one
    while (!sampled && nvoldone < nvoltot)
    {
	if (nsrcdone_cell < nvol(current_cell))
	{
	    rcon->set_num(vol_rnnum(current_cell) + nsrcdone_cell);
	    source_particle = get_evol(delta_t);
	    sampled = true;
	    nsrcdone_cell++;
	    nvoldone++;
	}
	else
	{
	    current_cell++;
	    nsrcdone_cell = 0;
	    if (current_cell > nvol.get_Mesh().num_cells())
	    {
		Check (nvoldone == nvoltot);
		current_cell = 1;
	    }
	}
    }

  // do all census particles, one-by-one
    while (!sampled && ncendone < ncentot)
    {
	    source_particle = get_census(delta_t);
	    sampled = true;
	    ncendone++;
    }

    Check (sampled);

    return source_particle;
}

//---------------------------------------------------------------------------//
// sample surface source particle

template<class MT, class PT>
SP<PT> IMC::Source<MT, PT>::get_ss(double delta_t)
{

  // face on which surface source resides
    int face = fss(current_cell);

  // get the random number object
    Sprng rand = rcon->get_rn();

  // sample location
    vector<double> r = nss.get_Mesh().sample_pos_on_face
	(current_cell, face, rand);

  // find inward normal, sample direction, and add
    vector<double> omega = nss.get_Mesh().get_normal_in(current_cell, face); 
    double costheta = sqrt(rand.ran());
    double phi = 2.0 * Global::pi * rand.ran();
    nss.get_Mesh().get_Coord().calc_omega(costheta, phi, omega);

    double ew = ew_ss(current_cell);
    int cell = current_cell;
    double fraction = 1.0;
    double time_left = rand.ran() * delta_t;

  // instantiate particle to return
    SP<PT> ss_particle = new PT(r, omega, ew, cell, rand, fraction, 
				time_left);

  // return the ss particle;
    return ss_particle;
}

//---------------------------------------------------------------------------//
// sample volume emission particle

template<class MT, class PT>
SP<PT> IMC::Source<MT, PT>::get_evol(double delta_t)
{
  // get the random number object
    Sprng rand = rcon->get_rn();

  // sample location using tilt
    double t4 = pow(material->get_T(current_cell), 4);
    vector<double> r = nvol.get_Mesh().
	sample_pos(current_cell, rand, t4_slope(current_cell), t4);

  // sample particle direction
    vector<double> omega = nvol.get_Mesh().get_Coord().
	sample_dir("isotropic", rand); 
 
    double ew = ew_vol(current_cell);
    int cell = current_cell;
    double fraction = 1.0;
    double time_left = rand.ran() * delta_t;

  // instantiate particle to return
    SP<PT> vol_particle = new PT(r, omega, ew, cell, rand, fraction, 
				 time_left); 

    return vol_particle;
}

//---------------------------------------------------------------------------//
// read census particle

template<class MT, class PT>
SP<PT> Source<MT, PT>::get_census(double delta_t)
{
    SP<PT> census_particle;

  // read a Census_Buffer from the census file
    SP<Particle_Buffer<PT>::Census_Buffer> cenpart;
    cenpart = buffer.read_census(census);
    Check (cenpart);

  // make a Particle from the Census buffer
    census_particle = new PT(cenpart->r, cenpart->omega, cenpart->ew,
			     cenpart->cell, cenpart->random,
			     cenpart->fraction, delta_t);

  // return the particle
    return census_particle;
}

//---------------------------------------------------------------------------//
// source diagnostic functions
//---------------------------------------------------------------------------//
// print out an ascii description of the source

template<class MT, class PT>
void Source<MT, PT>::print(ostream &out) const
{
    out << "*** SOURCE DATA ***" << endl;
    out << "-------------------" << endl;

  // numbers of each type of particle
    out << setw(25) << setiosflags(ios::right) << "Census Particles: "
	<< setw(10) << ncentot << endl;
    out << setw(25) << setiosflags(ios::right) << "Volume Particles: "
	<< setw(10) << nvoltot << endl;
    out << setw(25) << setiosflags(ios::right) << "Surface Particles: "
	<< setw(10) << nsstot << endl;

  // lets look at the number of particles in each cell
    out << endl << " ** Sources **" << endl;
    out << setw(10) << " " << setw(15) << setiosflags(ios::right)
	<< "Volume" << setw(5) << " " << setw(15) << setiosflags(ios::right)
	<< "Surface" << endl;
    out << setw(10) << setiosflags(ios::right) << "Cell"
	<< setw(10) << setiosflags(ios::right) << "Number"
	<< setw(10) << setiosflags(ios::right) << "Stream ID"
	<< setw(10) << setiosflags(ios::right) << "Number"
	<< setw(10) << setiosflags(ios::right) << "Stream ID" << endl;
    for (int i = 1; i <= nvol.get_Mesh().num_cells(); i++)
	out << setw(10) << i << setw(10) << nvol(i) 
	    << setw(10) << vol_rnnum(i)  << setw(10) << nss(i)
	    << setw(10) << ss_rnnum(i)   << endl;

  // lets look at the energy weight in each cell
    out << endl << " ** Source Energy ** " << endl;
    out << setw(10) << setiosflags(ios::right) << "Cell" 
	<< setw(15) << setiosflags(ios::right) << "Volume ew" 
	<< setw(15) << setiosflags(ios::right) << "Surface ew" << endl;
    out.precision(3);
    out << setiosflags(ios::fixed);
    for (int i = 1; i <= nvol.get_Mesh().num_cells(); i++)
	out << setw(10) << i << setw(15) << ew_vol(i) << setw(15)
	    << ew_ss(i) << endl;
}
	
CSPACE

//---------------------------------------------------------------------------//
//                              end of Source.cc
//---------------------------------------------------------------------------//
