//----------------------------------*-C++-*----------------------------------//
// Source.t.hh
// Thomas M. Evans
// Thu May 14 08:45:49 1998
//---------------------------------------------------------------------------//
// @> Source class implementation file
//---------------------------------------------------------------------------//

#include "Source.hh"
#include "Global.hh"
#include "ds++/Assert.hh"
#include <iomanip>
#include <cmath>
#include <vector>

namespace rtt_imc 
{

// Draco functions
using rtt_rng::Sprng;
using rtt_rng::Rnd_Control;
using rtt_mc::global::pi;
using dsxx::SP;

// STL functions
using std::ios;
using std::cos;
using std::sin;
using std::sqrt;
using std::pow;
using std::endl;
using std::setiosflags;
using std::setw;
using std::string;
using std::ostream;
using std::vector;

//---------------------------------------------------------------------------//
// constructor
//---------------------------------------------------------------------------//
// standard constructor for the Source class

template<class MT, class PT>
Source<MT, PT>::Source(ccsf_int &vol_rnnum_, 
		       ccsf_int &nvol_,
		       ccsf_double &ew_vol_,
		       ccsf_int &ss_rnnum_, 
		       ccsf_int &nss_, 
		       ccsf_int &fss_,
		       ccsf_double &ew_ss_,
		       SP_Census census_,
		       string ssd, 
		       int nvoltot_, 
		       int nsstot_,
		       SP_Rnd_Control rcon_, 
		       SP_Mat_State mat_state,
		       SP_Mesh_Op operations) 
    : vol_rnnum(vol_rnnum_), 
      nvol(nvol_), 
      ew_vol(ew_vol_), 
      ss_rnnum(ss_rnnum_), 
      nss(nss_), 
      fss(fss_), 
      ew_ss(ew_ss_),
      ss_dist(ssd), 
      census(census_), 
      nvoltot(nvoltot_), 
      nsstot(nsstot_),
      ncentot(census->size()), 
      rcon(rcon_),
      material(mat_state), 
      mesh_op(operations)
{
    // some assertions
    Check (vol_rnnum.get_Mesh() == nvol.get_Mesh());
    Check (nvol.get_Mesh()      == ss_rnnum.get_Mesh());
    Check (nvol.get_Mesh()      == nss.get_Mesh());
    Check (nvol.get_Mesh()      == ew_vol.get_Mesh());
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
	    if (nssdone == nsstot)
	    {
		current_cell  = 1;
		nsrcdone_cell = 0;
	    }
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
	    if (nvoldone == nvoltot)
	    {
		current_cell  = 1;
		nsrcdone_cell = 0;
	    }
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

    Ensure (sampled);
    Ensure (source_particle->status());
    
    // return the particle
    return source_particle;
}

//---------------------------------------------------------------------------//
// sample surface source particle

template<class MT, class PT>
SP<PT> Source<MT, PT>::get_ss(double delta_t)
{
    // draco directives
    using rtt_mc::global::pi;

    // face on which surface source resides
    int face = fss(current_cell);

    // get the random number object
    Sprng rand = rcon->get_rn();

    // sample location
    vector<double> r = nss.get_Mesh().sample_pos_on_face
	(current_cell, face, rand);

    // find inward normal, sample direction, and add

    // normal distributed surface source
    vector<double> omega = nss.get_Mesh().get_normal_in(current_cell, face); 
    
    // cosine distribution of surface source about normal
    if (ss_dist == "cosine")
    {
	double costheta = sqrt(rand.ran());
	double phi = 2.0 * pi * rand.ran();
	nss.get_Mesh().get_Coord().calc_omega(costheta, phi, omega);
    }

    // complete description of surface source particle    
    double ew = ew_ss(current_cell);
    int cell = current_cell;
    double fraction = 1.0;
    double time_left = rand.ran() * delta_t;

    // instantiate particle to return
    SP<PT> ss_particle(new PT(r, omega, ew, cell, rand, fraction, 
			      time_left));

    // return the ss particle;
    return ss_particle;
}

//---------------------------------------------------------------------------//
// sample volume emission particle

template<class MT, class PT>
SP<PT> Source<MT, PT>::get_evol(double delta_t)
{
    // get the random number object
    Sprng rand = rcon->get_rn();

    // sample location using tilt
    double t4 = pow(material->get_T(current_cell), 4);
    vector<double> r = mesh_op->sample_pos_tilt(current_cell, t4, rand);

    // sample particle direction
    vector<double> omega = nvol.get_Mesh().get_Coord().
	sample_dir("isotropic", rand); 
 
    double ew        = ew_vol(current_cell);
    int cell         = current_cell;
    double fraction  = 1.0;
    double time_left = rand.ran() * delta_t;

    // instantiate particle to return
    SP<PT> vol_particle(new PT(r, omega, ew, cell, rand, fraction, 
			       time_left)); 

    return vol_particle;
}

//---------------------------------------------------------------------------//
// read census particle

template<class MT, class PT>
SP<PT> Source<MT, PT>::get_census(double delta_t)
{
    Require (census->size() > 0);

    // get the census particle from the Census buffer
    SP<PT> census_particle = census->top();
    census_particle->set_time_left(delta_t);
    census_particle->set_descriptor("from_census");

    // remove the census particle from the bank
    census->pop();

    // make sure the particle is active
    census_particle->reset_status();

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
    out << endl;
    out << ">>> SOURCE DATA <<<" << endl;
    out << "===================" << endl;

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
    out.setf(ios::scientific, ios::floatfield);
    for (int i = 1; i <= nvol.get_Mesh().num_cells(); i++)
	out << setw(10) << i << setw(15) << ew_vol(i) << setw(15)
	    << ew_ss(i) << endl;
}
	
} // end namespace rtt_imc

//---------------------------------------------------------------------------//
//                              end of Source.t.hh
//---------------------------------------------------------------------------//
