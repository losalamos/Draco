//----------------------------------*-C++-*----------------------------------//
// Source.cc
// Thomas M. Evans
// Thu May 14 08:45:49 1998
//---------------------------------------------------------------------------//
// @> Source class implementation file
//---------------------------------------------------------------------------//

#include "imc/Source.hh"
#include "ds++/Assert.hh"
#include <iomanip>

IMCSPACE

// STL functions
using std::ios;
using Global::pi;
using std::cos;
using std::sin;
using std::sqrt;


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
		       typename MT::CCSF_int &ss_rnnum_, 
		       typename MT::CCSF_int &nss_, 
		       string title,
		       int nvoltot_, int nsstot_, int ncentot_,
		       SP<Rnd_Control> rcon_, 
		       const Particle_Buffer<PT> &buffer_) 
    : vol_rnnum(vol_rnnum_), nvol(nvol_), ss_rnnum(ss_rnnum_), nss(nss_),
      census(title.c_str(), ios::in), nvoltot(nvoltot_), nsstot(nsstot_),
      ncentot(ncentot_), rcon(rcon_), buffer(buffer_)
{
  // some assertions
    Check (vol_rnnum.get_Mesh() == nvol.get_Mesh());
    Check (nvol.get_Mesh()      == ss_rnnum.get_Mesh());
    Check (nvol.get_Mesh()      == nss.get_Mesh());

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
// sample surface source particle

// template<class MT, class PT>
// SP<PT> IMC::Source<MT, PT>::get_ss()
// {

//   // face on which surface source resides
//     int face = fss(current_cell);

//   // sample location
//     vector <double> r = mesh.sample_pos_on_face("uniform", current_cell,
// 						face, Sprng &rand);

//   // find inward normal, sample direction, and add
//     vector <double> omega = mesh::get_in_normal(current_cell, face);
//     double costheta = sqrt(rand.ran());
//     double phi = 2.0 * Global::pi * rand.ran();
//     mesh.get_Coord().calc_omega(costheta, phi, omega);

//     double ew = ew_ss(current_cell);
//     int cell = current_cell;
//     double fraction = 1.0;
//     double time_left = rand.ran() * delta_t;

//   // instantiate particle to return
//     SP<Particle<MT> > ss_particle(r, omega, ew, cell, rand, 
// 				  fraction, time_left);

//     return ss_particle;
// }

//---------------------------------------------------------------------------//
// sample volume emission particle

// template<class MT, class PT>
// SP<PT> IMC::Source<MT, PT>::get_evol()
// {

//   // sample location
//     vector <double> r = mesh.sample_pos("uniform", current_cell, 
// 					Sprng &rand);

//   // sample particle direction
//     vector<double> omega = mesh.get_Coord().sample_dir("isotropic", rand);

//     double ew = ew_vol(current_cell);
//     int cell = current_cell;
//     double fraction = 1.0;
//     double time_left = rand.ran() * delta_t;

//   // instantiate particle to return
//     SP<Particle<MT> > vol_particle(r, omega, ew, cell, rand, 
// 				  fraction, time_left);

//     return vol_particle;
// }

//---------------------------------------------------------------------------//
// read census particle

// template<class MT, class PT>
// SP<PT> IMC::Source<MT, PT>::get_census()
// {

//   // read census particle
//     vector <double> r;
//     vector <double> omega;
//     double ew;
//     int cell;
//     double fraction;

//   // time remaining is the entire time step
//     double time_left = delta_t;

//   // instantiate particle to return
//     SP<Particle<MT> > cen_particle(r, omega, ew, cell, rand, 
// 				   fraction, time_left);

//     return census_particle;
// }

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
    out << setw(10) << " " << setw(13) << setiosflags(ios::right)
	<< "Volume" << setw(7) << " " << setw(12) << setiosflags(ios::right)
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
}
	
CSPACE

//---------------------------------------------------------------------------//
//                              end of Source.cc
//---------------------------------------------------------------------------//
