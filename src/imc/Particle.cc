//----------------------------------*-C++-*----------------------------------//
// Particle.cc
// Thomas M. Evans
// Fri Jan 30 17:04:24 1998
//---------------------------------------------------------------------------//
// @> source file for Particle class
//---------------------------------------------------------------------------//

#include "imctest/Particle.hh"
#include "imctest/Global.hh"
#include <iostream>
#include <iomanip>
#include <cmath>

IMCSPACE

using std::cout;
using std::endl;
using std::setw;
using std::ios;
using std::setiosflags;

//---------------------------------------------------------------------------//
// constructors
//---------------------------------------------------------------------------//
// default constructor, explicit definition of Particle required to
// guarantee that the Particle is defined by a Mesh;
// vectors handle their own dynamic memory
template<class MT>
Particle<MT>::Particle(const MT &mesh, int seed, double weight_, 
		       double energy_)
    :weight(weight_), energy(energy_), r(mesh.Coord().Get_dim(), 0.0), 
     omega(mesh.Coord().Get_sdim(), 0.0), cell(0), alive(true), random(seed)
{
  // the syntax for initialization of the vector is equivalent to:
  // vector<double> r(dimension, 0.0),
  // which means, "make a vector of type double of dimension
  // initialized with all zeroes--the vectors must be
  // initialized by 0.0 because if the arguments match the iterator
  // constructor will be called instead, see Stroustrup pg. 450
}

// need to write assignment operators and copy constructors to take care of
// the random number

//---------------------------------------------------------------------------//
// public transport member functions
//---------------------------------------------------------------------------//
// calculate source from a given point
template<class MT>
void Particle<MT>::Source(vector<double> &r_, vector<double> &omega_,
			  const MT &mesh)
{
  // initial location
    r = r_;
    omega = omega_;

  // sample initial direction (isotropic)
  // mesh.Coord().Set_omega(omega, random);

  // find particle cell
    cell = mesh.Get_cell(r);
}

// transport particle through mesh
template<class MT>
void Particle<MT>::Transport(const MT &mesh, const Opacity<MT> &xs)
{
  // explicit calls from standard library namespace
    using std::log;
  // explicit calls from Global namespace
    using Global::pi;

  // transport loop, ended when transport = false
    while (alive)
    {
      // dist-to-collision and dist-to-boundary definitions
        double d_collision, d_boundary;
      // cell face definition
        int face = 0;
        
      // sample distance-to-collision
        d_collision = -log(random.ran()) / xs.Sigma(cell);

      // get distance-to-boundary and cell face
        d_boundary  = mesh.Get_db(r, omega, cell, face);

      // for 3D transport in 2D meshes, do a transform, for 3D transform
      // simply returns value
	double d_travel = mesh.Coord().Transform(d_boundary, omega);

      // do diagnostic print
	Print();
	cout << "Boundary:  " << d_boundary << endl;
	cout << "Collision: " << d_collision << endl;
	cout << "Travel:    " << d_travel << endl;

      // streaming
        if (d_collision <= d_boundary)
        {
            Stream(d_collision);
            alive = Collide(mesh, xs);
        }
        else
        {
            Stream(d_travel);
            cell  = mesh.Next_cell(cell, face);
	    alive = cell;
        }
    }
  
  // end-particle diagnostic print
    Print();
}

//---------------------------------------------------------------------------//
// private transport member functions
//---------------------------------------------------------------------------//
// stream particle a distance d
template<class MT>
void Particle<MT>::Stream(double distance)
{
  // calculate new location
    for (int i = 0; i <= r.size()-1; i++)
        r[i] = r[i] + distance * omega[i];
}

// calculate everything about a collision
template<class MT>
bool Particle<MT>::Collide(const MT &mesh, const Opacity<MT> &xs)
{   
    bool status;

  // determine absorption or collision
    if (random.ran() <= xs.Sigma(cell) / xs.Sigma(cell))
        status = false;
    else
    {
        status = true;
        
      // calculate theta and phi (isotropic)
        double costheta, phi;
        costheta = 1 - 2 * random.ran();
        phi      = 2 * Global::pi * random.ran();

      // get new direction cosines
        mesh.Coord().Calc_omega(costheta, phi, omega);
    }

    return status;
}

//---------------------------------------------------------------------------//
// public diagnostic member functions
//---------------------------------------------------------------------------//
template<class MT>
void Particle<MT>::Print()
{
  // print particulars of the particle based on its status
    if (alive == true)
	Print_alive();
    else
	Print_dead();
}

template<class MT>
void Particle<MT>::Print_alive()
{
  // print active particle (alive = true)
    cout << "Particle is alive." << endl;

  // coordinates
    cout << setw(15) << setiosflags(ios::right) << "Coordinates: ";
    for (int i = 0; i < r.size(); i++)
	cout << setw(8) << r[i] << " ";
    cout << endl;
    
  // direction
    cout << setw(15) << setiosflags(ios::right) << "Direction: ";
    for (int i = 0; i < omega.size(); i++)
	cout << setw(8) << omega[i] << " ";
    cout << endl;

  // cell
    cout << setw(15) << setiosflags(ios::right) << "Cell: " << setw(8) 
	 << cell << endl;

  // energy
    cout << setw(15) << setiosflags(ios::right) << "Energy: " << setw(8) 
	 << energy << endl;
	
  // weight
    cout << setw(15) << setiosflags(ios::right) << "Energy: " << setw(8) 
	 << weight << endl;

    cout << endl;
}

template<class MT>
void Particle<MT>::Print_dead()
{
  // print dead particle (alive = false)
    cout << "Particle is dead." << endl;

  // coordinates
    cout << setw(20) << setiosflags(ios::right) << " Last Coordinates: ";
    for (int i = 0; i < r.size(); i++)
	cout << setw(8) << r[i] << " ";
    cout << endl;
    
  // direction
    cout << setw(20) << setiosflags(ios::right) << " Last Direction: ";
    for (int i = 0; i < omega.size(); i++)
	cout << setw(8) << omega[i] << " ";
    cout << endl;
}

//---------------------------------------------------------------------------//
// overloaded operators
//---------------------------------------------------------------------------//

CSPACE

//---------------------------------------------------------------------------//
//                              end of Particle.cc
//---------------------------------------------------------------------------//
