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
using Global::Dot;

//===========================================================================//
// class Particle<MT>
//===========================================================================//

//---------------------------------------------------------------------------//
// constructors
//---------------------------------------------------------------------------//
// defined inline

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
void Particle<MT>::Transport(const MT &mesh, const Opacity<MT> &xs, 
			     SP<Diagnostic> diagnostic)
{
  // explicit calls from standard library namespace
    using std::log;
  // explicit calls from Global namespace
    using Global::pi;

  // initialize diagnostics
    if (diagnostic)
    {
	diagnostic->Header();
	diagnostic->Print(*this);
    }
  
  // !!! BEGIN TRANSPORT LOOP !!!

  // transport loop, ended when alive = false
    while (alive)
    {
      // dist-to-collision and dist-to-boundary definitions
        double d_collision, d_boundary;
      // cell face definition
        int face = 0;
        
      // sample distance-to-collision
        d_collision = -log(random.Ran()) / xs.Sigma(cell);

      // get distance-to-boundary and cell face
        d_boundary  = mesh.Get_db(r, omega, cell, face);

      // detailed diagnostics
	if (diagnostic)
	    if (diagnostic->Detail())
	    {
		diagnostic->Print_dist(d_collision, d_boundary, cell);
		diagnostic->Print_xs(xs, cell);
	    }

      // streaming
        if (d_collision <= d_boundary)
        {
            Stream(d_collision);
            alive = Collide(mesh, xs);
        }
	else
        {
	  // stream to cell boundary and find next cell
            Stream(d_boundary);
	    alive = Surface(mesh, face);
	}

      // do diagnostic print
	if (diagnostic)
	    diagnostic->Print(*this);
    } 

  // !!! END OF TRANSPORT LOOP !!!
}

//---------------------------------------------------------------------------//
// private transport member functions
//---------------------------------------------------------------------------//
// calculate everything about a collision
template<class MT>
bool Particle<MT>::Collide(const MT &mesh, const Opacity<MT> &xs)
{   
  // status from collision
    bool status;

  // determine absorption or collision
    if (random.Ran() <= xs.Sigma(cell) / xs.Sigma(cell))
    {
	descriptor = "absorption";
        status = false;
    }
    else
    {
        status = true;
        
      // calculate theta and phi (isotropic)
        double costheta, phi;
        costheta = 1 - 2 * random.Ran();
        phi      = 2 * Global::pi * random.Ran();

      // get new direction cosines
        mesh.Coord().Calc_omega(costheta, phi, omega);
    }

  // return outcome of the event
    return status;
}

//---------------------------------------------------------------------------//
// do surface crossings
template<class MT>
bool Particle<MT>::Surface(const MT &mesh, int face)
{
    using Global::Dot;

  // status from surface crossing
    bool status;

  // determine the next cell
    int next_cell = mesh.Next_cell(cell, face);

  // determine descriptor and outcome of this event

    if (next_cell == cell)
    {
      // reflection
	descriptor = "reflection";
	vector<double> normal = mesh.Get_normal(cell, face);
	double factor = Dot(omega, normal);
	for (int i = 0; i < mesh.Coord().Get_sdim(); i++)
	    omega[i] -= 2 * factor * normal[i];
	cell = next_cell;
    }
    else if (next_cell == 0)
    {
      // escape
	descriptor = "escape";
	cell = next_cell;
    }
    else 
    {
      // continue streaming
	descriptor = "stream";
	cell = next_cell;
    }

  // return outcome of the event
    status = cell;
    return status;
}

//---------------------------------------------------------------------------//
// overloaded operators
//---------------------------------------------------------------------------//

//===========================================================================//
// class Particle<MT>::Diagnostic
//===========================================================================//

//---------------------------------------------------------------------------//
// public diagnostic member functions
//---------------------------------------------------------------------------//
template<class MT>
void Particle<MT>::Diagnostic::Print(const Particle<MT> &particle) const
{
  // set output precision
    output.precision(3);
    output << setiosflags(ios::fixed);

  // print particulars of the particle based on its status
    if (particle.alive == true)
	Print_alive(particle);
    else
	Print_dead(particle);
}

template<class MT>
void Particle<MT>::Diagnostic::Print_alive(const Particle<MT> &particle) const 
{
  // print active particle (alive = true)
    output << " -- Particle is alive -- " << endl;

  // event
    output << setw(20) << setiosflags(ios::right) << "Event: " 
	   << setw(12) << particle.descriptor.c_str() << endl;
    
  // coordinates
    output << setw(20) << setiosflags(ios::right) << "Coordinates: ";
    for (int i = 0; i < particle.r.size(); i++)
	output << setw(12) << particle.r[i] << " ";
    output << endl;
    
  // direction
    output << setw(20) << setiosflags(ios::right) << "Direction: ";
    for (int i = 0; i < particle.omega.size(); i++)
	output << setw(12) << particle.omega[i] << " ";
    output << endl;
    
  // cell
    output << setw(20) << setiosflags(ios::right) << "Cell: " << setw(12) 
	   << particle.cell << endl;
    
  // energy
    output << setw(20) << setiosflags(ios::right) << "Energy: " << setw(12) 
	   << particle.energy << endl;
    
  // weight
    output << setw(20) << setiosflags(ios::right) << "Weight: " << setw(12) 
	   << particle.weight << endl;
    
    output << endl;
}

template<class MT>
void Particle<MT>::Diagnostic::Print_dead(const Particle<MT> &particle) const
{
  // print dead particle (alive = false)
    output << " -- Particle is dead -- " << endl;

  // event
    output << setw(20) << setiosflags(ios::right) << "Event: " 
	   << setw(12) << particle.descriptor.c_str() << endl;
    
  // coordinates
    output << setw(20) << setiosflags(ios::right) << " Last Coordinates: ";
    for (int i = 0; i < particle.r.size(); i++)
	output << setw(12) << particle.r[i] << " ";
    output << endl;
    
  // direction
    output << setw(20) << setiosflags(ios::right) << " Last Direction: ";
    for (int i = 0; i < particle.omega.size(); i++)
	output << setw(12) << particle.omega[i] << " ";
    output << endl;

  // cell
    output << setw(20) << setiosflags(ios::right) << "Last Cell: " 
	   << setw(12) << particle.cell << endl;
    
  // energy
    output << setw(20) << setiosflags(ios::right) << "Last Energy: "
	   << setw(12) << particle.energy << endl;
    
  // weight
    output << setw(20) << setiosflags(ios::right) << "Last Weight: " 
	   << setw(12) << particle.weight << endl;

    output << endl;
}

template<class MT>
void Particle<MT>::Diagnostic::Print_dist(double d_col, double d_bnd, 
					  int cell) const
{
  // do detailed diagnostic print of particle event distances
    output << setw(20) << setiosflags(ios::right) << "Present cell: "
	   << setw(12) << cell << endl;
    output << setw(20) << setiosflags(ios::right) << "Dist-collision: "
	   << setw(12) << d_col << endl;
    output << setw(20) << setiosflags(ios::right) << "Dist-boundary: "
	   << setw(12) << d_bnd << endl;
}

template<class MT>
void Particle<MT>::Diagnostic::Print_xs(const Opacity<MT> &xs,
					int cell) const
{
  // do detailed diagnostic print of particle event cross sections
    output << setw(20) << setiosflags(ios::right) << "Opacity: " 
	   << setw(12) << xs.Sigma(cell) << endl;
}

CSPACE

//---------------------------------------------------------------------------//
//                              end of Particle.cc
//---------------------------------------------------------------------------//
