//----------------------------------*-C++-*----------------------------------//
// AMR_Layout.hh
// Thomas M. Evans
// Thu May 28 15:05:44 1998
//---------------------------------------------------------------------------//
// @> AMR_Layout implementation file
//---------------------------------------------------------------------------//

#ifndef __imc_AMR_Layout_hh__
#define __imc_AMR_Layout_hh__

//===========================================================================//
// class AMR_Layout - 
//
// Purpose : provide levelized layout necessary for AMR meshes
//
// revision history:
// -----------------
//  0) original
// 
//===========================================================================//

#include "Names.hh"
#include <iostream>
#include <vector>

IMCSPACE

using std::vector;
using std::ostream;

class AMR_Layout
{
private:
  // cell-face-cell data in the form of coarse face-refined face
    vector<vector<int> > face;
    
  // number of cells in this mesh
    int nc;

public:
  // inline default constructor
    inline AMR_Layout(int = 0, int = 6);

  // get size member functions
    int num_cells() const { return nc; }
    inline int num_faces(int, int) const;

  // diagnostic functions
    void print(ostream &, int) const;

  // overloaded subscripting operators for assignment and retrieval
    inline int operator()(int, int, int) const;
    inline int& operator()(int, int, int);

  // overloaded operators for equality
    inline bool operator==(const AMR_Layout &) const;
    bool operator!=(const AMR_Layout &rhs) const { return !(*this == rhs); }
};

//---------------------------------------------------------------------------//
// overloaded operators for AMR_Layout
//---------------------------------------------------------------------------//
// overload operator for stream output

ostream& operator<<(ostream &, const AMR_Layout &);

//---------------------------------------------------------------------------//
// overload equality (==) operator for design-by-contract

inline bool AMR_Layout::operator==(const AMR_Layout &rhs) const
{
  // if the data is equal, the AMR_Layouts are equal
    if (face == rhs.face && nc == rhs.nc)
	return true;

  // if we haven't returned then the AMR_Layouts aren't equal
    return false;
}

//---------------------------------------------------------------------------//
// INLINE functions
//---------------------------------------------------------------------------//
// constructor

AMR_Layout::AMR_Layout(int ncells, int coarse_faces)
    : face(ncells * coarse_faces), nc(ncells) {}

CSPACE

#endif                          // __imc_AMR_Layout_hh__

//---------------------------------------------------------------------------//
//                              end of imc/AMR_Layout.hh
//---------------------------------------------------------------------------//
