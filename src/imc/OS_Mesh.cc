//----------------------------------*-C++-*----------------------------------//
// OS_Mesh.cc
// Thomas M. Evans
// Tue Feb  3 16:50:13 1998
//---------------------------------------------------------------------------//
// @> OS_Mesh class implementation file
//---------------------------------------------------------------------------//

#include "OS_Mesh.hh"
#include "Global.hh"
#include <iostream>
#include <iomanip>

IMCSPACE

//---------------------------------------------------------------------------//
// constructors
//---------------------------------------------------------------------------//
// defined inline

//---------------------------------------------------------------------------//
// member functions
//---------------------------------------------------------------------------//
int OS_Mesh::getCell(vector<double> &r) const
{
  // used variables
    int dim = coord->getDim();
    int return_cell = 1;
    int subcells = 1;
    
  // binary search of cells

    for (int i = 0; i < dim; i++)
    {
	int low_index  = 0;
	int high_index = sur[i].size() - 1;
	int index;
	while ( (high_index - low_index) != 1)
	{
	    index = (high_index + low_index) / 2;
	    if (r[i] < sur[i][index])
		high_index = index;
	    else
		low_index  = index;
	}
	return_cell += subcells * (high_index - 1);
	subcells    *= sur[i].size() - 1;
    }  
    
  // return cell index
    return return_cell;
}

double OS_Mesh::getDb(vector<double> &r, vector<double> &omega,
		      int cell, int &face) const
{
    using std::vector;
    using std::min_element;
    using Global::huge;
    
  // calculate distance to the vec(r) boundaries

  // boundary distances over each coordinate direction
    vector<double> dim_dist_boundary(coord->getDim(), 0.0);
    
  // loop to get the distances to boundary in each coordinate direction
    for (int i = 0; i < coord->getDim(); i++)
    {
	if (omega[i] == 0.0)
	    dim_dist_boundary[i] = Global::huge;
	else if (omega[i] > 0.0)
	    dim_dist_boundary[i] = ((pos[i][cell-1] + dim[i][cell-1]/2.0) - 
				    r[i]) / omega[i];
	else
	    dim_dist_boundary[i] = ((pos[i][cell-1] - dim[i][cell-1]/2.0) -
				    r[i]) / omega[i];
    }

  // calculate the distance to boundary
    vector<double>::iterator itor = min_element(dim_dist_boundary.begin(),
						dim_dist_boundary.end());
    double dist_boundary = *itor;

  // calculate the face that the boundary is on
    int index = itor - dim_dist_boundary.begin();
    if (omega[index] < 0.0)
	face = 1 + 2 * index;
    else
	face = 2 + 2 * index;

  // return the distance-to-boundary
    return dist_boundary;
}

void OS_Mesh::Print(int cell) const
{
    using std::cout;
    using std::endl;
    using std::setw;

  // print out content info for 1 cell
    cout << "+++++++++++++++" << endl;
    cout << "Cell        : " << cell << endl;
    cout << "---------------" << endl;
    cout << "Dimensions  : " << endl;
    if (coord->getDim() == 2)
    {
	cout << " x  : " << pos[0][cell-1] << endl;
	cout << " y  : " << pos[1][cell-1] << endl;
    	cout << " dx : " << dim[0][cell-1] << endl;
	cout << " dy : " << dim[1][cell-1] << endl;
    }
    else
    {
	cout << " x  : " << pos[0][cell-1] << endl;
	cout << " y  : " << pos[1][cell-1] << endl;
	cout << " z  : " << pos[2][cell-1] << endl;
    	cout << " dx : " << dim[0][cell-1] << endl;
	cout << " dy : " << dim[1][cell-1] << endl;
	cout << " dz : " << dim[2][cell-1] << endl;
    }	
    cout << "---------------" << endl;
    cout << "Layout      :  " << endl;
    layout.print(cell);
    cout << "+++++++++++++++" << endl;
}

CSPACE

//---------------------------------------------------------------------------//
//                              end of OS_Mesh.cc
//---------------------------------------------------------------------------//
