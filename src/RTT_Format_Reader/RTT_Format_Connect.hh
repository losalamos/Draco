//----------------------------------*-C++-*----------------------------------//
// Connectivity.hh
// B.T. Adams
// 7 June 00
/*! 
 * \file   RTT_Format_Reader/Connectivity.hh
 * \author B.T. Adams
 * \date   Wed Jun 7 10:33:26 2000
 * \brief  Header file for RTT_Format_Reader/Connectivity class.
 */
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __RTT_Format_Reader_Connectivity_hh__
#define __RTT_Format_Reader_Connectivity_hh__

#include "Dims.hh"
#include "CellDefs.hh"
#include "Nodes.hh"
#include "Sides.hh"
#include "Cells.hh"
#include <fstream>
#include <iostream>
#include <vector>
#include <set>
#include <map>

namespace rtt_RTT_Format_Reader
{
/*!
 * \brief Determines the mesh connectivity from the input mesh file data.
 */
class Connectivity
{
    // typedefs
    typedef std::set<int> set_int;
    typedef std::vector<int> vector_int;
    typedef std::vector<std::vector<int> > vector_vector_int;
    typedef std::vector<double> vector_dbl;

    const Dims & dims;
    const CellDefs & cellDefs;
    const Cells & cells;
    const Sides & sides;
    const Nodes & nodes;
    std::vector<std::vector<std::vector<int> > > adjCell;
    std::multimap<int, int> bndryFaces;
    std::multimap<int, vector_int > Side_to_Cell_Face;

  public:
    Connectivity(const Dims & dims_, const CellDefs & cellDefs_,
		 const Cells & cells_, const Sides & sides_, 
		 const Nodes & nodes_);
    ~Connectivity() {}

  private:
    void calcAdjacentCells();
  public:
/*!
 * \brief Returns the number of cells adjacent to the specified cell face.
 * \param cell Cell number.
 * \param face Face number.
 * \return The number of adjacent cells.
 */
    int get_adjCell_size(int cell, int face) const
    { return adjCell[cell][face].size(); }
/*!
 * \brief Returns the number of the cell adjacent to the specified cell, face,
 *        and optional adjacent cell index.
 * \param cell Cell number.
 * \param face Face number.
 * \param adjcell Adjacent cell number (defaults to 0).
 * \return The adjacent cell number.
 */
    int get_adjCell(int cell, int face, int adjcell = 0) const
    { return adjCell[cell][face][adjcell]; }
/*!
 * \brief Returns the number of boundary faces (i.e., faces that are either
 *        on the outer boundary of the problem geometry or a connection between
 *        cells with different refinement levels in an AMR mesh) with the 
 *        specified face number
 * \param face Face number.
 * \return The number of boundary faces.
 */
    int get_bndryFaces_count(int face) const { return bndryFaces.count(face); }

    set_int get_bndryCells(int face) const;
    bool check_bndryFace(int cell, int face) const;
    int get_Cell_from_Side(int side) const;
    int get_Cell_Face_from_Side(int side) const;
};

} // end namespace rtt_RTT_Format_Reader

#endif                          // __RTT_Format_Reader_Connect_hh__


//---------------------------------------------------------------------------//
//                end of RTT_Format_Reader/RTT_Format_Connect.hh
//---------------------------------------------------------------------------//
