//----------------------------------*-C++-*----------------------------------//
/*! 
 * \file   RTT_Format_Reader/Release.hh
 * \author B.T. Adams
 * \date   Wed June 7 10:33:26 2000
 * \brief  Header file for RTT_Format_Reader library release function.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __RTT_Format_Reader_Release_hh__
#define __RTT_Format_Reader_Release_hh__

//===========================================================================//
// namespace version - 
//
// Purpose : Return the version of RTT_Format_Reader; this can be used to get 
///          exact version information in codes that use RTT_Format_Reader
// 
//===========================================================================//

#include <string>

/*!
 * \brief Namespace to contain the RTT_Format_Reader utilities.
 *
 * Provides namespace protection for the Draco RTT_Format_Reader utilities.
 *
 *\sa The RTT_Format_Reader class constructor automatically instantiates and 
 *    executes the readMesh member function used to parse the mesh data. 
 *    Accessor functions are provided for all of the remaining member classes 
 *    to allow data retrieval. The \ref rtt_mesh_reader_overview page presents
 *    a summary of the capabilities provided by the namespace.
 */
namespace rtt_RTT_Format_Reader 
{
/*!
 * \brief  Gets the release number for the RTT_Format_Reader package. 
 * \return release number as a string in the form "RTT_Format_Reader-\#_\#_\#"
 */
const std::string release();

}  // end of rtt_RTT_Format_Reader namespace

#endif                          // __RTT_Format_Reader_Release_hh__

/*!
 * \page rtt_mesh_reader_overview Overview of the RTT_Format_Reader class
 *
 * \version 1_0_0
 *
 * <h3> Introduction </h3>
 * The RTT_Format_Reader class consists of member functions that are used to 
 * parse a mesh file in the \ref rtt_format_defined and to access the data.
 *
 * <h3> Intended Usage </h3>
 * The RTT_Format_Reader class constructor automatically parses the specified 
 * input file via a call to the private member functions readMesh. The mesh 
 * data can then be accessed using the public member accessor functions. The 
 * RTT_Format_Reader class contains several data members that are classes
 * corresponding to the organization of the data blocks in the \ref 
 * rtt_format_defined, with the addition of two member data classes:
 * <ul>
 *  <li> rtt_RTT_Format_Reader::Header
 *  <li> rtt_RTT_Format_Reader::Dims (dimensions)
 *  <li> rtt_RTT_Format_Reader::Flags (member data class of NodeFlags, 
 *                                     SideFlags, and CellFlags)
 *  <li> rtt_RTT_Format_Reader::NodeFlags
 *  <li> rtt_RTT_Format_Reader::SideFlags
 *  <li> rtt_RTT_Format_Reader::CellFlags
 *  <li> rtt_RTT_Format_Reader::NodeDataIDs
 *  <li> rtt_RTT_Format_Reader::SideDataIDs
 *  <li> rtt_RTT_Format_Reader::CellDataIDs
 *  <li> rtt_RTT_Format_Reader::CellDef (member data class of CellDefs)
 *  <li> rtt_RTT_Format_Reader::CellDefs (cell definitions)
 *  <li> rtt_RTT_Format_Reader::Nodes
 *  <li> rtt_RTT_Format_Reader::Sides
 *  <li> rtt_RTT_Format_Reader::Cells
 *  <li> rtt_RTT_Format_Reader::NodeData
 *  <li> rtt_RTT_Format_Reader::SideData
 *  <li> rtt_RTT_Format_Reader::CellData
 *  <li> rtt_RTT_Format_Reader::Connectivity
 * </ul> 
 * These classes provide a convenient grouping of the mesh data, and the 
 * RTT_Format_Reader public accessor member functions reflect the name of the 
 * associated class. Alternatively, the provided RTT_Mesh_Reader class is a 
 * derived type of the DRACO meshReaders package which utilizes the 
 * RTT_Format_Reader directly. A standard mesh reader interface has been 
 * specified for the meshReaders package and, thus, future packages should
 * incorporate this interface rather than the RTT_Format_Reader interface.
 *
 */

/*!
 * \page rtt_format_defined RTT Format File Structure
 * The following example "mesh" documents the format of the RTT file and 
 * explains the associated nomenclature. A graphical depiction of the \ref 
 * rtt_stdcell is  provided via the links.
 *
 * \include RTT_Format.defined
 */

/*!
 * \page rtt_stdcell RTT Format ICEM/DDN Cell Definitions
 * The RTT_Formatt_Reader side set numbering output by ICEM/DDN is depicted on
 * this page. Note that the "right hand rule" is used to return the direction
 * of the outward-directed normal when the nodes are traversed in the order 
 * that is specified in the side set node ordering. The RTT_Format cell 
 * definitions do not assume any particular orientation of the sides relative
 * to the problem coordinate system.
 *
 * <center>
 *   <table>
 *     <tr>
 *       <td align=center valign=center>
 *         <img src="../../draco/src/RTT_Format_Reader/doc/stdcell.jpg"> 
 *       </td>
 *     </tr>
 *   </table>
 * </center>
 *
 */

//---------------------------------------------------------------------------//
//                              end of RTT_Format_Reader/Release.hh
//---------------------------------------------------------------------------//
