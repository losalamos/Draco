//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   meshReaders/Element_Definition.hh
 * \author John McGhee
 * \date   Fri Feb 25 10:03:18 2000
 * \brief  Header file for the RTT Element_Definition class.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __meshReaders_Element_Definition_hh__
#define __meshReaders_Element_Definition_hh__

#include <vector>
#include <string>
#include <stdexcept>

namespace rtt_meshReaders
{
 
//===========================================================================//
/*!
 * \class Element_Definition
 *
 * \brief Provides some descriptive information on the standard mesh 
 *        elements used in the RTT meshReader class.
 *
 * A few high points, trying not to wax eloquent. It was originally desired
 * to create a simple class the would concisely, unambiguously, and 
 * completely describe any mesh element that could be conceived. While
 * this may be a laudable goal, it appears to be harder than it appears. Perhaps
 * we could get some help on this from some computational geometry experts
 * at some time. In the mean time here is my 80% solution.
 *
 *  First, we will
 * reduce the scope from any element to just the ones currently described
 * in the CGNS standard. Remember that it is only necessary to describe
 * the problem "geometry" with these elements, not the solution, or any
 * other field on the mesh, so this may not be as much of a restriction
 * as it first appears. This set consists of  18 elements including most of the
 * commonly used ones. Moreover, remember we currently have no means of
 * generating a mesh with weird custom elements. Any mesh anyone in the
 * group has ever run on can be expressed with just six of these elements.
 *
 * Second, we will not try to design a completely general element
 * description, but will settle for providing a limited set of services that can
 * be used to discover a lot of things about the elements in the CGNS sub-set,
 * but may not necessarily be a universal, complete, and unambiguous 
 * description. 
 * The ultimate authority on the element descriptions are the figures and
 * text found in the CGNS SIDS-Additions manual.
 *
 * The description implemented herein utilizes a hierarchical approach. 3D 
 * elements
 * are described as assemblies of 2D elements, which are composed of 1D
 * elements, which are themselves composed of nodes. For examples a hexahedra
 * is described in terms of its quadrilateral faces, which are described in 
 * terms of line edge elements, which are then described in terms of their
 * constituent nodes.  This approach appears to be adequate for
 * the subset of elements under consideration herein, but it is not clear
 * that this will suffice in the general case.
 *
 * Utilities are provided to inquired about the type of a face (i.e. quad,
 * triangle, etc....) as well as the nodes that compose the face.
 *
 * In addition to face types, there is a concept of "node-location" within
 * the element. Thus all nodes are given a location (i.e. "CORNER", "EDGE",
 * etc....) to aide in the description of the element. Again this appears 
 * to be adequate for the sub-set of elements under consideration herein
 * but may not be adequate in a more general case.
 *
 * It is hoped that the node-location, together with the data available
 * through recursively descending through element faces, edges, and nodes 
 * provides an adequate amount of information for our present needs. However,
 * it is difficult to show that this description is complete and unambiguous.
 *
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//


class Element_Definition 
{

    // NESTED CLASSES AND TYPEDEFS

  public:

    /*!
     * \brief Describes the location of a node within an element.
     *
     * For the purposes of this enumeration, the terms have the
     * following meaning: A "corner" node terminates one or more edges
     *  of an element.
     * The term "edge" describes a node that lies on the interior of
     * a 1D element. The term "face" describes a node that lies on the
     * interior of a 2D element. Finally the term "cell" connotates
     * a node that lies on the interior of a 3D element.
     *
     * All elements will
     * always have corner nodes. In addition, all elements may have
     * edge nodes. Two and three dimensional elements may also have face
     * nodes, and finally, three-dimensional elements may have cell
     * nodes. Under these definitions, note that a node's location is
     * unchanged in an element and all its sub-elements. i.e. the corner
     * nodes of a quadrilateral are also corner nodes in the line elements
     * which form the edges of the quadrilateral.
     *
     */
    enum Node_Location {
	CORNER,
	EDGE,
	FACE,
	CELL
    };
    
    /*!
     * \brief Standard element names. 
     *
     * These names and the elements that they
     * represent are the same as those defined in the CGNS SIDS Manual.
     *
     */
    enum Element_Type {
	NODE,
	BAR_2,
	BAR_3,
	TRI_3,
	TRI_6,
	QUAD_4,
	QUAD_8,
	QUAD_9,
	TETRA_4,
	TETRA_10,
	PYRA_5,
	PYRA_14,
	PENTA_6,
	PENTA_15,
	PENTA_18,
	HEXA_8,
	HEXA_20,
	HEXA_27
    };

  private:
    
    // DATA

    std::string name;
    Element_Type type;
    int dimension;
    int number_of_nodes;
    int number_of_sides;
    std::vector<Element_Definition> elem_defs;
    std::vector<int> side_type;
    std::vector<std::vector<int> > side_nodes;
    std::vector<Node_Location> node_loc;

  public:

    // CREATORS
    
    /*!
     * \brief Constructor for the Element_Definition class.
     *
     * \param name_ The element type to be constructed.
     */
    Element_Definition(const Element_Type &name_);

    // Defaulted
    // Element_Definition(const Element_Definition &rhs)

    /*!
     * \brief Destructor for the Element_Definition class.
     */
    ~Element_Definition();
    //
    // Empty
    //

    // MANIPULATORS
    
    // Defaulted
    // Element_Definition& operator=(const Element_Definition &rhs);

    // ACCESSORS

    /*!
     * \brief Returns the name of an element.
     * \return The element name.
     */
    std::string get_name() const
    {
	return name;
    }
    /*! 
     * \brief Returns the total number of nodes in an element.
     * \return Total number of nodes in an element.
     */
    int get_number_of_nodes() const
    {
	return number_of_nodes;
    }
    /*!
     * \brief Returns the dimension of an element. i.e. nodes return 0,
     * lines return 1, quads return 2, hexahedra return 3.
     *
     * \return The element dimension (0, 1, 2, or 3).
     */
    int get_dimension() const
    {
	return dimension;
    }
    /*!
     * \brief Returns the number of sides on an element.
     *
     * \return The number of n-1 dimensional entities that compose an n
     *         dimensional element. i.e. nodes return 0, lines return 2,
     *        quads return 4, hexahedra return 6.
     */
    int get_nsides() const
    {
	return number_of_sides;
    }
    /*!
     * \brief Returns the location of a node within the element.
     *
     * \param node_number the node number for which a location is
     *  desired. Node numbers must be in the range [0:number_of_nodes).
     *
     * \return The location of the node. See the 
     * Element_Definition::Node_Location enumeration
     * for additional discussion on node locations.
     *
     */
    Node_Location get_node_loc(const int node_number) const
    {
       if (node_number >= 0 && node_number < number_of_nodes) 
	   return node_loc[node_number];
       else
	   throw std::runtime_error("Node index out of range!");
    }

    /*!
     * \brief Returns the type (i.e. quad, tri, etc.) of a specified
     * element side.
     *
     * \return Returns a valid element definition that
     *  describes a element side. Can be queried using any of 
     *  the accessors provided in the Element_Definition class.

     * \param side_number Side number for which a type is desired.
     *   Side numbers are in the range [0:number_of_sides).
     *
     * Note that there is no valid side number for a "NODE" element.
     * "Side" in the context of this method means the
     * (n-1) dimensional element that composes a n dimensional
     * element.
     */
    Element_Definition get_side_type(const int side_number) const
    {
	if ( side_number >= 0 && side_number < side_type.size() )
	{
	    return elem_defs[side_type[side_number] ];
	}
	else 
	{
	    throw std::runtime_error("Side index out of range!");
	}
    }


    /*!
     * \brief Returns a vector of node numbers that are associated
     * with a particular element side.
     *
     * \param side_number The number of the element side for which
     *        the nodes are desired. Side numbers are in the 
     *         range [0:number_of_sides)
     *
     * \return A vector of the nodes associated with the side.
     *
     * "Side" in the context of this method means the
     * (n-1) dimensional element that composes a (n) dimensional
     * element. For example, on a hexahedra, a side is a quadrilateral,
     * whereas, on a quadrilateral a side is a line element.
     * The returned order of the side nodes is significant. The side-node
     * numbers are returned in the following order based on node
     * location: (corners, edges, faces, cells). For sides which are faces
     * of 3D elements, the vector cross product of the vector from
     * (side-node1 to side-node2) with the vector from (side-node1 to
     * side- node3) results in a vector that is oriented outward from
     * the parent element.
     * Equivalently, the side
     * corner-nodes are listed sequentially in a counter-clockwise
     * direction when viewed from outside the element. Both corner 
     * and edge nodes
     * are returned in a sequential order as one progresses around
     * a side. Moreover, the corner and edge nodes are returned so that
     * edge-node1 lies between corner-node1 and corner-node2, etc., etc.
     *
     * For sides which are edges of 2D elements, the vector cross 
     * product of the vector from
     * (side-node1 to side-node2) with a vector pointing towards the
     * observer results in a vector that is oriented outward from the 
     * parent element.
     *
     * Note that there is no valid side number for a "NODE" element.
     *
     */
    std::vector<int> get_side_nodes(const int side_number) const
    {
	if (side_number >= 0 && side_number < side_nodes.size() )
	{
	    return side_nodes[side_number];
	}
	else
	{
	    throw std::runtime_error("Side index out of range!");
	}
    }

    /*!
     * \brief Performs some simple sanity checks on the private data
     *        of the Element_Description class.
     */
    bool invariant_satisfied() const;

  private:
    
    // IMPLEMENTATION

    void construct_node();
    
    void construct_bar();
    
    void construct_tri();
    
    void construct_quad();
    
    void construct_tetra();

    void construct_pyra();

    void construct_penta();

    void construct_hexa();

};

} // end namespace rtt_meshReaders

#endif                          // __meshReaders_Element_Definition_hh__

//---------------------------------------------------------------------------//
//                              end of meshReaders/Element_Definition.hh
//---------------------------------------------------------------------------//
