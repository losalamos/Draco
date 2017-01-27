//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   RTT_Format_Reader/test/TestRTTPolyhedron.cc
 * \author James S. Warsa
 * \date   Fri Aug 19 12:31:47 2016
 * \brief  Testing the POLYHEDRON Element_Definition
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//

#include <sstream>

#include "ds++/Release.hh"
#include "ds++/SP.hh"
#include "ds++/ScalarUnitTest.hh"

#include "RTT_Format_Reader/RTT_Mesh_Reader.hh"

using rtt_RTT_Format_Reader::RTT_Mesh_Reader;
using rtt_mesh_element::Element_Definition;

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_RTT_Format_Reader;

class is_cell {
public:
  //! Construct the predicate for a specified dimensionality.
  is_cell(unsigned const dimensionality) : dimensionality(dimensionality) {}

  //! Returns \c true if the specified type index has the predicate
  //! dimensionality; \c false otherwise.
  bool operator()(Element_Definition::Element_Type const type) {
    return (static_cast<unsigned>(Element_Definition(type).get_dimension()) ==
            dimensionality);
  }

private:
  unsigned dimensionality;
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void test_polyhedron(rtt_dsxx::UnitTest &ut) {

  vector<string> filenames = {"rttquad5.mesh", "rttquad9.mesh", "rttquad.mesh"};

  for (unsigned i = 0; i < filenames.size(); ++i) {
    string filename(filenames[i]);
    SP<RTT_Mesh_Reader> mesh(new RTT_Mesh_Reader(filename));

    ostringstream m;
    m << "Read mesh file " << filename << std::endl;
    ut.passes(m.str());

    unsigned const ndim = mesh->get_dims_ndim();
    if (ndim != 2) {
      FAILMSG("Unexpected dimension.");
    }

    vector<SP<Element_Definition>> const element_defs(mesh->get_element_defs());
    for (size_t j = 0; j < element_defs.size(); ++j) {
      cout << "Element definition for element " << j << endl;
      element_defs[j]->print(cout);
    }
  }

  //-------------------------------------------------------------------------------------------------------------------------------//

  // Read polyhedron mesh file - this is the mesh that is of most interest in
  // this test
  {
    vector<string> filenames = {"rttpolyhedron.mesh", "rttpolyhedron.2.mesh",
                                "rttpolyhedron.2o.mesh", "rttpolyhedron.3.mesh",
                                "rttpolyhedron.4.mesh"};

    for (unsigned i = 0; i < filenames.size(); ++i) {
      string filename(filenames[i]);
      SP<RTT_Mesh_Reader> mesh(new RTT_Mesh_Reader(filename));

      ostringstream m;
      m << "Read mesh file " << filename << std::endl;
      ut.passes(m.str());

      // Investigate and report on the mesh

      // The element types begins with side types, followed by cell types.
      // We can distinguish these by their dimensionality. Cell types have
      // the full dimensionality of the mesh; we assume side types have one
      // less than the full dimensionality of the mesh.

      unsigned const ndim = mesh->get_dims_ndim();
      vector<Element_Definition::Element_Type> const element_types(
          mesh->get_element_types());
      vector<vector<int>> const element_nodes(mesh->get_element_nodes());
      map<string, set<int>> element_sets(mesh->get_element_sets());
      vector<SP<Element_Definition>> const element_defs(
          mesh->get_element_defs());

      if (ndim != 3) {
        FAILMSG("Unexpected dimension.");
      } else {
        PASSMSG("Correct dimension.");
      }

      for (size_t j = 0; j < element_defs.size(); ++j) {
        cout << "Element definition for element " << j << endl;
        element_defs[j]->print(cout);
      }

      unsigned const ncells =
          count_if(element_types.begin(), element_types.end(), is_cell(ndim));
      unsigned const mcells = mesh->get_dims_ncells();
      if (ncells != mcells) {
        FAILMSG("Unexpected number of cells.");
      } else {
        PASSMSG("Correct number of cells.");
        std::cout << " There are " << ncells << " cells in the mesh"
                  << std::endl;
      }

      unsigned const nsides =
          static_cast<unsigned>(element_types.size()) - ncells;
      unsigned const msides = mesh->get_dims_nsides();
      if (nsides != msides) {
        FAILMSG("Unexpected number of sides.");
      } else {
        PASSMSG("Correct number of sides.");
        std::cout << " There are " << nsides << " sides in the mesh"
                  << std::endl;
      }

      for (map<string, set<int>>::const_iterator it = element_sets.begin();
           it != element_sets.end(); ++it) {
        // See if this is a side or cell flag.

        unsigned const representative_element = *it->second.begin();

        Element_Definition::Element_Type const type_index =
            static_cast<Element_Definition::Element_Type>(
                element_types[representative_element]);

        SP<Element_Definition const> type(new Element_Definition(type_index));

        if (type->get_dimension() == ndim) {
          std::cout << " Elements with flags " + it->first
                    << " are cell elements " << std::endl;
        } else {
          // It has to be a side flag.
          if (type->get_dimension() + 1 != ndim) {
            throw invalid_argument(
                " Elements with flags " + it->first +
                " have the wrong dimension for a side element");
          } else {
            std::cout << " Elements with flags " + it->first
                      << " are side elements" << std::endl;
          }
        }
      }

      unsigned ncorner = 0;
      for (unsigned et = nsides; et < element_types.size(); et++) {
        unsigned const c = et - nsides;
        SP<Element_Definition> cell_def = element_defs[c];
        unsigned const number_of_nodes = cell_def->get_number_of_nodes();
        ncorner += number_of_nodes;
      }

      if (i == 0 && ncorner != 13) {
        FAILMSG("Number of corners should be 13 for mesh 1");
      } else if (i == 1 && ncorner != 17) {
        FAILMSG("Number of corners should be 17 for mesh 2");
      } else {
        PASSMSG("Number of corners is correct");
      }

    } // End polyhedron mesh tests
  }
}

//---------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    test_polyhedron(ut);
  }
  UT_EPILOG(ut);
}

//---------------------------------------------------------------------------//
// end of TestRTTPolyhedron.cc
//---------------------------------------------------------------------------//
