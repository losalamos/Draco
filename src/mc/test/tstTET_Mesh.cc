//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/test/tstTET_Mesh.cc
 * \author H. Grady Hughes
 * \date   Thu Jan 27 16:34:07 MST 2000
 * \brief  Tests for tetrahedral meshes and ThreeVectors.
 */
//---------------------------------------------------------------------------//

// revision history:
// -----------------
//  0) original   : Committed 2000-01-28
//  1) 2000-02-07 : Added many tests of TET_Mesh functions, for a "hand-coded"
//                  test mesh.
//  2) 2000-02-09 : Extension to include instantiation from a generic interface
//                  class, and tests based on a particular example of such an
//                  interface class, contained in TET_test_1.hh.
//  3) 2000-04-10 : Rewritten to be consistent with new meshReader classes.

#include "../TET_Mesh.hh"
#include "../TET_Builder.hh"
#include "../Layout.hh"
#include "../XYZCoord_sys.hh"
#include "../Release.hh"
#include "TET_test_1.hh"
#include "c4/global.hh"
#include "rng/Sprng.hh"
#include "ds++/SP.hh"
#include "meshReaders/RTT_Format.hh"

#include <iomanip>
#include <vector>
#include <cmath>
#include <string>

using namespace std;

using rtt_mc::Coord_sys;
using rtt_mc::XYZCoord_sys;
using rtt_mc::Layout;
using rtt_mc::TET_Mesh;
using rtt_mc::TET_Builder;
using rtt_mc::ThreeVector;
using rtt_mc::sample_in_triangle;
using rtt_rng::Sprng;
using rtt_dsxx::SP;
using rtt_mc_test::TET_test_1;
using rtt_meshReaders::RTT_Format;

//! Typedef for scalar field of ThreeVectors.
typedef std::vector<ThreeVector> SF_THREEVECTOR;

//! Typedef for scalar field of doubles.
typedef std::vector<double> SF_DOUBLE;

//! Typedef for vector field of integers.
typedef std::vector< std::vector<int> > VF_INT;

double MID_epsilon = 0.0001;
int seed = 493875348;
int num = 5;

bool passed = true;
#define ITFAILS passed = rtt_mc_test::fail(__LINE__);

//! Mesh proxy class
class Mesh_Proxy
{
  private:
    SP<TET_Mesh> mesh;
  public:
    Mesh_Proxy(SP<TET_Mesh> m) : mesh(m) {}
    const TET_Mesh& get_Mesh() const { return *mesh; }
};

// ThreeVector Tests
void Test_ThreeVector()
{
    // Test sampling in a triangle.

    int *idr = init_sprng(0, num, seed, 1);
    Sprng random(idr, 0);

    int numm = 10000;
    double factor = 36.0/static_cast<double>(numm);

    ThreeVector A(100.0,10.0,0.0);
    ThreeVector B(112.0,10.0,0.0);
    ThreeVector C(106.0,16.0,0.0);

    double count[6][12];
    for (int y = 0; y < 6 ; y++)
        for (int x = 0; x < 12 ; x++)
            count[y][x] = 0.;

    double x_cent = 0.0;
    double y_cent = 0.0;

    for (int i = 0 ; i < numm ; i++)
    {
        ThreeVector R = sample_in_triangle(A,B,C,random);

        x_cent += R.get_x();
        y_cent += R.get_y();

        if (R.get_x() <= 100.0 || R.get_x() >= 112.0)       ITFAILS;
        if (R.get_y() <= 10.0  || R.get_y() >= 16.0)        ITFAILS;
        if (R.get_z() != 0.0)                               ITFAILS;

        int xx = static_cast<int>(R.get_x() - 100.);
        int yy = static_cast<int>(R.get_y() - 10.);
        count[yy][xx] += 1.;
    }
    x_cent /= static_cast<double>(numm);
    y_cent /= static_cast<double>(numm);

    if (fabs(x_cent - 105.99714295) > MID_epsilon)  ITFAILS;
    if (fabs(y_cent - 12.017748255) > MID_epsilon)  ITFAILS;

    for (int yy = 5 ; yy >= 0 ; yy--)
        for (int xx = 0; xx < 12; xx++)
            count[yy][xx] *= factor;

    if (count[5][0] != 0.0)                        ITFAILS;
    if (count[5][1] != 0.0)                        ITFAILS;
    if (count[5][2] != 0.0)                        ITFAILS;
    if (count[5][3] != 0.0)                        ITFAILS;
    if (count[5][4] != 0.0)                        ITFAILS;
    if (fabs(count[5][5] - 0.5544) > MID_epsilon)  ITFAILS;
    if (fabs(count[5][6] - 0.5040) > MID_epsilon)  ITFAILS;
    if (count[5][7] != 0.0)                        ITFAILS;
    if (count[5][8] != 0.0)                        ITFAILS;
    if (count[5][9] != 0.0)                        ITFAILS;
    if (count[5][10] != 0.0)                       ITFAILS;
    if (count[5][11] != 0.0)                       ITFAILS;
    if (count[4][0] != 0.0)                        ITFAILS;
    if (count[4][1] != 0.0)                        ITFAILS;
    if (count[4][2] != 0.0)                        ITFAILS;
    if (count[4][3] != 0.0)                        ITFAILS;
    if (fabs(count[4][4] - 0.4968) > MID_epsilon)  ITFAILS;
    if (fabs(count[4][5] - 1.0476) > MID_epsilon)  ITFAILS;
    if (fabs(count[4][6] - 1.0368) > MID_epsilon)  ITFAILS;
    if (fabs(count[4][7] - 0.5040) > MID_epsilon)  ITFAILS;
    if (count[4][8] != 0.0)                        ITFAILS;
    if (count[4][9] != 0.0)                        ITFAILS;
    if (count[4][10] != 0.0)                       ITFAILS;
    if (count[4][11] != 0.0)                       ITFAILS;
    if (count[3][0] != 0.0)                        ITFAILS;
    if (count[3][1] != 0.0)                        ITFAILS;
    if (count[3][2] != 0.0)                        ITFAILS;
    if (fabs(count[3][3] - 0.5220) > MID_epsilon)  ITFAILS;
    if (fabs(count[3][4] - 0.9792) > MID_epsilon)  ITFAILS;
    if (fabs(count[3][5] - 0.9396) > MID_epsilon)  ITFAILS;
    if (fabs(count[3][6] - 0.9504) > MID_epsilon)  ITFAILS;
    if (fabs(count[3][7] - 1.0260) > MID_epsilon)  ITFAILS;
    if (fabs(count[3][8] - 0.5184) > MID_epsilon)  ITFAILS;
    if (count[3][9] != 0.0)                        ITFAILS;
    if (count[3][10] != 0.0)                       ITFAILS;
    if (count[3][11] != 0.0)                       ITFAILS;
    if (count[2][0] != 0.0)                        ITFAILS;
    if (count[2][1] != 0.0)                        ITFAILS;
    if (fabs(count[2][2] - 0.4932) > MID_epsilon)  ITFAILS;
    if (fabs(count[2][3] - 0.9972) > MID_epsilon)  ITFAILS;
    if (fabs(count[2][4] - 1.0800) > MID_epsilon)  ITFAILS;
    if (fabs(count[2][5] - 1.0872) > MID_epsilon)  ITFAILS;
    if (fabs(count[2][6] - 1.0404) > MID_epsilon)  ITFAILS;
    if (fabs(count[2][7] - 0.9324) > MID_epsilon)  ITFAILS;
    if (fabs(count[2][8] - 1.0368) > MID_epsilon)  ITFAILS;
    if (fabs(count[2][9] - 0.5508) > MID_epsilon)  ITFAILS;
    if (count[2][10] != 0.0)                       ITFAILS;
    if (count[2][11] != 0.0)                       ITFAILS;
    if (count[1][0] != 0.0)                        ITFAILS;
    if (fabs(count[1][1] - 0.4752) > MID_epsilon)  ITFAILS;
    if (fabs(count[1][2] - 1.0116) > MID_epsilon)  ITFAILS;
    if (fabs(count[1][3] - 1.0044) > MID_epsilon)  ITFAILS;
    if (fabs(count[1][4] - 0.9684) > MID_epsilon)  ITFAILS;
    if (fabs(count[1][5] - 1.0188) > MID_epsilon)  ITFAILS;
    if (fabs(count[1][6] - 0.9324) > MID_epsilon)  ITFAILS;
    if (fabs(count[1][7] - 1.0080) > MID_epsilon)  ITFAILS;
    if (fabs(count[1][8] - 1.0044) > MID_epsilon)  ITFAILS;
    if (fabs(count[1][9] - 0.8424) > MID_epsilon)  ITFAILS;
    if (fabs(count[1][10] - 0.5076) > MID_epsilon) ITFAILS;
    if (count[1][11] != 0.0)                       ITFAILS;
    if (fabs(count[0][0] - 0.5184) > MID_epsilon)  ITFAILS;
    if (fabs(count[0][1] - 0.9216) > MID_epsilon)  ITFAILS;
    if (fabs(count[0][2] - 0.9180) > MID_epsilon)  ITFAILS;
    if (fabs(count[0][3] - 1.0872) > MID_epsilon)  ITFAILS;
    if (fabs(count[0][4] - 1.0080) > MID_epsilon)  ITFAILS;
    if (fabs(count[0][5] - 0.9900) > MID_epsilon)  ITFAILS;
    if (fabs(count[0][6] - 0.9540) > MID_epsilon)  ITFAILS;
    if (fabs(count[0][7] - 1.0188) > MID_epsilon)  ITFAILS;
    if (fabs(count[0][8] - 1.0044) > MID_epsilon)  ITFAILS;
    if (fabs(count[0][9] - 1.0548) > MID_epsilon)  ITFAILS;
    if (fabs(count[0][10] - 0.9936) > MID_epsilon) ITFAILS;
    if (fabs(count[0][11] - 0.4608) > MID_epsilon) ITFAILS;

//  Beginning:  a way to see these results ordered on the page.
//  cout.setf(ios::fixed);
//  cout.precision(5);
//
//  for (int yy = 5 ; yy >= 0 ; yy--) {
//      for (int xx = 0; xx < 12; xx++)
//          cout << setw(10) << count[yy][xx];
//      cout << endl;
//  }
//  cout << endl << endl;
//  End:  a way to see these results ordered on the page.

    // End of test sampling in a triangle.

}   // end Test_ThreeVector()

// TET_Mesh Tests
void Test_TET()
{
    double TET_epsilon = 0.0000000001;

    // First, test a hand-built mesh, without using an interface.
    // This mesh is a 2-tet pyramid.

    SP<Coord_sys> coor(new XYZCoord_sys());
    if (!coor)                     ITFAILS;

    SF_THREEVECTOR vertex_vec;
    vertex_vec.push_back(ThreeVector(0.0,0.0,0.0));
    vertex_vec.push_back(ThreeVector(1.0,0.0,0.0));
    vertex_vec.push_back(ThreeVector(0.0,1.0,0.0));
    vertex_vec.push_back(ThreeVector(1.0,1.0,0.0));
    vertex_vec.push_back(ThreeVector(0.5,0.5,1.0));
    if (vertex_vec.size() != 5)            ITFAILS;

    VF_INT cells_ver(2);

    cells_ver[0].push_back(0);
    cells_ver[0].push_back(1);
    cells_ver[0].push_back(2);
    cells_ver[0].push_back(4);

    cells_ver[1].push_back(2);
    cells_ver[1].push_back(1);
    cells_ver[1].push_back(3);
    cells_ver[1].push_back(4);

    if (cells_ver.size() != 2)      ITFAILS;
    if (cells_ver[0].size() != 4)   ITFAILS;
    if (cells_ver[1].size() != 4)   ITFAILS;

    Layout layo;
    layo.set_size(2);
    layo.set_size(1,4);
    layo.set_size(2,4);

    for (int c = 1 ; c <=2 ; c++)
        for (int f = 1 ; f <= 4 ; f++)
            layo(c,f) = 0;

    layo(1,1) = 2;
    layo(2,3) = 1;

    SP<TET_Mesh> mesh_ptr_H(new TET_Mesh(coor, layo, vertex_vec, cells_ver));

    if ( !mesh_ptr_H )                                                 ITFAILS;
    if ( !mesh_ptr_H->full_Mesh() )                                    ITFAILS;

    if ( mesh_ptr_H->next_cell(1,1) != 2 )                             ITFAILS;
    if ( mesh_ptr_H->next_cell(1,2) != 0 )                             ITFAILS;
    if ( mesh_ptr_H->next_cell(1,3) != 0 )                             ITFAILS;
    if ( mesh_ptr_H->next_cell(1,4) != 0 )                             ITFAILS;

    if ( mesh_ptr_H->next_cell(2,1) != 0 )                             ITFAILS;
    if ( mesh_ptr_H->next_cell(2,2) != 0 )                             ITFAILS;
    if ( mesh_ptr_H->next_cell(2,3) != 1 )                             ITFAILS;
    if ( mesh_ptr_H->next_cell(2,4) != 0 )                             ITFAILS;

    if (fabs(1.0 - 6.0*mesh_ptr_H->volume(1)) > TET_epsilon)           ITFAILS;
    if (fabs(1.0 - 6.0*mesh_ptr_H->volume(2)) > TET_epsilon)           ITFAILS;

    if (fabs(0.707106781187-mesh_ptr_H->face_area(1,1)) > TET_epsilon) ITFAILS;
    if (fabs(0.559016994375-mesh_ptr_H->face_area(1,2)) > TET_epsilon) ITFAILS;
    if (fabs(0.559016994375-mesh_ptr_H->face_area(1,3)) > TET_epsilon) ITFAILS;
    if (fabs(0.5 - mesh_ptr_H->face_area(1,4)) > TET_epsilon)          ITFAILS;

    if (fabs(0.559016994375-mesh_ptr_H->face_area(2,1)) > TET_epsilon) ITFAILS;
    if (fabs(0.559016994375-mesh_ptr_H->face_area(1,2)) > TET_epsilon) ITFAILS;
    if (fabs(0.707106781187-mesh_ptr_H->face_area(2,3)) > TET_epsilon) ITFAILS;
    if (fabs(0.5 - mesh_ptr_H->face_area(1,4)) > TET_epsilon)          ITFAILS;

    // Later, add some vector position and direction tests here.

    // Test interface ===> TET_Builder instantiation with hand-coded interface.
    // This mesh should be the same 2-tet pyramid as in mesh_ptr_H.

    SP<TET_test_1> interface(new TET_test_1());
    if (!interface)                                   ITFAILS;

    TET_Builder builder(interface);

    SP<TET_Mesh> mesh_ptr_0 = builder.build_Mesh();
    if (!mesh_ptr_0)                                  ITFAILS;

    // The pointers themselves should not be equal...
    if (mesh_ptr_H == mesh_ptr_0)                     ITFAILS;

    // ... but, with luck, the meshes will be equal.
    if (*mesh_ptr_H != *mesh_ptr_0)                   ITFAILS;

    vector< vector<double> > m_coords = mesh_ptr_0->get_point_coord();

    if (fabs(m_coords[0][0] - 0.0) > TET_epsilon)     ITFAILS;
    if (fabs(m_coords[0][1] - 0.0) > TET_epsilon)     ITFAILS;
    if (fabs(m_coords[0][2] - 0.0) > TET_epsilon)     ITFAILS;
    if (fabs(m_coords[1][0] - 1.0) > TET_epsilon)     ITFAILS;
    if (fabs(m_coords[1][1] - 0.0) > TET_epsilon)     ITFAILS;
    if (fabs(m_coords[1][2] - 0.0) > TET_epsilon)     ITFAILS;
    if (fabs(m_coords[2][0] - 0.0) > TET_epsilon)     ITFAILS;
    if (fabs(m_coords[2][1] - 1.0) > TET_epsilon)     ITFAILS;
    if (fabs(m_coords[2][2] - 0.0) > TET_epsilon)     ITFAILS;
    if (fabs(m_coords[3][0] - 1.0) > TET_epsilon)     ITFAILS;
    if (fabs(m_coords[3][1] - 1.0) > TET_epsilon)     ITFAILS;
    if (fabs(m_coords[3][2] - 0.0) > TET_epsilon)     ITFAILS;
    if (fabs(m_coords[4][0] - 0.5) > TET_epsilon)     ITFAILS;
    if (fabs(m_coords[4][1] - 0.5) > TET_epsilon)     ITFAILS;
    if (fabs(m_coords[4][2] - 1.0) > TET_epsilon)     ITFAILS;

    // Test interface ===> TET_Builder instantiation with RTT_Format class.
    // This mesh should also be the same 2-tet pyramid as in mesh_ptr_H.

    SP<RTT_Format> reader_1(new RTT_Format("TET_RTT_1"));
    if (!reader_1)                                       ITFAILS;

    TET_Builder read_build_1(reader_1);

    SP<TET_Mesh> mesh_ptr_1 = read_build_1.build_Mesh();
    if (!mesh_ptr_1)                                     ITFAILS;

    // Again, the pointers themselves should not be equal...
    if (mesh_ptr_H == mesh_ptr_1)                        ITFAILS;
    if (mesh_ptr_0 == mesh_ptr_1)                        ITFAILS;

    // ... but the meshes should be equal.
    if (*mesh_ptr_H != *mesh_ptr_1)                      ITFAILS;

    // Test sampling in a tethedron.

    int *idr = init_sprng(0, num, seed, 1);
    Sprng random(idr, 0);

    int numm = 100000;
    double factor = 500.0/(3.0*static_cast<double>(numm));

    double count[10][10][10];
    for (int x = 0; x < 10 ; x++)
       for (int y = 0; y < 10 ; y++)
           for (int z = 0; z < 10 ; z++)
               count[x][y][z] = 0.;

    double x_cent = 0.0;
    double y_cent = 0.0;
    double z_cent = 0.0;

    for (int i = 0 ; i < numm ; i++)
    {
        SF_DOUBLE R = mesh_ptr_1->sample_pos(1,random);

        x_cent += R[0];
        y_cent += R[1];
        z_cent += R[2];

        int xx = static_cast<int>(10.0*R[0]);
        int yy = static_cast<int>(10.0*R[1]);
        int zz = static_cast<int>(10.0*R[2]);

        if (xx < 0 || xx > 9)                           ITFAILS;
        if (yy < 0 || yy > 9)                           ITFAILS;
        if (zz < 0 || zz > 9)                           ITFAILS;

        count[xx][yy][zz] += 1.;
    }
    x_cent /= static_cast<double>(numm);
    y_cent /= static_cast<double>(numm);
    z_cent /= static_cast<double>(numm);

    if (fabs(x_cent - 0.375357) > MID_epsilon)          ITFAILS;
    if (fabs(y_cent - 0.374783) > MID_epsilon)          ITFAILS;
    if (fabs(z_cent - 0.249045) > MID_epsilon)          ITFAILS;

    for (int xx = 0; xx < 10 ; xx++)
        for (int yy = 0; yy < 10 ; yy++)
            for (int zz = 0; zz < 10 ; zz++)
                count[xx][yy][zz] *= factor;

    if (fabs(count[1][6][0] - 1.01500) > MID_epsilon)   ITFAILS;
    if (fabs(count[0][8][1] - 0.23333) > MID_epsilon)   ITFAILS;
    if (fabs(count[5][4][2] - 0.54000) > MID_epsilon)   ITFAILS;
    if (fabs(count[8][1][3] - 0.04167) > MID_epsilon)   ITFAILS;
    if (fabs(count[2][7][4] - 0.31167) > MID_epsilon)   ITFAILS;
    if (fabs(count[5][3][5] - 0.94333) > MID_epsilon)   ITFAILS;
    if (fabs(count[4][4][6] - 1.00333) > MID_epsilon)   ITFAILS;
    if (fabs(count[3][6][7] - 0.02667) > MID_epsilon)   ITFAILS;
    if (fabs(count[4][4][8] - 0.50833) > MID_epsilon)   ITFAILS;
    if (fabs(count[4][4][9] - 0.11167) > MID_epsilon)   ITFAILS;

    if (count[1][9][0] != 0.0)                          ITFAILS;
    if (count[4][6][3] != 0.0)                          ITFAILS;
    if (count[5][2][7] != 0.0)                          ITFAILS;
    if (count[3][4][9] != 0.0)                          ITFAILS;



//  Beginning:  a way to see these results ordered on the page.
    cout.setf(ios::fixed);
    cout.precision(5);
  
    for (int zz = 0 ; zz < 10 ; zz++) {
       for (int yy = 9 ; yy >= 0 ; yy--) {
           for (int xx = 0; xx < 10; xx++)
               cout << setw(10) << count[xx][yy][zz];
           cout << endl;
       }
       cout << endl << endl;
    }
//  End:  a way to see these results ordered on the page.


    // End of test sampling in a tethedron.

    // Test interface ===> TET_Builder instantiation with RTT_Format class.
    // This mesh should be the 96-tet cube from Todd Wareing.

    SP<RTT_Format> reader_2(new RTT_Format("TET_RTT_2"));
    if (!reader_2)                                       ITFAILS;

    TET_Builder read_build_2(reader_2);

    SP<TET_Mesh> mesh_ptr_2 = read_build_2.build_Mesh();
    if (!mesh_ptr_2)                                     ITFAILS;

    // Again, the pointers themselves should not be equal...
    if (mesh_ptr_H == mesh_ptr_2)                        ITFAILS;
    if (mesh_ptr_0 == mesh_ptr_2)                        ITFAILS;
    if (mesh_ptr_1 == mesh_ptr_2)                        ITFAILS;

    // ... and now the meshes should not be equal either.
    if (*mesh_ptr_H == *mesh_ptr_2)                      ITFAILS;
    if (*mesh_ptr_0 == *mesh_ptr_2)                      ITFAILS;
    if (*mesh_ptr_1 == *mesh_ptr_2)                      ITFAILS;

    // Test interface ===> TET_Builder instantiation with RTT_Format class.
    // This mesh should be the one-tet case from the definition file.

    SP<RTT_Format> reader_3(new RTT_Format("TET_RTT_3"));
    if (!reader_3)                                       ITFAILS;

    TET_Builder read_build_3(reader_3);

    SP<TET_Mesh> mesh_ptr_3 = read_build_3.build_Mesh();
    if (!mesh_ptr_3)                                     ITFAILS;

    // Again, the pointers themselves should not be equal...
    if (mesh_ptr_H == mesh_ptr_3)                        ITFAILS;
    if (mesh_ptr_0 == mesh_ptr_3)                        ITFAILS;
    if (mesh_ptr_1 == mesh_ptr_3)                        ITFAILS;
    if (mesh_ptr_2 == mesh_ptr_3)                        ITFAILS;

    // ... and the meshes should not be equal.
    if (*mesh_ptr_H == *mesh_ptr_3)                      ITFAILS;
    if (*mesh_ptr_0 == *mesh_ptr_3)                      ITFAILS;
    if (*mesh_ptr_1 == *mesh_ptr_3)                      ITFAILS;
    if (*mesh_ptr_2 == *mesh_ptr_3)                      ITFAILS;

}   // end Test_TET()

int main(int argc, char *argv[])
{
 try {
    C4::Init(argc, argv);

    // this is a serial test
    if (C4::node())
    {
	C4::Finalize();
	return 0;
    }

    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	{
	    cout << argv[0] << ": version " << rtt_mc::release() << endl;
	    C4::Finalize();
	    return 0;
	}

    // ThreeVector tests
    Test_ThreeVector();

    // TET_Mesh tests
    Test_TET();

    // status of test
    cout << endl;
    cout <<     "************************************" << endl;
    if (passed)
    {
        cout << "**** TET_Mesh Self Test: PASSED ****" << endl;
    }
    cout <<     "************************************" << endl;
    cout << endl;

    cout << "Done testing TET_Mesh." << endl;

    C4::Finalize();

 }
 catch (rtt_dsxx::assertion &yucch) {
    cerr << "\a" << yucch.what() << endl;
    return 1;
 }
 catch (string &s) {
    cerr << "\aERROR:  " << s << endl;
    return 1;
 }
 catch (...) {
    cerr << "\aUnrecognized error." << endl;
    return 1;
 }

}   // end main(int, char *[])

//---------------------------------------------------------------------------//
//                              end of tstTET_Mesh.cc
//---------------------------------------------------------------------------//
