//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   plot2D/test/tstPlot2D.cc
 * \author Rob Lowrie
 * \date   In the past.
 * \brief  Plot2D test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <fstream>
#include <string>

#include "../Plot2D.hh"
#include "../Release.hh"

using std::cin;
using std::cout;
using std::endl;
using std::string;
using std::ofstream;
using rtt_plot2D::Plot2D;

void pause();
void tstPlot2D();
int main(int argc, char *argv[]);

//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
void pause()
{
    cout << "Press RETURN to continue: ";
    cin.get();
}
//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
void
tstPlot2D()
{
    string blockName("tmp.block");

    // Set batch to false if you want to call up the GUI
    bool batch = true;

    // do first plot
    
    ofstream block(blockName.c_str());

    const int n = 10;
    for ( int i = 0; i < n; i++ ) {
	block << i
	      << " " << i
	      << " " << i * i
	      << endl;
    }

    block.close();

    Plot2D p(2, "tstPlot2D.par", batch);

    p.readBlock(blockName);
    p.setTitles("plot0", "subtitle0", 0);
    p.setTitles("plot1", "subtitle1", 1);
    p.setAxesLabels("x", "y0", 0);
    p.setAxesLabels("x", "y1", 1);

    Plot2D::SetProps prop;
    prop.line.color = 2;
    p.setProps(0, 0, prop);

    if ( ! batch ) {
	pause();
    }

    p.save("plot1.agr");

    // second plot

    block.open(blockName.c_str());

    for ( int i = 0; i < n; i++ ) {
	block << i
	      << " " << i * i
	      << " " << i * i * i
	      << endl;
    }

    block.close();

    p.killAllSets();
    p.readBlock(blockName);

    if ( ! batch ) {
	pause();
    }

    p.save("plot2.agr");

    // third plot

    p.arrange(0, 1);

    if ( ! batch ) {
	pause();
    }

    p.save("plot3.agr");

    // fourth plot

    p.close();
    p.open(1, "tstPlot2D.par", batch);
    p.readBlock(blockName, 0);
    p.setTitles("Same Data, One Graph", "subtitle");

    if ( ! batch ) {
	pause();
    }

    p.save("plot4.agr");
}
//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
int
main(int argc, char *argv[])
{
    // version tag
    for (int arg = 1; arg < argc; arg++) {
        if (string(argv[arg]) == "--version") {
            cout << argv[0] << ": version " << rtt_plot2D::release() << endl;
            return 0;
        }
    }

    cout << endl;
    cout << "**********************************************" << endl;
 
    try {
        // tests
        tstPlot2D();
 
        // run python diff scrips
        system("python ./tstPlot2D_Diff.py");
    }
    catch(rtt_dsxx::assertion &ass) {
        cout << "Assertion failure on " << ass.what() << endl;
	cout << "Better luck next time!" << endl;
        return 1;
    }
 
    // status of test
    cout << "********* Plot2D Self Test: PASSED ***********" << endl;
    cout << "**********************************************" << endl;
    cout << endl;
 
    cout << "Done testing Plot2D." << endl;
}

//---------------------------------------------------------------------------//
// end of tstPlot2D.cc
//---------------------------------------------------------------------------//
