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
#include <sstream>
#include <string>
#include <cmath>
#include <cstdlib>

#include "../Plot2D.hh"
#include "ds++/Release.hh"

using std::cin;
using std::cout;
using std::endl;
using std::string;
using std::ofstream;
using rtt_plot2D::Plot2D;
using rtt_plot2D::SetProps;

namespace test
{
void pause();
}

void tstPlot2D(const bool batch);
int main(int argc, char *argv[]);

//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
namespace test
{

void pause()
{
    cout << "Press RETURN to continue: ";
    cin.get();
}

}
//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
void
tstPlot2D(const bool batch)
{
    string paramFile("tstPlot2D.par");

    const int numGraphs = 3; // number of graphs in plot
    const int n = 10; // number of points in each set

    //////////////////////////////////////////////////////////////////////
    // Plot 1.
    //////////////////////////////////////////////////////////////////////

    // Output the block data file that Grace will read

    string blockName("tmp.block");
    ofstream block(blockName.c_str());
    
    for ( int i = 0; i < n; i++ ) {
	block << i;
	for ( int j = 0; j < numGraphs; j++ ) {
	    block << " " << std::pow(double(i), j + 1);
	}
	block << endl;
    }

    block.close();

    // Generate the plot

    Plot2D p(numGraphs, paramFile, batch);

    p.readBlock(blockName);

    // Set titles and axis labels

    for ( int j = 0; j < numGraphs; j++ ) {
	std::ostringstream title, subtitle, ylabel;
	title << "title " << j;
	subtitle << "subtitle " << j;
	ylabel << "y" << j;
	p.setTitles(title.str(), subtitle.str(), j);
	p.setAxesLabels("x", ylabel.str(), j);
    }

    // Changes some set properties; most were set in the param
    // file.

    SetProps prop;
    prop.line.color = rtt_plot2D::COLOR_RED;
    p.setProps(0, 0, prop); // of graph 0, set 0

    if ( ! batch ) {
	test::pause();
    }

    p.save("plot1.agr");

    //////////////////////////////////////////////////////////////////////
    // Plot 2.  Add a set to each graph of the previous plot.
    //////////////////////////////////////////////////////////////////////

    // Generate new data.  Use a different temporary filename,
    // because Grace may not have read the data from plot 1 yet.

    blockName = "tmp2.block";
    block.open(blockName.c_str());

    for ( int i = 0; i < n; i++ ) {
	block << i;
	for ( int j = 0; j < numGraphs; j++ ) {
	    block << " " << i + std::pow(double(i), j + 2);
	}
	block << endl;
    }

    block.close();

    // Add the data to the graphs

    p.killAllSets();
    p.readBlock(blockName);

    if ( ! batch ) {
	test::pause();
    }

    p.save("plot2.agr");

    //////////////////////////////////////////////////////////////////////
    // Plot 3.  Rearrange graph matrix.
    //////////////////////////////////////////////////////////////////////

    p.arrange(2, 2);

    if ( ! batch ) {
	test::pause();
    }

    p.save("plot3.agr");

    //////////////////////////////////////////////////////////////////////
    // Plot 4.  Put all of the data of the previous plot
    // into one graph.
    //////////////////////////////////////////////////////////////////////

    p.close();
    p.open(1, "", batch); // 1 graph, no param file specified

    // Must specify the graph number (0) to put all sets into one
    // graph.
    p.readBlock(blockName, 0);
    
    p.setTitles("Same Data, One Graph", "subtitle");

    if ( ! batch ) {
	test::pause();
    }

    p.save("plot4.agr");
}
//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
int
main(int argc, char *argv[])
{
    bool batch = true;
    
    // version tag
    for (int arg = 1; arg < argc; arg++) {
        if (string(argv[arg]) == "--version") {
            cout << argv[0] << ": version " << rtt_dsxx::release() << endl;
            return 0;
        }
        else if (string(argv[arg]) == "--gui") {
	    batch = false;
	}
    }

    cout << endl;
    cout << "**********************************************" << endl;
 
    try {
        // tests
        if ( Plot2D::is_supported() ) {
            tstPlot2D(batch);
 
            // run python diff scrips
            if (std::system("python ./tstPlot2D_Diff.py"))
                throw rtt_dsxx::assertion ("Python script failed.");
        }
        else {
            cout << "Unsupported test: pass\n";
        }
    }
    catch(rtt_dsxx::assertion &ass) {
        cout << "Assertion: " << ass.what() << endl;
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
