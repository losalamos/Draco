//----------------------------------*-C++-*----------------------------------//
/*!
  \file   Plot2D.cc
  \author lowrie
  \date   2002-04-12
  \brief  Implementation for Plot2D.
*/
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Plot2D.hh"

#include <iostream>
#include <fstream>
#include <cmath>
#include <cctype>

#include <grace_np.h>

using namespace rtt_plot2D;

//---------------------------------------------------------------------------//
/*!
  \brief Default constructor.

  Must then call open() to actually generate plots.
*/
//---------------------------------------------------------------------------//
Plot2D::
Plot2D()
{
}
//---------------------------------------------------------------------------//
/*!
  \brief Constructor that calls open(); see open() for arguments.
*/
//---------------------------------------------------------------------------//
Plot2D::
Plot2D(const int numGraphs,
       const std::string &paramFile,
       const bool batch)
{
    open(numGraphs, paramFile, batch);
}
//---------------------------------------------------------------------------//
/*!
  \brief The destructor.

  Closes the communication pipe if it is open.
*/
//---------------------------------------------------------------------------//
Plot2D::
~Plot2D()
{
    close();
}
//---------------------------------------------------------------------------//
/*!
  \brief Opens the Grace communication pipe.

  \param numGraphs Number of graphs to use in the graph matrix.

  \param paramFile The Grace parameter file (*.par) to load.  If "",
  no parameter file is loaded.

  \param batch If true, do not open the Grace GUI window.  Plots may
  be generated and then saved with the save() function.
*/
//---------------------------------------------------------------------------//
void
Plot2D::
open(const int numGraphs,
     const std::string &paramFile,
     const bool batch)
{
    Require(numGraphs > 0);
    Require(! GraceIsOpen());

    d_numGraphs = numGraphs;
    d_batch = batch;
    d_numSets.resize(numGraphs);
    
    for ( int i = 0; i < d_numGraphs; i++ ) {
	d_numSets[i] = 0;
    }

    // Open the grace pipe
    
    const int bufferSize = 4096;
    int openStatus;

    if ( d_batch ) {
	openStatus = GraceOpenVA("gracebat", bufferSize, "-noprint", NULL);
    }
    else {
	openStatus = GraceOpen(bufferSize);
    }

    Insist(openStatus != -1, "Error opening grace.");

    if ( ! paramFile.empty() ) {
	// Grab the param file...
	GracePrintf("getp \"%s\"", paramFile.c_str());
    }

    // Set up graph matrix

    arrange(0, 0);
}
//---------------------------------------------------------------------------//
/*!
  \brief Arranges the number of graphs into the specified matrix.

  \param numRows Number of rows to use in the graph matrix.  If less than
  1, automatically computed.

  \param numCols Number of columns to use in the graph matrix.  If less than
  1, automatically computed.

  Specifying both numRows and numCols less than 1 means both are computed
  automatically.
*/
//---------------------------------------------------------------------------//
void
Plot2D::
arrange(const int numRows,
	const int numCols)
{
    Require(GraceIsOpen());

    if ( numRows > 0 ) {
	// number of rows specified
	d_numRows = numRows;
	if ( numCols > 0 ) {
	    // .. as are number of columns
	    d_numCols = numCols;
	    Require(d_numCols * d_numRows >= d_numGraphs);
	}
	else {
	    // .. compute number of columns
	    d_numCols = d_numGraphs / d_numRows;
	    if ( d_numCols * d_numRows != d_numGraphs ) {
		++d_numCols;
	    }
	}
    }
    else if ( numCols > 0 ) {
	// only number of columns specified
	d_numCols = numCols;
	d_numRows = d_numGraphs / d_numCols;
	if ( d_numCols * d_numRows != d_numGraphs ) {
	    ++d_numRows;
	}
    }
    else {
	// neither number of columns or rows was specified
	d_numRows = int(std::sqrt(double(d_numGraphs)));
	d_numCols = d_numGraphs / d_numRows;
    
	if ( d_numCols * d_numRows != d_numGraphs ) {
	    ++d_numCols;
	}
    }

    GracePrintf("arrange(%d, %d, 0.1, 0.3, 0.2)",
		d_numRows, d_numCols);
    
    // turn off unused graphs
    
    for ( int i = d_numGraphs; i < d_numCols * d_numRows; i++ ) {
	GracePrintf("focus g%d", graphNum(i, true));
	GracePrintf("frame off");
	GracePrintf("xaxis off");
	GracePrintf("yaxis off");
    }
    
    redraw();
}
//---------------------------------------------------------------------------//
/*!
  \brief Closes the Grace communication pipe.

  All sets and properties are erased.  This function is called
  by the destructor.
*/
//---------------------------------------------------------------------------//
void
Plot2D::
close()
{
    if ( GraceIsOpen() ) {
	if ( GraceClose() == -1 ) {
	    std::cerr << "WARNING: Error closing xmgrace." << std::endl;
	}
    }
}
//---------------------------------------------------------------------------//
/*!
  \brief Saves the current plot in a Grace project file.

  \param filename The file name to use.
*/
//---------------------------------------------------------------------------//
void
Plot2D::
save(const std::string filename)
{
    Require(GraceIsOpen());
    
    if ( ! filename.empty() ) {
	GracePrintf("saveall \"%s\"", filename.c_str());
    }
}
//---------------------------------------------------------------------------//
/*!
  \brief Sets the graph title and subtitle.

  \param title The title.

  \param subTitle The subtitle.

  \param iG The graph number to apply the titles to.
*/
//---------------------------------------------------------------------------//
void
Plot2D::
setTitles(const std::string title,
	  const std::string subTitle,
	  const int iG)
{
    Require(GraceIsOpen());
    
    GracePrintf("focus g%d", graphNum(iG));
    GracePrintf("title \"%s\"", title.c_str());
    GracePrintf("subtitle \"%s\"", subTitle.c_str());
    redraw();
}
//---------------------------------------------------------------------------//
/*!
  \brief Sets the axes labels.

  \param xLabel The label for the x-axis.

  \param yLabel The label for the y-axis.

  \param iG The graph number.

  \param charSize The character size in [0,1].
*/
//---------------------------------------------------------------------------//
void
Plot2D::
setAxesLabels(const std::string xLabel,
	      const std::string yLabel,
	      const int iG,
	      const double charSize)
{
    Require(GraceIsOpen());

    GracePrintf("focus g%d", graphNum(iG));

    GracePrintf("xaxis label \"%s\"", xLabel.c_str());
    if ( xLabel.size() > 0 ) {
	GracePrintf("xaxis label char size %f", charSize);
    }

    GracePrintf("yaxis label \"%s\"", yLabel.c_str());
    if ( yLabel.size() > 0 ) {
	GracePrintf("yaxis label layout perp");
	GracePrintf("yaxis label char size %f", charSize);
    }

    GracePrintf("xaxis ticklabel char size %f", charSize);
    GracePrintf("yaxis ticklabel char size %f", charSize);

    redraw();
}
//---------------------------------------------------------------------------//
/*!
  \brief Kills all sets from graphs.

  Note that the set parameters are saved, by telling Grace to "saveall".
  This way we don't have to reload the parameter file, or require the user to
  call setProps(), after new sets are read.
*/
//---------------------------------------------------------------------------//
void
Plot2D::
killAllSets()
{
    Require(GraceIsOpen());
    
    for ( int iG = 0; iG < d_numGraphs; iG++ ) {

	for ( int j = 0; j < d_numSets[iG]; j++ ) {
	    GracePrintf("kill g%d.s%d saveall", graphNum(iG), j);
	}

	d_numSets[iG] = 0;
    }
}
//---------------------------------------------------------------------------//
/*!
  \brief Sends a command to the Grace communications pipe.

  \param command The command string to be sent.

  This allows full access to the Grace capabilities, albeit at a
  crude level.  See the Grace documentation for GracePrintf.
*/
//---------------------------------------------------------------------------//
void
Plot2D::
rawCom(const std::string command)
{
    Require(GraceIsOpen());
    GracePrintf(command.c_str());
}
//---------------------------------------------------------------------------//
/*!
  \brief Reads block data from file, one set per graph.

  \param blockFilename The name of the block data file.

  One set of data is added to each graph.  If one wants to replace the sets
  that are currently plotted, call killAllSets() before calling this
  function.

  The datafile must be columns in the format

  x y(1) y(2) .... y(numGraphs)

  where (x, y(N)) is the set added to graph number N.
*/
//---------------------------------------------------------------------------//
void
Plot2D::
readBlock(const std::string blockFilename)
{
    Require(GraceIsOpen());
    Require(numColumnsInFile(blockFilename) == d_numGraphs + 1);

    GracePrintf("read block \"%s\"", blockFilename.c_str());
    
    for ( int iG = 0; iG < d_numGraphs; iG++ ) {
	GracePrintf("focus g%d", graphNum(iG));
	GracePrintf("block xy \"1:%d\"", iG + 2);

	++d_numSets[iG];
    }

    redraw();
}
//---------------------------------------------------------------------------//
/*!
  \brief Reads block data from file, all sets into one graph.

  \param blockFilename The name of the block data file.

  \param iG The graph number to add the sets to.

  The sets are added to the graph.  If one wants to replace the sets that are
  currently plotted, call killAllSets() before calling this function.
  
  The datafile must be columns in the format

  x y(1) y(2) .... y(nSets)

  where each pair (x, y(N)) is a set added to graph number \a iG.
*/
//---------------------------------------------------------------------------//
void
Plot2D::
readBlock(const std::string blockFilename,
	  const int iG)
{
    Require(GraceIsOpen());
    Require(iG >= 0 && iG < d_numGraphs);

    const int nSets = numColumnsInFile(blockFilename) - 1;

    GracePrintf("read block \"%s\"", blockFilename.c_str());
    GracePrintf("focus g%d", graphNum(iG));
    
    for ( int i = 0; i < nSets; i++ ) {
	GracePrintf("block xy \"1:%d\"", i + 2);
	++d_numSets[iG];
    }

    redraw();
}
//---------------------------------------------------------------------------//
/*!
  \brief Redraws the graph.

  \param autoscale If true, autoscale all of the graphs.
*/
//---------------------------------------------------------------------------//
void
Plot2D::
redraw(const bool autoscale)
{
    Require(GraceIsOpen());

    if ( autoscale ) {
	for ( int iG = 0; iG < d_numGraphs; iG++ ) {
	    GracePrintf("focus g%d", graphNum(iG));
	    GracePrintf("autoscale");
	}
    }

    GracePrintf("redraw");
    GraceFlush();
}
//---------------------------------------------------------------------------//
/*!
  \brief Changes the set properties.

  \param iG The graph number for the set.

  \param iSet The set number.  The set is NOT required to exist,
  so that setProps() may be called before any data is read.

  \param p The set properties.
*/
//---------------------------------------------------------------------------//
void
Plot2D::
setProps(const int iG,
	 const int iSet,
	 const SetProps &p)
{
    Require(GraceIsOpen());
    Require(iSet >= 0);

    GracePrintf("focus g%d", graphNum(iG));

    GracePrintf("s%d line type %d", iSet, p.line.type);
    if ( p.line.type > 0 ) {
	GracePrintf("s%d line color %d", iSet, p.line.color);
	GracePrintf("s%d line linewidth %f", iSet, p.line.width);
    }
    
    GracePrintf("s%d symbol %d", iSet, p.symbol.type);
    if ( p.symbol.type > 0 ) {
	GracePrintf("s%d symbol size %f", iSet, p.symbol.size);
	GracePrintf("s%d symbol color %d", iSet, p.symbol.color);
	GracePrintf("s%d symbol linewidth %f", iSet, p.symbol.width);
	GracePrintf("s%d symbol fill pattern %d",
		    iSet, p.symbol.fillPattern);
	GracePrintf("s%d symbol fill color %d", iSet, p.symbol.fillColor);
    }
    
    if ( p.legend.size() > 0 ) {
	GracePrintf("s%d legend \"%s\"", iSet, p.legend.c_str());
    }

    redraw();
}
//---------------------------------------------------------------------------//
/*!
  \brief Changes the set properties for all graphs.

  \param iSet The set number.  The set is NOT required to exist
  in all graphs.
*/
//---------------------------------------------------------------------------//
void
Plot2D::
setProps(const int iSet,
	 const SetProps &p)
{
    for ( int iG = 0; iG < d_numGraphs; iG++ ) {
	setProps(iG, iSet, p);
    }
}
//---------------------------------------------------------------------------//
/*!
  \brief Determines the Grace graph number for the given graph number.
  
  \param iG The graph number used by Plot2D.

  \param allowVacant Allow access to vacant graph locations.
  
  \returns The graph number used by Grace.

  For a 3x3 layout, Grace lays out the graph numbers as

  0 1 2

  3 4 5

  6 7 8

  To follow this, we could just return iG.  This code allows more
  general layouts and does some checking.
*/
//---------------------------------------------------------------------------//
int
Plot2D::
graphNum(const int iG,
	 const bool allowVacant) const
{
    Require(iG >= 0);
    Require((allowVacant && iG < d_numRows * d_numCols)
	    || iG < d_numGraphs);
    
    int iRow = iG / d_numCols;
    int iCol = iG - d_numCols * iRow;

    Ensure(iCol < d_numCols);
    Ensure(iRow < d_numRows);

    return iRow * d_numCols + iCol;
}
//---------------------------------------------------------------------------//
/*!
  \brief Computes the number of columns of data in a file.
  
  \param filename The name of the file to parse.
  
  \returns The number of columns.
*/
//---------------------------------------------------------------------------//
int
Plot2D::
numColumnsInFile(const std::string filename) const
{
    std::ifstream f(filename.c_str());

    std::string buf;
    std::getline(f, buf);

    f.close();

    int n = 0; // return value
    bool whitespace = true;

    for ( int i = 0; i < buf.size(); i++ ) {
	if ( std::isspace(buf[i]) ) {
	    whitespace = true;
	}
	else {
	    if ( whitespace ) {
		++n;
	    }
	    whitespace = false;
	}
    }

    return n;
}

//---------------------------------------------------------------------------//
// end of Plot2D.cc
//---------------------------------------------------------------------------//
