//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   Plot2D.cc
 * \author lowrie
 * \date   2002-04-12
 * \brief  Implementation for Plot2D.
 * \note   Copyright (C) 2016-2020 Triad National Security, LLC.
 *         All rights reserved. */
//----------------------------------------------------------------------------//

#include "Plot2D.hh"
#include "plot2D_grace.h"
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>

namespace rtt_plot2D {

//----------------------------------------------------------------------------//
//! Constructor that calls open(); see open() for arguments.
Plot2D::Plot2D(const int numGraphs, const std::string &paramFile,
               const bool batch)
    : d_autoscale(0), d_numGraphs(numGraphs), d_numRows(0), d_numCols(0),
      d_batch(batch), d_graceVersion(), d_numSets(std::vector<int>()),
      d_setsBeenRead(std::vector<bool>()) {
  // is_supported must be checked for all constructors.
  Insist(is_supported(), "Plot2D unsupported on this platform!");

  d_graceVersion = graceVersion();

  open(numGraphs, paramFile, batch);
}

//----------------------------------------------------------------------------//
/*!
 * \brief The destructor.
 *
 * Closes the communication pipe if it is open.
 */
//----------------------------------------------------------------------------//
Plot2D::~Plot2D() { close(); }

//----------------------------------------------------------------------------//
/*!
 * \brief Gets grace's version number.
 *
 * Unfortunately, grace doesn't store the version number in the include file, so
 * 'xmgrace -version' is used.  Ugh.
 *
 * Note: popen and pclose are in the global namespace, at least linux.
 *
 * \return The version number.
 */
//----------------------------------------------------------------------------//
Plot2D::VersionNumber Plot2D::graceVersion() {
  // Grab the version string from a pipe.
  std::FILE *pipe = popen("gracebat -version", "r");

  Insist(pipe, "plot2D::graceVersion: Unable to get grace version!");

  const int bufferSize = 20;
  char buffer[bufferSize];
  bool fgets_ok = true;

  // get the blank line
  fgets_ok = fgets_ok && std::fgets(buffer, bufferSize, pipe) != nullptr;
  // get the version string
  fgets_ok = fgets_ok && std::fgets(buffer, bufferSize, pipe) != nullptr;
  Insist(fgets_ok, "Unable to retrieve Grace's version number.");

  pclose(pipe);

  std::string s(buffer);

  // Assume that the version string is of the form 'Grace-x.x.x'

  Insist(s.size() > 11, "plot2D::graceVersion: Version string too small");

  s = s.substr(6); // skip over 'Grace-'

  VersionNumber n;

  for (int i = 0; i < 3; i++) {
    std::string sn(s);

    if (i < 2) {
      size_t dot = s.find(".");
      Insist(dot != std::string::npos,
             "plot2D::graceVersion: Version string wrong format.");
      sn = s.substr(0, dot);
      s = s.substr(dot + 1);
    }

    std::istringstream is(sn);

    is >> n.v[i];
  }

  return n;
}

//----------------------------------------------------------------------------//
/*!
 * \brief Opens the Grace communication pipe.
 *
 * \param numGraphs Number of graphs to use in the graph matrix.
 * \param paramFile The Grace parameter file (*.par) to load.  If "", no
 *                  parameter file is loaded, and autoscaling is turned on on
 *                  reads.  If a parameter file is specificed, autoscaling is
 *                  turned off on reads.
 * \param batch If true, do not open the Grace GUI window.  Plots may be
 *                  generated and then saved with the save() function.
 */
//----------------------------------------------------------------------------//
void Plot2D::open(const int numGraphs, const std::string &paramFile,
                  const bool batch) {
  Require(numGraphs > 0);
  Require(!GraceIsOpen());

  d_numGraphs = numGraphs;
  d_batch = batch;
  d_numSets.resize(numGraphs);
  d_setsBeenRead.resize(numGraphs);

  for (int i = 0; i < d_numGraphs; i++) {
    d_numSets[i] = 0;
    d_setsBeenRead[i] = false;
  }

  // Open the grace pipe

  const int bufferSize = 4096;
  int openStatus;

  if (d_batch) {
    // Add -nosafe and -noask options for versions 5.1.8 and later.

    // GraceOpenVA takes a char * as its first arg, when it should take a const
    // char *.  So to avoid the justified warnings/errors, jump through some
    // hoops:
    char *exe = new char[9];
    std::strcpy(exe, "gracebat");

    if (d_graceVersion.v[0] >= 5 && d_graceVersion.v[1] >= 1 &&
        d_graceVersion.v[2] >= 8) {
      openStatus =
          GraceOpenVA(exe, bufferSize, "-noprint", "-nosafe", "-noask", NULL);
    } else {
      openStatus = GraceOpenVA(exe, bufferSize, "-noprint", NULL);
    }

    delete[] exe; // end jumping through hoops
  } else {
    // Still want "safe" for non-batch mode?
    openStatus = GraceOpen(bufferSize);
  }

  Insist(openStatus != -1, "Error opening grace.");

  if (paramFile.empty()) {
    d_autoscale = AUTOSCALE_ON;
  } else {
    d_autoscale = AUTOSCALE_OFF;
    // Grab the param file...
    GracePrintf("getp \"%s\"", paramFile.c_str());
  }

  // Set up graph matrix
  arrange(0, 0);
}
//----------------------------------------------------------------------------//
/*!
 *  \brief Arranges the number of graphs into the specified matrix.
 *
 * \param numRows Number of rows to use in the graph matrix.  If less than 1,
 *                automatically computed.
 *  \param numCols Number of columns to use in the graph matrix.  If less than
 *                 1, automatically computed.
 *
 *  Specifying both numRows and numCols less than 1 means both are computed
 *  automatically.
 */
//----------------------------------------------------------------------------//
void Plot2D::arrange(const int numRows, const int numCols) {
  Require(GraceIsOpen());

  if (numRows > 0) {
    // number of rows specified
    d_numRows = numRows;
    if (numCols > 0) {
      // .. as are number of columns
      d_numCols = numCols;
      Require(d_numCols * d_numRows >= d_numGraphs);
    } else {
      // .. compute number of columns
      d_numCols = d_numGraphs / d_numRows;
      if (d_numCols * d_numRows != d_numGraphs) {
        ++d_numCols;
      }
    }
  } else if (numCols > 0) {
    // only number of columns specified
    d_numCols = numCols;
    d_numRows = d_numGraphs / d_numCols;
    if (d_numCols * d_numRows != d_numGraphs) {
      ++d_numRows;
    }
  } else {
    // neither number of columns or rows was specified
    d_numRows = int(std::sqrt(double(d_numGraphs)));
    d_numCols = d_numGraphs / d_numRows;

    if (d_numCols * d_numRows != d_numGraphs) {
      ++d_numCols;
    }
  }

  GracePrintf("arrange(%d, %d, 0.1, 0.3, 0.2)", d_numRows, d_numCols);

  // turn off unused graphs

  for (int i = d_numGraphs; i < d_numCols * d_numRows; i++) {
    GracePrintf("focus g%d", graphNum(i, true));
    GracePrintf("frame off");
    GracePrintf("xaxis off");
    GracePrintf("yaxis off");
  }

  redraw();
}

//----------------------------------------------------------------------------//
/*!
 *  \brief Closes the Grace communication pipe.
 *
 *  All sets and properties are erased.  This function is called by the
 *  destructor.
*/
//----------------------------------------------------------------------------//
void Plot2D::close() {
  if (GraceIsOpen()) {
    if (GraceClose() == -1) {
      std::cerr << "WARNING: Error closing xmgrace." << std::endl;
    }
  }
}

//----------------------------------------------------------------------------//
/*!
 *  \brief Saves the current plot in a Grace project file.
 *
 *  \param filename The file name to use.
*/
//----------------------------------------------------------------------------//
void Plot2D::save(const std::string filename) {
  Require(GraceIsOpen());

  if (!filename.empty()) {
    GracePrintf("saveall \"%s\"", filename.c_str());
  }
}

//----------------------------------------------------------------------------//
/*!
 *  \brief Sets the graph title and subtitle.
 *
 *  \param title The title.
 *  \param subTitle The subtitle.
 *  \param iG The graph number to apply the titles to.
 */
//----------------------------------------------------------------------------//
void Plot2D::setTitles(const std::string title, const std::string subTitle,
                       const int iG) {
  Require(GraceIsOpen());

  GracePrintf("focus g%d", graphNum(iG));
  GracePrintf("title \"%s\"", title.c_str());
  GracePrintf("subtitle \"%s\"", subTitle.c_str());
  redraw();
}

//----------------------------------------------------------------------------//
/*!
 *  \brief Sets the axes labels.
 *
 *  \param xLabel The label for the x-axis.
 *  \param yLabel The label for the y-axis.
 *  \param iG The graph number.
 *  \param charSize The character size in [0,1].
*/
//----------------------------------------------------------------------------//
void Plot2D::setAxesLabels(const std::string xLabel, const std::string yLabel,
                           const int iG, const double charSize) {
  Require(GraceIsOpen());

  GracePrintf("focus g%d", graphNum(iG));

  GracePrintf("xaxis label \"%s\"", xLabel.c_str());
  if (xLabel.size() > 0) {
    GracePrintf("xaxis label char size %f", charSize);
  }

  GracePrintf("yaxis label \"%s\"", yLabel.c_str());
  if (yLabel.size() > 0) {
    GracePrintf("yaxis label layout perp");
    GracePrintf("yaxis label char size %f", charSize);
  }

  GracePrintf("xaxis ticklabel char size %f", charSize);
  GracePrintf("yaxis ticklabel char size %f", charSize);

  redraw();
}

//----------------------------------------------------------------------------//
/*!
 *  \brief Kills all sets from graphs.
 *
 *  Note that the set parameters are saved, by telling Grace to "saveall".  This
 *  way we don't have to reload the parameter file, or require the user to call
 *  setProps(), after new sets are read.
 */
//----------------------------------------------------------------------------//
void Plot2D::killAllSets() {
  Require(GraceIsOpen());

  for (int iG = 0; iG < d_numGraphs; iG++) {

    for (int j = 0; j < d_numSets[iG]; j++) {
      GracePrintf("kill g%d.s%d saveall", graphNum(iG), j);
    }

    d_numSets[iG] = 0;
  }
}

//----------------------------------------------------------------------------//
/*!
 *  \brief Reads block data from file, one set per graph.
 *
 *  \param blockFilename The name of the block data file.
 *
 *  One set of data is added to each graph.  If one wants to replace the sets
 *  that are currently plotted, call killAllSets() before calling this function.
 *
 *  The datafile must be columns in the format
 *
 *  x y(1) y(2) .... y(numGraphs)
 *
 *  where (x, y(N)) is the set added to graph number N.
 */
//----------------------------------------------------------------------------//
void Plot2D::readBlock(const std::string blockFilename) {
  Require(GraceIsOpen());
  Require(numColumnsInFile(blockFilename) == d_numGraphs + 1);

  GracePrintf("read block \"%s\"", blockFilename.c_str());

  for (int iG = 0; iG < d_numGraphs; iG++) {
    GracePrintf("focus g%d", graphNum(iG));
    setAutoscale(iG);
    GracePrintf("block xy \"1:%d\"", iG + 2);

    ++d_numSets[iG];
    d_setsBeenRead[iG] = true;
  }

  redraw();
}

//----------------------------------------------------------------------------//
/*!
 *  \brief Reads block data from file, all sets into one graph.
 *
 *  \param blockFilename The name of the block data file.
 *  \param iG The graph number to add the sets to.
 *
 *  The sets are added to the graph.  If one wants to replace the sets that are
 *  currently plotted, call killAllSets() before calling this function.
 *
 *  The datafile must be columns in the format
 *
 *  x y(1) y(2) .... y(nSets)
 *
 *  where each pair (x, y(N)) is a set added to graph number \a iG.
 */
//----------------------------------------------------------------------------//
void Plot2D::readBlock(const std::string blockFilename, const int iG) {
  Require(GraceIsOpen());
  Require(iG >= 0 && iG < d_numGraphs);

  const int nSets = numColumnsInFile(blockFilename) - 1;

  GracePrintf("autoscale onread none");
  GracePrintf("read block \"%s\"", blockFilename.c_str());
  GracePrintf("focus g%d", graphNum(iG));
  setAutoscale(iG);

  for (int i = 0; i < nSets; i++) {
    GracePrintf("block xy \"1:%d\"", i + 2);
    ++d_numSets[iG];
  }

  d_setsBeenRead[iG] = true;

  redraw();
}

//----------------------------------------------------------------------------//
/*!
 *  \brief Redraws the graph.
 *
 *  \param autoscale If true, autoscale all of the graphs.
*/
//----------------------------------------------------------------------------//
void Plot2D::redraw(const bool autoscale) {
  Require(GraceIsOpen());

  if (autoscale) {
    for (int iG = 0; iG < d_numGraphs; iG++) {
      GracePrintf("focus g%d", graphNum(iG));
      GracePrintf("autoscale");
    }
  }

  GracePrintf("redraw");
  GraceFlush();
}

//----------------------------------------------------------------------------//
/*!
 *  \brief Changes the set properties.
 *
 *  \param iG The graph number for the set.
 *  \param iSet The set number.  The set is NOT required to exist, so that
 *              setProps() may be called before any data is read.
 *  \param p The set properties.
*/
//----------------------------------------------------------------------------//
void Plot2D::setProps(const int iG, const int iSet, const SetProps &p) {
  Require(GraceIsOpen());
  Require(iSet >= 0);

  GracePrintf("focus g%d", graphNum(iG));

  GracePrintf("s%d line linestyle %d", iSet, p.line.style);
  if (p.line.style != LineProps::STYLE_NONE) {
    GracePrintf("s%d line color %d", iSet, p.line.color);
    GracePrintf("s%d line linewidth %f", iSet, p.line.width);
  }

  GracePrintf("s%d symbol %d", iSet, p.symbol.shape);
  if (p.symbol.shape != SymbolProps::SHAPE_NONE) {
    GracePrintf("s%d symbol size %f", iSet, p.symbol.size);
    GracePrintf("s%d symbol color %d", iSet, p.symbol.color);
    GracePrintf("s%d symbol linewidth %f", iSet, p.symbol.width);
    GracePrintf("s%d symbol fill pattern %d", iSet, p.symbol.fillPattern);
    GracePrintf("s%d symbol fill color %d", iSet, p.symbol.fillColor);
  }

  if (p.legend.size() > 0) {
    GracePrintf("s%d legend \"%s\"", iSet, p.legend.c_str());
  }

  redraw();
}

//----------------------------------------------------------------------------//
/*!
 *  \brief Determines the Grace graph number for the given graph number.
 *
 *  \param iG The graph number used by Plot2D.
 *  \param allowVacant Allow access to vacant graph locations.
 *
 *  \return The graph number used by Grace.
 *
 *  For a 3x3 layout, Grace lays out the graph numbers as
 *
 *  0 1 2
 *
 *  3 4 5
 *
 *  6 7 8
 *
 *  To follow this, we could just return iG.  This code allows more
 *  general layouts and does some checking.
 */
//----------------------------------------------------------------------------//
int Plot2D::graphNum(const int iG, const bool allowVacant) const {
  Require(iG >= 0);
  Require((allowVacant && iG < d_numRows * d_numCols) || iG < d_numGraphs);

  int iRow = iG / d_numCols;
  int iCol = iG - d_numCols * iRow;

  Ensure(iCol < d_numCols);
  Ensure(iRow < d_numRows);

  return iRow * d_numCols + iCol;
}

//----------------------------------------------------------------------------//
/*!
 *  \brief Computes the number of columns of data in a file.
 *
 *  \param filename The name of the file to parse.
 *
 *  \return The number of columns.
*/
//----------------------------------------------------------------------------//
int Plot2D::numColumnsInFile(const std::string filename) const {
  using std::isspace;

  std::ifstream f(filename.c_str());

  std::string buf;
  std::getline(f, buf);

  f.close();

  int n = 0; // return value
  bool whitespace = true;

  for (size_t i = 0; i < buf.size(); i++) {
    if (isspace(buf[i])) {
      whitespace = true;
    } else {
      if (whitespace) {
        ++n;
      }
      whitespace = false;
    }
  }

  return n;
}

//----------------------------------------------------------------------------//
/*!
 *  \brief Sets the current autoscaling mode for a graph.
 *
 *  \param iG Graph number.
*/
//----------------------------------------------------------------------------//
void Plot2D::setAutoscale(const int iG) {
  switch (d_autoscale) {
  case AUTOSCALE_ON:
    GracePrintf("autoscale onread xyaxes");
    break;
  case AUTOSCALE_OFF:
    GracePrintf("autoscale onread none");
    break;
  case AUTOSCALE_FIRSTREAD:
    if (d_setsBeenRead[iG]) {
      GracePrintf("autoscale onread none");
    } else {
      GracePrintf("autoscale onread xyaxes");
    }
    break;
  }
}

} // end namespace rtt_plot2D

//----------------------------------------------------------------------------//
// end of Plot2D.cc
//----------------------------------------------------------------------------//
