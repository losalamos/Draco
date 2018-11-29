//----------------------------------*-C++-*----------------------------------//
/*!
  \file   Plot2D.hh
  \author lowrie
  \date   2002-04-12
  \brief  Header for Plot2D.
  \note   Copyright (C) 2016-2018 Los Alamos National Security, LLC.
          All rights reserved.*/
//---------------------------------------------------------------------------//

#ifndef rtt_plot2D_Plot2D_hh
#define rtt_plot2D_Plot2D_hh

#include "SetProps.hh"
#include <vector>

namespace rtt_plot2D {

//===========================================================================//
/*!

  \class Plot2D

  \brief Generates 2D line plots via the Grace plotting program.

  The purpose of this interface is to allow one to generate diagnostic
  graphs from inside a running code. The interface is NOT designed to provide
  access to every available capability of Grace.  This documentation assumes
  a basic understanding of Grace.

  Plot2D uses similar notation as Grace.  A plot is a collection of graphs.
  A graph is a coordinate axes, with a (optional) title, and any number of
  sets (y vs. x data) graphed with respect to the axes.  A given set is
  always associated with a graph.

  Typically, Plot2D opens a Grace window and arranges the specified
  number of graphs in a matrix.  For example, by default 8 graphs are
  arranged in 3 rows and columns, numbered as

  0 1 2

  3 4 5

  6 7

  The plot window stays open until close() is called, or the destructor
  is called.  As long as the window is open, one has full access to the
  Grace menus.

  Plot2D provides a limited capability for setting titles, axes labels, and set
  properties.  For more detailed settings, one has two alternatives:

  1) The rawCom() function.  See the Grace documentation for GracePrintf.

  2) Make changes while the plot window is up using the Grace menus.  Once
  these changes are made, a parameter file (*.par) file may be saved under
  Grace (Plot -> Save Parameters).  The Plot2D constructor, or open(), takes
  the parameter file as an optional argument for future instances.  Changes
  to the graph matrix layout are NOT saved in the parameter file, nor are any
  sets that were entered manually, such as via the Grace spreadsheet.

  Sets are loaded into Grace through Grace's block data files.  See the
  readBlock() member function.  Although crude, according to Grace's author,
  this is the most efficient way of passing large amounts of data to Grace.
  The major drawback is that because Grace is running asynchronously, two (or
  more) successive calls to readBlock() may erase the data on the first call
  before Grace has a chance to read it.  Aside from putting in a delay
  between calls, there is no way around this without modifying Grace.  On the
  other hand, Grace WILL plot the data from the last call to readBlock()
  correctly, but older data sets may appear corrupted.

  Finally, future hackers should not consider using the "point" manner
  of communicating data to Grace, such as

  GracePrintf("g0.s0 point %f, %f", x, y)

  This is horribly inefficient and will slow Grace to a standstill unless
  the amount of data is small.  Ideally, Grace should be modified to
  allow something like "block data start" and "block data end" commands.
*/
//===========================================================================//
class Plot2D {
  // TYPES

  // Mode for autoscaling on reads
  enum AutoscaleMode { AUTOSCALE_ON, AUTOSCALE_OFF, AUTOSCALE_FIRSTREAD };

  // Version number (used for Grace's version).
  // Assumes format v[0].v[1].v[2]
  struct VersionNumber {
    static const int indices = 3;
    int v[indices];

    VersionNumber() {
      for (int i = 0; i < indices; i++) {
        v[i] = 0;
      }
    }
  };

  // DATA

  // current autoscale mode
  int d_autoscale;

  // number of graphs to be plotted
  int d_numGraphs;

  // number of rows in plot matrix
  int d_numRows;

  // number of columns in plot matrix
  int d_numCols;

  // true if in batch mode
  bool d_batch;

  // grace's version number
  VersionNumber d_graceVersion;

  // current number of sets in each graph.  Needed so we can
  // kill sets.
  std::vector<int> d_numSets;

  // If true for a particular graph, then sets have been read into
  // the particular graph.
  std::vector<bool> d_setsBeenRead;

public:
  // CREATORS

  // default constructor; window left unopened
  // Plot2D();

  // another constructor; window is opened
  explicit Plot2D(const int numGraphs, const std::string &paramFile = "",
                  const bool batch = false);

  // destructor
  ~Plot2D();

  // MANIPULATORS

  //! Returns true if platform is supported.
  static bool is_supported();

  // arranges the graph matrix
  void arrange(const int numRows = 0, const int numCols = 0);

  // closes Grace window
  void close();

  // deletes all data sets
  void killAllSets();

  // opens Grace window
  void open(const int numGraphs = 1, const std::string &paramFile = "",
            const bool batch = false);

  // sends Grace a command
  // void rawCom(const std::string command);

  // reads block data from file, one set per graph
  void readBlock(const std::string blockFilename);

  // reads block data from file, all sets into one graph
  void readBlock(const std::string blockFilename, const int iG);

  // redraws the graphs
  void redraw(const bool autoscale = false);

  // saves plot into Grace project file
  void save(const std::string filename);

  // sets the titles
  void setTitles(const std::string title, const std::string subTitle,
                 const int iG = 0);

  // sets the axes labels
  void setAxesLabels(const std::string xLabel, const std::string yLabel,
                     const int iG, const double charSize = 1.0);

  // Turns on autoscale when reading sets.
  // void autoscaleOnRead();

  // Turns off autoscale when reading sets.
  //void noAutoscaleOnRead();

  // Turns on autoscaling for first set read into a graph.
  // void autoscaleOnFirstRead();

  // sets the properties for a data set
  void setProps(const int iG, const int iSet, const SetProps &setProps);

  // sets the properties for a data set in all graphs
  // void setProps(const int iSet,
  //     	  const SetProps &setProps);

  // ACCESSORS

  /// Returns true if running in batch mode.
  bool batch() const { return d_batch; }

private:
  // IMPLEMENTATION

  // Returns grace's version number
  VersionNumber graceVersion();

  // Returns the Grace graph number for the given graph number.
  int graphNum(const int iEqn, const bool allowVacant = false) const;

  // Returns the number of columns of data in file
  int numColumnsInFile(const std::string filename) const;

  // sets the current autoscaling mode for a graph
  void setAutoscale(const int iG);

  // We don't implement these because it's not clear how
  // the grace pipe would be copied.  Moreover, these
  // ops are not well defined conceptually.

  /// NOT IMPLEMENTED
  Plot2D(const Plot2D &from);

  /// NOT IMPLEMENTED
  Plot2D &operator=(const Plot2D &from);
};

} // namespace rtt_plot2D

#endif // rtt_plot2D_Plot2D_hh

//---------------------------------------------------------------------------//
// end of plot2D/Plot2D.hh
//---------------------------------------------------------------------------//
