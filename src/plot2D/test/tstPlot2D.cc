//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   plot2D/test/tstPlot2D.cc
 * \author Rob Lowrie
 * \brief  Plot2D test.
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/XGetopt.hh"
#include "plot2D/Plot2D.hh"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

using rtt_plot2D::Plot2D;
using rtt_plot2D::SetProps;
using std::cin;
using std::cout;
using std::endl;
using std::ofstream;
using std::string;

namespace test {
void pause();
}

void tstPlot2D(const bool batch);
int main(int argc, char *argv[]);
bool pass(true);

//---------------------------------------------------------------------------//
// Helper function
//---------------------------------------------------------------------------//
namespace test {

void pause() {
  cout << "Press RETURN to continue: ";
  cin.get();
}
} // namespace test

//---------------------------------------------------------------------------//
// The main test
//---------------------------------------------------------------------------//
void tstPlot2D(const bool batch) {
  string paramFile("tstPlot2D.par");

  const int numGraphs = 3; // number of graphs in plot
  const int n = 10;        // number of points in each set

  //////////////////////////////////////////////////////////////////////
  // Plot 1.
  //////////////////////////////////////////////////////////////////////

  // Output the block data file that Grace will read

  string blockName("tmp.block");
  ofstream block(blockName.c_str());

  for (int i = 0; i < n; i++) {
    block << i;
    for (int j = 0; j < numGraphs; j++)
      block << " " << std::pow(double(i), j + 1);
    block << endl;
  }

  block.close();

  // Generate the plot

  Plot2D p(numGraphs, paramFile, batch);

  p.readBlock(blockName);

  // Set titles and axis labels

  for (int j = 0; j < numGraphs; j++) {
    std::ostringstream title, subtitle, ylabel;
    title << "title " << j;
    subtitle << "subtitle " << j;
    ylabel << "y" << j;
    p.setTitles(title.str(), subtitle.str(), j);
    p.setAxesLabels("x", ylabel.str(), j);
  }

  // Changes some set properties; most were set in the param file.

  SetProps prop;
  prop.line.color = rtt_plot2D::COLOR_RED;
  p.setProps(0, 0, prop); // of graph 0, set 0

  if (!batch)
    test::pause();

  p.save("plot1.agr");

  //////////////////////////////////////////////////////////////////////
  // Plot 2.  Add a set to each graph of the previous plot.
  //////////////////////////////////////////////////////////////////////

  // Generate new data.  Use a different temporary filename, because Grace may
  // not have read the data from plot 1 yet.

  blockName = "tmp2.block";
  block.open(blockName.c_str());

  for (int i = 0; i < n; i++) {
    block << i;
    for (int j = 0; j < numGraphs; j++)
      block << " " << i + std::pow(double(i), j + 2);
    block << endl;
  }

  block.close();

  // Add the data to the graphs

  p.killAllSets();
  p.readBlock(blockName);

  if (!batch)
    test::pause();

  p.save("plot2.agr");

  //////////////////////////////////////////////////////////////////////
  // Plot 3.  Rearrange graph matrix.
  //////////////////////////////////////////////////////////////////////

  p.arrange(2, 2);

  if (!batch)
    test::pause();

  p.save("plot3.agr");

  //////////////////////////////////////////////////////////////////////
  // Plot 4.  Put all of the data of the previous plot into one graph.
  //////////////////////////////////////////////////////////////////////

  p.close();
  p.open(1, "", batch); // 1 graph, no param file specified

  // Must specify the graph number (0) to put all sets into one
  // graph.
  p.readBlock(blockName, 0);

  p.setTitles("Same Data, One Graph", "subtitle");

  if (!batch)
    test::pause();

  p.save("plot4.agr");
}

//---------------------------------------------------------------------------//
/* Compare the generated plot file to gold standard.
 *
 * Note: The XM Grace plot file has this format:
 *
 * Line(s)  Description
 * -------  -------------
 * 1-2      Header
 * 3        XM-GRACE Version
 * 4-<A>    Meta data (fonts, colors, etc)
 *          This block ends with the line "@timestamp def "..."
 *          The length of this block varies by machine and XM-Grace version.
 * <A>-EOF  Plot data
 *
 * The top portion of these files are often mismatched due to different
 * timestamps, grace versions, fonts, colors, etc. Don't compare these.
 *
 * We need to find the line number ofr the timestamp and only compare text
 * that appears after this line.
 */
//----------------------------------------------------------------------------//
void checkOutputFiles(std::string const &filename) {
  using std::string;

  Check(!filename.empty());

  cout << "\nChecking contents of generated plot file \"" << filename
       << ".agr\"..." << endl;

  string benchfile(string("bench_") + filename);

  // Look at the output file and locate the line number of the timestamp.
  int end_of_metadata_block = 1; // ndselect starts counting at 1.
  // open the file for reading
  std::ifstream rawfile(filename + string(".agr"), std::ifstream::in);
  Insist(rawfile.good(), (string("Could not open file ") + filename).c_str());
  string line; // store each line of the file in this temporary.
  while (getline(rawfile, line)) {
    ++end_of_metadata_block;
    if (line.find("@timestamp def") != string::npos)
      break;
  }
  rawfile.close();

  // Trim the meta-data text from the top of the generated files and then
  // compare the remaining content with gold files.

  std::vector<int> begin = {end_of_metadata_block};
  std::vector<int> end = {-1};

  Check(begin.size() == end.size());
  for (size_t i = 1; i < begin.size(); ++i) {
    Check(begin[i] >= end[i - 1]);
    if (end[i] > 0)
      Check(end[i] >= begin[i]);
  }

  // break the file into portions that do not include these lines and then run
  // numdiff on each portion. Don't compare the first block as this contains
  // font names that can differ between machines.
  for (size_t i = 0; i < begin.size(); ++i) {
    std::ostringstream ndselect_file;
    std::ostringstream ndselect_bench;
    std::ostringstream numdiff_cmd;

    cout << "\nGenerating part " << i << " (lines ";

    // output file
    ndselect_file << "ndselect ";
    if (begin[i] >= 0) {
      ndselect_file << "-b " << begin[i];
      cout << begin[i];
    }
    cout << " - ";
    if (end[i] >= 0) {
      ndselect_file << " -e " << end[i];
      cout << end[i];
    }
    cout << ") for comparison..." << std::endl;
    ndselect_file << " -o " << filename << ".part" << i << " " << filename
                  << ".agr";
    cout << ndselect_file.str() << endl;
    int ret = system(ndselect_file.str().c_str());
    if (ret != 0)
      pass = false;

    // benchmark file (meta data section already pruned)
    ndselect_bench << "ndselect ";
    if (i == 0)
      ndselect_bench << "-b 1";
    else
      ndselect_bench << "-b " << begin[i];
    if (end[i] >= 0)
      ndselect_bench << " -e " << end[i];
    ndselect_bench << " -o "
                   << "bench_" << filename << ".part" << i << " "
                   << "bench_" << filename << ".agr";
    cout << ndselect_bench.str() << endl;
    ret = system(ndselect_bench.str().c_str());
    if (ret != 0)
      pass = false;

    // Do the diff on each file part via numdiff...
    numdiff_cmd << "numdiff " << filename << ".part" << i << " bench_"
                << filename << ".part" << i;
    cout << "\n" << numdiff_cmd.str() << endl;
    ret = system(numdiff_cmd.str().c_str());
    if (ret != 0)
      pass = false;
  }

  return;
}

//---------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  bool batch = true;

  // Process known command line arguments:
  rtt_dsxx::XGetopt::csmap long_options;
  long_options['v'] = "version";
  long_options['g'] = "gui";
  long_options['h'] = "help";
  std::map<char, std::string> help_strings;
  help_strings['h'] = "print this message.";
  help_strings['v'] = "print version information and exit.";
  rtt_dsxx::XGetopt program_options(argc, argv, long_options, help_strings);

  int c(0);
  while ((c = program_options()) != -1) {
    switch (c) {
    case 'v': // --version
      cout << argv[0] << ": version " << rtt_dsxx::release() << endl;
      break;

    case 'h':
      std::cout << program_options.display_help("tstPlot2D") << std::endl;
      return 0;
      break;

    case 'g': // --gui
      batch = false;
      break;

    default:
      break; // nothing to do.
    }
  }

  cout << "\n**********************************************\n";

  try {
    // tests
    if (Plot2D::is_supported()) {
      tstPlot2D(batch);

      // Grace doesn't apper to flush the data to disk immediately.  We need to
      // wait a bit before comparing the output to the gold standards.
      std::cout << "\nWaiting for Grace to finish writing files...\n"
                << std::endl;
      bool system_call_ok = system("sleep 5") == 0;
      Insist(system_call_ok, "system call returned an error.");

      // Check the created files by comparing to a gold standard.
      checkOutputFiles(std::string("plot1"));
      checkOutputFiles(std::string("plot2"));
      checkOutputFiles(std::string("plot3"));
      checkOutputFiles(std::string("plot4"));
    } else {
      cout << "Unsupported test: pass\n";
    }
  } catch (rtt_dsxx::assertion const &error) {
    cout << "Assertion: " << error.what() << "\nBetter luck next time!\n\n"
         << endl;
    pass = false;
  }

  // status of test
  cout << "\n**********************************************\n";
  if (pass)
    cout << "********* Plot2D Self Test: PASSED ***********\n";
  else
    cout << "********* Plot2D Self Test: FAILED ***********\n";
  cout << "**********************************************\n\n"
       << "Done testing Plot2D.\n"
       << endl;
  return 0;
}

//---------------------------------------------------------------------------//
// end of tstPlot2D.cc
//---------------------------------------------------------------------------//
