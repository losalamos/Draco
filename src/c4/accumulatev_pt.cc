//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/accumulatev_pt.cc
 * \author Kelly Thompson
 * \date   Wednesday, Oct 31, 2012, 14:26 pm
 * \brief  C4 MPI accumulate instantiations.
 */
//---------------------------------------------------------------------------//
// $Id: gatherv_pt.cc 6273 2011-11-22 19:58:32Z kellyt $
//---------------------------------------------------------------------------//

// #include <c4/config.h>
// // #include "C4_Functions.hh"
// // #include "C4_Req.hh"
// #include "accumulatev.t.hh"
// #include "Functors.hh"

// namespace rtt_c4
// {
// using std::vector;

// //---------------------------------------------------------------------------//
// // EXPLICIT INSTANTIATIONS OF NON-BLOCKING SEND/RECEIVE
// //---------------------------------------------------------------------------//

// template
// void accumulatev<int>( int* localBegin, int* localEnd,
//                        int init, rtt_c4::max<int>() );
    

// template
// void indeterminate_gatherv<unsigned>(
//     std::vector<unsigned> &outgoing_data,
//     std::vector<std::vector<unsigned> > &incoming_data);

// template
// void indeterminate_gatherv<char>(
//     std::vector<char> &outgoing_data,
//     std::vector<std::vector<char> > &incoming_data);

//---------------------------------------------------------------------------//

// template
// void determinate_gatherv<unsigned>(
//     std::vector<unsigned> &outgoing_data,
//     std::vector<std::vector<unsigned> > &incoming_data);

// template
// void determinate_gatherv<int>(
//     std::vector<int> &outgoing_data,
//     std::vector<std::vector<int> > &incoming_data);

// template
// void determinate_gatherv<double>(
//     std::vector<double> &outgoing_data,
//     std::vector<std::vector<double> > &incoming_data);

// template
// void determinate_gatherv<char>(
//     std::vector<char> &outgoing_data,
//     std::vector<std::vector<char> > &incoming_data);

// } // end namespace rtt_c4


//---------------------------------------------------------------------------//
//                              end of C4_MPI_gatherv_pt.cc
//---------------------------------------------------------------------------//
