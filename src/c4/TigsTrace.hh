//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/TigsTrace.hh
 * \author Kelly Thompson, Bob Webster, Thomas M. Evans
 * \date   Fri Sep 30 12:53:25 2005
 * \brief  TigsTrace class definitions.
 * \note   Copyright (C) 2016 Los Alamos National Security, LLC.
 *         All rights reserved.
 *
 * These comm classes are used to perform gather/scatter operations.  They
 * were originally developed by Bob Webster for the Antero Project.  They have
 * been imported here and refactored for use in wedgehog census-mapping
 * operations.  They are used with permission of the original author.
 */
//---------------------------------------------------------------------------//
// $Id: Tigs.hh 6056 2012-06-19 19:05:27Z kellyt $
//---------------------------------------------------------------------------//

#ifndef c4_TigsTrace_hh
#define c4_TigsTrace_hh

#include "ds++/config.h"
#include <map>
#include <vector>

namespace rtt_c4
{

//===========================================================================//
/*!
 * \class TigsTrace
 *
 * This class allows users to perform gather/scatter operations on
 * user-supplied traces.
 */
//===========================================================================//

class DLL_PUBLIC_c4 TigsTrace
{
  public:
    // Useful typedefs.
    typedef std::map<int, std::vector<int> > TigsComm_map;
    typedef std::vector<std::vector<int> >   VVec_int;

    // Enumerations.
    enum GS_Tag
    {
        ANY_TAG    = 0,
        OK         = 1,
        FATAL      = 2,
        WARN       = 3,
        TRACE_INIT = 4
    };

  private:
    // >>> DATA

    // These variables store the "local" sizes of the range and domain of
    // the trace map.
    unsigned const onProcDomain;
    unsigned const onProcRange;

    // These variables store the indirection information for the
    // serial, or single processor case.
    std::map<int, std::vector<int> > IM;
    std::vector<int> IMV;
    std::vector<unsigned> counts;

    // These variables store the indirection information for the parallel,
    // multiprocessor case. The Domain of the trace map is refered to as
    // the Iside, and the Range of the map is the Jside. Where the map
    // should be thought of as Jmap:I->J
    //
    // The connects vectors store the processor to which a message must be
    // communicated to fill or send the data for the gather/scatter
    //
    // the Indirect vectors store the indirection vector to load/unload the
    // communciation buffer with data from the arrays provided to the
    // gather/scatter call.
    //
    // #####Connects[s]    is the s-th processor to communicate with
    // (send/receive one buffer)
    //
    // #####Indirect[s] is a vector to [un]load the buffer to be
    // communicated with the above processor
    //
    // BmapList is an indirection array to [un]pack the scatterList result
    // following communicationint a CSR like structure. countsList is
    // The count of items in the CSR like list for each range location
    // (on the current processor)
    //
    unsigned          IsideBufferSize;
    unsigned          JsideBufferSize;
    std::vector<int>  IsideConnects;
    VVec_int          IsideIndirect;
    std::vector<int>  JsideConnects;
    VVec_int          JsideIndirect;
    std::vector<int>  BmapList;
    std::vector<int>  countsList;

  public:
    // Constructor.
    TigsTrace(std::vector<int> const & JMap, unsigned const J);

    // >>> GATHER OPERATIONS

    // Copy data from the range to the domain.
    template < typename iterA, typename iterB >
    void gather( iterB Bfirst, iterB Blast,
                 iterA Afirst, iterA Alast );

    // >>> SCATTER OPERATIONS

    // Copy data defined on the domain (I) to the range (J).
    template < typename iterA , typename iterC, typename iterB >
    void scatterList( iterA Afirst,   iterA Alast,
                      iterC Cntfirst, iterC Cntlast,
                      iterB Bfirst,   iterB Blast ) const;

    // Copy with reduction of data from the domain to the range.
    template < typename iterA , typename iterB , typename BinaryOp >
    void scatter( iterA Afirst, iterA Alast,
                  iterB Bfirst, iterB Blast,
                  BinaryOp op );

    // >>> QUERIES

    /*!
     * \brief A routine to determine allocation size for the scatterList
     * output.
     *
     * The scatterList data size may not be obvious from the map input to
     * construct a trace. This method returns the resulting data container
     * size.
     *
     * \return The size of the container that holds scatterList data
     */
    int getListSize() const { return JsideBufferSize; }
};

} // end namespace rtt_c4

//---------------------------------------------------------------------------//
// TEMPLATE MEMBERS
//---------------------------------------------------------------------------//

#include "TigsTrace.i.hh"

#endif // c4_TigsTrace_hh

//---------------------------------------------------------------------------//
// end of c4/Tigs.hh
//---------------------------------------------------------------------------//
