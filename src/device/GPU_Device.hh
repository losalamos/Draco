//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   device/GPU_Device.hh
 * \author Kelly (KT) Thompson
 * \brief  Define class GPU_Device
 * \note   Copyright (C) 2011 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef device_GPU_Device_hh
#define device_GPU_Device_hh

#include <vector>
#include <iostream>
#include <string>

// All this garbage suppresses warnings found in "cuda.h".
// http://wiki.services.openoffice.org/wiki/Writing_warning-free_code#When_all_else_fails
#if defined __GNUC__
#pragma GCC system_header
#elif defined __SUNPRO_CC
#pragma disable_warn
#elif defined _MSC_VER
#pragma warning(push, 1)
#endif
#include <cuda.h>
#if defined __SUNPRO_CC
#pragma enable_warn
#elif defined _MSC_VER
#pragma warning(pop)
#endif

namespace rtt_device
{

//===========================================================================//
/*!
 * \class GPU_Device
 * \brief
 *
 * Long description or discussion goes here.  Information about Doxygen
 * commands can be found at http://www.doxygen.org.
 *
 * \sa GPU_Device.cc for detailed descriptions.
 *
 * \par Code Sample:
 * \code
 *     cout << "Hello, world." << endl;
 * \endcode
 */
/*! 
 * \example device/test/tstGPU_Device.cc
 *
 * Test of GPU_Device.
 */
//===========================================================================//

class GPU_Device 
{
  public:

    // NESTED CLASSES AND TYPEDEFS

    // CREATORS
    
    //! Default constructors.
    GPU_Device(void);
    // GPU_Device( int /*argc*/, char */*argv*/[] );

    //! Copy constructor (the long doxygen description is in the .cc file).
    // GPU_Device(const GPU_Device &rhs);

    //! Destructor.
    // ~GPU_Device();

    // MANIPULATORS
    
    //! Assignment operator for GPU_Device.
    // GPU_Device& operator=(const GPU_Device &rhs);

    // ACCESSORS

    //! maximum number of threads per block
    int maxThreadsPerBlock(int devId=0) const {
        return deviceProperties[devId].maxThreadsPerBlock; }
    //! maxThreadsDim[3] is the maximum sizes of each dimension of a block.
    int maxThreadsDim(int dim, int devId=0) const {
        return deviceProperties[devId].maxThreadsDim[dim]; }
    //! maxGridSize[3] is the maximum sizes of each dimension of a grid;
    int maxGridSize(int dim, int devId=0) const {
        return deviceProperties[devId].maxGridSize[dim]; }
    //! the total amount of shared memory available per block in bytes
    int sharedMemPerBlock(int devId=0) const {
        return deviceProperties[devId].sharedMemPerBlock; }
    //! the total amount of constant memory available on the device in bytes;
    int totalConstantMemory(int devId=0) const {
        return deviceProperties[devId].totalConstantMemory; }
    //! the warp size
    int SIMDWidth(int devId=0) const {
        return deviceProperties[devId].SIMDWidth; }
    /*! the maximum pitch allowed by the memory copy functions that involve
     *  memory regions allocated through cuMemAllocPitch() */
    int memPitch(int devId=0) const {
        return deviceProperties[devId].memPitch; }
    //! the total number of registers available per block
    int regsPerBlock(int devId=0) const {
        return deviceProperties[devId].regsPerBlock; }
    //! the clock frequency in kilohertz
    int clockRate(int devId=0) const {
        return deviceProperties[devId].clockRate; }
    /*! the alignment requirement; texture base addresses that are aligned to
     *  textureAlign bytes do not need an offset applied to texture fetches */
    int textureAlign(int devId=0) const {
        return deviceProperties[devId].textureAlign; }

    // SERVICES
    void printDeviceSummary(std::ostream & out = std::cout) const;

    // STATICS

  protected:

    // IMPLEMENTATION

  private:

    // NESTED CLASSES AND TYPEDEFS

    // IMPLEMENTATION

    // DATA
    int deviceCount;
    std::vector< std::vector<size_t> > computeCapability;
    std::vector< std::string > deviceName;
    std::vector< CUdevprop_st > deviceProperties;

};

} // end namespace rtt_device

#endif // device_GPU_Device_hh

//---------------------------------------------------------------------------//
//              end of device/GPU_Device.hh
//---------------------------------------------------------------------------//
