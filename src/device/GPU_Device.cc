//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   device/GPU_Device.cc
 * \author Kelly (KT) Thompson
 * \date   Thu Oct 20 15:28:48 2011
 * \brief
 * \note   Copyright (C) 2016-2020 Triad National Security, LLC.
 *         All rights reserved. */
//----------------------------------------------------------------------------//

#include "GPU_Device.hh"
#include <iostream>
#include <sstream>

namespace rtt_device {

//----------------------------------------------------------------------------//
// Constructor
//----------------------------------------------------------------------------//
/*!
 * \brief Constructor
 *
 * Create a GPU_Device object.
 * - initialize the CUDA library.
 * - Set device and context handles.
 * - Query the devices for features.
 */
GPU_Device::GPU_Device(void)
    : deviceCount(0), computeCapability(), deviceName() {

  // Get a device count, determine compute capability
  cudaError_t err = cudaGetDeviceCount(&deviceCount);
  checkForCudaError(err);
  Insist(deviceCount > 0, "No GPU devices found!");

  // Collect information about each GPU device found
  computeCapability.resize(deviceCount);

  for (int device = 0; device < deviceCount; device++) {
    int cudaDevice;
    err = cudaSetDevice(device);
    checkForCudaError(err);
    err = cudaGetDevice(&cudaDevice);
    checkForCudaError(err);
    Check(cudaDevice == device);

    // Compute capability revision
    int major = 0;
    int minor = 0;
    cudaDeviceGetAttribute(&major, cudaDevAttrComputeCapabilityMajor, device);
    cudaDeviceGetAttribute(&minor, cudaDevAttrComputeCapabilityMinor, device);
    computeCapability[device].push_back(major);
    computeCapability[device].push_back(minor);

    // Device name
    cudaDeviceProp device_properties;
    err = cudaGetDeviceProperties(&device_properties, device);
    checkForCudaError(err);
    deviceName.push_back(std::string(device_properties.name));

    // Query and archive device properties.
    {
      int tmp(0);
      m_maxthreadsperblock.push_back(device_properties.maxThreadsPerBlock);
      m_maxthreadsdim.push_back(
          std::array<int, 3>{device_properties.maxThreadsDim[0],
                             device_properties.maxThreadsDim[1],
                             device_properties.maxThreadsDim[2]});
      m_maxgridsize.push_back(std::array<int, 3>{
          device_properties.maxGridSize[0], device_properties.maxGridSize[1],
          device_properties.maxGridSize[2]});
      m_sharedmemperblock.push_back(device_properties.sharedMemPerBlock);

      err =
          cudaDeviceGetAttribute(&tmp, cudaDevAttrTotalConstantMemory, device);
      checkForCudaError(err);
      m_totalconstantmemory.push_back(tmp);

      m_simdwidth.push_back(device_properties.warpSize);
      m_mempitch.push_back(device_properties.memPitch);

      err =
          cudaDeviceGetAttribute(&tmp, cudaDevAttrMaxRegistersPerBlock, device);
      checkForCudaError(err);
      m_regsperblock.push_back(tmp);

      err = cudaDeviceGetAttribute(&tmp, cudaDevAttrClockRate, device);
      checkForCudaError(err);
      m_clockrate.push_back(tmp);

      err = cudaDeviceGetAttribute(&tmp, cudaDevAttrTextureAlignment, device);
      checkForCudaError(err);
      m_texturealign.push_back(tmp);
    }
  }
}

GPU_Device::~GPU_Device() {}

//----------------------------------------------------------------------------//
// Print a summary of all GPU devices found
//----------------------------------------------------------------------------//
/*!
 * \brief Print a summary of device features
 * \arg idevice device index
 * \arg out Send output to this ostream (cout by default).
 */
void GPU_Device::printDeviceSummary(int const idevice,
                                    std::ostream &out) const {
  out << "Device: " << idevice
      << "\n   Name               : " << deviceName[idevice]
      << "\n   Compute capability : " << computeCapability[idevice][0] << "."
      << computeCapability[idevice][1]
      << "\n   maxThreadsPerBlock : " << maxThreadsPerBlock(idevice)
      << "\n   maxThreadsDim      : " << maxThreadsDim(0, idevice) << " x "
      << maxThreadsDim(1, idevice) << " x " << maxThreadsDim(2, idevice)
      << "\n   maxGridSize        : " << maxGridSize(0, idevice) << " x "
      << maxGridSize(1, idevice) << " x " << maxGridSize(2, idevice)
      << "\n   sharedMemPerBlock  : " << sharedMemPerBlock(idevice)
      << "\n   totalConstantMemory: " << totalConstantMemory(idevice)
      << "\n   SIMDWidth          : " << SIMDWidth(idevice)
      << "\n   memPitch           : " << memPitch(idevice)
      << "\n   regsPerBlock       : " << regsPerBlock(idevice)
      << "\n   clockRate          : " << clockRate(idevice)
      << "\n   textureAlign       : " << textureAlign(idevice) << "\n"
      << std::endl;
  return;
}

#ifdef DBC
//----------------------------------------------------------------------------//
/*!
 * \brief Convert a CUDA return enum value into a descriptive string.
 *
 * \param errorCode CUDA enum return value
 * \return descriptive string associated with
 *
 * For optimized builds with DRACO_DBC_LEVEL=0, this function will be empty
 * and any decent compiler will optimize this call away.
 */
void GPU_Device::checkForCudaError(cudaError_t const errorCode) {
  std::ostringstream msg;
  msg << "A CUDA call returned the error: \"" << getErrorMessage(errorCode)
      << "\"";
  Insist(errorCode == cudaSuccess, msg.str());
}

#else
//----------------------------------------------------------------------------//
/*!
 * \brief Convert a CUDA return enum value into a descriptive string.
 * \return descriptive string associated with
 */
void GPU_Device::checkForCudaError(cudaError_t const) { /* empty */
}
#endif

//----------------------------------------------------------------------------//
/*!
 * \brief Return a text string that corresponds to a CUDA error enum.
 */
std::string GPU_Device::getErrorMessage(cudaError_t const err) {
  auto raw_message = cudaGetErrorString(err);
  std::string message(raw_message);

  /*
  switch (err) {
  case cudaSuccess:
    message = std::string("No errors.");
    break;
  case cudaErrorInvalidValue:
    message = std::string("Invalid value.");
    break;
  case cudaErrorMemoryAllocation:
    message = std::string("Out of memory.");
    break;
  case cudaErrorInitializationError:
    message = std::string("Driver not initialized.");
    break;
  case cudaErrorDeviceUninitilialized:
    message = std::string("Driver deinitialized.");
    break;
  case cudaErrorNoDevice:
    message = std::string("No CUDA-capable device available.");
    break;
  case cudaErrorInvalidDevice:
    message = std::string("Invalid device.");
    break;
  case cudaErrorInvalidKernelImage:
    message = std::string("Invalid kernel image.");
    break;
  case cudaErrorIncompatibleDriverContext:
    message = std::string("Invalid context.");
    break;
  case cudaErrorMapBufferObjectFailed:
    message = std::string("Map failed.");
    break;
  case cudaErrorUnmapBufferObjectFailed:
    message = std::string("Unmap failed.");
    break;
  case cudaErrorArrayIsMapped:
    message = std::string("Array is mapped.");
    break;
  case cudaErrorAlreadyMapped:
    message = std::string("Already mapped.");
    break;
  case cudaErrorAlreadyAcquired:
    message = std::string("Already acquired.");
    break;
  case cudaErrorNotMapped:
    message = std::string("Not mapped.");
    break;
  case cudaErrorInvalidSource:
    message = std::string("Invalid source.");
    break;
  case cudaErrorFileNotFound:
    message = std::string("File not found.");
    break;
  case cudaErrorInvalidResourceHandle:
    message = std::string("Invalid handle.");
    break;
  case cudaErrorNotReady:
    message = std::string("CUDA not ready.");
    break;
  case cudaErrorLaunchFailure:
    message = std::string("Launch failed.");
    break;
  case cudaErrorLaunchOutOfResources:
    message = std::string("Launch exceeded resources.");
    break;
  case cudaErrorLaunchTimeout:
    message = std::string("Launch exceeded timeout.");
    break;
  case cudaErrorLaunchIncompatibleTexturing:
    message = std::string("Launch with incompatible texturing.");
    break;
  case cudaErrorUnknown:
    message = std::string("Unknown error. ");
    break;
  default:
    // cudaErrorProfilerDisabled
    // cudaErrorNotMappedAsArray
    // cudaErrorNotMappedAsPointer
    // cudaErrorECCUncorrectable
    // cudaErrorUnsupportedLimit
    // cudaErrorSharedObjectSymbolNotFound
    // cudaErrorSharedObjectInitFailed
    // cudaErrorOperatingSystem
    // cudaErrorPeerAccessAlreadyEnabled
    // cudaErrorPeerAccessUnsupported
    // cudaErrorPeerAccessNotEnabled
    // cudaErrorContextIsDestroyed
    message = std::string("Unknown error. ");
    break;
  }
  */
  return message;
}

//----------------------------------------------------------------------------//
/*!
 * \brief Wrap the cudaMemAlloc funtion to include error checking
 *
 * \param nbytes number of bytes to allocate (e.g.: len*sizeof(double) ).
 * \return GPU device pointer to allocated memory.
 */

void *GPU_Device::MemAlloc(unsigned const nbytes) {
  void *ptr;
  cudaError_t err = cudaMalloc(&ptr, nbytes);
  checkForCudaError(err);
  return ptr;
}

void GPU_Device::MemcpyHtoD(void *ptr, void const *loc, unsigned nbytes) {
  cudaError_t err = cudaMemcpy(ptr, loc, nbytes, cudaMemcpyHostToDevice);
  checkForCudaError(err);
  return;
}

void GPU_Device::MemcpyDtoH(void *loc, void *ptr, unsigned nbytes) {
  cudaError_t err = cudaMemcpy(loc, ptr, nbytes, cudaMemcpyDeviceToHost);
  checkForCudaError(err);
  return;
}

void GPU_Device::MemFree(void *ptr) {
  cudaError_t err = cudaFree(ptr);
  checkForCudaError(err);
  return;
}

} // end namespace rtt_device

//----------------------------------------------------------------------------//
// end of GPU_Device.cc
//----------------------------------------------------------------------------//
