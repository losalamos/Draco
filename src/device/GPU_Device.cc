//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   device/GPU_Device.cc
 * \author Kelly (KT) Thompson
 * \date   Thu Oct 20 15:28:48 2011
 * \brief
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "GPU_Device.hh"
#include <sstream>

namespace rtt_device {

//---------------------------------------------------------------------------//
// Constructor
//---------------------------------------------------------------------------//
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
  // Initialize the library
  cudaError_enum err = cuInit(0); // currently must be 0.
  checkForCudaError(err);

  // Get a device count, determine compute capability
  err = cuDeviceGetCount(&deviceCount);
  checkForCudaError(err);
  Insist(deviceCount > 0, "No GPU devices found!");

  // Collect information about each GPU device found
  computeCapability.resize(deviceCount);

  for (int device = 0; device < deviceCount; device++) {
    CUdevice cuDevice;
    err = cuDeviceGet(&cuDevice, device);
    checkForCudaError(err);

    // Compute capability revision
    int major = 0;
    int minor = 0;
    cuDeviceGetAttribute(&major, CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR,
                         cuDevice);
    cuDeviceGetAttribute(&minor, CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MINOR,
                         cuDevice);
    computeCapability[device].push_back(major);
    computeCapability[device].push_back(minor);

    // Device name
    char name[200];
    err = cuDeviceGetName(name, 200, cuDevice);
    checkForCudaError(err);
    deviceName.push_back(std::string(name));

    // Query and archive device properties.
    {
      int tmp1(0), tmp2(0), tmp3(0);
      err = cuDeviceGetAttribute(
          &tmp1, CU_DEVICE_ATTRIBUTE_MAX_THREADS_PER_BLOCK, cuDevice);
      checkForCudaError(err);
      m_maxthreadsperblock.push_back(tmp1);

      err = cuDeviceGetAttribute(&tmp1, CU_DEVICE_ATTRIBUTE_MAX_BLOCK_DIM_X,
                                 cuDevice);
      checkForCudaError(err);
      err = cuDeviceGetAttribute(&tmp2, CU_DEVICE_ATTRIBUTE_MAX_BLOCK_DIM_Y,
                                 cuDevice);
      checkForCudaError(err);
      err = cuDeviceGetAttribute(&tmp3, CU_DEVICE_ATTRIBUTE_MAX_BLOCK_DIM_Z,
                                 cuDevice);
      checkForCudaError(err);
      m_maxthreadsdim.push_back(std::array<int, 3>{tmp1, tmp2, tmp3});

      err = cuDeviceGetAttribute(&tmp1, CU_DEVICE_ATTRIBUTE_MAX_GRID_DIM_X,
                                 cuDevice);
      checkForCudaError(err);
      err = cuDeviceGetAttribute(&tmp2, CU_DEVICE_ATTRIBUTE_MAX_GRID_DIM_Y,
                                 cuDevice);
      checkForCudaError(err);
      err = cuDeviceGetAttribute(&tmp3, CU_DEVICE_ATTRIBUTE_MAX_GRID_DIM_Z,
                                 cuDevice);
      checkForCudaError(err);
      m_maxgridsize.push_back(std::array<int, 3>{tmp1, tmp2, tmp3});

      err = cuDeviceGetAttribute(
          &tmp1, CU_DEVICE_ATTRIBUTE_MAX_SHARED_MEMORY_PER_BLOCK, cuDevice);
      checkForCudaError(err);
      m_sharedmemperblock.push_back(tmp1);

      err = cuDeviceGetAttribute(
          &tmp1, CU_DEVICE_ATTRIBUTE_TOTAL_CONSTANT_MEMORY, cuDevice);
      checkForCudaError(err);
      m_totalconstantmemory.push_back(tmp1);

      err =
          cuDeviceGetAttribute(&tmp1, CU_DEVICE_ATTRIBUTE_WARP_SIZE, cuDevice);
      checkForCudaError(err);
      m_simdwidth.push_back(tmp1);

      err =
          cuDeviceGetAttribute(&tmp1, CU_DEVICE_ATTRIBUTE_MAX_PITCH, cuDevice);
      checkForCudaError(err);
      m_mempitch.push_back(tmp1);

      err = cuDeviceGetAttribute(
          &tmp1, CU_DEVICE_ATTRIBUTE_MAX_REGISTERS_PER_BLOCK, cuDevice);
      checkForCudaError(err);
      m_regsperblock.push_back(tmp1);

      err =
          cuDeviceGetAttribute(&tmp1, CU_DEVICE_ATTRIBUTE_CLOCK_RATE, cuDevice);
      checkForCudaError(err);
      m_clockrate.push_back(tmp1);

      err = cuDeviceGetAttribute(&tmp1, CU_DEVICE_ATTRIBUTE_TEXTURE_ALIGNMENT,
                                 cuDevice);
      checkForCudaError(err);
      m_texturealign.push_back(tmp1);
    }

    // CUdevprop_st properties;
    // err = cuDeviceGetProperties(&properties, cuDevice);
    // checkForCudaError(err);
    // deviceProperties.push_back(properties);
  }

  // Save the handle and context for each device
  device_handle.resize(2);
  context.resize(2);
  for (int device = 0; device < deviceCount; device++) {
    // Only initialize if compute capability >= 2.0
    if (computeCapability[device][0] >= 2) {
      // Save the handle for each device
      err = cuDeviceGet(&device_handle[device], device);
      checkForCudaError(err);

      // Save the handle for each context
      err = cuCtxCreate(&context[device], device, device_handle[device]);
      checkForCudaError(err);
    }
  }
}

/*!
 * \brief destructor
 *
 * Free the device context and unload any modules.
 */
GPU_Device::~GPU_Device() {
  // Free reserved contexts:
  for (int device = 0; device < deviceCount; device++) {
    // Only initialize if compute capability >= 2.0
    if (computeCapability[device][0] >= 2) {
      cudaError_enum err = cuCtxDestroy(context[device]);
      checkForCudaError(err);
    }
  }
}

//---------------------------------------------------------------------------//
// Print a summary of all GPU devices found
//---------------------------------------------------------------------------//
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
//---------------------------------------------------------------------------//
/*!
 * \brief Convert a CUDA return enum value into a descriptive string.
 *
 * \param errorCode CUDA enum return value
 * \return descriptive string associated with
 *
 * For optimized builds with DRACO_DBC_LEVEL=0, this function will be empty
 * and any decent compiler will optimize this call away.
 */
void GPU_Device::checkForCudaError(cudaError_enum const errorCode) {
  std::ostringstream msg;
  msg << "A CUDA call returned the error: \"" << getErrorMessage(errorCode)
      << "\"";
  Insist(errorCode == CUDA_SUCCESS, msg.str());
}

#else
//---------------------------------------------------------------------------//
/*!
 * \brief Convert a CUDA return enum value into a descriptive string.
 * \return descriptive string associated with
 */
void GPU_Device::checkForCudaError(cudaError_enum const) { /* empty */
}
#endif

//---------------------------------------------------------------------------//
/*!
 * \brief Return a text string that corresponds to a CUDA error enum.
 */
std::string GPU_Device::getErrorMessage(cudaError_enum const err) {
  std::string message;
  switch (err) {
  case CUDA_SUCCESS:
    message = std::string("No errors.");
    break;
  case CUDA_ERROR_INVALID_VALUE:
    message = std::string("Invalid value.");
    break;
  case CUDA_ERROR_OUT_OF_MEMORY:
    message = std::string("Out of memory.");
    break;
  case CUDA_ERROR_NOT_INITIALIZED:
    message = std::string("Driver not initialized.");
    break;
  case CUDA_ERROR_DEINITIALIZED:
    message = std::string("Driver deinitialized.");
    break;
  case CUDA_ERROR_NO_DEVICE:
    message = std::string("No CUDA-capable device available.");
    break;
  case CUDA_ERROR_INVALID_DEVICE:
    message = std::string("Invalid device.");
    break;
  case CUDA_ERROR_INVALID_IMAGE:
    message = std::string("Invalid kernel image.");
    break;
  case CUDA_ERROR_INVALID_CONTEXT:
    message = std::string("Invalid context.");
    break;
  case CUDA_ERROR_CONTEXT_ALREADY_CURRENT:
    message = std::string("Context already current.");
    break;
  case CUDA_ERROR_MAP_FAILED:
    message = std::string("Map failed.");
    break;
  case CUDA_ERROR_UNMAP_FAILED:
    message = std::string("Unmap failed.");
    break;
  case CUDA_ERROR_ARRAY_IS_MAPPED:
    message = std::string("Array is mapped.");
    break;
  case CUDA_ERROR_ALREADY_MAPPED:
    message = std::string("Already mapped.");
    break;
  case CUDA_ERROR_NO_BINARY_FOR_GPU:
    message = std::string("No binary for GPU.");
    break;
  case CUDA_ERROR_ALREADY_ACQUIRED:
    message = std::string("Already acquired.");
    break;
  case CUDA_ERROR_NOT_MAPPED:
    message = std::string("Not mapped.");
    break;
  case CUDA_ERROR_INVALID_SOURCE:
    message = std::string("Invalid source.");
    break;
  case CUDA_ERROR_FILE_NOT_FOUND:
    message = std::string("File not found.");
    break;
  case CUDA_ERROR_INVALID_HANDLE:
    message = std::string("Invalid handle.");
    break;
  case CUDA_ERROR_NOT_FOUND:
    message = std::string("Not found.");
    break;
  case CUDA_ERROR_NOT_READY:
    message = std::string("CUDA not ready.");
    break;
  case CUDA_ERROR_LAUNCH_FAILED:
    message = std::string("Launch failed.");
    break;
  case CUDA_ERROR_LAUNCH_OUT_OF_RESOURCES:
    message = std::string("Launch exceeded resources.");
    break;
  case CUDA_ERROR_LAUNCH_TIMEOUT:
    message = std::string("Launch exceeded timeout.");
    break;
  case CUDA_ERROR_LAUNCH_INCOMPATIBLE_TEXTURING:
    message = std::string("Launch with incompatible texturing.");
    break;
  case CUDA_ERROR_UNKNOWN:
    message = std::string("Unknown error. ");
    break;
  default:
    // CUDA_ERROR_PROFILER_DISABLED
    // CUDA_ERROR_PROFILER_NOT_INITIALIZED
    // CUDA_ERROR_PROFILER_ALREADY_STARTED
    // CUDA_ERROR_PROFILER_ALREADY_STOPPED
    // CUDA_ERROR_NOT_MAPPED_AS_ARRAY
    // CUDA_ERROR_NOT_MAPPED_AS_POINTER
    // CUDA_ERROR_ECC_UNCORRECTABLE
    // CUDA_ERROR_UNSUPPORTED_LIMIT
    // CUDA_ERROR_CONTEXT_ALREADY_IN_USE
    // CUDA_ERROR_SHARED_OBJECT_SYMBOL_NOT_FOUND
    // CUDA_ERROR_SHARED_OBJECT_INIT_FAILED
    // CUDA_ERROR_OPERATING_SYSTEM
    // CUDA_ERROR_PEER_ACCESS_ALREADY_ENABLED
    // CUDA_ERROR_PEER_ACCESS_NOT_ENABLED
    // CUDA_ERROR_PEER_MEMORY_ALREADY_REGISTERED
    // CUDA_ERROR_PEER_MEMORY_NOT_REGISTERED
    // CUDA_ERROR_PRIMARY_CONTEXT_ACTIVE
    // CUDA_ERROR_CONTEXT_IS_DESTROYED
    message = std::string("Unknown error. ");
    break;
  }
  return message;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Wrap the cuMemAlloc funtion to include error checking
 *
 * \param nbytes number of bytes to allocate (e.g.: len*sizeof(double) ).
 * \return GPU device pointer to allocated memory.
 */

CUdeviceptr GPU_Device::MemAlloc(unsigned const nbytes) {
  CUdeviceptr ptr;
  cudaError_enum err = cuMemAlloc(&ptr, nbytes);
  checkForCudaError(err);
  return ptr;
}

void GPU_Device::MemcpyHtoD(CUdeviceptr ptr, void const *loc, unsigned nbytes) {
  cudaError_enum err = cuMemcpyHtoD(ptr, loc, nbytes);
  checkForCudaError(err);
  return;
}

void GPU_Device::MemcpyDtoH(void *loc, CUdeviceptr ptr, unsigned nbytes) {
  cudaError_enum err = cuMemcpyDtoH(loc, ptr, nbytes);
  checkForCudaError(err);
  return;
}

void GPU_Device::MemFree(CUdeviceptr ptr) {
  cudaError_enum err = cuMemFree(ptr);
  checkForCudaError(err);
  return;
}

} // end namespace rtt_device

//---------------------------------------------------------------------------//
// end of GPU_Device.cc
//---------------------------------------------------------------------------//
