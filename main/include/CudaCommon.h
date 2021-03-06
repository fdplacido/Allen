#pragma once

#include <stdexcept>
#include <iostream>
#include "LoggerCommon.h"

#ifdef CPU

#include <cmath>
#include <cstring>

// -----------
// CPU support
// -----------

using std::signbit;

#define __host__
#define __device__
#define __shared__
#define __global__
#define __constant__
#define __syncthreads()
#define __launch_bounds__(_i)
#define cudaError_t int
#define cudaEvent_t int
#define cudaStream_t int
#define cudaSuccess 0
#define cudaErrorMemoryAllocation 2
#define half_t short
#define __popcll __builtin_popcountll
#define cudaEventBlockingSync 0x01

enum cudaMemcpyKind {
  cudaMemcpyHostToHost,
  cudaMemcpyHostToDevice,
  cudaMemcpyDeviceToHost,
  cudaMemcpyDeviceToDevice,
  cudaMemcpyDefault
};

struct float3 {
  float x;
  float y;
  float z;
};

struct float2 {
  float x;
  float y;
};

struct dim3 {
  unsigned int x = 1;
  unsigned int y = 1;
  unsigned int z = 1;

  dim3() = default;
  dim3(const dim3&) = default;

  dim3(const unsigned int& x);
  dim3(const unsigned int& x, const unsigned int& y);
  dim3(const unsigned int& x, const unsigned int& y, const unsigned int& z);
};

struct GridDimensions {
  unsigned int x;
  unsigned int y;
  unsigned int z;
};

struct BlockIndices {
  unsigned int x;
  unsigned int y;
  unsigned int z;
};

struct BlockDimensions {
  unsigned int x = 1;
  unsigned int y = 1;
  unsigned int z = 1;
};

struct ThreadIndices {
  unsigned int x = 0;
  unsigned int y = 0;
  unsigned int z = 0;
};

extern thread_local GridDimensions gridDim;
extern thread_local BlockIndices blockIdx;
constexpr BlockDimensions blockDim {1, 1, 1};
constexpr ThreadIndices threadIdx {0, 0, 0};

cudaError_t cudaMalloc(void** devPtr, size_t size);
cudaError_t cudaMallocHost(void** ptr, size_t size);
cudaError_t cudaMemcpy(void* dst, const void* src, size_t count, enum cudaMemcpyKind kind);
cudaError_t cudaMemcpyAsync(void* dst, const void* src, size_t count, enum cudaMemcpyKind kind, cudaStream_t stream);
cudaError_t cudaMemset(void* devPtr, int value, size_t count);
cudaError_t cudaMemsetAsync(void* devPtr, int value, size_t count, cudaStream_t stream);
cudaError_t cudaPeekAtLastError();
cudaError_t cudaEventCreate(cudaEvent_t* event);
cudaError_t cudaEventCreateWithFlags(cudaEvent_t* event, int flags);
cudaError_t cudaEventSynchronize(cudaEvent_t event);
cudaError_t cudaEventRecord(cudaEvent_t event, cudaStream_t stream);
cudaError_t cudaFreeHost(void* ptr);
cudaError_t cudaDeviceReset();
cudaError_t cudaStreamCreate(cudaStream_t* pStream);
cudaError_t cudaMemcpyToSymbol(
  void* symbol,
  const void* src,
  size_t count,
  size_t offset = 0,
  enum cudaMemcpyKind kind = cudaMemcpyDefault);

template<class T>
cudaError_t cudaMemcpyToSymbol(
  const T& symbol,
  const void* src,
  size_t count,
  size_t offset = 0,
  enum cudaMemcpyKind = cudaMemcpyHostToDevice)
{
  std::memcpy(reinterpret_cast<void*>(((char*) &symbol) + offset), src, count);
  return 0;
}

template<class T>
cudaError_t cudaMemcpyFromSymbol(
  void* dst,
  const T& symbol,
  size_t count,
  size_t offset = 0,
  enum cudaMemcpyKind = cudaMemcpyHostToDevice)
{
  std::memcpy(dst, reinterpret_cast<void*>(((char*) &symbol) + offset), count);
  return 0;
}

template<class T, class S>
T atomicAdd(T* address, S val)
{
  const T old = *address;
  *address += val;
  return old;
}

template<class T, class S>
T atomicOr(T* address, S val)
{
  const T old = *address;
  *address |= val;
  return old;
}

template<class T>
T max(const T& a, const T& b)
{
  return std::max(a, b);
}

template<class T>
T min(const T& a, const T& b)
{
  return std::min(a, b);
}

unsigned int atomicInc(unsigned int* address, unsigned int val);

half_t __float2half(float value);

#define cudaCheck(stmt)                                    \
  {                                                        \
    cudaError_t err = stmt;                                \
    if (err != cudaSuccess) {                              \
      std::cerr << "Failed to run " << #stmt << std::endl; \
      throw std::invalid_argument("cudaCheck failed");     \
    }                                                      \
  }

#define cudaCheckKernelCall(stmt, kernel_name)                      \
  {                                                                 \
    cudaError_t err = stmt;                                         \
    if (err != cudaSuccess) {                                       \
      std::cerr << "Failed to invoke " << kernel_name << std::endl; \
      throw std::invalid_argument("cudaCheckKernelCall failed");    \
    }                                                               \
  }

namespace Configuration {
  extern uint verbosity_level;
}

#elif defined(HIP)

// ---------------
// Support for HIP
// ---------------

#if defined(__HCC__) || defined(__HIP__) || defined(__NVCC__) || defined(__CUDACC__)
#include <hip/hip_runtime.h>
#else
#define __HIP_PLATFORM_HCC__
#include <hip/hip_runtime_api.h>
#endif

// Support for CUDA to HIP conversion
#define cudaMalloc hipMalloc
#define cudaMallocHost hipHostMalloc
#define cudaMemcpy hipMemcpy
#define cudaMemcpyAsync hipMemcpyAsync
#define cudaMemset hipMemset
#define cudaMemsetAsync hipMemsetAsync
#define cudaPeekAtLastError hipPeekAtLastError
#define cudaEventCreate hipEventCreate
#define cudaEventCreateWithFlags hipEventCreateWithFlags
#define cudaEventSynchronize hipEventSynchronize
#define cudaEventRecord hipEventRecord
#define cudaFreeHost hipHostFree
#define cudaDeviceReset hipDeviceReset
#define cudaStreamCreate hipStreamCreate
#define cudaMemGetInfo hipMemGetInfo
#define cudaGetDeviceCount hipGetDeviceCount
#define cudaSetDevice hipSetDevice
#define cudaGetDeviceProperties hipGetDeviceProperties
#define cudaDeviceProp hipDeviceProp_t

#define cudaError_t hipError_t
#define cudaEvent_t hipEvent_t
#define cudaStream_t hipStream_t
#define cudaSuccess hipSuccess
#define cudaErrorMemoryAllocation hipErrorMemoryAllocation
#define cudaEventBlockingSync hipEventBlockingSync

#define cudaMemcpyHostToHost hipMemcpyHostToHost
#define cudaMemcpyHostToDevice hipMemcpyHostToDevice
#define cudaMemcpyDeviceToHost hipMemcpyDeviceToHost
#define cudaMemcpyDeviceToDevice hipMemcpyDeviceToDevice
#define cudaMemcpyDefault hipMemcpyDefault

#define cudaCheck(stmt)                                                                                           \
  {                                                                                                               \
    hipError_t err = stmt;                                                                                        \
    if (err != hipSuccess) {                                                                                      \
      fprintf(                                                                                                    \
        stderr, "Failed to run %s\n%s (%d) at %s: %d\n", #stmt, hipGetErrorString(err), err, __FILE__, __LINE__); \
      throw std::invalid_argument("cudaCheck failed");                                                            \
    }                                                                                                             \
  }

#define cudaCheckKernelCall(stmt, kernel_name)                   \
  {                                                              \
    cudaError_t err = stmt;                                      \
    if (err != cudaSuccess) {                                    \
      fprintf(                                                   \
        stderr,                                                  \
        "Failed to invoke %s\n%s (%d) at %s: %d\n",              \
        kernel_name.c_str(),                                     \
        hipGetErrorString(err),                                  \
        err,                                                     \
        __FILE__,                                                \
        __LINE__);                                               \
      throw std::invalid_argument("cudaCheckKernelCall failed"); \
    }                                                            \
  }

#define half_t short

__device__ __host__ half_t __float2half(float value);

namespace Configuration {
  extern __constant__ uint verbosity_level;
}

#else

// ------------
// CUDA support
// ------------

#include "cuda_runtime.h"
#include <mma.h>
#define half_t half

/**
 * @brief Macro to check cuda calls.
 */
#define cudaCheck(stmt)                                    \
  {                                                        \
    cudaError_t err = stmt;                                \
    if (err != cudaSuccess) {                              \
      std::cerr << "Failed to run " << #stmt << std::endl; \
      std::cerr << cudaGetErrorString(err) << std::endl;   \
      throw std::invalid_argument("cudaCheck failed");     \
    }                                                      \
  }

#define cudaCheckKernelCall(stmt, kernel_name)                      \
  {                                                                 \
    cudaError_t err = stmt;                                         \
    if (err != cudaSuccess) {                                       \
      std::cerr << "Failed to invoke " << kernel_name << std::endl; \
      std::cerr << cudaGetErrorString(err) << std::endl;            \
      throw std::invalid_argument("cudaCheckKernelCall failed");    \
    }                                                               \
  }

namespace Configuration {
  extern __constant__ uint verbosity_level;
}

#endif

// Replacement for gsl::span in CUDA code
namespace cuda {
  template<class T>
  struct span {
    T* __ptr;
    size_t __size;

    __device__ __host__ T* data() const { return __ptr; }
    __device__ __host__ size_t size() const { return __size; }
    __device__ __host__ T& operator[](int i) { return __ptr[i]; }
    __device__ __host__ const T operator[](int i) const { return __ptr[i]; }
  };
} // namespace cuda

/**
 * @brief Cross architecture for statement.
 * @details It can be used to iterate with variable _TYPE _I from 0 through _END.
 */
#ifdef __CUDA_ARCH__
#define FOR_STATEMENT(_TYPE, _I, _END) for (_TYPE _I = threadIdx.x; _I < _END; _I += blockDim.x)
#else
#define FOR_STATEMENT(_TYPE, _I, _END) for (_TYPE _I = 0; _I < _END; ++_I)
#endif

/**
 * @brief Macro to avoid warnings on Release builds with variables used by asserts.
 */
#define _unused(x) ((void) (x))

void print_gpu_memory_consumption();

std::tuple<bool, std::string> set_device(int cuda_device, size_t stream_id);
