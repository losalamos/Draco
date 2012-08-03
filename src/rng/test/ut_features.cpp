/*
Copyright 2010-2011, D. E. Shaw Research.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

* Redistributions of source code must retain the above copyright
  notice, this list of conditions, and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions, and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

* Neither the name of D. E. Shaw Research nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
// This "unit test" is basically a test of the completeness
// of compilerfeatures.hpp.  Each of the pp-symbols in compilerfeatures.hpp
// is supposed to have a definition.  We check them all, and
// in some cases, emit some appropriate code to check that
// they reflect reality.
#include <assert.h>
#include <Random123/features/compilerfeatures.h>

#ifndef R123_USE_X86INTRIN_H
#error "No  R123_USE_X86INTRIN_H"
#endif
#if R123_USE_X86INTRIN_H
#include <x86intrin.h>
#endif

#ifndef R123_USE_IA32INTRIN_H
#error "No  R123_USE_IA32INTRIN_H"
#endif
#if R123_USE_IA32INTRIN_H
#include <ia32intrin.h>
#endif

#ifndef R123_USE_EMMINTRIN_H
#error "No  R123_USE_EMMINTRIN_H"
#endif
#if R123_USE_EMMINTRIN_H
#include <emmintrin.h>
#endif

#ifndef R123_USE_SMMINTRIN_H
#error "No  R123_USE_SMMINTRIN_H"
#endif
#if R123_USE_SMMINTRIN_H
#include <smmintrin.h>
#endif

#ifndef R123_USE_WMMINTRIN_H
#error "No  R123_USE_WMMINTRIN_H"
#endif
#if R123_USE_WMMINTRIN_H
#include <wmmintrin.h>
#endif

#ifndef R123_USE_INTRIN_H
#error "No  R123_USE_INTRIN_H"
#endif
#if R123_USE_INTRIN_H
#include <intrin.h>
#endif

#ifndef R123_USE_SSE
#error "No  R123_USE_SSE"
#endif
#if R123_USE_SSE
#include <Random123/features/sse.h>
__m128i mm;
#endif

#ifndef R123_CUDA_DEVICE
#error "No  R123_CUDA_DEVICE"
#endif
R123_CUDA_DEVICE void cuda_device_func(){}

#ifndef R123_USE_CXX0X
#error NO R123_USE_CXX0X
#endif
#if R123_USE_CXX0X
// We *could* check for a lot of things here.
// We do check for:
//   static_assert
//   <tuple>
//   auto
//   decltype
static_assert(1, "always true");
#include <tuple>
bool cxx0x(){
    auto t = std::make_tuple(3, 5, 9);
    decltype(t) s;
    return std::get<0>(t) == std::get<0>(s);
}
#endif

#ifndef R123_FORCE_INLINE
#error "No  R123_FORCE_INLINE"
#endif
inline R123_FORCE_INLINE(int forcibly_inlined(int i));
inline int forcibly_inlined(int i){ return i+1;}

#ifndef R123_USE_AES_NI
#error "No  R123_USE_AES_NI"
#endif
#if R123_USE_AES_NI
__m128i aes(__m128i in){
    if( haveAESNI() )
        return _mm_aesenc_si128(in, in);
    else
        return _mm_setzero_si128();
}
#endif

#ifndef R123_USE_SSE4_2
#error "No  R123_USE_SSE4_2"
#endif
#if R123_USE_SSE4_2
__m128i sse42(__m128i in){
    return _mm_cmpgt_epi64(in, in);
}
#endif

#ifndef R123_USE_SSE4_1
#error "No  R123_USE_SSE4_1"
#endif
#if R123_USE_SSE4_1
int sse41(__m128i in){
    return _mm_testz_si128(in, in);
}
#endif

#ifndef R123_USE_STD_RANDOM
#error "No  R123_USE_STD_RANDOM"
#endif
#if !R123_USE_STD_RANDOM
#else
#include <random>
std::seed_seq ss;
#endif

#ifndef R123_USE_AES_OPENSSL
#error "No  R123_USE_AES_OPENSSL"
#endif
#if R123_USE_AES_OPENSSL
#include <openssl/aes.h>
#endif

#ifndef R123_USE_GNU_UINT128
#error "No  R123_USE_GNU_UINT128"
#endif
#if R123_USE_GNU_UINT128
__uint128_t u128;
#endif

#ifndef R123_USE_ASM_GNU
#error "No  R123_USE_ASM_GNU"
#endif
#if R123_USE_ASM_GNU
R123_STATIC_INLINE int anotherAESNI(){
    unsigned int eax, ebx, ecx, edx;
    __asm__ __volatile__ ("cpuid": "=a" (eax), "=b" (ebx), "=c" (ecx), "=d" (edx) :
                      "a" (1));
    return (ecx>>25) & 1;
}
#endif

#ifndef R123_USE_CPUID_MSVC
#error "No  R123_USE_CPUID_MSVC"
#endif
#if R123_USE_CPUID_MSVC
int chkcpuid(){
    int CPUInfo[4];
    __cpuid(CPUInfo, 1);
    return CPUInfo[2]&(1<<25);
}
#endif

#ifndef R123_USE_MULHILO32_ASM
#error "No  R123_USE_MULHILO32_ASM"
#endif

#ifndef R123_USE_MULHILO64_ASM
#error "No  R123_USE_MULHILO64_ASM"
#endif

#ifndef R123_USE_MULHILO64_MSVC_INTRIN
#error "No  R123_USE_MULHILO_MSVC_INTRIN"
#endif
#if R123_USE_MULHILO64_MSVC_INTRIN
#include <cstdint>
void msvc64mul(){
    uint64_t a=1000000000000000000;
    uint64_t b=a;
    uint64_t h, l;
    l = _umul128(a, b, &h);
    assert( l == a*b);
    assert( h == 54210108624275221ULL );
}
#endif

#ifndef R123_USE_MULHILO64_CUDA_INTRIN
#error "No  R123_USE_MULHILO64_CUDA_INTRIN"
#endif

#ifndef R123_USE_MULHILO64_OPENCL_INTRIN
#error "No  R123_USE_MULHILO64_OPENCL_INTRIN"
#endif

#ifndef R123_64BIT
#error "No R123_64BIT"
#else
void xx() {
    uint64_t a = R123_64BIT(0x1234567890abcdef);
    assert ( (a >> 60) == 0x1 );
}
#endif

#ifndef R123_USE_PHILOX_64BIT
#error "No  R123_USE_PHILOX_64BIT"
#endif

#ifndef R123_ASSERT
#error "No  R123_ASSERT"
#else
void chkassert(){
    R123_ASSERT(1);
}
#endif

#ifndef R123_STATIC_ASSERT
#error "No  R123_STATIC_ASSERT"
#else
R123_STATIC_ASSERT(1, "looks true to me");
void chkstaticassert(){
    R123_STATIC_ASSERT(1, "it's ok inside a function too");
}
#endif

#ifndef R123_USE_U01_DOUBLE
#error "No R123_USE_U01_DOUBLE"
#else
#if R123_USE_U01_DOUBLE
double foodbl = (1./(4294967296.*4294967296.));
#endif
#endif

int main(int , char **){return 0;}
