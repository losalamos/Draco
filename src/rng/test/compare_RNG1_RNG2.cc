// compare_RNG1_RNG2.cc
// T. M. Kelley
// Aug 15, 2018
// (c) Copyright 2018 LANSLLC, all rights reserved

#include "c4/Timer.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "rng/Counter_RNG.hh"
#include "rng/Counter_RNG2.hh"
#include <vector>

#define NUM_ELEMS 1000000

double out_1[NUM_ELEMS];
double out_2[NUM_ELEMS];

uint64_t const time_step = 18436;
uint64_t const seed = 0xFEED;
uint64_t const stream_num = 30001;

using namespace rtt_rng;

void run_rng_1(){
  uint64_t data[4];
  data[0] = 0;
  data[1] = stream_num;
  data[2] = seed;
  data[3] = time_step;
  Counter_RNG rng(data, data+4);
  rtt_c4::Timer timer;
  timer.start();
  for(size_t i = 0; i < NUM_ELEMS; ++i){
    out_1[i] = rng.ran();
  }
  timer.stop();
  std::cout << "RNG 1: Generating " << NUM_ELEMS << " randoms took "
    << timer.wall_clock() << " seconds.\n";
  return;
}

void run_rng_2(){
  Counter_RNG2 rng(seed, stream_num, time_step);
  rtt_c4::Timer timer;
  timer.start();
  for(size_t i = 0; i < NUM_ELEMS; ++i){
    out_2[i] = rng.ran();
  }
  timer.stop();
  std::cout << "RNG 2: Generating " << NUM_ELEMS << " randoms took "
    << timer.wall_clock() << " seconds.\n";
  return;
}

void compare_rngs(rtt_dsxx::UnitTest &ut){
  bool ok = true;
  size_t num_fails = 0;
  constexpr size_t max_fails = 10;
  size_t fails[max_fails];
  // zero fails array
  for(size_t i = 0; i < max_fails; ++i){
    fails[i] = 0;
  }
  for(size_t i = 0; i < NUM_ELEMS; ++i){
    bool this_ok = (out_1 == out_2);
    if(!this_ok){
      if(num_fails < max_fails){
        fails[num_fails] = i;
      }
      num_fails++;
    }
    this_ok = ok && this_ok;
  }
  if(!ok){
    std::cout << "Arrays were not the same, at least " << num_fails
      << " failed to matchl eading indices: \n";
    for(size_t i = 0; i < NUM_ELEMS; ++i){
      std::cout << fails[i] << "; ";
      std::cout << "\n";
    }

  }
  else{
    std::cout << "Arrays agreed\n";
  }
  if (ut.numFails == 0){
    PASSMSG("compare old & new RNGs passed");
  }

  return;
}

void init_arrays(){
  for(size_t i = 0; i < NUM_ELEMS; ++i){
    out_1[i] = 0.0;
    out_2[i] = 0.0;
  }
  return;
}

int main(int argc, char **argv){
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try{
    init_arrays();
    run_rng_1();
    run_rng_2();
    compare_rngs(ut);
  }
  UT_EPILOG(ut);
  return 0;
}

// End of file
