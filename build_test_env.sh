#!/bin/bash
# PUT THIS SCRIPT OUTSIDE THE milkywayathome_client directory
# After this, to run the whole testing suite, go to the newly created test_env folder and type 'make test' to run a specific test say 'ctest -R __test_name_here__' with the -VV option on if you want verbose output (otherwise you get no output)

# Path to milkywayathome_client directory 
cd "$(dirname "$0")"
PathToMilkyWayAtHomeClientDirectory="$(pwd)"
# Note, you need to specify if you plan to run the GPU testing suite. 'true' for true, 'false' for false 
includeGPUtesting=false

# In that case, you would need to first run this with includeGPUtesting to false to test everything except the GPU
# Then, set to true to run the rest of test with "make test" when in the test_env directory and stop after they have completed. Note, the GPU tests are very quick when compared to the other tests, so this should take around 30 minutes - 1 hour at MOST even on older GPUs.

if $includeGPUtesting; then
  rm $PathToMilkyWayAtHomeClientDirectory/nbody/tests/CMakeLists.txt
  cp $PathToMilkyWayAtHomeClientDirectory/nbody/tests/GPUTesting/makeGPUfile.txt $PathToMilkyWayAtHomeClientDirectory/nbody/tests/CMakeLists.txt
  
  cp $PathToMilkyWayAtHomeClientDirectory/nbody/tests/GPUTesting/gpu_test.c $PathToMilkyWayAtHomeClientDirectory/nbody/tests/gpu_test.c
  cp $PathToMilkyWayAtHomeClientDirectory/nbody/tests/GPUTesting/gpuTest.lua $PathToMilkyWayAtHomeClientDirectory/nbody/tests/gpuTest.lua
  cp $PathToMilkyWayAtHomeClientDirectory/nbody/tests/GPUTesting/gpu_test_basicmodels.c $PathToMilkyWayAtHomeClientDirectory/nbody/tests/gpu_test_basicmodels.c
  cp $PathToMilkyWayAtHomeClientDirectory/nbody/tests/GPUTesting/gpu_test_checkpoint.lua $PathToMilkyWayAtHomeClientDirectory/nbody/tests/gpu_test_checkpoint.lua
  cp $PathToMilkyWayAtHomeClientDirectory/nbody/tests/GPUTesting/gpu_test_morebasicmodels.c $PathToMilkyWayAtHomeClientDirectory/nbody/tests/gpu_test_morebasicmodels.c 
  cp $PathToMilkyWayAtHomeClientDirectory/nbody/tests/GPUTesting/gpu_test_moremodels.c $PathToMilkyWayAtHomeClientDirectory/nbody/tests/gpu_test_moremodels.c
else
  rm $PathToMilkyWayAtHomeClientDirectory/nbody/tests/CMakeLists.txt
  cp $PathToMilkyWayAtHomeClientDirectory/nbody/tests/GPUTesting/makeCPUfile.txt $PathToMilkyWayAtHomeClientDirectory/nbody/tests/CMakeLists.txt
  
  rm $PathToMilkyWayAtHomeClientDirectory/nbody/tests/gpu_test.c
  rm $PathToMilkyWayAtHomeClientDirectory/nbody/tests/gpuTest.lua
  rm $PathToMilkyWayAtHomeClientDirectory/nbody/tests/gpu_test_basicmodels.c
  rm $PathToMilkyWayAtHomeClientDirectory/nbody/tests/gpu_test_checkpoint.lua
  rm $PathToMilkyWayAtHomeClientDirectory/nbody/tests/gpu_test_morebasicmodels.c
  rm $PathToMilkyWayAtHomeClientDirectory/nbody/tests/gpu_test_moremodels.c
fi

  cd $PathToMilkyWayAtHomeClientDirectory/
  rm -rf test_env
  mkdir test_env
  cd test_env

if $includeGPUtesting; then 
  cmake -DNBODY_STATIC=OFF -DDOUBLEPREC=ON -DSEPARATION=OFF -DNBODY_GL=OFF -DBOINC_APPLICATION=OFF -DNBODY_OPENCL=ON $PathToMilkyWayAtHomeClientDirectory

  make gpu_sanity 
  make gpu_checkpoint
  make gpu_basicModelsPart1
  make gpu_basicModelsPart2
  make gpu_advanceModels
fi

  cmake -DNBODY_DEV_OPTIONS=ON -DNBODY_STATIC=OFF -DDOUBLEPREC=ON -DSEPARATION=OFF -DNBODY_GL=OFF -DBOINC_APPLICATION=OFF -DNBODY_OPENMP=ON -DNBODY_OPENCL=OFF $PathToMilkyWayAtHomeClientDirectory

  make blind_test
  make sanity
  make test_primitives
  make test_data
  make test_barriers
  make test_queue
  make emd_test
  make bessel_test
  make poisson_test
  make virial_test
  make mixeddwarf_test
  make nbody_test_driver
  make average_bins_test
  make stability_test
  
  make all

