#!/bin/bash
# PUT THIS SCRIPT OUTSIDE THE milkywayathome_client directory
# After this, to run the whole testing suite, go to the newly created test_env folder and type 'make test' to run a specific test say 'ctest -R __test_name_here__' with the -VV option on if you want verbose output (otherwise you get no output)

# Specify YOUR OWN path here, in the folder which contains the milkywayathome_client directory
PathToMilkyWayAtHomeClientFolder='/home/dylansheils/Desktop/pushToLMCBranch/'
# Note, you need to specify if you plan to run the GPU testing suite. 'true' for true, 'false' for false 
includeGPUtesting=false

# In that case, you would need to first run this with includeGPUtesting to false to test everything except the GPU
# Then, set to true to run the rest of test with "make test" when in the test_env directory and stop after they have completed. Note, the GPU tests are very quick when compared to the other tests, so this should take around 30 minutes - 1 hour at MOST even on older GPUs.

if $includeGPUtesting; then
  rm $PathToMilkyWayAtHomeClientFolder/milkywayathome_client/nbody/tests/CMakeLists.txt
  cp $PathToMilkyWayAtHomeClientFolder/milkywayathome_client/nbody/tests/GPUTesting/makeGPUfile.txt $PathToMilkyWayAtHomeClientFolder/milkywayathome_client/nbody/tests/CMakeLists.txt
  
  cp $PathToMilkyWayAtHomeClientFolder/milkywayathome_client/nbody/tests/GPUTesting/gpu_test.c $PathToMilkyWayAtHomeClientFolder/milkywayathome_client/nbody/tests/gpu_test.c
  cp $PathToMilkyWayAtHomeClientFolder/milkywayathome_client/nbody/tests/GPUTesting/gpuTest.lua $PathToMilkyWayAtHomeClientFolder/milkywayathome_client/nbody/tests/gpuTest.lua
  cp $PathToMilkyWayAtHomeClientFolder/milkywayathome_client/nbody/tests/GPUTesting/gpu_test_basicmodels.c $PathToMilkyWayAtHomeClientFolder/milkywayathome_client/nbody/tests/gpu_test_basicmodels.c
  cp $PathToMilkyWayAtHomeClientFolder/milkywayathome_client/nbody/tests/GPUTesting/gpu_test_checkpoint.lua $PathToMilkyWayAtHomeClientFolder/milkywayathome_client/nbody/tests/gpu_test_checkpoint.lua
  cp $PathToMilkyWayAtHomeClientFolder/milkywayathome_client/nbody/tests/GPUTesting/gpu_test_morebasicmodels.c $PathToMilkyWayAtHomeClientFolder/milkywayathome_client/nbody/tests/gpu_test_morebasicmodels.c 
  cp $PathToMilkyWayAtHomeClientFolder/milkywayathome_client/nbody/tests/GPUTesting/gpu_test_moremodels.c $PathToMilkyWayAtHomeClientFolder/milkywayathome_client/nbody/tests/gpu_test_moremodels.c
  
  cd $PathToMilkyWayAtHomeClientFolder/
  rm -r test_env
  mkdir test_env
  cd test_env
  
  cmake -DNBODY_STATIC=OFF -DDOUBLEPREC=ON -DSEPARATION=OFF -DNBODY_GL=OFF -DBOINC_APPLICATION=OFF -DNBODY_OPENCL=ON $PathToMilkyWayAtHomeClientFolder/milkywayathome_client

  make gpu_sanity 
  make gpu_checkpoint
  make gpu_basicModelsPart1
  make gpu_basicModelsPart2
  make gpu_advanceModels

  cmake -DNBODY_STATIC=OFF -DDOUBLEPREC=ON -DSEPARATION=OFF -DNBODY_GL=OFF -DBOINC_APPLICATION=OFF -DNBODY_OPENMP=ON -DNBODY_OPENCL=OFF $PathToMilkyWayAtHomeClientFolder/milkywayathome_client

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
  make all
else
  rm $PathToMilkyWayAtHomeClientFolder/milkywayathome_client/nbody/tests/CMakeLists.txt
  cp $PathToMilkyWayAtHomeClientFolder/milkywayathome_client/nbody/tests/GPUTesting/makeCPUfile.txt $PathToMilkyWayAtHomeClientFolder/milkywayathome_client/nbody/tests/CMakeLists.txt
  
  rm $PathToMilkyWayAtHomeClientFolder/milkywayathome_client/nbody/tests/gpu_test.c
  rm $PathToMilkyWayAtHomeClientFolder/milkywayathome_client/nbody/tests/gpuTest.lua
  rm $PathToMilkyWayAtHomeClientFolder/milkywayathome_client/nbody/tests/gpu_test_basicmodels.c
  rm $PathToMilkyWayAtHomeClientFolder/milkywayathome_client/nbody/tests/gpu_test_checkpoint.lua
  rm $PathToMilkyWayAtHomeClientFolder/milkywayathome_client/nbody/tests/gpu_test_morebasicmodels.c
  rm $PathToMilkyWayAtHomeClientFolder/milkywayathome_client/nbody/tests/gpu_test_moremodels.c
  
  cd $PathToMilkyWayAtHomeClientFolder/
  rm -r test_env
  mkdir test_env
  cd test_env

  cmake -DNBODY_STATIC=OFF -DDOUBLEPREC=ON -DSEPARATION=OFF -DNBODY_GL=OFF -DBOINC_APPLICATION=OFF -DNBODY_OPENMP=ON -DNBODY_OPENCL=OFF $PathToMilkyWayAtHomeClientFolder/milkywayathome_client

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
  make all
fi
