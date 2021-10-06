rm -r test_env
mkdir test_env
cd test_env

cmake -DNBODY_STATIC=OFF -DDOUBLEPREC=ON -DSEPARATION=OFF -DNBODY_GL=OFF -DBOINC_APPLICATION=OFF -DNBODY_OPENCL=ON /home/dylansheils/Desktop/MilkyWayGithub/milkywayathome_client

make gpu_sanity 
make gpu_checkpoint
make gpu_basicModelsPart1
make gpu_basicModelsPart2
make gpu_advanceModels

cmake -DNBODY_STATIC=OFF -DDOUBLEPREC=ON -DSEPARATION=OFF -DNBODY_GL=OFF -DBOINC_APPLICATION=OFF -DNBODY_OPENMP=ON -DNBODY_OPENCL=OFF /home/dylansheils/Desktop/MilkyWayGithub/milkywayathome_client

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
