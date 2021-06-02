rm -r test_env
mkdir test_env
cd test_env

cmake -DNBODY_STATIC=OFF -DDOUBLEPREC=OFF -DSEPARATION=OFF -DNBODY_GL=OFF -DBOINC_APPLICATION=OFF -DCMAKE_BUILD_TYPE=Debug -DNBODY_OPENCL=ON /home/dylansheils/Desktop/Milkyway/milkywayathome_client/

make gpu_test

cmake -DNBODY_STATIC=OFF -DDOUBLEPREC=ON -DNBODY_DEV_OPTIONS=ON -DEBUG=ON -DSEPARATION=OFF -DNBODY_GL=OFF -DBOINC_APPLICATION=OFF -DCMAKE_BUILD_TYPE=Debug -DNBODY_OPENCL=OFF /home/dylansheils/Desktop/Milkyway/milkywayathome_client/

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
