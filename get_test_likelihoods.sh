#Should be in same folder as get_test_likelihoods.py when run

cd test_env/

ctest -R $1__$2_test -VV > $3
