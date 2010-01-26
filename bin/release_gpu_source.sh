echo "creating files for milkyway_gpu app release "$1

mkdir milkyway
cd milkyway
mkdir license
mkdir astronomy
mkdir astronomy_gpu
mkdir searches
mkdir util
mkdir bin
mkdir evaluation
mkdir bin/boinc/
mkdir bin/test_files/

cp ../../license/GPL.txt ./license/
cp ../../astronomy/boinc_astronomy.C ./astronomy/
cp ../../astronomy/atSurveyGeometry.c ./astronomy/
cp ../../astronomy/atSurveyGeometry.h ./astronomy/
cp ../../astronomy/numericalIntegration.c ./astronomy/
cp ../../astronomy/numericalIntegration.h ./astronomy/
cp ../../astronomy/parameters.c ./astronomy/
cp ../../astronomy/parameters.h ./astronomy/
cp ../../astronomy/probability.c ./astronomy/
cp ../../astronomy/probability.h ./astronomy/
cp ../../astronomy/stCoords.c ./astronomy/
cp ../../astronomy/stCoords.h ./astronomy/
cp ../../astronomy/stCnum.c ./astronomy/
cp ../../astronomy/stCnum.h ./astronomy/
cp ../../astronomy/stMath.c ./astronomy/
cp ../../astronomy/stMath.h ./astronomy/
cp ../../astronomy/stVector.c ./astronomy/
cp ../../astronomy/stVector.h ./astronomy/
cp ../../astronomy/star_points.c ./astronomy/
cp ../../astronomy/star_points.h ./astronomy/
cp ../../astronomy/evaluation_optimized.c ./astronomy/
cp ../../astronomy/evaluation_optimized.h ./astronomy/
cp ../../astronomy/evaluation_state.c ./astronomy/
cp ../../astronomy/evaluation_state.h ./astronomy/

cp ../../astronomy_gpu/evaluation_gpu6.cu ./astronomy_gpu/
cp ../../astronomy_gpu/evaluation_gpu6_double.cu ./astronomy_gpu/
cp ../../astronomy_gpu/evaluation_gpu6_float.cu ./astronomy_gpu/
cp ../../astronomy_gpu/evaluation_gpu6_likelihood.cu ./astronomy_gpu/
cp ../../astronomy_gpu/evaluation_gpu6_macros.cu ./astronomy_gpu/
cp ../../astronomy_gpu/evaluation_gpu.h ./astronomy_gpu/
cp ../../astronomy_gpu/coords.h ./astronomy_gpu/
cp ../../astronomy_gpu/cpu_coords.h ./astronomy_gpu/
cp ../../astronomy_gpu/r_constants.h ./astronomy_gpu/
cp ../../astronomy_gpu/pi_constants.h ./astronomy_gpu/
cp ../../astronomy_gpu/gauss_legendre.h ./astronomy_gpu/
cp ../../astronomy_gpu/cpu_integrals.h ./astronomy_gpu/
cp ../../astronomy_gpu/coords.c ./astronomy_gpu/
cp ../../astronomy_gpu/gauss_legendre.c ./astronomy_gpu/
cp ../../astronomy_gpu/cpu_coords.c ./astronomy_gpu/
cp ../../astronomy_gpu/cpu_integrals.c ./astronomy_gpu/

cp ../../evaluation/simple_evaluator.c ./evaluation/
cp ../../evaluation/simple_evaluator.h ./evaluation/
cp ../../evaluation/evaluator.h ./evaluation/

cp ../../searches/hessian.c ./searches/
cp ../../searches/hessian.h ./searches/
cp ../../searches/gradient.c ./searches/
cp ../../searches/gradient.h ./searches/
cp ../../searches/line_search.c ./searches/
cp ../../searches/line_search.h ./searches/
cp ../../searches/newton_method.c ./searches/
cp ../../searches/newton_method.h ./searches/
cp ../../searches/regression.c ./searches/
cp ../../searches/regression.h ./searches/
cp ../../searches/result.c ./searches/
cp ../../searches/result.h ./searches/
cp ../../searches/search_parameters.c ./searches/
cp ../../searches/search_parameters.h ./searches/

cp ../../util/matrix.c ./util/
cp ../../util/matrix.h ./util/
cp ../../util/io_util.c ./util/
cp ../../util/io_util.h ./util/
cp ../../util/settings.c ./util/
cp ../../util/settings.h ./util/

cp ../Makefile ./bin/
cp ../milkyway.vcproj ./bin/
cp ../milkyway_vc90.sln ./bin/
cp ../Cuda.Rules ./bin/
cp ../Cuda2.3.Rules ./bin/
cp ../test_boinc/set_parameters.sh ./bin/test_files/set_parameters.sh
cp ../test_boinc/set_parameters.bat ./bin/test_files/set_parameters.bat
cp ../test_boinc/search_parameters-20.txt ./bin/test_files/search_parameters-20.txt
cp ../test_boinc/search_parameters-21.txt ./bin/test_files/search_parameters-21.txt
cp ../test_boinc/search_parameters-79.txt ./bin/test_files/search_parameters-79.txt
cp ../test_boinc/search_parameters-82.txt ./bin/test_files/search_parameters-82.txt
cp ../test_boinc/search_parameters-86.txt ./bin/test_files/search_parameters-86.txt
cp ../test_boinc/astronomy_parameters-20.txt ./bin/test_files/astronomy_parameters-20.txt
cp ../test_boinc/astronomy_parameters-21.txt ./bin/test_files/astronomy_parameters-21.txt
cp ../test_boinc/astronomy_parameters-79.txt ./bin/test_files/astronomy_parameters-79.txt
cp ../test_boinc/astronomy_parameters-82.txt ./bin/test_files/astronomy_parameters-82.txt
cp ../test_boinc/astronomy_parameters-86.txt ./bin/test_files/astronomy_parameters-86.txt
cp ../test_boinc/stars-20.txt ./bin/test_files/stars-20.txt
cp ../test_boinc/stars-21.txt ./bin/test_files/stars-21.txt
cp ../test_boinc/stars-79.txt ./bin/test_files/stars-79.txt
cp ../test_boinc/stars-82.txt ./bin/test_files/stars-82.txt
cp ../test_boinc/stars-86.txt ./bin/test_files/stars-86.txt

cp ../boinc/* ./bin/boinc/

cd ..
tar cvzf mw_gpu_v$1.tar.gz milkyway
zip -r mw_gpu_v$1.zip milkyway
#rm -r milkyway
