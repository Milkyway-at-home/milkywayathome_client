echo "creating files for milkyway3 app release "$1

mkdir milkyway
cd milkyway
mkdir license
mkdir astronomy
mkdir evaluation
mkdir searches
mkdir util
mkdir bin
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

cp ../../searches/search_parameters.c ./searches/
cp ../../searches/search_parameters.h ./searches/
cp ../../searches/result.c ./searches/
cp ../../searches/result.h ./searches/

cp ../../evaluation/evaluator.h ./evaluation/
cp ../../evaluation/simple_evaluator.c ./evaluation/
cp ../../evaluation/simple_evaluator.h ./evaluation/

cp ../../util/matrix.c ./util/
cp ../../util/matrix.h ./util/
cp ../../util/io_util.c ./util/
cp ../../util/io_util.h ./util/
cp ../../util/settings.c ./util/
cp ../../util/settings.h ./util/

cp ../Makefile ./bin/
cp ../test_boinc/astronomy_parameters-*.txt ./bin/test_files/
cp ../test_boinc/stars-*.txt ./bin/test_files/
cp ../test_boinc/*.sh ./bin/test_files/

cd ..
tar cvzf mw3_v$1.tar milkyway
zip -r mw3_v$1.zip milkyway
rm -r milkyway
