echo "creating files for milkyway app release "$1

mkdir milkyway
cd milkyway
mkdir license
mkdir astronomy
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

cp ../../util/io_util.c ./util/
cp ../../util/io_util.h ./util/
cp ../../util/settings.c ./util/
cp ../../util/settings.h ./util/

cp ../make.osx ./bin/
cp ../make.linux ./bin/
cp ../test_boinc/set_parameters.sh ./bin/test_files/test_application.sh
cp ../test_boinc/astronomy_parameters-11.txt ./bin/test_files/astronomy_parameters-11.txt
cp ../test_boinc/astronomy_parameters-12.txt ./bin/test_files/astronomy_parameters-12.txt
cp ../test_boinc/astronomy_parameters-20.txt ./bin/test_files/astronomy_parameters-20.txt
cp ../test_boinc/astronomy_parameters-21.txt ./bin/test_files/astronomy_parameters-21.txt
cp ../test_boinc/astronomy_parameters-79.txt ./bin/test_files/astronomy_parameters-79.txt
cp ../test_boinc/astronomy_parameters-82.txt ./bin/test_files/astronomy_parameters-82.txt
cp ../test_boinc/astronomy_parameters-86.txt ./bin/test_files/astronomy_parameters-86.txt
cp ../test_boinc/stars-11.txt ./bin/test_files/stars-11.txt
cp ../test_boinc/stars-12.txt ./bin/test_files/stars-12.txt
cp ../test_boinc/stars-20.txt ./bin/test_files/stars-20.txt
cp ../test_boinc/stars-21.txt ./bin/test_files/stars-21.txt
cp ../test_boinc/stars-79.txt ./bin/test_files/stars-79.txt
cp ../test_boinc/stars-82.txt ./bin/test_files/stars-82.txt
cp ../test_boinc/stars-86.txt ./bin/test_files/stars-86.txt

cd ..
tar cvzf mw_v$1.tar milkyway
zip -r mw_v$1.zip milkyway
rm -r milkyway
