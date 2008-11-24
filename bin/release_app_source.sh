echo "creating files for milkyway app relase "$1

mkdir milkyway
cd milkyway
mkdir astronomy
mkdir searches
mkdir util
mkdir bin
mkdir bin/test_files/

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
cp ../stars.txt ./bin/test_files/stars.txt
cp ../search_parameters.txt ./bin/test_files/search_parameters.txt
cp ../parameters.txt ./bin/test_files/astronomy_parameters.txt
cp ../parameters-cut-large.txt ./bin/test_files/astronomy_parameters-cut-large.txt
cp ../parameters-cut-medium.txt ./bin/test_files/astronomy_parameters-cut-medium.txt
cp ../parameters-nocut-small.txt ./bin/test_files/astronomy_parameters-nocut-small.txt
cp ../parameters-unconvolved-small.txt ./bin/test_files/astronomy_parameters-unconvolved-small.txt

cd ..
tar cvzf milkyway_release_$1.tar milkyway
zip -r milkyway_release_$1.zip milkyway
rm -r milkyway
