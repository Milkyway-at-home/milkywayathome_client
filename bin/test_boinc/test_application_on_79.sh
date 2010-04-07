rm -f boinc_finish_called
rm -f init_data.xml
rm -f boinc_lockfile
rm -f astronomy_checkpoint
rm -f out
rm -f stderr.txt
touch out
touch stderr.txt

cp stars-79.txt stars.txt
cp astronomy_parameters-79-small.txt astronomy_parameters.txt

./$1 -np 8 -p 0.34217373320392042 25.9517910846623 -2.1709414738826602 38.272511356953906 30.225190442596112 2.2149060013372885 0.32316169064291655 2.7740244716285285 
echo "finished test on small parameters"

cp astronomy_parameters-79.txt astronomy_parameters.txt
./$1 -np 8 -p 0.34217373320392042 25.9517910846623 -2.1709414738826602 38.272511356953906 30.225190442596112 2.2149060013372885 0.32316169064291655 2.7740244716285285 

