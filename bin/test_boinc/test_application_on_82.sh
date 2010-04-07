rm -f boinc_finish_called
rm -f init_data.xml
rm -f boinc_lockfile
rm -f astronomy_checkpoint
rm -f out
rm -f stderr.txt
touch out
touch stderr.txt

cp stars-82.txt stars.txt
cp astronomy_parameters-82-small.txt astronomy_parameters.txt

./$1 -np 8 -p 0.40587961154742185 17.529961843393409 -1.8575145272144837 29.360893891378243 31.228263575178566 -1.551741065334 .064096152599308373 2.55428209912781
echo "finished test on small parameters"

cp astronomy_parameters-82.txt astronomy_parameters.txt
./$1 -np 8 -p 0.40587961154742185 17.529961843393409 -1.8575145272144837 29.360893891378243 31.228263575178566 -1.551741065334 .064096152599308373 2.55428209912781
