rm -f boinc_finish_called
rm -f init_data.xml
rm -f boinc_lockfile
rm -f astronomy_checkpoint
rm -f out
rm -f stderr.txt
touch out
touch stderr.txt

cp stars-11.txt stars.txt
cp astronomy_parameters-11-small.txt astronomy_parameters.txt

./$1 -np 20 -p 0.571713 12.312119 -3.305187 148.010257 22.453902 0.42035 -0.468858 0.760579 -1.361644 177.884238 23.882892 1.210639 -1.611974 8.534378 -1.361644 177.884238 10.882892 1.210639 -1.611974 8.534378 

echo "finished test on small parameters"

cp astronomy_parameters-11.txt astronomy_parameters.txt
./$1 -np 20 -p 0.571713 12.312119 -3.305187 148.010257 22.453902 0.42035 -0.468858 0.760579 -1.361644 177.884238 23.882892 1.210639 -1.611974 8.534378 -1.361644 177.884238 10.882892 1.210639 -1.611974 8.534378 
