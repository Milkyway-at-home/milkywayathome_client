rm -f boinc_finish_called
rm -f init_data.xml
rm -f astronomy_checkpoint
rm -f out
rm -f stderr.txt
touch out

cp stars-$1.txt stars.txt
cp astronomy_parameters-$1.txt astronomy_parameters.txt
cp search_parameters-$1.txt search_parameters.txt
