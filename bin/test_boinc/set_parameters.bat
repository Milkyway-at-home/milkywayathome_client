del boinc_finish_called
del init_data.xml
del boinc_lockfile
del astronomy_checkpoint
del out
del stderr.txt
echo "" > out
echo "" > stderr.txt

copy stars-%1.txt stars.txt
copy astronomy_parameters-%1.txt astronomy_parameters.txt
copy search_parameters-%1.txt search_parameters.txt