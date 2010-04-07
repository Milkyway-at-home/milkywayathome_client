rm -f boinc_finish_called
rm -f init_data.xml
rm -f boinc_lockfile
rm -f astronomy_checkpoint
rm -f out
rm -f stderr.txt
touch out
touch stderr.txt

cp stars-86.txt stars.txt
cp astronomy_parameters-86-small.txt astronomy_parameters.txt

./$1 -np 8 -p 0.73317163557524425 14.657212876628332 -1.7054653473950408 16.911711745343633 28.077212666463502 -1.2032908515814611 3.5273606439247287 2.2248214505875008 
echo "finished test on small parameters"

cp astronomy_parameters-86.txt astronomy_parameters.txt
./$1 -np 8 -p 0.73317163557524425 14.657212876628332 -1.7054653473950408 16.911711745343633 28.077212666463502 -1.2032908515814611 3.5273606439247287 2.2248214505875008 
