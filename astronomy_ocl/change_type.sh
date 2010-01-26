#!/bin/bash
for f in *.cpp
do 
    mv $f $f.temp
    sed "s/$1/$2/g" $f.temp > $f
    rm -f $f.temp
done

for f in *.h
do 
    mv $f $f.temp
    sed "s/$1/$2/g" $f.temp > $f
    rm -f $f.temp
done

for f in *.cl
do 
    mv $f $f.temp
    sed "s/$1/$2/g" $f.temp > $f
    rm -f $f.temp
done