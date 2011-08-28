#!/bin/sh

#
# Take AMD generated binary, strip extra sections and replace IL
# section with custom IL
#


if [ $# -ne 3 ]
then
    echo "Usage: `basename $0` <input binary> <IL file> <output binary>"
    exit 1
fi

input=${1}  # AMD generated ELF file
kernel=${2} # text file of new IL
output=${3} # New binary file to write


# How to find this:
# readelf -aW kernel.bin | grep Machine
#
machineCode=0x3f8


# Remove extra sections (.source and .llvmir). Also remove .amdil since it seems to reject --add-section with it already existing. objcopy also seems to reject it ("Unable to recognise the format") without specifying the machine stuff
objcopy -I elf32-i386 -O elf32-i386 --alt-machine-code=${machineCode} -R ".amdil" -R ".source" -R ".llvmir" --add-section ".amdil"=${kernel} ${input} ${output}



