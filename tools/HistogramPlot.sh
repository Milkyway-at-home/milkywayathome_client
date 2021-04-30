#!/bin/sh
#
# Take an arbitrary number of histograms produced from nbody runs
# and plot them on a single histogram
#

gnuplotPath=gnuplot

if [ $# -eq 0 ]; then
    echo "Need argument"
    exit 1
fi

# 0 = number 1 .. n
# 1 = filename
titleType=1

cmd="set xlabel \"λ\"; "
cmd="${cmd} plot "

n=0
for i in "$@"
do
    n=$((${n} + 1))

    if [ ${titleType} -eq 1 ]; then
        cmd="${cmd} \"${i}\"  using 2:3:4 title \"${i}\" with boxerrorbars"
    else
        cmd="${cmd} \"${i}\"  using 2:3:4 title \"${n}\" with boxerrorbars"
    fi

    if [ ${n} -ne $# ]; then
        cmd="${cmd}, "
    fi
done

echo ${cmd}
${gnuplotPath} -p -e "${cmd}"




