#!/usr/bin/env python2.7
#
#
# Script to run batches of simulation while varying a few parameters
# and keeping things somewhat organized
#
# Can run the same simulations on multiple devices with the --devices argument
#    --devices {device index 1} {device index 2}...
#    Use -1 for using the non-OpenCL version
#
# Can run the same simulation with multiple seeds with --seeds
#   --seeds {seed1} {seed2} ...
#   If this is omitted, a single time-seeded simulation will be run for each device
#
# --help should give more information on other arguments. Other
# required arguments are the input files and output directory. Output
# will be in subdirectories named for the type of device
#
# The simulation will run a simulation for each seed on each device.
# The tasks for different devices will run at the same time
#
#
#
# Requires PyOpenCL and decorator (apparently)
#
#

import os
import subprocess
import pyopencl as cl
import argparse
import time
from threading import Thread


nbodyBin = None


def maybeMkdir(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)

def maybeMakeAllDir(outputDir):
    maybeMkdir(os.path.join(outputDir, "histograms"))
    maybeMkdir(os.path.join(outputDir, "stderrs"))
    maybeMkdir(os.path.join(outputDir, "outputs"))
    maybeMkdir(os.path.join(outputDir, "checkpoints"))




def getArguments():
    parser = argparse.ArgumentParser(description='Run batch of nbody')
    parser.add_argument('--outdir', '-o',
                        type=str,
                        dest='baseOutputDirectory',
                        required=True,
                        help='base output directory')

    parser.add_argument('--input-file', '-f',
                        type=str,
                        dest='inputFile',
                        required=True,
                        help='input file for run')

    parser.add_argument('--histogram', '-i', # -h taken by the help
                        type=str,
                        dest='histogram',
                        required=True,
                        help='input histogram file for run')

    parser.add_argument('--binary', '-b',
                        type=str,
                        default='./milkyway_nbody',
                        help='path to milkyway_nbody')

    parser.add_argument('--platform', '-p',
                        type=int,
                        default=0,
                        help='platform id')

    parser.add_argument('--devices', '-d',
                        type=int,
                        nargs='+',
                        default=None,
                        help='device ids (< 0 for plain CPU)')

    parser.add_argument('--seeds', '-e',
                        type=int,
                        nargs='+',
                        default=None,
                        help='seeds')

    return parser.parse_known_args()


def realMain(outDir, inputFile, inputHistogram, seed=None, platform=0, device=None, extraArgs=[]):
    def getDeviceName(device):
        # No OpenCL if no device
        if device == None:
            return "cpu"
        else:
            return cl.get_platforms()[platform].get_devices()[device].name

    deviceName = getDeviceName(device)
    outputDir = os.path.join(outDir, deviceName)
    maybeMakeAllDir(outputDir)

    def makeRun(seed):
        runName = "run_" + str(seed)

        histOutFileName = os.path.join(outputDir, "histograms", runName)
        outFileName = os.path.join(outputDir, "outputs", runName)
        stderrFileName = os.path.join(outputDir, "stderrs", runName)
        checkpointFileName = os.path.join(outputDir, "checkpoints", runName)

        nbodyCmd = [nbodyBin,
                    "--output-cartesian",
                    "--timing",
                    "--no-clean-checkpoint",
                   # "--progress", # Bad things seem to happen with multiple at the same time

                    "--seed",           str(seed),

                    "--input-file",     inputFile,
                    "--histogram-file", inputHistogram,

                    "--checkpoint",     checkpointFileName,
                    "--output-file",    outFileName,
                    "--histoout-file",  histOutFileName
                    ]

        if platform:
            nbodyCmd.extend(["--platform", str(platform)])

        if device:
            nbodyCmd.extend(["--device",  str(device)])
        else:
            nbodyCmd.append("--disable-opencl")


        nbodyCmd.extend(extraArgs)

        fd = open(stderrFileName, "w+")
        proc = subprocess.Popen(nbodyCmd, stderr = fd)
        fd.close()
        proc.wait()
        print "Device '%s' finished run of seed %d\n" % (deviceName, seed)
        return proc.returncode


    return makeRun(seed)

def main():
    global nbodyBin

    args, extraArgs = getArguments()
    nbodyBin = args.binary

    if not args.seeds:
        seeds = [int(time.time())]
    else:
        seeds = args.seeds


    def runDevice(device):
        for seed in seeds:
            proc = realMain(outDir         = args.baseOutputDirectory,
                            inputFile      = args.inputFile,
                            inputHistogram = args.histogram,
                            seed           = seed,
                            device         = device,
                            platform       = args.platform,
                            extraArgs      = extraArgs)


    # Launch a thread for each device.
    for dev in args.devices:
        # Make -1 device map to using normal CPU one
        if dev < 0:
            dev = None
        t = Thread(target = runDevice, args = (dev, ))
        t.start()


if __name__ == '__main__':
    main()

