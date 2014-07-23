import sys
import subprocess

# Simple script to start builds using jenkins CI
# No dependecy checking, etc

# Usage
def usage():
    print "Usage: python build.py (win,linux,mac) (32|64)"

# Simple wrapper function
def execute(args):
    proc = subprocess.Popen(args, shell=False)
    proc.communicate()

# Check correct number of args
if len(sys.argv) != 3:
    usage()
    exit(1)

# Make sure os and architecture are valid
os  = sys.argv[1]
arch = sys.argv[2]
assert os in ["win", "linux", "mac"], "ERROR: Unknown OS " + os
assert arch in ["32", "64"], "ERROR: Unknown arch " + bit 

if os == "linux":
    execute(["cmake", "."])
    execute(["make"])

if os == "win":
    execute(["cmake", ".", "-G", "MinGW Makefiles"])
    execute(["mingw32-make"])
