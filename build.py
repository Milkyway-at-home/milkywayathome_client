import sys
import subprocess

# Simple script to start builds using jenkins CI
# No dependecy checking, etc

# Usage
def usage():
    print("Usage: python build.py (win,linux,mac) (32|64)")

# Simple wrapper function
def execute(args):
    proc = subprocess.Popen(args, shell=False)
    proc.communicate()
    if proc.returncode != 0:
        print("ERROR: \"" + " ".join(args) + "\" exited with a non-zero return code")
        exit(1)

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
    if arch == "64":
        execute(["cmake", ".", "-DCMAKE_FIND_ROOT_PATH=/srv/chroot/hardy_amd64", "-DCMAKE_C_FLAGS=-m64", "-DCMAKE_CXX_FLAGS=-m64"])
    if arch == "32":
        execute(["cmake", ".", "-DBUILD_32=ON", "-DCMAKE_FIND_ROOT_PATH=/srv/chroot/hardy_i386", "-DCMAKE_C_FLAGS=-m32", "-DCMAKE_CXX_FLAGS=-m32", "-DCMAKE_AR=/usr/bin/ar"])
        
    execute(["make"])

if os == "win":
    if arch == "64":
        execute(["cmake", ".", "-G", "MinGW Makefiles", "-DSEPARATION=OFF"])
    if arch == "32":
        execute(["cmake", ".", "-G", "MinGW Makefiles", "-DBUILD_32=ON","-DSEPARATION=OFF"])
    execute(["mingw32-make"])
