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

# CMake flags used for all platforms
cmake_shared_flags = ["-DBOINC_RELEASE_NAMES=ON", "-DCMAKE_FIND_ROOT_PATH_MODE_PROGRAM=NEVER"]

# CMake flags used for windows
cmake_windows_flags = ["-DNBODY_STATIC=ON"]

# Linux
if os == "linux":

    if arch == "64":
        execute(["cmake", ".", "-DCMAKE_FIND_ROOT_PATH=/srv/chroot/hardy_amd64"] + cmake_shared_flags)

    if arch == "32":
        execute(["cmake", ".", "-DBUILD_32=ON", "-DCMAKE_FIND_ROOT_PATH=/srv/chroot/hardy_i386"] + cmake_shared_flags)
    
    execute(["make", "clean"])
    execute(["make"])

# Windows
if os == "win":

    if arch == "64":
        execute(["cmake", ".", "-G", "MinGW Makefiles", "-DSEPARATION=OFF"] + cmake_shared_flags + cmake_windows_flags)

    if arch == "32":
        execute(["cmake", ".", "-G", "MinGW Makefiles", "-DBUILD_32=ON","-DSEPARATION=OFF"] + cmake_shared_flags + cmake_windows_flags)
    
    execute(["mingw32-make", "clean"])
    execute(["mingw32-make"])
