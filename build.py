import sys
import subprocess

# Simple script to start builds using jenkins CI
# No dependecy checking, etc

# Usage
def usage():
    print("Usage: python build.py (win,linux,mac) (32|64) (ON|OFF)")

# Simple wrapper function
def execute(args):
    proc = subprocess.Popen(args, shell=False)
    proc.communicate()
    if proc.returncode != 0:
        print("ERROR: \"" + " ".join(args) + "\" exited with a non-zero return code")
        exit(1)

# Check correct number of args
if len(sys.argv) != 4:
    usage()
    exit(1)

# Make sure os and architecture are valid
os  = sys.argv[1]
arch = sys.argv[2]
nbody_openmp_sep_opencl = sys.argv[3]
assert os in ["win", "linux", "mac"], "ERROR: Unknown OS " + os
assert arch in ["32", "64"], "ERROR: Unknown arch " + bit 
assert nbody_openmp_sep_opencl in ["ON", "OFF"], "ERROR: Set NBODY_OPENMP to ON or OFF"

# CMake flags used for all platforms
cmake_shared_flags = ["-DBOINC_RELEASE_NAMES=ON", "-DCMAKE_FIND_ROOT_PATH_MODE_PROGRAM=NEVER", "-DNBODY_OPENMP=" + nbody_openmp_sep_opencl, "-DSEPARATION_OPENCL=" + nbody_openmp_sep_opencl]

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
        execute(["cmake", ".", "-G", "MinGW Makefiles", "-DBUILD_32=ON","-DSEPARATION=ON"] + cmake_shared_flags + cmake_windows_flags)
    
    execute(["mingw32-make", "clean"])
    execute(["mingw32-make"])

# Mac OS X
if os == "mac":

    if arch == "64":
        execute(["/opt/local/bin/cmake", "."] + cmake_shared_flags)

    if arch == "32":
        print("ERROR: Mac OSX 32 bit not supported")

    execute(["/usr/bin/make", "clean"])
    execute(["/usr/bin/make"])

