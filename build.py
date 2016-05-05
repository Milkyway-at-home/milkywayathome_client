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
cmake_shared_flags = ["-DBOINC_RELEASE_NAMES=ON", "-DSEPARATION=OFF", "-DCMAKE_FIND_ROOT_PATH_MODE_PROGRAM=NEVER", "-DNBODY_OPENMP=" + nbody_openmp_sep_opencl, "-DSEPARATION_OPENCL=" + nbody_openmp_sep_opencl]

# CMake flags used for windows
cmake_static_flag = ["-DNBODY_STATIC=ON",  "-DBOINC_APPLICATION=ON",  "-DCMAKE_BUILD_TYPE=Release"]

# Linux
if os == "linux":

    if arch == "64":
        execute(["cmake", ".", "-DCMAKE_FIND_ROOT_PATH=/srv/chroot/hardy_amd64", "-DOPENCL_LIBRARIES=/srv/chroot/hardy_amd64/usr/lib/libOpenCL.so", "-DOPENCL_INCLUDE_DIRS=/srv/chroot/hardy_amd64/usr/local/cuda/include/"] + cmake_shared_flags + cmake_static_flag)

    if arch == "32":
        print("ERROR: Linux 32 bit not supported")
        #execute(["cmake", ".", "-DBUILD_32=ON", "-DCMAKE_FIND_ROOT_PATH=/srv/chroot/hardy_i386", "-DOPENCL_LIBRARIES=/srv/chroot/hardy_i386/usr/lib/libOpenCL.so", "-DOPENCL_INCLUDE_DIRS=/srv/chroot/hardy_i386/usr/local/cuda/include/"] + cmake_shared_flags + cmake_static_flag)
    
    execute(["make", "clean"])
    execute(["make"])

# Windows
if os == "win":

    if arch == "64":
        execute(["cmake", ".", "-G", "MinGW Makefiles"] + cmake_shared_flags + cmake_static_flag)

    if arch == "32":
        print("ERROR: Windows 32 bit not supported")
        #execute(["cmake", ".", "-G", "MinGW Makefiles", "-DBUILD_32=ON"] + cmake_shared_flags + cmake_static_flag)
    
    execute(["mingw32-make", "clean"])
    execute(["mingw32-make"])

# Mac OS X
if os == "mac":

    if arch == "64":
        execute(["/opt/local/bin/cmake", "."] + cmake_shared_flags + cmake_static_flag)

    if arch == "32":
        print("ERROR: Mac OSX 32 bit not supported")

    execute(["/usr/bin/make", "clean"])
    execute(["/usr/bin/make"])

