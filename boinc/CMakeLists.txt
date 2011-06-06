
cmake_minimum_required(VERSION 2.6)
project(boinc)

if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_BUILD_TYPE MinSizeRel)
endif()

include(CheckIncludeFiles)
include(CheckLibraryExists)
include(CheckFunctionExists)
include(IgnoreWarningFlags)
include(CheckTypeSize)
include(CheckStructHasMember)

set(BOINC_MAJOR_VERSION 6)
set(BOINC_MINOR_VERSION 13)
set(BOINC_RELEASE 0)

set(BOINC_PRERELEASE TRUE)
set(BOINC_BUGREPORT "")

# Check headers
check_include_files(alloca.h HAVE_ALLOCA_H)
check_include_files(arpa/inet.h HAVE_ARPA_INET_H)
check_include_files(dirent.h HAVE_DIRENT_H)
check_include_files(dlfcn.h HAVE_DLFCN_H)
check_include_files(errno.h HAVE_ERRNO_H)
check_include_files(fcntl.h HAVE_FCNTL_H)

check_include_files(GLUT/glut.h HAVE_GLUT_GLUT_H)
check_include_files(glut.h HAVE_GLUT_H)
check_include_files(glu.h HAVE_GLU_H)
check_include_files(GL/glaux.h HAVE_GLAUX_H)
check_include_files(GL/glut.h HAVE_GL_GLUT_H)
check_include_files(GL/glu.h HAVE_GL_GLU_H)
check_include_files(GL/gl.h HAVE_GL_GL_H)
check_include_files(gl.h HAVE_GL_H)
check_include_files(OpenGL/glaux.h HAVE_OPENGL_GLAUX_H)
check_include_files(OpenGL/glut.h HAVE_OPENGL_GLUT_H)
check_include_files(OpenGL/glu.h HAVE_OPENGL_GLU_H)
check_include_files(OpenGL/gl.h HAVE_OPENGL_GL_H)

check_include_files(ieeefp.h HAVE_IEEEFP_H)
check_include_files(inttypes.h HAVE_INTTYPES_H)
check_include_files(malloc.h HAVE_MALLOC_H)
check_include_files(netdb.h HAVE_NETDB_H)
check_include_files(netinet/ether.h HAVE_NETINET_ETHER_H)
check_include_files(netinet/if_ether.h HAVE_NETINET_IF_ETHER_H)
check_include_files(netinet/in.h HAVE_NETINET_IN_H)
check_include_files(netinet/tcp.h HAVE_NETINET_TCP_H)
check_include_files(net/if.h HAVE_NET_IF_H)
check_include_files(stdint.h HAVE_STDINT_H)
check_include_files(stdlib.h HAVE_STDLIB_H)
check_include_files(procfs.h HAVE_PROCFS_H)
check_include_files(signal.h HAVE_SIGNAL_H)

check_include_files(sys/fcntl.h HAVE_SYS_FCNTL_H)
check_include_files(sys/file.h HAVE_SYS_FILE_H)
check_include_files(sys/ioctl.h HAVE_SYS_IOCTL_H)
check_include_files(sys/ipc.h HAVE_SYS_IPC_H)
check_include_files(sys/mount.h HAVE_SYS_MOUNT_H)
check_include_files(sys/msg.h HAVE_SYS_MSG_H)
check_include_files(sys/param.h HAVE_SYS_PARAM_H)
check_include_files(sys/resource.h HAVE_SYS_RESOURCE_H)
check_include_files(sys/shm.h HAVE_SYS_SHM_H)
check_include_files(sys/socket.h HAVE_SYS_SOCKET_H)
check_include_files(sys/sockio.h HAVE_SYS_SOCKIO_H)
check_include_files(sys/statfs.h HAVE_SYS_STATFS_H)
check_include_files(sys/statvfs.h HAVE_SYS_STATVFS_H)
check_include_files(sys/stat.h HAVE_SYS_STAT_H)
check_include_files(sys/swap.h HAVE_SYS_SWAP_H)
check_include_files(sys/sysctl.h HAVE_SYS_SYSCTL_H)
check_include_files(sys/systeminfo.h HAVE_SYS_SYSTEMINFO_H)
check_include_files(sys/time.h HAVE_SYS_TIME_H)
check_include_files(sys/types.h HAVE_SYS_TYPES_H)
check_include_files(sys/un.h HAVE_SYS_UN_H)
check_include_files(sys/utsname.h HAVE_SYS_UTSNAME_H)
check_include_files(sys/vmmeter.h HAVE_SYS_VMMETER_H)
check_include_files(sys/wait.h HAVE_SYS_WAIT_H)

check_include_files(unistd.h HAVE_UNISTD_H)
check_include_files(utmp.h HAVE_UTMP_H)
check_include_files(windows.h HAVE_WINDOWS_H)

# Check functions
check_function_exists(alloca HAVE_ALLOCA)
check_function_exists(daemon HAVE_DAEMON)
check_function_exists(ether_ntoa HAVE_ETHER_NTOA)
check_function_exists(strcasecmp HAVE_STRCASECMP)
check_function_exists(strcasestr HAVE_STRCASESTR)
check_function_exists(strdup HAVE_STRDUP)
check_function_exists(strdupa HAVE_STRDUPA)
check_function_exists(strlcat HAVE_STRLCAT)
check_function_exists(strlcpy HAVE_STRLCPY)
check_include_files(strings.h HAVE_STRINGS_H)
check_include_files(string.h HAVE_STRING_H)
check_function_exists(res_init HAVE_RES_INIT)
check_function_exists(sched_setscheduler HAVE_SCHED_SETSCHEDULER)
check_function_exists(setenv HAVE_SETENV)
check_function_exists(setpriority HAVE_SETPRIORITY)
check_function_exists(setutent HAVE_SETUTENT)
check_function_exists(getutent HAVE_GETUTENT)
check_function_exists(sigaction HAVE_SIGACTION)

# Checking for types
if(APPLE)
  set(CMAKE_EXTRA_INCLUDE_FILES sys/socket.h)
  check_type_size(socklen_t SOCKLEN_T)
  set(CMAKE_EXTRA_INCLUDE_FILES)
else()
  check_type_size(socklen_t SOCKLEN_T)
endif()

check_type_size("struct ether_addr" STRUCT_ETHER_ADDR)
check_type_size("struct ifconf" STRUCT_IFCONF)
check_type_size("struct ifreq" STRUCT_IFREQ)
check_type_size("struct lifconf" STRUCT_LIFCONF)
check_struct_has_member("struct tm" "tm_zone" "time.h" HAVE_STRUCT_TM_TM_ZONE)

# Check libraries
find_package(Threads)
set(HAVE_PTHREAD ${CMAKE_USE_PTHREADS_INIT})
check_library_exists("xss" "XScreenSaverRegister" "/usr/lib" HAVE_XSS)


# CHECKME: Is there a better way to check this?
if(EXISTS "/proc/self/psinfo")
  set(HAVE__PROC_SELF_PSINFO TRUE)
endif()

if(EXISTS "/proc/self/stat")
  set(HAVE__PROC_SELF_STAT TRUE)
endif()


################################################################################

set(boinc_api_src boinc/api/boinc_api.cpp)
set(boinc_api_headers boinc/api/boinc_api.h)


set(boinc_lib_headers boinc/version.h
                      boinc/lib/app_ipc.h
                      boinc/lib/boinc_win.h
	                  boinc/lib/url.h
	                  boinc/lib/common_defs.h
	                  boinc/lib/diagnostics.h
	                  boinc/lib/diagnostics_win.h
	                  boinc/lib/filesys.h
	                  boinc/lib/hostinfo.h
	                  boinc/lib/proxy_info.h
	                  boinc/lib/prefs.h
	                  boinc/lib/miofile.h
	                  boinc/lib/mfile.h
	                  boinc/lib/parse.h
	                  boinc/lib/util.h
	                  boinc/lib/coproc.h
	                  boinc/lib/cal_boinc.h
	                  boinc/lib/cl_boinc.h
	                  boinc/svn_version.h
	                  boinc/win_build/config.h
	                  boinc/lib/str_util.h)

set(boinc_lib_src boinc/lib/util.cpp
                  boinc/lib/app_ipc.cpp
                  boinc/lib/diagnostics.cpp
                  boinc/lib/filesys.cpp
	              boinc/lib/hostinfo.cpp
                  boinc/lib/md5.c
                  boinc/lib/md5_file.cpp
                  boinc/lib/mem_usage.cpp
                  boinc/lib/mfile.cpp
                  boinc/lib/miofile.cpp
                  boinc/lib/parse.cpp
                  boinc/lib/prefs.cpp
	              boinc/lib/proxy_info.cpp
                  boinc/lib/str_util.cpp
                  boinc/lib/shmem.cpp
                  boinc/lib/base64.cpp
                  boinc/lib/url.cpp
                  boinc/lib/coproc.cpp)

if(WIN32)
  list(APPEND boinc_lib_src boinc/lib/win_util.cpp
                            boinc/lib/diagnostics_win.cpp
                            boinc/lib/stackwalker_win.cpp)
endif()


set(boinc_graphics_src boinc/api/graphics2.cpp boinc/api/graphics2_util.cpp)
set(boinc_graphics_headers boinc/api/graphics2.h)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/boinc/config.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/version.h.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/boinc/version.h)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/boinc)

if(APPLE)
  set(MAC_LIB_SRC_DIR boinc/lib/mac)
  include_directories(${MAC_LIB_SRC_DIR})
  list(APPEND boinc_lib_src ${MAC_LIB_SRC_DIR}/mac_backtrace.cpp
                            ${MAC_LIB_SRC_DIR}/QBacktrace.c
                            ${MAC_LIB_SRC_DIR}/QCrashReport.c
                            ${MAC_LIB_SRC_DIR}/QMachOImage.c
                            ${MAC_LIB_SRC_DIR}/QMachOImageList.c
                            ${MAC_LIB_SRC_DIR}/QSymbols.c
                            ${MAC_LIB_SRC_DIR}/QTaskMemory.c)

  list(APPEND boinc_headers ${MAC_LIB_SRC_DIR}/dyld_gdb.h
                            ${MAC_LIB_SRC_DIR}/mac_backtrace.h
                            ${MAC_LIB_SRC_DIR}/QBacktrace.h
                            ${MAC_LIB_SRC_DIR}/QCrashReport.h
                            ${MAC_LIB_SRC_DIR}/QMachOImage.h
                            ${MAC_LIB_SRC_DIR}/QMachOImageList.h
                            ${MAC_LIB_SRC_DIR}/QSymbols.h
                            ${MAC_LIB_SRC_DIR}/QTaskMemory.h)
endif()

add_library(boinc STATIC ${boinc_lib_src} ${boinc_lib_headers})
add_library(boinc_api STATIC ${boinc_api_src} ${boinc_api_headers})
add_library(boinc_graphics2 STATIC ${boinc_graphics_src} ${boinc_graphics_headers})

set(BOINC_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/boinc/lib" CACHE INTERNAL "BOINC headers")
set(BOINC_LIBRARIES boinc boinc_api boinc_graphics2 CACHE INTERNAL "BOINC libraries")

add_custom_command(OUTPUT ${PROJECT_SOURCE_DIR}/svn_version.h
                   COMMAND ${PROJECT_SOURCE_DIR}/boinc/generate_svn_version.sh
                   COMMENT "Creates/updates header SVN Revision Header")

install(TARGETS boinc boinc_api boinc_graphics2
        ARCHIVE       DESTINATION lib
        PUBLIC_HEADER DESTINATION include/boinc)

#export(TARGETS boinc boinc_api boinc_graphics2 FILE boinc-exports.cmake)


