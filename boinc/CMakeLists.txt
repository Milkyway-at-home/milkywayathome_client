
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
include(CheckSSE2)

set(BOINC_MAJOR_VERSION 6)
set(BOINC_MINOR_VERSION 13)
set(BOINC_RELEASE 0)

set(BOINC_PRERELEASE TRUE)
set(BOINC_BUGREPORT "")

find_package(OpenSSL)


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
check_library_exists("Xss" "XScreenSaverRegister" "/usr/lib" HAVE_XSS)


# CHECKME: Is there a better way to check this?
if(EXISTS "/proc/self/psinfo")
  set(HAVE__PROC_SELF_PSINFO TRUE)
endif()

if(EXISTS "/proc/self/stat")
  set(HAVE__PROC_SELF_STAT TRUE)
endif()

execute_process(COMMAND ${PROJECT_SOURCE_DIR}/config.guess
                 OUTPUT_VARIABLE CANONICAL_SYSTEM_NAME
                 OUTPUT_STRIP_TRAILING_WHITESPACE)

# Generate configuration headers
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/boinc/config.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/version.h.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/boinc/version.h)

add_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/boinc/svn_version.h
                     COMMAND ${CMAKE_COMMAND}
                       ARGS -DHEADER_BASENAME:STRING="svn_version"
                            -P ${MILKYWAYATHOME_CLIENT_CMAKE_MODULES}/FindGitVersion.cmake
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/boinc/
                     COMMENT "Creates/updates header SVN Revision Header")

################################################################################
# Libraries

set(BOINC_LIB_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/boinc/lib)
set(BOINC_API_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/boinc/api)
set(BOINC_CLIENT_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/boinc/client)


include_directories(${BOINC_API_SRC_DIR})
include_directories(${BOINC_LIB_SRC_DIR})

# FIXME: Things in libraries should not depend on things from here but they do
include_directories(${BOINC_CLIENT_SRC_DIR})

set(boinc_api_src ${BOINC_API_SRC_DIR}/boinc_api.cpp)
set(boinc_api_headers ${BOINC_API_SRC_DIR}/boinc_api.h)



set(boinc_lib_headers boinc/version.h
                      ${BOINC_LIB_SRC_DIR}/app_ipc.h
                      ${BOINC_LIB_SRC_DIR}/boinc_win.h
                      ${BOINC_LIB_SRC_DIR}/url.h
                      ${BOINC_LIB_SRC_DIR}/common_defs.h
                      ${BOINC_LIB_SRC_DIR}/diagnostics.h
                      ${BOINC_LIB_SRC_DIR}/diagnostics_win.h
                      ${BOINC_LIB_SRC_DIR}/filesys.h
                      ${BOINC_LIB_SRC_DIR}/hostinfo.h
                      ${BOINC_LIB_SRC_DIR}/md5_file.h
                      ${BOINC_LIB_SRC_DIR}/proxy_info.h
                      ${BOINC_LIB_SRC_DIR}/prefs.h
                      ${BOINC_LIB_SRC_DIR}/miofile.h
                      ${BOINC_LIB_SRC_DIR}/mfile.h
                      ${BOINC_LIB_SRC_DIR}/parse.h
                      ${BOINC_LIB_SRC_DIR}/util.h
                      ${BOINC_LIB_SRC_DIR}/coproc.h
                      ${BOINC_LIB_SRC_DIR}/cal_boinc.h
                      ${BOINC_LIB_SRC_DIR}/cl_boinc.h
                      boinc/svn_version.h
                      boinc/win_build/config.h
                      ${BOINC_LIB_SRC_DIR}/str_util.h
                      ${BOINC_LIB_SRC_DIR}/procinfo.h
                      ${BOINC_LIB_SRC_DIR}/proc_control.h)

set(boinc_lib_src ${BOINC_LIB_SRC_DIR}/util.cpp
                  ${BOINC_LIB_SRC_DIR}/app_ipc.cpp
                  ${BOINC_LIB_SRC_DIR}/diagnostics.cpp
                  ${BOINC_LIB_SRC_DIR}/filesys.cpp
                  ${BOINC_LIB_SRC_DIR}/hostinfo.cpp
                  ${BOINC_LIB_SRC_DIR}/md5.c
                  ${BOINC_LIB_SRC_DIR}/md5_file.cpp
                  ${BOINC_LIB_SRC_DIR}/mem_usage.cpp
                  ${BOINC_LIB_SRC_DIR}/mfile.cpp
                  ${BOINC_LIB_SRC_DIR}/miofile.cpp
                  ${BOINC_LIB_SRC_DIR}/parse.cpp
                  ${BOINC_LIB_SRC_DIR}/prefs.cpp
                  ${BOINC_LIB_SRC_DIR}/proxy_info.cpp
                  ${BOINC_LIB_SRC_DIR}/str_util.cpp
                  ${BOINC_LIB_SRC_DIR}/shmem.cpp
                  ${BOINC_LIB_SRC_DIR}/base64.cpp
                  ${BOINC_LIB_SRC_DIR}/url.cpp
                  ${BOINC_LIB_SRC_DIR}/coproc.cpp
                  ${BOINC_LIB_SRC_DIR}/procinfo.cpp
                  ${BOINC_LIB_SRC_DIR}/proc_control.cpp)

if(WIN32)
  list(APPEND boinc_lib_src ${BOINC_LIB_SRC_DIR}/win_util.cpp
                            ${BOINC_LIB_SRC_DIR}/diagnostics_win.cpp
                            ${BOINC_LIB_SRC_DIR}/stackwalker_win.cpp
                            ${BOINC_LIB_SRC_DIR}/procinfo_win.cpp)
endif()

if(UNIX)
  list(APPEND boinc_lib_src ${BOINC_LIB_SRC_DIR}/procinfo_unix.cpp)
endif()

if(APPLE)
  list(APPEND boinc_lib_src ${BOINC_LIB_SRC_DIR}/procinfo_mac.cpp)
endif()


set(boinc_graphics_src ${BOINC_API_SRC_DIR}/graphics2.cpp
                       ${BOINC_API_SRC_DIR}/graphics2_util.cpp)
set(boinc_graphics_headers ${BOINC_API_SRC_DIR}/graphics2.h)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/boinc)

if(APPLE)
  set(MAC_LIB_SRC_DIR ${BOINC_LIB_SRC_DIR}/mac)
  include_directories(${MAC_LIB_SRC_DIR})
  list(APPEND boinc_lib_src ${MAC_LIB_SRC_DIR}/mac_backtrace.cpp
                            ${MAC_LIB_SRC_DIR}/QBacktrace.c
                            ${MAC_LIB_SRC_DIR}/QCrashReport.c
                            ${MAC_LIB_SRC_DIR}/QMachOImage.c
                            ${MAC_LIB_SRC_DIR}/QMachOImageList.c
                            ${MAC_LIB_SRC_DIR}/QSymbols.c
                            ${MAC_LIB_SRC_DIR}/QTaskMemory.c)

  list(APPEND boinc_lib_headers ${MAC_LIB_SRC_DIR}/dyld_gdb.h
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

maybe_disable_ssen(boinc_api boinc boinc_graphics2)



# install(TARGETS boinc boinc_api boinc_graphics2
#         ARCHIVE       DESTINATION lib
#         PUBLIC_HEADER DESTINATION include/boinc)

# install(FILES ${boinc_lib_headers} ${boinc_api_headers} ${boinc_graphics_headers}
#         DESTINATION include/boinc)

if(OPENSSL_FOUND)
  include_directories(${OPENSSL_INCLUDE_DIR})
  add_library(boinc_crypt STATIC ${BOINC_LIB_SRC_DIR}/crypt.cpp ${BOINC_LIB_SRC_DIR}/crypt.h)
  set_target_properties(boinc_crypt PROPERTIES PUBLIC_HEADER ${BOINC_LIB_SRC_DIR}/crypt.h)
  # install(TARGETS boinc_crypt
  #           ARCHIVE       DESTINATION lib
  #           PUBLIC_HEADER DESTINATION include/boinc)
else()
  message(WARNING "OpenSSL not found. Not building libboinc_crypt")
endif()

set(BOINC_INCLUDE_DIR "${BOINC_API_SRC_DIR};${BOINC_LIB_SRC_DIR}" CACHE INTERNAL "BOINC headers" FORCE)
set(BOINC_LIBRARIES boinc_graphics2 boinc_api boinc CACHE INTERNAL "BOINC libraries" FORCE)
set(BOINC_GRAPHICS_LIBRARY boinc_graphics2 CACHE INTERNAL "BOINC graphics library" FORCE)
set(BOINC_GRAPHICS_FOUND TRUE CACHE INTERNAL "Built BOINC libraries" FORCE)


#export(TARGETS boinc boinc_api boinc_graphics2 FILE boinc-exports.cmake)
#export(PACKAGE BOINC)

# These random utilities do not belong in 'lib'
# add_executable(md5_test ${BOINC_LIB_SRC_DIR}/md5_test.cpp)
# target_link_libraries(md5_test boinc)

# # add_executable(shmem_test ${BOINC_LIB_SRC_DIR}/shmem_test.cpp)
# # target_link_libraries(shmem_test boinc)


# set(msg_test_src ${BOINC_LIB_SRC_DIR}/msg_test.cpp
#                  ${BOINC_LIB_SRC_DIR}/msg_queue.cpp)
# set(msg_test_headers ${BOINC_LIB_SRC_DIR}/msg_queue.h)
# add_executable(msg_test ${msg_test_src} ${msg_test_headers})
# target_link_libraries(msg_test boinc)

# include_directories(${OPENSSL_INCLUDE_DIR})
# add_executable(crypt_prog ${BOINC_LIB_SRC_DIR}/crypt_prog.cpp)
# target_link_libraries(crypt_prog boinc boinc_crypt ${OPENSSL_LIBRARIES})



################################################################################
# Client

set(boinc_client_src ${BOINC_CLIENT_SRC_DIR}/acct_mgr.cpp
                     ${BOINC_CLIENT_SRC_DIR}/acct_setup.cpp
                     ${BOINC_CLIENT_SRC_DIR}/app.cpp
                     ${BOINC_CLIENT_SRC_DIR}/app_control.cpp
                     ${BOINC_CLIENT_SRC_DIR}/app_graphics.cpp
                     ${BOINC_CLIENT_SRC_DIR}/app_start.cpp
                     # auto_update.cpp # Not in original Makefile.am
                     ${BOINC_CLIENT_SRC_DIR}/check_state.cpp
                     ${BOINC_CLIENT_SRC_DIR}/client_msgs.cpp
                     ${BOINC_CLIENT_SRC_DIR}/client_state.cpp
                     ${BOINC_CLIENT_SRC_DIR}/client_types.cpp
                     ${BOINC_CLIENT_SRC_DIR}/coproc_detect.cpp
                     ${BOINC_CLIENT_SRC_DIR}/cpu_sched.cpp
                     ${BOINC_CLIENT_SRC_DIR}/cs_account.cpp
                     ${BOINC_CLIENT_SRC_DIR}/cs_apps.cpp
                     ${BOINC_CLIENT_SRC_DIR}/cs_benchmark.cpp
                     ${BOINC_CLIENT_SRC_DIR}/cs_cmdline.cpp
                     ${BOINC_CLIENT_SRC_DIR}/cs_files.cpp
                     ${BOINC_CLIENT_SRC_DIR}/cs_notice.cpp
                     ${BOINC_CLIENT_SRC_DIR}/cs_platforms.cpp
                     ${BOINC_CLIENT_SRC_DIR}/cs_prefs.cpp
                     ${BOINC_CLIENT_SRC_DIR}/cs_proxy.cpp
                     ${BOINC_CLIENT_SRC_DIR}/cs_scheduler.cpp
                     ${BOINC_CLIENT_SRC_DIR}/cs_statefile.cpp
                     ${BOINC_CLIENT_SRC_DIR}/cs_trickle.cpp
                     ${BOINC_CLIENT_SRC_DIR}/current_version.cpp
                     ${BOINC_CLIENT_SRC_DIR}/dhrystone.cpp
                     ${BOINC_CLIENT_SRC_DIR}/dhrystone2.cpp
                     ${BOINC_CLIENT_SRC_DIR}/file_names.cpp
                     ${BOINC_CLIENT_SRC_DIR}/file_xfer.cpp
                     ${BOINC_CLIENT_SRC_DIR}/gui_http.cpp
                     ${BOINC_CLIENT_SRC_DIR}/gui_rpc_server.cpp
                     ${BOINC_CLIENT_SRC_DIR}/gui_rpc_server_ops.cpp
                     ${BOINC_CLIENT_SRC_DIR}/hostinfo_network.cpp
                     ${BOINC_CLIENT_SRC_DIR}/hostinfo_unix.cpp
                     ${BOINC_CLIENT_SRC_DIR}/http_curl.cpp
                     ${BOINC_CLIENT_SRC_DIR}/log_flags.cpp
                     ${BOINC_CLIENT_SRC_DIR}/main.cpp
                     ${BOINC_CLIENT_SRC_DIR}/net_stats.cpp
                     ${BOINC_CLIENT_SRC_DIR}/pers_file_xfer.cpp
                     ${BOINC_CLIENT_SRC_DIR}/rr_sim.cpp
                     ${BOINC_CLIENT_SRC_DIR}/sandbox.cpp
                     ${BOINC_CLIENT_SRC_DIR}/scheduler_op.cpp
                     ${BOINC_CLIENT_SRC_DIR}/time_stats.cpp
                     ${BOINC_CLIENT_SRC_DIR}/whetstone.cpp
                     ${BOINC_CLIENT_SRC_DIR}/work_fetch.cpp)

# app_stats_mac.cpp

set(boinc_client_headers ${BOINC_CLIENT_SRC_DIR}/acct_mgr.h
                         ${BOINC_CLIENT_SRC_DIR}/client_types.h
                         ${BOINC_CLIENT_SRC_DIR}/dhrystone.h
                         ${BOINC_CLIENT_SRC_DIR}/http_curl.h
                         ${BOINC_CLIENT_SRC_DIR}/sandbox.h
                         ${BOINC_CLIENT_SRC_DIR}/acct_setup.h
                         ${BOINC_CLIENT_SRC_DIR}/cpp.h
                         ${BOINC_CLIENT_SRC_DIR}/file_names.h
                         ${BOINC_CLIENT_SRC_DIR}/log_flags.h
                         ${BOINC_CLIENT_SRC_DIR}/scheduler_op.h
                         ${BOINC_CLIENT_SRC_DIR}/app.h
                         ${BOINC_CLIENT_SRC_DIR}/cpu_benchmark.h
                         ${BOINC_CLIENT_SRC_DIR}/file_xfer.h
                         ${BOINC_CLIENT_SRC_DIR}/main.h
                         ${BOINC_CLIENT_SRC_DIR}/sim.h
                         ${BOINC_CLIENT_SRC_DIR}/auto_update.h
                         ${BOINC_CLIENT_SRC_DIR}/cs_notice.h
                         ${BOINC_CLIENT_SRC_DIR}/gui_http.h
                         ${BOINC_CLIENT_SRC_DIR}/net_stats.h
                         ${BOINC_CLIENT_SRC_DIR}/sysmon_win.h
                         ${BOINC_CLIENT_SRC_DIR}/client_msgs.h
                         ${BOINC_CLIENT_SRC_DIR}/cs_proxy.h
                         ${BOINC_CLIENT_SRC_DIR}/gui_rpc_server.h
                         ${BOINC_CLIENT_SRC_DIR}/pers_file_xfer.h
                         ${BOINC_CLIENT_SRC_DIR}/time_stats.h
                         ${BOINC_CLIENT_SRC_DIR}/client_state.h
                         ${BOINC_CLIENT_SRC_DIR}/current_version.h
                         ${BOINC_CLIENT_SRC_DIR}/hostinfo_network.h
                         ${BOINC_CLIENT_SRC_DIR}/rr_sim.h
                         ${BOINC_CLIENT_SRC_DIR}/work_fetch.h)


set(switcher_src ${BOINC_CLIENT_SRC_DIR}/switcher.cpp)

# add_executable(boinc_client ${boinc_client_src} ${boinc_client_headers})
# target_link_libraries(boinc_client boinc)

# install(TARGETS boinc_client
#           RUNTIME DESTINATION bin)

