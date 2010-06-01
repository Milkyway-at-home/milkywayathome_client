
find_path(POPT_INCLUDE_DIR popt.h)

find_library(POPT_LIBRARY popt)

if(POPT_INCLUDE_DIR AND POPT_LIBRARY)
   set(POPT_FOUND TRUE)
endif(POPT_INCLUDE_DIR AND POPT_LIBRARY)


if(POPT_FOUND)
   if(NOT Popt_FIND_QUIETLY)
      message(STATUS "Found POPT Library: ${POPT_LIBRARY}")
   endif(NOT Popt_FIND_QUIETLY)
else(POPT_FOUND)
   if(Popt_FIND_REQUIRED)
      message(FATAL_ERROR "Could not find POPT Library")
   endif(Popt_FIND_REQUIRED)
endif(POPT_FOUND)

