
find_path(BOINC_INCLUDE_DIR boinc_api.h /usr/local/include/boinc )

find_library(BOINC_LIBRARY boinc)
find_library(BOINC_API_LIBRARY boinc_api)

if(BOINC_LIBRARY AND BOINC_API_LIBRARY)
    set(BOINC_LIBRARIES ${BOINC_LIBRARY} ${BOINC_API_LIBRARY})
endif(BOINC_LIBRARY AND BOINC_API_LIBRARY)

if(BOINC_INCLUDE_DIR AND BOINC_LIBRARIES)
   set(BOINC_FOUND TRUE)
endif(BOINC_INCLUDE_DIR AND BOINC_LIBRARIES)


if(BOINC_FOUND)
   if(NOT Boinc_FIND_QUIETLY)
      message(STATUS "Found BOINC Libraries: ${BOINC_LIBRARIES}")
   endif(NOT Boinc_FIND_QUIETLY)
else(BOINC_FOUND)
   if(Boinc_FIND_REQUIRED)
      message(FATAL_ERROR "Could not find BOINC Libraries")
   endif(Boinc_FIND_REQUIRED)
endif(BOINC_FOUND)

