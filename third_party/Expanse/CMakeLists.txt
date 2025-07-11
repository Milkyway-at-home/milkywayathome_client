# File: third_party/Expanse/CMakeLists.txt

# Manually specify the source files for the library.
set(EXPANSE_SOURCES
    src/bfe.c
    src/basis.c
    src/simulation.c
)

# Define the library as STATIC
add_library(Expanse STATIC ${EXPANSE_SOURCES})

# Specify all necessary include directories for this library.
target_include_directories(Expanse PUBLIC
    # Path to its own public headers
    "${CMAKE_CURRENT_SOURCE_DIR}/include"

    # Path to your vendored headers from MW@H
    "${CMAKE_CURRENT_SOURCE_DIR}/src/vendor"

    # Absolute path to the main milkyway library's headers
    "${CMAKE_SOURCE_DIR}/milkyway/include"
)