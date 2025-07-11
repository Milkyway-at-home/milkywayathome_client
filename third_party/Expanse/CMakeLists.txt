# File: third_party/Expanse/CMakeLists.txt

# Automatically find all source files in the src/ directory
file(GLOB_RECURSE EXPANSE_SOURCES "src/*.c")

# Define your library as a STATIC library called 'Expanse'
add_library(Expanse STATIC ${EXPANSE_SOURCES})

# Specify the include directories your library needs.
# This makes its headers available to the rest of the MW@H project.
target_include_directories(Expanse PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
)
