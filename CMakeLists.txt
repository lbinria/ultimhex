cmake_minimum_required(VERSION 3.12...3.26)

project(ultimhex)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# set(
#   GEOGRAM_SOURCE_DIR "${CMAKE_SOURCE_DIR}/_deps/geogram-src/"
#   CACHE PATH "full path to the Geogram (or Vorpaline) installation"
# )

# Get config
get_property(is_multi_config GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(NOT is_multi_config AND NOT (CMAKE_BUILD_TYPE OR DEFINED ENV{CMAKE_BUILD_TYPE}))
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Release default")
endif()

# Try to find OpenMP package (if installed)
find_package(OpenMP COMPONENTS CXX)

# Display warning and error levels
if (MSVC)
    # warning level 4 (and all warnings as errors, /WX)
    add_compile_options(/W4)
else()
    # lots of warnings and all warnings as errors
    add_compile_options(-Wall -Wextra -pedantic)
endif()


# Get development version of ultimaille directly from git repository 
include(FetchContent)
FetchContent_Declare(
  ultimaille
  GIT_REPOSITORY https://github.com/ssloy/ultimaille
  GIT_TAG        v2.1.0
)
FetchContent_MakeAvailable(ultimaille)

# Add Geogram
add_subdirectory(ext/geogram)
# Add Nicostuff
add_subdirectory(ext/nicostuff)

# Add executable binary
file(GLOB SOURCES *.h *.cpp tools/*.h tools/*.cpp)
add_executable(ultimhex ${SOURCES})

# Add ultimaille directory to get ultimaille CMakeLists
include_directories(${ultimaille_SOURCE_DIR})

# Include ext
target_include_directories(ultimhex PRIVATE 
  ext/nicostuff/algo/framework
  ext/nicostuff/algo
  ext/nicostuff
  ext
)



# Link binary with ultimaille
target_link_libraries(ultimhex PRIVATE ${CMAKE_DL_LIBS} ultimaille geogram geogram_gfx algo_nicostuff m $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>)