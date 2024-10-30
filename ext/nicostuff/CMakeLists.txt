cmake_minimum_required(VERSION 3.0...3.9)



project(nicostuff)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_RPATH_USE_ORIGIN TRUE)

if (CMAKE_COMPILER_IS_GNUCC)
    set(CMAKE_CXX_FLAGS "-lstdc++fs") 
endif()
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

find_package(OpenMP COMPONENTS CXX)

if (WIN32)
	add_definitions(-DWIN32)
    add_definitions(-DNICOSTUFF_FRAMEWORK)
if(MSVC)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()
endif()

# Get development version of ultimaille directly from git repository 
include(FetchContent)
FetchContent_Declare(
  ultimaille
  GIT_REPOSITORY https://github.com/ssloy/ultimaille
  GIT_TAG        master
)
FetchContent_MakeAvailable(ultimaille) 

# Add ultimaille directory to get ultimaille CMakeLists
include_directories(${ultimaille_SOURCE_DIR})

include_directories(. algo)
add_subdirectory(algo)


file(GLOB source_files exe/*.cpp exe/*.h)
# Loop through each source file
foreach(source_file IN LISTS source_files)
    # Get the name of the file without the extension (e.g. 'a' from src/a.cpp'
    get_filename_component(target_name ${source_file} NAME_WE)
    
	    if (MSVC)
        add_compile_options(/W4 /openmp:llvm   /wd26451)
    endif()
	# Create an executable with the above name, building the above source
    add_executable("${target_name}" "${source_file}" "${SOURCES}")
	target_link_libraries("${target_name}" ultimaille  algo_nicostuff OpenNL_psm $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>)
	set_target_properties( "${target_name}" PROPERTIES FOLDER "NICO/exe")
endforeach()



file(GLOB source_files_UTILS exe/UTILS/*.cpp exe/UTILS/*.h)
foreach(source_file IN LISTS source_files_UTILS)
    get_filename_component(target_name ${source_file} NAME_WE)
    if (MSVC)
        add_compile_options(/W4 /openmp:llvm   /wd26451)
    endif()
	add_executable("${target_name}" "${source_file}" "${SOURCES}")
	target_link_libraries("${target_name}" ultimaille  algo_nicostuff OpenNL_psm $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>)
	set_target_properties( "${target_name}" PROPERTIES FOLDER "NICO/UTILS")
endforeach()

file(GLOB source_files_OLD exe/OLD/*.cpp exe/OLD/*.h)
foreach(source_file IN LISTS source_files_OLD)
    get_filename_component(target_name ${source_file} NAME_WE)
    if (MSVC)
        add_compile_options(/W4 /openmp:llvm   /wd26451)
    endif()
	add_executable("${target_name}" "${source_file}" "${SOURCES}")
	target_link_libraries("${target_name}" ultimaille  algo_nicostuff OpenNL_psm $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>)
	set_target_properties( "${target_name}" PROPERTIES FOLDER "NICO/OLD")
endforeach()

file(GLOB source_files_API exe/API/*.cpp exe/API/*.h)
foreach(source_file IN LISTS source_files_API)
    get_filename_component(target_name ${source_file} NAME_WE)
    if (MSVC)
        add_compile_options(/W4 /openmp:llvm   /wd26451)
    endif()
	add_executable("${target_name}" "${source_file}" "${SOURCES}")
	target_link_libraries("${target_name}" ultimaille  algo_nicostuff OpenNL_psm $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>)
	set_target_properties( "${target_name}" PROPERTIES FOLDER "NICO/API")
endforeach()





