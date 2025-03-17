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
	set( CMAKE_CXX_FLAGS_DEBUG_INIT "/GS /W4 /openmp:llvm /wd26451 /fsanitize=address" )
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



file(GLOB source_files_utils exe/utils/*.cpp exe/utils/*.h)
foreach(source_file IN LISTS source_files_UTILS)
    get_filename_component(target_name ${source_file} NAME_WE)
    if (MSVC)
        add_compile_options(/W4 /openmp:llvm   /wd26451)
    endif()
	add_executable("${target_name}" "${source_file}" "${SOURCES}")
	target_link_libraries("${target_name}" ultimaille  algo_nicostuff OpenNL_psm $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>)
	set_target_properties( "${target_name}" PROPERTIES FOLDER "NICO/utils")
endforeach()

file(GLOB source_files_old exe/old/*.cpp exe/old/*.h)
foreach(source_file IN LISTS source_files_old)
    get_filename_component(target_name ${source_file} NAME_WE)
    if (MSVC)
        add_compile_options(/W4 /openmp:llvm   /wd26451)
    endif()
	add_executable("${target_name}" "${source_file}" "${SOURCES}")
	target_link_libraries("${target_name}" ultimaille  algo_nicostuff OpenNL_psm $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>)
	set_target_properties( "${target_name}" PROPERTIES FOLDER "NICO/old")
endforeach()

file(GLOB source_files_api exe/api/*.cpp exe/api/*.h)
foreach(source_file IN LISTS source_files_api)
    get_filename_component(target_name ${source_file} NAME_WE)
    if (MSVC)
        add_compile_options(/W4 /openmp:llvm   /wd26451)
    endif()
	add_executable("${target_name}" "${source_file}" "${SOURCES}")
	target_link_libraries("${target_name}" ultimaille  algo_nicostuff OpenNL_psm $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>)
	set_target_properties( "${target_name}" PROPERTIES FOLDER "NICO/api")
endforeach()

file(GLOB source_files_michelin exe/michelin/*.cpp exe/michelin/*.h)
foreach(source_file IN LISTS source_files_michelin)
    get_filename_component(target_name ${source_file} NAME_WE)
    if (MSVC)
        add_compile_options(/W4 /openmp:llvm   /wd26451)
    endif()
	add_executable("${target_name}" "${source_file}" "${SOURCES}")
	target_link_libraries("${target_name}" ultimaille  algo_nicostuff OpenNL_psm $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>)
	set_target_properties( "${target_name}" PROPERTIES FOLDER "NICO/michelin")
endforeach()


file(GLOB source_files_polygraph exe/polygraph/*.cpp exe/polygraph/*.h)
foreach(source_file IN LISTS source_files_POLYGRAPH)
    get_filename_component(target_name ${source_file} NAME_WE)
    if (MSVC)
        add_compile_options(/W4 /openmp:llvm   /wd26451)
    endif()
	add_executable("${target_name}" "${source_file}" "${SOURCES}")
	target_link_libraries("${target_name}" ultimaille  algo_nicostuff OpenNL_psm $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>)
	set_target_properties( "${target_name}" PROPERTIES FOLDER "NICO/polygraph")
endforeach()




