add_subdirectory(third_party/tetgen)

include_directories(algo)
file(GLOB SOURCES *.cpp *.h framework/*.cpp framework/*.h surface/*.cpp surface/*.h third_party/*.cpp third_party/*.h volume/*.cpp volume/*.h fullhex/*.cpp fullhex/*.h )
# if (WIN32)
	set(BUILD_SHARED_LIBS OFF)
# endif()
if (MSVC)
    add_compile_options(/W4 /openmp:llvm)
endif()
add_library(algo_nicostuff ${SOURCES})
target_link_libraries(algo_nicostuff ultimaille OpenNL_psm tet  OpenNL_psm $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>)

set_target_properties( algo_nicostuff PROPERTIES FOLDER "NICO")