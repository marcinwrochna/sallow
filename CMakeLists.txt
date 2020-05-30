cmake_minimum_required(VERSION 3.13)
project (Sallow
		 VERSION 1.0
		 DESCRIPTION "An algorithm for finding treedepth decompositions.")

# Options that affect libraries as well.
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -O3 -march=native -DNDEBUG")

add_subdirectory(libs/GSL)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)  # Generate compile_commands.json for IDEs.

if (MSVC)
	set(W_FLAGS /W4)
else()
	set(W_FLAGS  -Wall -Wextra -Wpedantic
	-Wcast-align
	# -Wconversion
	#	-Wno-shorten-64-to-32 -Wno-sign-conversion -Wno-implicit-int-float-conversion # not supported by GCC
	-Wfloat-equal
	-Wformat=2
	# -Wmissing-prototypes # triggers the "extern C" signal handler, which can't be static.
	# -Wmissing-declarations # same in GCC
	-Woverlength-strings
	-Wshadow
	-Wunreachable-code
	-Wno-sign-compare
	-Wno-unused-function)
endif()

add_executable(sallow)
target_sources(sallow PRIVATE src/environment.h
                              src/flowCutter.h src/flowCutter.cpp
                              src/graph.h
                              src/graphUtils.h src/graphUtils.cpp
                              src/greedy.h src/greedy.cpp
                              src/io.h src/io.cpp
                              src/main.cpp
                              src/mainAlgo.h src/mainAlgo.cpp
                              src/minHeap.h
                              src/simpleCutter.h src/simpleCutter.cpp
                              src/tinyGraphDB.h src/tinyGraphDB.cpp
							  src/utils.h
							  src/queue.h)
target_include_directories(sallow PUBLIC "libs/GSL/include/")
target_link_libraries(sallow PUBLIC GSL)
set_property(TARGET sallow PROPERTY CXX_STANDARD 17)

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
	target_compile_options(sallow PRIVATE -O0 -fsanitize=address -fno-omit-frame-pointer -fno-optimize-sibling-calls ${W_FLAGS})
	target_link_options(sallow PRIVATE -fsanitize=address -fno-omit-frame-pointer -fno-optimize-sibling-calls)
	# `target_link_options` requires CMake >= 3.13, otherwise just add this to link_libraries instead
	# To enable LeakSanitizer: set -fsanitize=address,leak here and update ASAN_OPTIONS in .vscode/launch.json.
elseif(CMAKE_BUILT_TYPE STREQUAL "RelWithDebInfo")
	target_compile_options(sallow PRIVATE -DNDEBUG -O3 -march=native -DWITHGPERFTOOLS -g -fno-omit-frame-pointer ${W_FLAGS})
	target_link_options(sallow PRIVATE -lprofiler -g)
	set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
else()
	# Options for submission of static binary to optil.io.
	target_compile_options(sallow PRIVATE -DNDEBUG -O3 -march=ivybridge ${W_FLAGS})
	target_link_options(sallow PRIVATE -static)
	set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()


