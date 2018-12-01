cmake_minimum_required (VERSION 3.8)

project (Darwin)

#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /bigobj")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++11 -O4 -mbmi2 -msse4.2 -mavx2")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O4 -mbmi2 -msse4.2 -mavx2")

find_package(ZLIB REQUIRED)
find_package (Threads)

# To find and use zlib
find_path(ZLIB_INCLUDE_DIR zlib.h)
include_directories(${ZLIB_INCLUDE_DIR})

# To find and use tbb
include(${TBB_ROOT}/cmake/TBBBuild.cmake)

tbb_build(TBB_ROOT ${TBB_ROOT} CONFIG_DIR TBB_DIR MAKE_ARGS tbb_cpf=1)

find_package(TBB REQUIRED tbbmalloc tbbmalloc_proxy tbb_preview)

add_executable(Darwin 
Darwin_apply.cpp
Darwin_apply.h
Darwin_reflection.h
Darwin_types.cpp
Darwin_types.h
DRAM.cpp
DRAM.h
Index.cpp
Index.h
ntcoding.cpp
ntcoding.h
ConfigFile.cpp
ConfigFile.h
Chameleon.cpp
Chameleon.h
seed_pos_table.cpp
seed_pos_table.h
Processor.cpp
Processor.h
graph.h
extender.cpp
filter.cpp
maf_printer.cpp
seeder.cpp
sender.cpp
main.cpp)

target_link_libraries(Darwin PRIVATE ${TBB_IMPORTED_TARGETS})
target_link_libraries(Darwin PRIVATE ZLIB::ZLIB)
target_link_libraries (Darwin PRIVATE ${CMAKE_THREAD_LIBS_INIT})
