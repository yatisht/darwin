cmake_minimum_required (VERSION 3.8)

project (Darwin)

#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /bigobj")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}-mbmi2 -msse4.2 -mavx2")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mbmi2 -msse4.2 -mavx2")

#find_package(bond REQUIRED)
find_package(TBB CONFIG REQUIRED)
find_package(ZLIB REQUIRED)
find_package (Threads)

# To find and use zlib
find_path(ZLIB_INCLUDE_DIR zlib.h)
include_directories(${ZLIB_INCLUDE_DIR})

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

target_link_libraries(Darwin PRIVATE TBB::tbb)
target_link_libraries(Darwin PRIVATE ZLIB::ZLIB)
target_link_libraries (Darwin PRIVATE ${CMAKE_THREAD_LIBS_INIT})

find_library(INIH_LIBRARY libinih)
link_libraries(${INIH_LIBRARY})

find_library(TBB_MALLOC_LIBRARY libtbbmalloc)
link_libraries(${TBB_MALLOC_LIBRARY})
