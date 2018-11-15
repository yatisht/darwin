cmake_minimum_required (VERSION 3.8)

project (Darwin)

#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /bigobj")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++11 -mbmi2 -msse4.2 -mavx2 -ltbb -ltbbmalloc")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -mbmi2 -msse4.2 -mavx2 -ltbb -ltbbmalloc")

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

target_link_libraries(Darwin PRIVATE TBB::tbb)
target_link_libraries(Darwin PRIVATE ZLIB::ZLIB)
target_link_libraries (Darwin PRIVATE ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries (Darwin PRIVATE tbbmalloc)

#find_library(INIH_LIBRARY libinih)
#link_libraries(${INIH_LIBRARY})

#set(CMAKE_PREFIX_PATH /usr/local/lib)
#find_library(TBB_MALLOC_LIBRARY libtbbmalloc)
#message (${TBB_MALLOC_LIBRARY})
#link_libraries(${TBB_MALLOC_LIBRARY})
#link_directories(/usr/local/lib)
