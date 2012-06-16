project (statisticaltests)

cmake_minimum_required (VERSION 2.8)
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules ${CMAKE_MODULE_PATH})

find_package(GSL REQUIRED)

include_directories(${GSL_INCLUDE_DIRS})

set (SRCS
    src/main.cpp
    src/StatisticalTests.cpp)

add_executable (statisticaltests ${SRCS})
target_link_libraries (statisticaltests ${GSL_LIBRARIES})
install (TARGETS statisticaltests RUNTIME DESTINATION bin)
