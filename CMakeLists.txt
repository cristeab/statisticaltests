project (statisticaltests)

cmake_minimum_required (VERSION 2.8)

set (SRCS
    src/main.cpp
    src/StatisticalTests.cpp)

add_executable (statisticaltests ${SRCS})
install (TARGETS statisticaltests RUNTIME DESTINATION bin)
