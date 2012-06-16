project (statisticaltests)

cmake_minimum_required (VERSION 2.8)
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules ${CMAKE_MODULE_PATH})

find_package(GSL REQUIRED)
find_package(Doxygen)

include_directories(${GSL_INCLUDE_DIRS})

set (SRCS
    src/main.cpp
    src/StatisticalTests.cpp)

add_executable (statisticaltests ${SRCS})
target_link_libraries (statisticaltests ${GSL_LIBRARIES})
install (TARGETS statisticaltests RUNTIME DESTINATION bin)

if(DOXYGEN_FOUND)
    add_custom_target(help ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile 
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen" VERBATIM
        )
endif()
