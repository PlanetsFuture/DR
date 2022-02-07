cmake_minimum_required(VERSION 3.2)
project(lalgebra C CXX)

set(SOURCE_FILES
  lu_decomposition.h
  lu_decomposition.cpp
  )

add_library(lalgebra SHARED STATIC ${SOURCE_FILES})

install(TARGETS lalgebra DESTINATION ${DRLIB_INSTALL_LIB_DIR})
install(FILES lu_decomposition.h DESTINATION ${DRLIB_INSTALL_INCLUDE_DIR})
