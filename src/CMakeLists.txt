cmake_minimum_required(VERSION 3.2)
project(drcomp)

#add_subdirectory(const)
add_subdirectory(lalgebra)

set(SOURCE_FILES
  main.cpp
  dr.h
  )

add_executable(drcomp ${SOURCE_FILES})
target_link_libraries(drcomp lalgebra)
install(TARGETS drcomp DESTINATION ${DRLIB_INSTALL_BIN_DIR})
