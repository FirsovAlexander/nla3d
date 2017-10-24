# Collect sources into the variable MATH_SOURCES with 
# having to explicitly list each header and source file
file (GLOB MATH_SOURCES
	  "*.cpp")

file (GLOB MATH_HEADERS
	  "*.h")

  add_library(math SHARED ${MATH_SOURCES} ${CMAKE_SOURCE_DIR}/src/lib/sys.cpp)

# Creates a folder "libraries" and adds target project (math.vcproj) under it
set_property(TARGET math PROPERTY FOLDER "libraries")
set_target_properties(math PROPERTIES COTIRE_PREFIX_HEADER_INCLUDE_PATH
  "${CMAKE_SOURCE_DIR}/site-src")
cotire(math)

install(TARGETS math
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(FILES ${MATH_HEADERS}
     DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/math/)