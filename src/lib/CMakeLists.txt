

add_subdirectory (math)

file (GLOB NLA3D_LIB_SOURCES "*.cpp" "materials/*.cpp" "elements/*.cpp")


file (GLOB NLA3D_LIB_HEADERS "*.h" )
file (GLOB NLA3D_LIB_MATERIALS_HEADERS "materials/*.h")
file (GLOB NLA3D_LIB_ELEMENTS_HEADERS "elements/*.h")

add_library(nla3d_lib SHARED ${NLA3D_LIB_SOURCES})

target_link_libraries(nla3d_lib math ${MKL_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

set_target_properties(nla3d_lib PROPERTIES COTIRE_PREFIX_HEADER_INCLUDE_PATH
  "${CMAKE_SOURCE_DIR}/site-src")
cotire(nla3d_lib)

install(TARGETS nla3d_lib
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(FILES ${NLA3D_LIB_HEADERS}
     DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(FILES ${NLA3D_LIB_MATERIALS_HEADERS}
     DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/materials)

install(FILES ${NLA3D_LIB_ELEMENTS_HEADERS}
     DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/elements)