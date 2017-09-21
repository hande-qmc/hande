add_subdirectory(ccmc)
add_subdirectory(lua_hande)

add_library(hande-shared SHARED
  $<TARGET_OBJECTS:lua_hande>
  $<TARGET_OBJECTS:ccmc>
  )
set_target_properties(hande-shared PROPERTIES SOVERSION 1 #${PROJECT_VERSION_MAJOR}
                                            OUTPUT_NAME "hande"
                                            EXPORT_NAME "hande")

                                          #install(TARGETS hande-shared
                                          #        EXPORT "${PROJECT_NAME}Targets-shared"
                                          #        RUNTIME DESTINATION ${CMAKE_INSTALL_LIBDIR}
                                          #        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

add_executable(hande.x hande.f90)
target_link_libraries(hande.x hande)
