find_package(Lua 5.3 EXACT REQUIRED)

list(APPEND LuaFortran-sources
  LuaFortran/dump_lua_fif_module.f90
  LuaFortran/flu_binding.f90
  LuaFortran/flu_kinds_module.f90
  LuaFortran/lua_fif.f90
  LuaFortran/lua_parameters.f90
  LuaFortran/wrap_lua_dump.c
  )

list(APPEND sources
  source/aot_err_module.f90
  source/aot_fun_declaration_module.f90
  source/aot_fun_module.f90
  source/aot_out_general_module.f90
  source/aot_out_module.f90
  source/aot_path_module.f90
  source/aot_references_module.f90
  source/aot_table_module.f90
  source/aot_table_ops_module.f90
  source/aot_top_module.f90
  source/aot_vector_module.f90
  source/aotus_module.f90
  source/dummy_extdouble_fun_module.f90
  source/dummy_extdouble_out_module.f90
  source/dummy_extdouble_table_module.f90
  source/dummy_extdouble_top_module.f90
  source/dummy_extdouble_vector_module.f90
  source/dummy_quadruple_fun_module.f90
  source/dummy_quadruple_out_module.f90
  source/dummy_quadruple_table_module.f90
  source/dummy_quadruple_top_module.f90
  source/dummy_quadruple_vector_module.f90
  )

add_library(aotus ${LuaFortran-sources} ${sources})
target_include_directories(aotus
  PRIVATE
    ${LUA_INCLUDE_DIR}
  )
target_link_libraries(aotus
  PUBLIC
    ${LUA_LIBRARIES}
  )
add_library(HANDE::aotus ALIAS aotus)
set_target_properties(aotus
  PROPERTIES
    POSITION_INDEPENDENT_CODE 1
  )
