find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

list(APPEND HANDE_sources
  # ccmc
  ccmc/ccmc.F90
  ccmc/ccmc_death_spawning.f90
  ccmc/ccmc_linked.f90
  ccmc/ccmc_multispawn.F90
  ccmc/ccmc_selection.f90
  ccmc/ccmc_utils.f90
  # core: types and general utilities
  core/basis_types.f90
  core/bit_utils.F90
  core/calc.F90
  core/molecular_integral_types.f90
  core/search.f90
  core/sort.f90
  core/symmetry_types.f90
  core/ueg_types.f90
  # determinants
  determinants/bloom_handler.F90
  determinants/determinant_enumeration.f90
  determinants/determinants.f90
  determinants/excit_gen_hub_k.f90
  determinants/excit_gen_mol.f90
  determinants/excit_gen_op_hub_k.f90
  determinants/excit_gen_op_mol.f90
  determinants/excit_gen_real_lattice.f90
  determinants/excit_gen_ringium.f90
  determinants/excit_gen_ueg.f90
  determinants/excit_gens.f90
  determinants/excitations.F90
  determinants/reference_determinant.f90
  # dmqmc
  dmqmc/dmqmc.f90
  dmqmc/dmqmc_data.f90
  dmqmc/dmqmc_estimators.F90
  dmqmc/dmqmc_initialise_dm.F90
  dmqmc/dmqmc_procedures.F90
  dmqmc/idmqmc.f90
  # fci
  fci/ct_fciqmc.F90
  fci/fci_lanczos.F90
  fci/fci_lapack.F90
  fci/fci_utils.F90
  fci/fciqmc.f90
  fci/ifciqmc.f90
  fci/semi_stoch.F90
  fci/simple_fciqmc.f90
  # hellmann_feynman_sampling
  hellmann_feynman_sampling/hfs_data.f90
  hellmann_feynman_sampling/hfs_fciqmc.f90
  # importance_sampling
  importance_sampling/importance_sampling.f90
  importance_sampling/importance_sampling_data.f90
  # qmc
  qmc/annihilation.f90
  qmc/death.f90
  qmc/energy_evaluation.F90
  qmc/particle_t_utils.f90
  qmc/spawning.F90
  qmc/spawn_data.F90
  qmc/qmc.F90
  qmc/qmc_common.F90
  qmc/qmc_data.f90
  qmc/qmc_io.f90
  # system
  system/basis.f90
  system/chem_pot.f90
  system/hamiltonian.f90
  system/hamiltonian_chung_landau.f90
  system/hamiltonian_data.f90
  system/hamiltonian_heisenberg.f90
  system/hamiltonian_hub_k.f90
  system/hamiltonian_hub_real.f90
  system/hamiltonian_molecular.f90
  system/hamiltonian_periodic_complex.f90
  system/hamiltonian_ringium.f90
  system/hamiltonian_ueg.f90
  system/hdf5_system.F90
  system/hubbard_k.f90
  system/kpoints.f90
  system/molecular_integrals.F90
  system/momentum_sym_read_in.f90
  system/momentum_symmetry.f90
  system/operators.F90
  system/pg_symmetry.f90
  system/read_in.F90
  system/read_in_symmetry.f90
  system/real_lattice.f90
  system/ringium.f90
  system/symmetry.f90
  system/system.f90
  system/ueg.f90
  # top_level
  top_level/calc_system_init.f90
  top_level/canonical_estimates.F90
  top_level/check_input.f90
  top_level/dealloc.F90
  top_level/hande_top_level.f90
  top_level/heisenberg_estimators.F90
  top_level/hilbert_space.F90
  top_level/interact.F90
  top_level/load_balancing.F90
  top_level/logging.f90
  top_level/lua_hande.F90
  top_level/lua_hande_calc.f90
  top_level/lua_hande_fns.f90
  top_level/lua_hande_systems.f90
  top_level/lua_hande_utils.F90
  top_level/non_blocking_comm.F90
  top_level/proc_pointers.f90
  top_level/restart_hdf5.F90
  top_level/restart_utils.F90
  top_level/stoch_utils.f90
  top_level/test.f90
  )

#add_subdirectory(ccmc)
#add_subdirectory(core)
#add_subdirectory(determinants)
#add_subdirectory(dmqmc)
#add_subdirectory(fci)
#add_subdirectory(hellmann_feynman_sampling)
#add_subdirectory(importance_sampling)
#add_subdirectory(qmc)
#add_subdirectory(system)
#add_subdirectory(top_level)
#
#add_library(hande-shared SHARED
#  $<TARGET_OBJECTS:ccmc>
#  $<TARGET_OBJECTS:core>
#  $<TARGET_OBJECTS:determinants>
#  $<TARGET_OBJECTS:dmqmc>
#  $<TARGET_OBJECTS:fci>
#  $<TARGET_OBJECTS:hellmann_feynman_sampling>
#  $<TARGET_OBJECTS:importance_sampling>
#  $<TARGET_OBJECTS:qmc>
#  $<TARGET_OBJECTS:system>
#  $<TARGET_OBJECTS:top_level>
#  )
add_library(hande-shared SHARED ${HANDE_sources}
            $<TARGET_OBJECTS:HANDE::cephes>
            $<TARGET_OBJECTS:HANDE::coreutils>
            $<TARGET_OBJECTS:HANDE::dSFMT>
            $<TARGET_OBJECTS:HANDE::external>
            )
add_dependencies(hande-shared
  HANDE::aotus
  HANDE::cephes
  HANDE::coreutils
  HANDE::dSFMT
  HANDE::external
  )
if(NOT USE_HDF5)
  target_compile_definitions(hande-shared
    PUBLIC
      DISABLE_HDF5
    )
endif()
if(NOT USE_TRLAN)
  target_compile_definitions(hande-shared
    PUBLIC
      DISABLE_LANCZOS
    )
endif()
set_target_properties(hande-shared PROPERTIES SOVERSION 1 #${PROJECT_VERSION_MAJOR}
                                            OUTPUT_NAME "hande"
                                            EXPORT_NAME "hande")

                                          #install(TARGETS hande-shared
                                          #        EXPORT "${PROJECT_NAME}Targets-shared"
                                          #        RUNTIME DESTINATION ${CMAKE_INSTALL_LIBDIR}
                                          #        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

add_executable(hande.x hande.f90)
target_link_libraries(hande.x
  hande-shared
  HANDE::aotus
  ${LUA_LIBRARY}
  ${BLAS_LIBRARIES}
  ${LAPACK_LIBRARIES}
  ${UUID_LIBRARIES}
  )
