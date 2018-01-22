add_library(hande-static STATIC
            $<TARGET_OBJECTS:HANDE::cephes>
            $<TARGET_OBJECTS:HANDE::coreutils>
            $<TARGET_OBJECTS:HANDE::dSFMT>
            $<TARGET_OBJECTS:HANDE::external>
            )

# Sources that depend on Lua
target_sources(hande-static
  PRIVATE
    interact.F90
    lua_hande.F90
    lua_hande_calc.f90
    lua_hande_calc_utils.F90
    lua_hande_fns.f90
    lua_hande_systems.f90
    lua_hande_utils.F90
    test.f90
)

target_sources(hande-static
  PRIVATE
    annihilation.f90
    basis.f90
    basis_types.f90
    bit_utils.F90
    blocking.F90
    bloom_handler.F90
    calc.F90
    calc_system_init.f90
    canonical_estimates.F90
    ccmc.f90
    ccmc_data.F90
    ccmc_death_spawning.f90
    ccmc_linked.f90
    ccmc_selection.f90
    ccmc_utils.F90
    check_input.f90
    chem_pot.f90
    core.f90
    ct_fciqmc.F90
    dealloc.F90
    death.f90
    determinant_enumeration.f90
    determinants.f90
    dmqmc.f90
    dmqmc_data.f90
    dmqmc_estimators.F90
    dmqmc_initialise_dm.F90
    dmqmc_procedures.F90
    energy_evaluation.F90
    excit_gen_hub_k.f90
    excit_gen_mol.f90
    excit_gen_op_hub_k.f90
    excit_gen_op_mol.f90
    excit_gen_real_lattice.f90
    excit_gen_ringium.f90
    excit_gen_ueg.f90
    excit_gens.f90
    excitations.F90
    fci_lanczos.F90
    fci_lapack.F90
    fci_utils.F90
    fciqmc.f90
    hamiltonian.f90
    hamiltonian_chung_landau.f90
    hamiltonian_data.f90
    hamiltonian_heisenberg.f90
    hamiltonian_hub_k.f90
    hamiltonian_hub_real.f90
    hamiltonian_molecular.f90
    hamiltonian_periodic_complex.f90
    hamiltonian_ringium.f90
    hamiltonian_ueg.f90
    hande_top_level.f90
    hdf5_system.F90
    heisenberg_estimators.F90
    hfs_data.f90
    hfs_fciqmc.f90
    hilbert_space.F90
    hubbard_k.f90
    idmqmc.f90
    ifciqmc.f90
    importance_sampling.f90
    importance_sampling_data.f90
    kpoints.f90
    load_balancing.F90
    logging.f90
    molecular_integral_types.f90
    molecular_integrals.F90
    momentum_sym_read_in.f90
    momentum_symmetry.f90
    non_blocking_comm.F90
    operators.F90
    particle_t_utils.f90
    pg_symmetry.f90
    proc_pointers.f90
    qmc.F90
    qmc_common.F90
    qmc_data.f90
    qmc_io.f90
    read_in.F90
    read_in_symmetry.f90
    real_lattice.f90
    reference_determinant.f90
    replica_rdm.F90
    restart_hdf5.F90
    restart_utils.F90
    ringium.f90
    search.f90
    semi_stoch.F90
    simple_fciqmc.f90
    sort.f90
    spawn_data.F90
    spawning.F90
    stoch_utils.f90
    symmetry.f90
    symmetry_types.f90
    system.f90
    ueg.f90
    ueg_types.f90
  )

target_compile_definitions(hande-static
  PRIVATE
    $<$<BOOL:${USE_INTRINSIC_POPCNT}>:USE_POPCNT>
    $<$<NOT:$<BOOL:${USE_HDF5}>>:DISABLE_HDF5>
    $<$<NOT:$<BOOL:${USE_LANCZOS}>>:DISABLE_LANCZOS>
    $<$<BOOL:${USE_MPI}>:PARALLEL>
    $<$<NOT:$<BOOL:${USE_ScaLAPACK}>>:DISABLE_SCALAPACK>
    $<$<BOOL:${USE_SINGLE_PRECISION}>:SINGLE_PRECISION>
    $<$<CONFIG:Debug>:DEBUG>
    )
target_compile_options(hande-static
  PRIVATE
    ${HANDE_Fortran_FLAGS}
    $<$<CONFIG:Debug>:${HANDE_Fortran_FLAGS_DEBUG}>
    $<$<CONFIG:Release>:${HANDE_Fortran_FLAGS_RELEASE}>
  )
target_include_directories(hande-static
  PUBLIC
    $<$<BOOL:${USE_HDF5}>:${HDF5_INCLUDE_DIRS}>
  )
target_link_libraries(hande-static
  PRIVATE
    HANDE::aotus
  PUBLIC
    "$<$<BOOL:${USE_LANCZOS}>:${TRLan_LIBRARIES}>"
    ${MATH_LIBS}
    "$<$<BOOL:${USE_MPI}>:${MPI_Fortran_LIBRARIES}>"
    "$<$<BOOL:${USE_UUID}>:${UUID_LIBRARIES}>"
    "$<$<BOOL:${USE_HDF5}>:${HDF5_Fortran_LIBRARIES}>"
  )

#set_target_properties(hande-shared PROPERTIES SOVERSION 1 #${PROJECT_VERSION_MAJOR}
#                                            OUTPUT_NAME "hande"
#                                            EXPORT_NAME "hande")

add_executable(hande.cmake.x core.f90)
target_link_libraries(hande.cmake.x
  PRIVATE
    hande-static
  )
set_target_properties(hande.cmake.x
  PROPERTIES
    POSITION_INDEPENDENT_CODE 1
  )

add_custom_command(TARGET hande.cmake.x
  POST_BUILD
  COMMAND
    ${CMAKE_COMMAND} -E copy_directory ${PROJECT_BINARY_DIR}/bin ${PROJECT_SOURCE_DIR}/bin
  COMMAND
    ${CMAKE_COMMAND} -E create_symlink ${PROJECT_SOURCE_DIR}/bin/hande.cmake.x ${PROJECT_SOURCE_DIR}/bin/hande.x
  COMMENT
    "Copying and symlinking executable hande.cmake.x to ${PROJECT_SOURCE_DIR}/bin/hande.x"
  )
