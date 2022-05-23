# Do not proceed further if BLAS and LAPACK were not found
if(NOT BLAS_FOUND OR NOT LAPACK_FOUND)
  message(FATAL_ERROR "HANDE links against BLAS/LAPACK and a valid implementation was not found!")
endif()

add_library(hande-static STATIC
  $<TARGET_OBJECTS:HANDE::aotus>
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

# All other sources
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
    check_input.F90
    chem_pot.f90
    core.f90
    ct_fciqmc.F90
    dealloc.F90
    death.f90
    determinant_data.f90
    determinant_decoders.f90
    determinant_enumeration.f90
    determinants.f90
    dmqmc.f90
    dmqmc_data.f90
    dmqmc_estimators.F90
    dmqmc_initialise_dm.F90
    dmqmc_procedures.F90
    energy_evaluation.F90
    excit_gen_heat_bath_mol.F90
    excit_gen_hub_k.f90
    excit_gen_mol.f90
    excit_gen_op_hub_k.f90
    excit_gen_op_mol.f90
    excit_gen_power_pitzer_mol.F90
    excit_gen_real_lattice.f90
    excit_gen_ringium.f90
    excit_gens.f90
    excit_gen_ueg.f90
    excit_gen_utils.f90
    excitations.F90
    fci_lapack.F90
    fci_davidson.F90
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
    state_histograms.F90
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
    mp1.F90
    non_blocking_comm.F90
    operators.F90
    particle_t_utils.f90
    pg_symmetry.f90
    proc_pointers.f90
    propagators.f90
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
    $<$<BOOL:${USE_MPI}>:PARALLEL>
    $<$<NOT:$<BOOL:${HAS_MPI3}>>:DISABLE_MPI3>
    $<$<NOT:$<BOOL:${USE_ScaLAPACK}>>:DISABLE_SCALAPACK>
    $<$<BOOL:${USE_SINGLE_PRECISION}>:SINGLE_PRECISION>
    $<$<CONFIG:Debug>:DEBUG>
    )
target_include_directories(hande-static
  PUBLIC
    "$<$<BOOL:${USE_HDF5}>:${HDF5_INCLUDE_DIRS}>"
    "${LUA_INCLUDE_DIR}"
  )
target_link_libraries(hande-static
  PUBLIC
    "${MATH_LIBS}"
    "${EXPLICIT_LIBS}"
    "$<$<BOOL:${USE_MPI}>:${MPI_Fortran_LIBRARIES}>"
    "$<$<BOOL:${USE_UUID}>:PkgConfig::UUID>"
    "$<$<BOOL:${USE_HDF5}>:${HDF5_Fortran_LIBRARIES}>"
    "${LUA_LIBRARIES}"
  )
set_target_properties(hande-static
  PROPERTIES
    OUTPUT_NAME "hande"
  )

add_executable(${HANDE_EXE_NAME} core.f90)
target_link_libraries(${HANDE_EXE_NAME}
  PRIVATE
    hande-static
  )
add_custom_command(TARGET ${HANDE_EXE_NAME}
  POST_BUILD
  COMMAND
    ${CMAKE_COMMAND} -E copy_directory ${PROJECT_BINARY_DIR}/bin ${PROJECT_SOURCE_DIR}/bin
  COMMAND
  ${CMAKE_COMMAND} -E create_symlink ${PROJECT_SOURCE_DIR}/bin/${HANDE_EXE_NAME} ${PROJECT_SOURCE_DIR}/bin/hande.x
  COMMENT
    "Copying executable ${HANDE_EXE_NAME} to ${PROJECT_SOURCE_DIR}/bin and symlinking it to ${PROJECT_SOURCE_DIR}/bin/hande.x"
  )
install(
  TARGETS
    hande-static ${HANDE_EXE_NAME}
  ARCHIVE
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
    COMPONENT lib
  RUNTIME
    DESTINATION ${CMAKE_INSTALL_BINDIR}
    COMPONENT bin
  LIBRARY
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
    COMPONENT lib
  )
