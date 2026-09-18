chaste_do_test_project(SbmlRefModels)

# Give each test a root of its own to avoid MPI_INIT races under `ctest -j`.
# Keep this short: Open MPI puts socket paths under it, capped near 108 chars.
set(SbmlRefModels_MPI_SESSION_ROOT "/tmp/chaste-sbml-mpi" CACHE PATH
    "Root for the per-test MPI session directories")
get_property(sbml_ref_models_tests DIRECTORY PROPERTY TESTS)
foreach(sbml_ref_models_test IN LISTS sbml_ref_models_tests)
    set(session_dir "${SbmlRefModels_MPI_SESSION_ROOT}/${sbml_ref_models_test}")
    file(MAKE_DIRECTORY "${session_dir}")
    # TMPDIR is a fallback that works for all MPI implementations.
    # For Open MPI, a site-wide openmpi-mca-params.conf may set tmpdir_base,
    # which outranks TMPDIR; the environment form set here outranks the file.
    set_property(TEST ${sbml_ref_models_test} APPEND PROPERTY ENVIRONMENT
        "TMPDIR=${session_dir}"                      # any implementation
        "OMPI_MCA_orte_tmpdir_base=${session_dir}"   # Open MPI 3 and 4 (ORTE)
        "PRTE_MCA_prte_tmpdir_base=${session_dir}"   # Open MPI 5 (PRRTE)
        )
endforeach()
