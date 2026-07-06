chaste_do_test_project(SbmlRefModels)

# Precompile the heavy Chaste headers shared by every generated case test (the CVODE solver and
# test helpers), to speed up rebuilds. FakePetscSetup and the per-test model header are left out.
get_property(_sbml_test_targets DIRECTORY PROPERTY BUILDSYSTEM_TARGETS)
foreach(_sbml_test_target ${_sbml_test_targets})
    get_target_property(_sbml_test_type ${_sbml_test_target} TYPE)
    if (_sbml_test_type STREQUAL "EXECUTABLE")
        target_precompile_headers(${_sbml_test_target} PRIVATE
            <cxxtest/TestSuite.h>
            <CvodeAdaptor.hpp>
            <OdeSolution.hpp>
            <SbmlTestHelpers.hpp>)
    endif()
endforeach()
