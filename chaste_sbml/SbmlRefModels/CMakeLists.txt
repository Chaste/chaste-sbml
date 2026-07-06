find_package(Chaste COMPONENTS cell_based)
chaste_do_project(SbmlRefModels)

# Precompile the heavy header every generated case model includes, to speed up rebuilds. The
# serialization export wrappers are deliberately excluded (they are order-sensitive and belong in
# each .cpp, not a shared PCH).
if (TARGET chaste_project_SbmlRefModels)
    target_precompile_headers(chaste_project_SbmlRefModels PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/src/AbstractSbmlOdeSystem.hpp")
endif()
