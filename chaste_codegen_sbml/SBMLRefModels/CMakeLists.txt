find_package(Chaste COMPONENTS cell_based)
chaste_do_project(SBMLRefModels)
target_compile_options(chaste_project_SBMLRefModels PRIVATE -Wno-unused-variable)
target_compile_options(chaste_project_SBMLRefModels PRIVATE -Wno-unused-but-set-variable)
