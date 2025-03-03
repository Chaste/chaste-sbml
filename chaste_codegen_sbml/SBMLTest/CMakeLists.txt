find_package(Chaste COMPONENTS crypt)
chaste_do_project(SBMLTest)
target_compile_options(chaste_project_SBMLTest PRIVATE -Wno-unused-variable)
target_compile_options(chaste_project_SBMLTest PRIVATE -Wno-unused-but-set-variable)
