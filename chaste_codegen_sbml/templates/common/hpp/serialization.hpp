namespace
{
  namespace serialization
  {
    /* Serialize information required to construct a {{class_name}}. */
    template <class Archive>
    inline void save_construct_data(
      Archive &ar,
      const {{class_name}} *t,
      const unsigned int file_version)
    {
      const std::vector<double> state_variables = t->rGetConstStateVariables();
      ar & state_variables;
    }

    /* De-serialize constructor parameters and intiialise a {{class_name}}. */
    template <class Archive>
    inline void load_construct_data(
      Archive &ar,
      {{class_name}} *t,
      const unsigned int file_version)
    {
      std::vector<double> state_variables;
      ar & state_variables;

      // Invoke inplace constructor to initialise instance
      ::new (t) {{class_name}}(state_variables);
    }
  }
}

/* Define SRN model using Wrappers. */
#include "SbmlSrnWrapperModel.hpp"
#include "SbmlSrnWrapperModel.cpp"

typedef SbmlSrnWrapperModel<{{class_name}}, {{num_state_vars}}> Goldbeter1991SrnModel;

// Declare identifiers for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT({{class_name}})
EXPORT_TEMPLATE_CLASS2(SbmlSrnWrapperModel, {{class_name}}, {{num_state_vars}})

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER({{class_name}})