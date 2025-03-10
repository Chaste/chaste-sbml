
#ifndef {{header_guard}}
#define {{header_guard}}

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/shared_ptr.hpp>

#include <cmath>
#include <iostream>

#include "AbstractOdeSystem.hpp"

class {{ode_system_name}} : public AbstractOdeSystem
{
private:

    /* Declare model compartments. */
    {{compartment_decls}}

    /* Declare model parameters. */
    {{parameter_decls}}

    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive &archive, const unsigned int version)
    {
        archive &boost::serialization::base_object<AbstractOdeSystem>(*this);
    }

public:

    /* Default constructor. */
    {{ode_system_name}}(std::vector<double> stateVariables = std::vector<double>());

    /* Destructor. */
    ~{{ode_system_name}}();

    /* Declare model functions. */
    {{function_decls}}

    void Init();

    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY);

};

namespace
{
namespace serialization
{
/* Serialize information required to construct a {{ode_system_name}}. */
template <class Archive>
inline void save_construct_data(
    Archive &ar, const {{ode_system_name}} *t, const unsigned int file_version)
{
    const std::vector<double> state_variables = t->rGetConstStateVariables();
    ar & state_variables;
}

/* De-serialize constructor parameters and intiialise a {{ode_system_name}}. */
template <class Archive>
inline void load_construct_data(
    Archive &ar, {{ode_system_name}} *t, const unsigned int file_version)
{
    std::vector<double> state_variables;
    ar & state_variables;

    // Invoke inplace constructor to initialise instance
    ::new (t) {{ode_system_name}}(state_variables);
}
}
} // namespace ...

/* Define SRN model using wrappers. */
#include "SbmlSrnWrapperModel.hpp"
#include "SbmlSrnWrapperModel.cpp"

typedef SbmlSrnWrapperModel<{{ode_system_name}}, {{num_state_vars}}> {{srn_name}};

/* Declare identifiers for the serializer */
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT({{ode_system_name}})
EXPORT_TEMPLATE_CLASS2(SbmlSrnWrapperModel, {{ode_system_name}}, {{num_state_vars}})

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER({{srn_name}})

#endif // {{header_guard}}
