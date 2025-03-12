#ifndef {{header_guard}}
#define {{header_guard}}

#include <cmath>
#include <iostream>

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/shared_ptr.hpp>

#include "AbstractOdeSystem.hpp"
#include "ChasteSerialization.hpp"

class {{ode_class_name}} : public AbstractOdeSystem
{
private:

{% if compartments %}
    /* Declare model compartments. */
{% for compartment in compartments %}
    double {{compartment["id"]}}; // {{compartment["varname"]}}
{% endfor %}
{% endif %}

{% if parameters %}
    /* Declare model parameters. */
{% for parameter in parameters %}
    double {{parameter["id"]}}; // {{parameter["varname"]}}
{% endfor %}
{% endif %}

    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive &archive, const unsigned int version)
    {
        archive &boost::serialization::base_object<AbstractOdeSystem>(*this);
    }

public:

    /* Default constructor. */
    {{ode_class_name}}(std::vector<double> stateVariables = std::vector<double>());

    /* Destructor. */
    ~{{ode_class_name}}();

{% if function_definitions %}
    /* Declare model functions. */
{% for fd in function_definitions %}
    double {{fd["id"]}}({{ fd["args"]|join(', ') }});
{% endfor %}
{% endif %}

    void Init();

    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY);

{% if num_events > 0 %}
    // Stopping event is required for Cell Cycle Models to divide.
    bool CalculateStoppingEvent(double time, const std::vector<double>& rY);

    // Checks if any events have been triggered and updates the system accordingly.
    void CheckAndUpdateEvents(double time, const std::vector<double>& rY);

    // Checks whether all events are satisifed.
    bool AreAllEventsSatisfied(double time, const std::vector<double>& rY);
{% endif %}
};

namespace
{
namespace serialization
{
/* Serialize information required to construct a {{ode_class_name}}. */
template <class Archive>
inline void save_construct_data(
    Archive &ar, const {{ode_class_name}} *t, const unsigned int file_version)
{
    const std::vector<double> state_variables = t->rGetConstStateVariables();
    ar & state_variables;
}

/* De-serialize constructor parameters and intiialise a {{ode_class_name}}. */
template <class Archive>
inline void load_construct_data(
    Archive &ar, {{ode_class_name}} *t, const unsigned int file_version)
{
    std::vector<double> state_variables;
    ar & state_variables;

    // Invoke inplace constructor to initialise instance
    ::new (t) {{ode_class_name}}(state_variables);
}
}
} // namespace ...

/* Define cell cycle model using wrappers. */
#include "{{wrapper_class_name}}.hpp"
#include "{{wrapper_class_name}}.cpp"

typedef {{wrapper_class_name}}<{{ode_class_name}}, {{num_state_vars}}> {{model_class_name}};

/* Declare identifiers for the serializer */
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT({{ode_class_name}})
EXPORT_TEMPLATE_CLASS2({{wrapper_class_name}}, {{ode_class_name}}, {{num_state_vars}})

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER({{model_class_name}})

#endif // {{header_guard}}
