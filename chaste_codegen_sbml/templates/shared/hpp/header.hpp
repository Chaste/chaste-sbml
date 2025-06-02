#ifndef {{ header_guard }}
#define {{ header_guard }}

#include <cmath>
#include <iostream>

#include <boost/serialization/base_object.hpp>

#include "AbstractOdeSystem.hpp"
#include "ChasteSerialization.hpp"

class {{ ode_class_name }} : public AbstractOdeSystem
{
private:
    // (De-)serialize {{ ode_class_name }}
    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(AbstractOdeSystem);
    }

{% if compartments %}
    // COMPARTMENTS:
{% for compartment in compartments %}
    double {{ compartment["id"] }};
{% endfor %}
{% endif %}

    // CONST PARAMETERS:
{% for parameter in parameters %}
{% if parameter["is_state_parameter"] is false() %}
    const double {{ parameter["id"] }} = {{ parameter["value"] }};
{% endif %}
{% endfor %}

    // STATE VARIABLES:
{% for sp in species %}
{% if sp["is_state_variable"] is true() %}
    double {{ sp["id"] }}; // {{ sp["name"] }}
{% endif %}
{% endfor %}

    // STATE PARAMETERS:
{% for sp in species %}
{% if sp["is_state_parameter"] is true() %}
    double {{ sp["id"] }};
{% endif %}
{% endfor %}

{% for parameter in parameters %}
{% if parameter["is_state_parameter"] is true() %}
    double {{ parameter["id"] }};
{% endif %}
{% endfor %}

{% if events %}
    // EVENTS:
    std::vector<bool> eventsSatisfied;
{% endif %}

public:
    {{ ode_class_name }}(std::vector<double> stateVariables = std::vector<double>());

    ~{{ ode_class_name }}();

    void Init();

    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY);

{% if events %}
    bool CalculateStoppingEvent(double time, const std::vector<double>& rY);
{% endif %}

    // FUNCTION DEFINITIONS:
{% for fd in function_definitions %}
    double {{ fd["id"] }}({{ fd["args"] }});
{% endfor %}
};

namespace
{
    namespace serialization
    {
        // Provide constructor for serializing {{ ode_class_name }}
        template <class Archive>
        inline void save_construct_data(Archive &ar, const {{ode_class_name}} * t, const unsigned int version)
        {
            // Save data required to construct instance
            const std::vector<double> state_variables = t->rGetConstStateVariables();
            ar << state_variables;
        }

        // Provide constructor for de-serializing {{ ode_class_name }}
        template <class Archive>
        inline void load_construct_data(Archive &ar, {{ode_class_name}} * t, const unsigned int version)
        {
            // Retrieve data from archive required to construct new instance
            std::vector<double> state_variables;
            ar >> state_variables;

            // Invoke inplace constructor to initialise instance
            ::new (t){{ode_class_name}}(state_variables);
        }
    } // namespace serialization
} // namespace ...

// Define {{ wrapper_class_name }} using wrappers
#include "{{ wrapper_class_name }}.hpp"
#include "{{ wrapper_class_name }}.cpp"

typedef {{ wrapper_class_name }}<{{ ode_class_name }}, {{ num_state_vars }}> {{ model_class_name }};

// Declare identifiers for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT({{ ode_class_name }})
EXPORT_TEMPLATE_CLASS2({{ wrapper_class_name }}, {{ ode_class_name }}, {{ num_state_vars }})

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER({{ model_class_name }})

#endif // {{ header_guard }}
