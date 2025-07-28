#ifndef {{ header_guard }}
#define {{ header_guard }}

#include <vector>

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

    // CONSTANT PARAMETERS
{% for param in constant_parameters %}
    const double {{ param["id"] }} = {{ param["value"] }}; // {{ param["label"] }}
{% endfor %}

    // STATE VARIABLES
{% for var in state_variables %}
    double {{ var["id"] }}; // {{ var["label"] }}
{% endfor %}

    // DERIVED QUANTITIES
{% for dq in derived_quantities %}
    double {{ dq["id"] }}; // {{ dq["label"] }}
{% endfor %}

    // VARIABLE PARAMETERS
{% for param in variable_parameters %}
    double {{ param["id"] }}; // {{ param["label"] }}
{% endfor %}

    // RULE-BASED PARAMETERS
{% for param in rule_based_parameters %}
    double {{ param["id"] }}; // {{ param["label"] }}
{% endfor %}

    // REACTIONS
{% for reaction in reactions %}
    double {{ reaction["id"] }}; // {{ reaction["label"] }}
{% endfor %}

{% if events %}
    // EVENTS
    std::vector<bool> mEventsSatisfied;
    bool mEventsInitialised;

    std::vector<bool> mStatesAdjusted;
    std::vector<double> mStatesAdjustedValues;
{% endif %}

public:
    {{ ode_class_name }}(std::vector<double> stateVariables = std::vector<double>());

    ~{{ ode_class_name }}();

    void AdjustOdeParameters(double time);

    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY);
    void ProcessRules(double time, const std::vector<double>& rY);

{% if derived_quantities %}
    std::vector<double> ComputeDerivedQuantities(double time, const std::vector<double> &rY);
{% endif %}

{% if events %}
    double CalculateRootFunction(double time, const std::vector<double>& rY);
    bool CalculateStoppingEvent(double time, const std::vector<double>& rY);
    double ProcessEvents(double time, const std::vector<double>& rY);
{% endif %}

    // FUNCTIONS
{% for func in functions %}
    inline double {{ func["id"] }}({{ func["args"] }});
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

typedef {{ wrapper_class_name }}<{{ ode_class_name }}, {{ state_variables|length }}> {{ model_class_name }};

// Declare identifiers for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT({{ ode_class_name }})
EXPORT_TEMPLATE_CLASS2({{ wrapper_class_name }}, {{ ode_class_name }}, {{ state_variables|length }})

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER({{ model_class_name }})

#endif // {{ header_guard }}
