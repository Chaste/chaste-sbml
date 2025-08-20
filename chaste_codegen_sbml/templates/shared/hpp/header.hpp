#ifndef {{ header_guard }}
#define {{ header_guard }}

#include <vector>

#include <boost/serialization/base_object.hpp>

#include "AbstractOdeSystem.hpp"
#include "ChasteSerialization.hpp"
#include "SbmlEventType.hpp"

class {{ ode_class_name }} : public AbstractOdeSystem
{
private:
    // (De-)Serialize {{ ode_class_name }}
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

    // EVENTS
    std::vector<bool> mEventSatisfied;
    std::vector<bool> mEventTriggered;
    std::vector<SbmlEventType> mEventType;
    std::vector<bool> mEventAdjustedStateVars;
    std::vector<double> mEventAdjustedStateValues;
    std::vector<bool> mEventAdjustedParameters;
    std::vector<double> mEventAdjustedParameterValues;

public:
    /** 
     * Default constructor
     * 
     * @param stateVariables Initial state variables (optional)
     */
    {{ ode_class_name }}(std::vector<double> stateVariables = std::vector<double>());

    /**
     * Copy constructor
     * 
     * @param rOdeSystem Reference to the original instance
     */
    {{ ode_class_name }}(const {{ ode_class_name }}& rOdeSystem);

    /**
     * Destructor
     */
    ~{{ ode_class_name }}();

    /**
     * Adjust parameters and state variables after a stopping event
     * 
     * @param time The current time
     */
    void AdjustParameters(double time);

    /**
     * Calculate whether the conditions to trigger an event have been met
     * (Used by CVODE solver to find exact stopping position)
     * 
     * @param time The current time
     * @param rY The current state variables
     * 
     * @return How close we are to the root of the stopping condition
     */
    double CalculateRootFunction(double time, const std::vector<double>& rY) override;

    /** 
     * Calculate whether the conditions to trigger an event have been met
     * 
     * @param time The current time
     * @param rY The current state variables
     * 
     * @return True if conditions for an event are met, false otherwise
     */
    bool CalculateStoppingEvent(double time, const std::vector<double>& rY) override;

    /**
     * Compute the derived quantities from the given system state.
     *
     * @param time  the time at which to compute the derived quantities
     * @param rY a vector of values for the state variables
     * 
     * @return a vector of derived quantities
     */
    std::vector<double> ComputeDerivedQuantities(double time, const std::vector<double> &rY);

    /**
     * Compute the RHS of the ODE system.
     * 
     * An ODE solver will call this function repeatedly to solve for y = [y1 ... yn].
     *
     * @param time the time used to evaluate the RHS.
     * @param rY an input solution vector used to evaluate the RHS.
     * @param rDY an output vector to be filled in with the resulting derivatives y' = [y1' ... yn'].
     */
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY) override;

    /**
     * Check if a specific type of event has occurred.
     *
     * @param eventType The type of event to check
     * 
     * @return True if the type of event has occurred, false otherwise
     */
    bool HasEventOccurred(SbmlEventType eventType);

    /**
     * Process the events in the model.
     * 
     * @param time The current time
     * @param rY The current state variables
     * 
     * @return How close we are to the time of the next event
     */
    double ProcessModelEvents(double time, const std::vector<double>& rY);

    /**
     * Reset the flags that indicate which events have been triggered.
     */
    void ResetEventsOccurred();

    /** 
     * Run the equations governing the model to update state.
     * 
     * @param time The current time
     * @param rY The current state variables
     */
    void RunModelRules(double time, const std::vector<double>& rY);

    // MODEL FUNCTIONS
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
