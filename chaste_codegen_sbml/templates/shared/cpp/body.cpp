#include <cmath>
#include <limits>
#include <vector>

#include "CellwiseOdeSystemInformation.hpp"
#include "SbmlEventType.hpp"
#include "SbmlMath.hpp"

#include "{{ model_hpp_file }}"

namespace sm = sbmlmath;

{{ ode_class_name }}::{{ ode_class_name }}(std::vector<double> stateVariables)
    : AbstractOdeSystem({{ state_variables|length }})
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<{{ ode_class_name }}>);

    // STATE VARIABLES
{% for var in state_variables %}
    {{ var["id"] }} = {{ var["initial_value"] }};
{% endfor %}

{% for var in state_variables %}
    SetDefaultInitialCondition({{ var["index"] }}, {{ var["id"] }});
{% endfor %}

    if (stateVariables.size() == {{ state_variables|length }})
    {
{% for var in state_variables %}
        {{ var["id"] }} = stateVariables[{{ var["index"] }}];
{% endfor %}
    }
    else if (stateVariables.size() != 0)
    {
        EXCEPTION("{{ ode_class_name }}: Expected {{ state_variables|length }} state variables, got " + std::to_string(stateVariables.size()));
    }

{% for var in state_variables %}
    mStateVariables.push_back({{ var["id"] }});
{% endfor %}

    // DERIVED QUANTITIES
{% for dq in derived_quantities %}
    {{ dq["id"] }} = 0.0;
{% endfor %}

    // VARIABLE PARAMETERS
{% for param in variable_parameters %}
    {{ param["id"] }} = {{ param["initial_value"] }};
{% endfor %}

{% for param in variable_parameters %}
    mParameters.push_back({{ param["id"] }});
{% endfor %}

    // RULE-BASED PARAMETERS
{% for var in rule_based_parameters %}
    {{ var['id'] }} = 0.0;
{% endfor %}

    // REACTIONS
{% for reaction in reactions %}
    {{ reaction["id"] }} = 0.0;
{% endfor %}

    // EVENTS
    mEventType.resize({{ events|length }}, SbmlEventType::UNKNOWN);

    // Uncomment lines below for events that should trigger cell division

{% for event in events %}
{% if event["type"] == EventType.CELL_DIVISION %}
    mEventType[{{ event["index"] }}] = SbmlEventType::CELL_DIVISION; // {{ event["label"] }}
{% else %}
    // mEventType[{{ event["index"] }}] = SbmlEventType::CELL_DIVISION; // {{ event["label"] }}
{% endif %}
{% endfor %}

    mEventSatisfied.resize({{ events|length }}, true); // Prevent events from triggering at the start
    mEventTriggered.resize({{ events|length }}, false);

    mEventAdjustedParameters.resize({{ variable_parameters|length }}, false);
    mEventAdjustedParameterValues.resize({{ variable_parameters|length }}, 0.0);

    mEventAdjustedStateVars.resize({{ state_variables|length }}, false);
    mEventAdjustedStateValues.resize({{ state_variables|length }}, 0.0);

    // Run model rules to complete state initialisation
    RunModelRules(0.0, mStateVariables);
}

{{ ode_class_name }}::{{ ode_class_name }}(const {{ ode_class_name }}& rOdeSystem)
        : {{ ode_class_name }}(rOdeSystem.mStateVariables)
{
    mEventSatisfied = rOdeSystem.mEventSatisfied;
    mEventTriggered = rOdeSystem.mEventTriggered;

    mEventAdjustedParameters = rOdeSystem.mEventAdjustedParameters;
    mEventAdjustedParameterValues = rOdeSystem.mEventAdjustedParameterValues;

    mEventAdjustedStateVars = rOdeSystem.mEventAdjustedStateVars;
    mEventAdjustedStateValues = rOdeSystem.mEventAdjustedStateValues;
}

{{ ode_class_name }}::~{{ ode_class_name }}()
{
}

void {{ ode_class_name }}::AdjustParameters(double time)
{
    for (unsigned i = 0; i < mEventAdjustedParameters.size(); ++i)
    {
        if (mEventAdjustedParameters[i])
        {
            SetParameter(i, mEventAdjustedParameterValues[i]);
        }
    }

    for (unsigned i = 0; i < mEventAdjustedStateVars.size(); ++i)
    {
        if (mEventAdjustedStateVars[i])
        {
            SetStateVariable(i, mEventAdjustedStateValues[i]);
            mEventAdjustedStateVars[i] = false;
        }
    }
}

double {{ ode_class_name }}::CalculateRootFunction(double time, const std::vector<double> &rY)
{
    return ProcessModelEvents(time, rY);
}

bool {{ ode_class_name }}::CalculateStoppingEvent(double time, const std::vector<double> &rY)
{
    return ProcessModelEvents(time, rY) == 0.0;
}

std::vector<double> {{ ode_class_name }}::ComputeDerivedQuantities(double time, const std::vector<double> &rY)
{
    RunModelRules(time, rY);

    std::vector<double> dqs;
{% for dq in derived_quantities %}
    dqs.push_back({{ dq["id"] }});
{% endfor %}
    return dqs;
}

void {{ ode_class_name }}::EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY)
{
    RunModelRules(time, rY);

{% for var in state_variables %}
    rDY[{{ var["index"] }}] = {{ var["rhs"] }}; // d[{{ var["id"] }}]/dt
{% endfor %}

    // TODO: Scale time appropriately
}

bool {{ ode_class_name }}::HasEventOccurred(SbmlEventType eventType)
{
    for (unsigned i = 0; i < mEventTriggered.size(); ++i)
    {
        if (mEventTriggered[i] && mEventType[i] == eventType)
        {
            return true;
        }
    }
    return false;
}

double {{ ode_class_name }}::ProcessModelEvents(double time, const std::vector<double> &rY)
{
    std::fill(std::begin(mEventAdjustedParameters), std::end(mEventAdjustedParameters), false);
    std::fill(std::begin(mEventAdjustedStateVars), std::end(mEventAdjustedStateVars), false);

    double min_dist = std::numeric_limits<double>::max();

{% for event in events %}
    //========================================
    // EVENT: {{ event["label"] }}
    //========================================
    {
        double event_dist = {{ event["distance"] }};

        // Avoid oscillation by ensuring event_dist is not close to 0 unless triggered
        if (std::abs(event_dist) < 1.0)
        {
            event_dist = 1.0;
        }

        // Update min_dist
        if (std::abs(event_dist) < std::abs(min_dist))
        {
            min_dist = event_dist;
        }

        // Process the event
        if ({{ event["trigger"] }})
        {
            if (!mEventSatisfied[{{ event["index"] }}])
            {
                // The condition is transitioning from false -> true: trigger the event
                mEventTriggered[{{ event["index"] }}] = true;
                event_dist = 0.0;
                min_dist = 0.0;

                // Adjust relevant state variables and parameters
{% for assignment in event["assignments"] %}
{% if ( assignment["type"] == VarType.STATE_VARIABLE ) %}
                // {{ assignment["lhs"] }} = {{ assignment["rhs"] }}
                mEventAdjustedStateVars[{{ assignment["index"] }}] = true;
                mEventAdjustedStateValues[{{ assignment["index"] }}] = {{ assignment["rhs"] }};

{% elif ( assignment["type"] == VarType.VARIABLE_PARAMETER ) %}
                // {{ assignment["lhs"] }} = {{ assignment["rhs"] }}
                mEventAdjustedParameters[{{ assignment["index"] }}] = true;
                mEventAdjustedParameterValues[{{ assignment["index"] }}] = {{ assignment["rhs"] }};

{% else %}
                {{ assignment["lhs"] }} = {{ assignment["rhs"] }}; {# TODO: does this case exist? #}

{% endif %}
{% endfor %} {# 'for assignment in event["assignments"]' #}

            }
            mEventSatisfied[{{ event["index"] }}] = true;
        }
        else
        {
            mEventSatisfied[{{ event["index"] }}] = false;
            mEventTriggered[{{ event["index"] }}] = false;
        }
    }

{% endfor %} {# 'for event in events' #}

    return min_dist; // Distance to closest event
}

void {{ ode_class_name }}::ResetEventsOccurred()
{
    std::fill(mEventTriggered.begin(), mEventTriggered.end(), false);
}

void {{ ode_class_name }}::RunModelRules(double time, const std::vector<double>& rY)
{
    // STATE VARIABLES
{% for var in state_variables %}
    {{ var["id"] }} = rY[{{ var["index"] }}];
{% endfor %}

    // VARIABLE PARAMETERS
{% for var in variable_parameters %}
    {{ var["id"] }} = GetParameter({{ var['index'] }});
{% endfor %}

    // ASSIGNMENT RULES
{% for rule in assignment_rules %}
    {{ rule["lhs"] }} = {{ rule["rhs"] }};
{% endfor %}

    // REACTIONS
{% for reaction in reactions %}
  {% if reaction["label"] %}
    // {{ reaction["label"] }}
  {% endif %}
  {% if reaction["parameters"] %}
    {{ reaction["id"] }} = 0.0;
    {
  {% for param in reaction["parameters"] %}
        double {{ param["id"] }} = {{ param["value"] }};
  {% endfor %}
        {{ reaction["id"] }} = {{ reaction["rhs"] }};
    }
  {% else %}
    {{ reaction["id"] }} = {{ reaction["rhs"] }};
  {% endif %}

{% endfor %}
}

// MODEL FUNCTIONS
{% for func in functions %}
inline double {{ ode_class_name }}::{{ func["id"] }}({{ func["args"] }})
{
    return {{ func["body"] }};
}

{% endfor %}

template <>
void CellwiseOdeSystemInformation<{{ ode_class_name }}>::Initialise()
{
    // STATE VARIABLES
{% for var in state_variables %}
    this->mVariableNames.push_back("{{ var['id'] }}");
    this->mVariableUnits.push_back("{{ var['units'] }}");
    this->mInitialConditions.push_back({{ var['initial_value'] }});

{% endfor %}

    // DERIVED QUANTITIES
{% for dq in derived_quantities %}
    this->mDerivedQuantityNames.push_back("{{ dq['id'] }}");
    this->mDerivedQuantityUnits.push_back("{{ dq['units'] }}");

{% endfor %}

    // PARAMETERS
{% for var in variable_parameters %}
    this->mParameterNames.push_back("{{ var['id'] }}");
    this->mParameterUnits.push_back("{{ var['units'] }}");

{% endfor %}
    this->mInitialised = true;
}

// Define {{ wrapper_class_name }} using wrappers
#include "{{ wrapper_class_name }}.hpp"
#include "{{ wrapper_class_name }}.cpp"

typedef {{ wrapper_class_name }}<{{ ode_class_name }}, {{ state_variables|length }}> {{ model_class_name }};

// Declare identifiers for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT({{ ode_class_name }})
EXPORT_TEMPLATE_CLASS2({{ wrapper_class_name }}, {{ ode_class_name }}, {{ state_variables|length }})

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER({{ model_class_name }})
