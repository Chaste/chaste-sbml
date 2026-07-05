#ifndef {{ ode_header_guard }}
#define {{ ode_header_guard }}

#include <vector>

#include <boost/serialization/base_object.hpp>

#include "AbstractSbmlOdeSystem.hpp"
#include "ChasteSerialization.hpp"
#include "SbmlEventType.hpp"

class {{ ode_class_name }} : public AbstractSbmlOdeSystem
{
private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Save / load {{ ode_class_name }} archive
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template <class Archive>
    void serialize(Archive &archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractSbmlOdeSystem>(*this);
    }

    // PARAMETERS
{% for param in parameters %}
    double {{ param["id"] }}; // {{ param["label"] }}
{% endfor %}

    // STATE VARIABLES
{% for var in state_variables %}
    double {{ var["id"] }}; // {{ var["label"] }}
{% endfor %}

{% for var in state_variables %}
   double {{ var["derivative_id"] }};
{% endfor %}

    // DERIVED QUANTITIES
{% for dq in derived_quantities %}
{% if dq["is_conversion"] is false() and dq["is_reaction"] is false() %}
    double {{ dq["id"] }}; // {{ dq["label"] }}
{% endif %}
{% endfor %}

    // STOICHIOMETRY VARIABLES
{% for var in stoichiometry_variables %}
    double {{ var["id"] }}; // {{ var["label"] }}
{% endfor %}

    // REACTIONS
{% for reaction in reactions %}
    double {{ reaction["id"] }}; // {{ reaction["label"] }}
{% endfor %}

    /**
     * Process the events in the model.
     * 
     * @param time The current time
     * @param rY The current state variables
     * 
     * @return How close we are to the time of the next event
     */
    double ProcessModelEvents(double time, const std::vector<double>& rY) override;

    // MODEL FUNCTIONS
{% for func in functions %}
    inline double {{ func["id"] }}({{ func["args"] }});
{% endfor %}

public:
    /** 
     * Default constructor
     */
    {{ ode_class_name }}();

    /**
     * Destructor
     */
    ~{{ ode_class_name }}();
    
    /**
     * Compute the derived quantities from the given system state.
     *
     * @param time  the time at which to compute the derived quantities
     * @param rY a vector of values for the state variables
     * 
     * @return a vector of derived quantities
     */
    std::vector<double> ComputeDerivedQuantities(double time, const std::vector<double> &rY) override;

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

    void Initialise(double time = 0.0);

    std::vector<double> RunModelEquations(double time, const std::vector<double>& rStateVariables);
};

// Register the ODE system with Boost serialization
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT({{ ode_class_name }})

#endif // {{ ode_header_guard }}
