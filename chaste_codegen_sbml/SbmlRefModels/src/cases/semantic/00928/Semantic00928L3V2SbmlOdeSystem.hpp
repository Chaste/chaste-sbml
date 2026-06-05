#ifndef SEMANTIC_00928_L3_V2_SBML_ODE_SYSTEM_HPP_
#define SEMANTIC_00928_L3_V2_SBML_ODE_SYSTEM_HPP_

#include <vector>

#include <boost/serialization/base_object.hpp>

#include "AbstractSbmlOdeSystem.hpp"
#include "ChasteSerialization.hpp"
#include "SbmlEventType.hpp"

class Semantic00928L3V2SbmlOdeSystem : public AbstractSbmlOdeSystem
{
private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Save / load Semantic00928L3V2SbmlOdeSystem archive
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive& boost::serialization::base_object<AbstractSbmlOdeSystem>(*this);
    }

    // PARAMETERS
    double k1; // k1

    // STATE VARIABLES
    double S1; // S1
    double S2; // S2

    double d_S1_dt;
    double d_S2_dt;

    // DERIVED QUANTITIES
    double C; //

    // STOICHIOMETRY VARIABLES

    // REACTIONS
    double reaction1; // reaction1

    /**
     * Process the events in the model.
     *
     * @param time The current time
     * @param rY The current state variables
     *
     * @return How close we are to the time of the next event
     */
    double ProcessModelEvents(double time, const std::vector<double>& rY) override;

    /**
     * Run the assignment rules to update state.
     *
     * @param time The current time
     * @param rY The current state variables
     */
    void RunAssignmentRules(double time) override;

    /**
     * Run the initial assignments to set initial state.
     *
     * @param time The current time
     */
    void RunInitialAssignments(double time) override;

    /**
     * Run the reactions to update state.
     *
     * @param time The current time
     * @param rY The current state variables
     */
    void RunReactions(double time) override;

    /**
     * Update variable parameters from current ODE system parameter settings.
     *
     * @param time The current time
     */
    void UpdateParameters(double time) override;

    /**
     * Update state variables from the given ODE system state.
     *
     * @param time The current time
     * @param rStateVariables The state variables to use
     */
    void UpdateStateVariables(double time, const std::vector<double>& rStateVariables) override;

    // MODEL FUNCTIONS

public:
    /**
     * Default constructor
     */
    Semantic00928L3V2SbmlOdeSystem();

    /**
     * Destructor
     */
    ~Semantic00928L3V2SbmlOdeSystem();

    /**
     * Compute the derived quantities from the given system state.
     *
     * @param time  the time at which to compute the derived quantities
     * @param rY a vector of values for the state variables
     *
     * @return a vector of derived quantities
     */
    std::vector<double> ComputeDerivedQuantities(double time, const std::vector<double>& rY) override;

    /**
     * Compute the RHS of the ODE system.
     *
     * An ODE solver will call this function repeatedly to solve for y = [y1 ... yn].
     *
     * @param time the time used to evaluate the RHS.
     * @param rY an input solution vector used to evaluate the RHS.
     * @param rDY an output vector to be filled in with the resulting derivatives y' = [y1' ... yn'].
     */
    void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY) override;

    void Initialise(double time = 0.0);

    std::vector<double> RunModelEquations(double time, const std::vector<double>& rStateVariables);
};

// Register the ODE system with Boost serialization
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(Semantic00928L3V2SbmlOdeSystem)

#endif // SEMANTIC_00928_L3_V2_SBML_ODE_SYSTEM_HPP_