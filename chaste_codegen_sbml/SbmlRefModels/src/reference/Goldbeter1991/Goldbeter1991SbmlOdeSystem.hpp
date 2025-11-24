#ifndef GOLDBETER_1991_SBML_ODE_SYSTEM_HPP_
#define GOLDBETER_1991_SBML_ODE_SYSTEM_HPP_

#include <vector>

#include <boost/serialization/base_object.hpp>

#include "AbstractSbmlOdeSystem.hpp"
#include "ChasteSerialization.hpp"
#include "SbmlEventType.hpp"

class Goldbeter1991SbmlOdeSystem : public AbstractSbmlOdeSystem
{
private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Save / load Goldbeter1991SbmlOdeSystem archive
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive& boost::serialization::base_object<AbstractSbmlOdeSystem>(*this);
    }

    // FIXED PARAMETERS
    const double VM1 = 3.0; // VM1
    const double VM3 = 1.0; // VM3
    const double Kc = 0.5;  // Kc

    // STATE VARIABLES
    double C; // Cyclin
    double M; // cdc_2_kinase
    double X; // Cyclin Protease

    // DERIVED QUANTITIES
    double V1; // V1
    double V3; // V3

    // VARIABLE PARAMETERS
    double cell; // cell

    // REACTIONS
    double reaction1; // creation of cyclin
    double reaction2; // default degradation of cyclin
    double reaction3; // cdc2 kinase triggered degration of cyclin
    double reaction4; // activation of cdc2 kinase
    double reaction5; // deactivation of cdc2 kinase
    double reaction6; // activation of cyclin protease
    double reaction7; // deactivation of cyclin protease

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
    Goldbeter1991SbmlOdeSystem();

    /**
     * Destructor
     */
    ~Goldbeter1991SbmlOdeSystem();

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
};

// Register the ODE system with Boost serialization
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(Goldbeter1991SbmlOdeSystem)

#endif // GOLDBETER_1991_SBML_ODE_SYSTEM_HPP_