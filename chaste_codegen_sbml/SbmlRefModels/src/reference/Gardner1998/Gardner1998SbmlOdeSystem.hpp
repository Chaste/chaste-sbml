#ifndef GARDNER_1998_SBML_ODE_SYSTEM_HPP_
#define GARDNER_1998_SBML_ODE_SYSTEM_HPP_

#include <vector>

#include <boost/serialization/base_object.hpp>

#include "AbstractSbmlOdeSystem.hpp"
#include "ChasteSerialization.hpp"
#include "SbmlEventType.hpp"

class Gardner1998SbmlOdeSystem : public AbstractSbmlOdeSystem
{
private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Save / load Gardner1998SbmlOdeSystem archive
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
    const double K6 = 0.3;   // K6
    const double V1p = 0.75; // V1p
    const double V3p = 0.3;  // V3p

    // STATE VARIABLES
    double C; // Convert C amount to concentration
    double X; // Convert X amount to concentration
    double M; // Convert M amount to concentration
    double Y; // Convert Y amount to concentration
    double Z; // Convert Z amount to concentration

    // DERIVED QUANTITIES
    double V1; // V1
    double V3; // V3

    // VARIABLE PARAMETERS
    double Cell; // Cell

    // REACTIONS
    double reaction1;  // creation of cyclin
    double reaction2;  // cdc2 kinase triggered degration of cyclin
    double reaction3;  // default degradation of cyclin
    double reaction4;  // activation of cdc2 kinase
    double reaction5;  // deactivation of cdc2 kinase
    double reaction6;  // activation of cyclin protease
    double reaction7;  // deactivation of cyclin protease
    double reaction8;  // reaction8
    double reaction9;  // reaction9
    double reaction10; // desinhibition of cyclin
    double reaction11; // degradation of inhibited cyclin
    double reaction12; // creation of cyclin inhibitor
    double reaction13; // degradation of cyclin inhibitor

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
    Gardner1998SbmlOdeSystem();

    /**
     * Destructor
     */
    ~Gardner1998SbmlOdeSystem();

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
CHASTE_CLASS_EXPORT(Gardner1998SbmlOdeSystem)

#endif // GARDNER_1998_SBML_ODE_SYSTEM_HPP_