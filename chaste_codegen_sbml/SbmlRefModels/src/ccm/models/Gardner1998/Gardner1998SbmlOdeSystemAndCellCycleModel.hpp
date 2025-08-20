#ifndef GARDNER1998SBMLODESYSTEMANDCELLCYCLEMODEL_HPP_
#define GARDNER1998SBMLODESYSTEMANDCELLCYCLEMODEL_HPP_

#include <vector>

#include <boost/serialization/base_object.hpp>

#include "AbstractOdeSystem.hpp"
#include "ChasteSerialization.hpp"
#include "SbmlEventType.hpp"

class Gardner1998SbmlOdeSystem : public AbstractOdeSystem
{
private:
    // (De-)Serialize Gardner1998SbmlOdeSystem
    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar& BOOST_SERIALIZATION_BASE_OBJECT_NVP(AbstractOdeSystem);
    }

    // CONSTANT PARAMETERS
    const double K6 = 0.3;   // K6
    const double V1p = 0.75; // V1p
    const double V3p = 0.3;  // V3p

    // STATE VARIABLES
    double C; // cyclin
    double X; // protease
    double M; // cdc2k
    double Y; // cyclin inhibitor
    double Z; // complex inhibitor-cyclin

    // DERIVED QUANTITIES

    // VARIABLE PARAMETERS
    double Cell; // Cell

    // RULE-BASED PARAMETERS
    double V1; // V1
    double V3; // V3

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
    Gardner1998SbmlOdeSystem(std::vector<double> stateVariables = std::vector<double>());

    /**
     * Copy constructor
     *
     * @param rOdeSystem Reference to the original instance
     */
    Gardner1998SbmlOdeSystem(const Gardner1998SbmlOdeSystem& rOdeSystem);

    /**
     * Destructor
     */
    ~Gardner1998SbmlOdeSystem();

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
    std::vector<double> ComputeDerivedQuantities(double time, const std::vector<double>& rY);

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
};

namespace
{
namespace serialization
{
    // Provide constructor for serializing Gardner1998SbmlOdeSystem
    template <class Archive>
    inline void save_construct_data(Archive& ar, const Gardner1998SbmlOdeSystem* t, const unsigned int version)
    {
        // Save data required to construct instance
        const std::vector<double> state_variables = t->rGetConstStateVariables();
        ar << state_variables;
    }

    // Provide constructor for de-serializing Gardner1998SbmlOdeSystem
    template <class Archive>
    inline void load_construct_data(Archive& ar, Gardner1998SbmlOdeSystem* t, const unsigned int version)
    {
        // Retrieve data from archive required to construct new instance
        std::vector<double> state_variables;
        ar >> state_variables;

        // Invoke inplace constructor to initialise instance
        ::new (t) Gardner1998SbmlOdeSystem(state_variables);
    }
} // namespace serialization
} // namespace

// Define SbmlCellCycleWrapperModel using wrappers
#include "SbmlCellCycleWrapperModel.cpp"
#include "SbmlCellCycleWrapperModel.hpp"

typedef SbmlCellCycleWrapperModel<Gardner1998SbmlOdeSystem, 5> Gardner1998SbmlCellCycleModel;

// Declare identifiers for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(Gardner1998SbmlOdeSystem)
EXPORT_TEMPLATE_CLASS2(SbmlCellCycleWrapperModel, Gardner1998SbmlOdeSystem, 5)

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(Gardner1998SbmlCellCycleModel)

#endif // GARDNER1998SBMLODESYSTEMANDCELLCYCLEMODEL_HPP_