#ifndef TAN2014SBMLODESYSTEMANDSRNMODEL_HPP_
#define TAN2014SBMLODESYSTEMANDSRNMODEL_HPP_

#include <vector>

#include <boost/serialization/base_object.hpp>

#include "AbstractOdeSystem.hpp"
#include "ChasteSerialization.hpp"
#include "SbmlEventType.hpp"

class Tan2014SbmlOdeSystem : public AbstractOdeSystem
{
private:
    // (De-)Serialize Tan2014SbmlOdeSystem
    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar& BOOST_SERIALIZATION_BASE_OBJECT_NVP(AbstractOdeSystem);
    }

    // CONSTANT PARAMETERS
    const double kdegradation = 0.0163;         //
    const double Bsyn = 1.306;                  //
    const double K_n_active_k = 17.16;          //
    const double kC_k1 = 1e-05;                 //
    const double kC_k2 = 0.000647;              //
    const double kN_k1 = 0.0001;                //
    const double kN_k2 = 0.00349;               //
    const double kdiffusion_k = 39.13;          //
    const double K_c_active_k = 4.5;            //
    const double ComplexTransitThreshold = 1.0; // ComplexTransitThreshold

    // STATE VARIABLES
    double bcat_cm;    //
    double ligand_cm;  //
    double complex_cm; //
    double bcat_nu;    //
    double ligand_nu;  //
    double complex_nu; //

    // DERIVED QUANTITIES
    double drag; // drag

    // VARIABLE PARAMETERS
    double compartment;     //
    double CytosolMembrane; //
    double nucleus;         //
    double wnt_level;       // wnt_level
    double gamma;           // gamma

    // RULE-BASED PARAMETERS

    // REACTIONS
    double Bsynthesis;   //
    double kDegradation; //
    double kC;           //
    double kN;           //
    double kdiffusion;   //
    double K_c_active;   //
    double K_n_active;   //

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
    Tan2014SbmlOdeSystem(std::vector<double> stateVariables = std::vector<double>());

    /**
     * Copy constructor
     *
     * @param rOdeSystem Reference to the original instance
     */
    Tan2014SbmlOdeSystem(const Tan2014SbmlOdeSystem& rOdeSystem);

    /**
     * Destructor
     */
    ~Tan2014SbmlOdeSystem();

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
    // Provide constructor for serializing Tan2014SbmlOdeSystem
    template <class Archive>
    inline void save_construct_data(Archive& ar, const Tan2014SbmlOdeSystem* t, const unsigned int version)
    {
        // Save data required to construct instance
        const std::vector<double> state_variables = t->rGetConstStateVariables();
        ar << state_variables;
    }

    // Provide constructor for de-serializing Tan2014SbmlOdeSystem
    template <class Archive>
    inline void load_construct_data(Archive& ar, Tan2014SbmlOdeSystem* t, const unsigned int version)
    {
        // Retrieve data from archive required to construct new instance
        std::vector<double> state_variables;
        ar >> state_variables;

        // Invoke inplace constructor to initialise instance
        ::new (t) Tan2014SbmlOdeSystem(state_variables);
    }
} // namespace serialization
} // namespace

// Define SbmlSrnWrapperModel using wrappers
#include "SbmlSrnWrapperModel.cpp"
#include "SbmlSrnWrapperModel.hpp"

typedef SbmlSrnWrapperModel<Tan2014SbmlOdeSystem, 6> Tan2014SbmlSrnModel;

// Declare identifiers for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(Tan2014SbmlOdeSystem)
EXPORT_TEMPLATE_CLASS2(SbmlSrnWrapperModel, Tan2014SbmlOdeSystem, 6)

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(Tan2014SbmlSrnModel)

#endif // TAN2014SBMLODESYSTEMANDSRNMODEL_HPP_