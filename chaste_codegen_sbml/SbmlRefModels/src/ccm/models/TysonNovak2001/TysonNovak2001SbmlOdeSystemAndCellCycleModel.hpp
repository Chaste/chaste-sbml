#ifndef TYSONNOVAK2001SBMLODESYSTEMANDCELLCYCLEMODEL_HPP_
#define TYSONNOVAK2001SBMLODESYSTEMANDCELLCYCLEMODEL_HPP_

#include <vector>

#include <boost/serialization/base_object.hpp>

#include "AbstractOdeSystem.hpp"
#include "ChasteSerialization.hpp"
#include "SbmlEventType.hpp"

class TysonNovak2001SbmlOdeSystem : public AbstractOdeSystem
{
private:
    // (De-)Serialize TysonNovak2001SbmlOdeSystem
    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar& BOOST_SERIALIZATION_BASE_OBJECT_NVP(AbstractOdeSystem);
    }

    // CONSTANT PARAMETERS
    const double k1 = 0.04;      // k1
    const double k2p = 0.04;     // k2p
    const double k2pp = 1.0;     // k2pp
    const double k2ppp = 1.0;    // k2ppp
    const double k3p = 1.0;      // k3p
    const double k3pp = 10.0;    // k3pp
    const double J3 = 0.04;      // J3
    const double k4 = 35.0;      // k4
    const double k5p = 0.005;    // k5p
    const double k5pp = 0.2;     // k5pp
    const double J5 = 0.3;       // J5
    const double k6 = 0.1;       // k6
    const double n = 4.0;        // n
    const double k7 = 1.0;       // k7
    const double J7 = 0.001;     // J7
    const double k8 = 0.5;       // k8
    const double J8 = 0.001;     // J8
    const double k9 = 0.1;       // k9
    const double k10 = 0.02;     // k10
    const double mu = 0.005;     // mu
    const double k11 = 1.0;      // k11
    const double k12p = 0.2;     // k12p
    const double k12pp = 50.0;   // k12pp
    const double mmax = 10.0;    // mmax
    const double k12ppp = 100.0; // k12ppp
    const double Keq = 1000.0;   // Keq
    const double k13 = 1.0;      // k13
    const double k14 = 1.0;      // k14
    const double k15p = 1.5;     // k15p
    const double k15pp = 0.05;   // k15pp
    const double k16p = 1.0;     // k16p
    const double k16pp = 3.0;    // k16pp
    const double J15 = 0.01;     // J15
    const double J16 = 0.01;     // J16
    const double k4p = 2.0;      // k4p
    const double J4 = 0.04;      // J4

    // STATE VARIABLES
    double CycBt;  // CycBt
    double Cdc20a; // Cdc20a
    double Cdh1;   // Cdh1
    double m;      // m
    double Cdc20t; // Cdc20t
    double IEP;    // IEP
    double CKIt;   // CKIt
    double SK;     // SK

    // DERIVED QUANTITIES
    double CycB;   // CycB
    double Trimer; // Trimer
    double Mad;    // Mad

    // VARIABLE PARAMETERS
    double cell; // cell

    // RULE-BASED PARAMETERS
    double TF; // TF

    // REACTIONS
    double CycBt_synthesis;           // CycBt synthesis
    double CycBdegradation;           // CycBt degradation
    double CycBdegradationviaCdh1;    // CycBt degradation via Cdh1
    double CycBtdegradationviaCdc20a; // CycBt degradation via Cdc20a
    double Cdh1synthesis;             // Cdh1 synthesis
    double Cdh1degradation;           // Cdh1 degradation
    double Cdc20tsynthesis;           // Cdc20t synthesis
    double Cdc20t_deg;                // Cdc20t degradation
    double Cdc20activation;           // Cdc20 activation
    double Cdc20ainhibition;          // Cdc20a inhibition
    double Cdc20adegradation;         // Cdc20a degradation
    double IEPsynthesis;              // IEP synthesis
    double IEPdegradation;            // IEP degradation
    double growth;                    // growth
    double CKItsynthesis;             // CKIt synthesis
    double CKIdegradation;            // CKIt degradation
    double CKItphosphorilationviaSK;  // CKIt phosphorilation via SK
    double eq_7;                      // CKIt Trimer sequestred
    double SKsynthesis;               // SK synthesis
    double SKdegradation;             // SK degradation

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
    TysonNovak2001SbmlOdeSystem(std::vector<double> stateVariables = std::vector<double>());

    /**
     * Copy constructor
     *
     * @param rOdeSystem Reference to the original instance
     */
    TysonNovak2001SbmlOdeSystem(const TysonNovak2001SbmlOdeSystem& rOdeSystem);

    /**
     * Destructor
     */
    ~TysonNovak2001SbmlOdeSystem();

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
    inline double GK(double A1, double A2, double A3, double A4);
};

namespace
{
namespace serialization
{
    // Provide constructor for serializing TysonNovak2001SbmlOdeSystem
    template <class Archive>
    inline void save_construct_data(Archive& ar, const TysonNovak2001SbmlOdeSystem* t, const unsigned int version)
    {
        // Save data required to construct instance
        const std::vector<double> state_variables = t->rGetConstStateVariables();
        ar << state_variables;
    }

    // Provide constructor for de-serializing TysonNovak2001SbmlOdeSystem
    template <class Archive>
    inline void load_construct_data(Archive& ar, TysonNovak2001SbmlOdeSystem* t, const unsigned int version)
    {
        // Retrieve data from archive required to construct new instance
        std::vector<double> state_variables;
        ar >> state_variables;

        // Invoke inplace constructor to initialise instance
        ::new (t) TysonNovak2001SbmlOdeSystem(state_variables);
    }
} // namespace serialization
} // namespace

// Define SbmlCellCycleWrapperModel using wrappers
#include "SbmlCellCycleWrapperModel.cpp"
#include "SbmlCellCycleWrapperModel.hpp"

typedef SbmlCellCycleWrapperModel<TysonNovak2001SbmlOdeSystem, 8> TysonNovak2001SbmlCellCycleModel;

// Declare identifiers for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(TysonNovak2001SbmlOdeSystem)
EXPORT_TEMPLATE_CLASS2(SbmlCellCycleWrapperModel, TysonNovak2001SbmlOdeSystem, 8)

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(TysonNovak2001SbmlCellCycleModel)

#endif // TYSONNOVAK2001SBMLODESYSTEMANDCELLCYCLEMODEL_HPP_