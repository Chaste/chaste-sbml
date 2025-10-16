#include <cmath>
#include <limits>
#include <vector>

#include "CellwiseOdeSystemInformation.hpp"
#include "SbmlEventType.hpp"
#include "SbmlMath.hpp"

#include "TysonNovak2001SbmlOdeSystem.hpp"

namespace sm = sbmlmath;

TysonNovak2001SbmlOdeSystem::TysonNovak2001SbmlOdeSystem()
        : AbstractSbmlOdeSystem(8, 1, 1)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<TysonNovak2001SbmlOdeSystem>);

    // STATE VARIABLES
    CycBt = 0.001;
    Cdc20a = 0.001;
    Cdh1 = 0.001;
    m = 0.5;
    Cdc20t = 0.001;
    IEP = 0.001;
    CKIt = 0.001;
    SK = 0.001;

    // DERIVED QUANTITIES
    CycB = 0.0;
    Trimer = 0.0;
    Mad = 0.0;
    TF = 0.0;

    // VARIABLE PARAMETERS
    cell = 1.0;

    // INITIAL ASSIGNMENTS

    // ODE SYSTEM INFORMATION
    SetDefaultInitialCondition(0, CycBt);
    SetDefaultInitialCondition(1, Cdc20a);
    SetDefaultInitialCondition(2, Cdh1);
    SetDefaultInitialCondition(3, m);
    SetDefaultInitialCondition(4, Cdc20t);
    SetDefaultInitialCondition(5, IEP);
    SetDefaultInitialCondition(6, CKIt);
    SetDefaultInitialCondition(7, SK);

    mStateVariables.push_back(CycBt);
    mStateVariables.push_back(Cdc20a);
    mStateVariables.push_back(Cdh1);
    mStateVariables.push_back(m);
    mStateVariables.push_back(Cdc20t);
    mStateVariables.push_back(IEP);
    mStateVariables.push_back(CKIt);
    mStateVariables.push_back(SK);

    mParameters.push_back(cell);

    // REACTIONS
    CycBt_synthesis = 0.0;
    CycBdegradation = 0.0;
    CycBdegradationviaCdh1 = 0.0;
    CycBtdegradationviaCdc20a = 0.0;
    Cdh1synthesis = 0.0;
    Cdh1degradation = 0.0;
    Cdc20tsynthesis = 0.0;
    Cdc20t_deg = 0.0;
    Cdc20activation = 0.0;
    Cdc20ainhibition = 0.0;
    Cdc20adegradation = 0.0;
    IEPsynthesis = 0.0;
    IEPdegradation = 0.0;
    growth = 0.0;
    CKItsynthesis = 0.0;
    CKIdegradation = 0.0;
    CKItphosphorilationviaSK = 0.0;
    eq_7 = 0.0;
    SKsynthesis = 0.0;
    SKdegradation = 0.0;

    // EVENTS
    mEventType.resize(1, SbmlEventType::UNKNOWN);

    // Uncomment lines below for events that should trigger cell division

    mEventType[0] = SbmlEventType::CELL_DIVISION; // Cell division

    mEventSatisfied.resize(1, true); // Prevent events from triggering at the start
    mEventTriggered.resize(1, false);

    mEventAdjustedParameters.resize(1, false);
    mEventAdjustedParameterValues.resize(1, 0.0);

    mEventAdjustedStateVars.resize(8, false);
    mEventAdjustedStateValues.resize(8, 0.0);

    // Run model rules to complete state initialisation
    RunModelRules(0.0, mStateVariables);
}

TysonNovak2001SbmlOdeSystem::~TysonNovak2001SbmlOdeSystem()
{
}

std::vector<double> TysonNovak2001SbmlOdeSystem::ComputeDerivedQuantities(double time, const std::vector<double>& rY)
{
    RunModelRules(time, rY);

    std::vector<double> dqs;
    dqs.push_back(CycB);
    dqs.push_back(Trimer);
    dqs.push_back(Mad);
    dqs.push_back(TF);
    return dqs;
}

void TysonNovak2001SbmlOdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    RunModelRules(time, rY);

    rDY[0] = (CycBt_synthesis - CycBdegradation - CycBdegradationviaCdh1 - CycBtdegradationviaCdc20a) / cell; // d[CycBt]/dt
    rDY[1] = (Cdc20activation - Cdc20ainhibition - Cdc20adegradation) / cell;                                 // d[Cdc20a]/dt
    rDY[2] = (Cdh1synthesis - Cdh1degradation) / cell;                                                        // d[Cdh1]/dt
    rDY[3] = growth / cell;                                                                                   // d[m]/dt
    rDY[4] = (Cdc20tsynthesis - Cdc20t_deg) / cell;                                                           // d[Cdc20t]/dt
    rDY[5] = (IEPsynthesis - IEPdegradation) / cell;                                                          // d[IEP]/dt
    rDY[6] = (CKItsynthesis - CKIdegradation - CKItphosphorilationviaSK - eq_7) / cell;                       // d[CKIt]/dt
    rDY[7] = (SKsynthesis - SKdegradation) / cell;                                                            // d[SK]/dt

    // TODO: Scale time appropriately
}

double TysonNovak2001SbmlOdeSystem::ProcessModelEvents(double time, const std::vector<double>& rY)
{
    std::fill(std::begin(mEventAdjustedParameters), std::end(mEventAdjustedParameters), false);
    std::fill(std::begin(mEventAdjustedStateVars), std::end(mEventAdjustedStateVars), false);

    double min_dist = std::numeric_limits<double>::max();

    //========================================
    // EVENT: Cell division
    //========================================
    {
        double event_dist = (0.1) - (CycB)-std::numeric_limits<double>::epsilon();

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
        if (sm::lt(CycB, 0.1))
        {
            if (!mEventSatisfied[0])
            {
                // The condition is transitioning from false -> true: trigger the event
                mEventTriggered[0] = true;
                event_dist = 0.0;
                min_dist = 0.0;

                // Adjust relevant state variables and parameters
                // m = m / 2.0
                mEventAdjustedStateVars[3] = true;
                mEventAdjustedStateValues[3] = m / 2.0;
            }
            mEventSatisfied[0] = true;
        }
        else
        {
            mEventSatisfied[0] = false;
            mEventTriggered[0] = false;
        }
    }

    return min_dist; // Distance to closest event
}

void TysonNovak2001SbmlOdeSystem::RunModelRules(double time, const std::vector<double>& rY)
{
    // STATE VARIABLES
    CycBt = rY[0];
    Cdc20a = rY[1];
    Cdh1 = rY[2];
    m = rY[3];
    Cdc20t = rY[4];
    IEP = rY[5];
    CKIt = rY[6];
    SK = rY[7];

    // VARIABLE PARAMETERS
    cell = GetParameter(0);

    // ASSIGNMENT RULES
    CycB = CycBt - 2.0 * CycBt * CKIt / (CycBt + CKIt + 1.0 / Keq + std::pow(std::pow(CycBt + CKIt + 1.0 / Keq, 2.0) - 4.0 * CycBt * CKIt, 1.0 / 2.0));
    Trimer = 2.0 * CycBt * CKIt / (CycBt + CKIt + 1.0 / Keq + std::pow(std::pow(CycBt + CKIt + 1.0 / Keq, 2.0) - 4.0 * CycBt * CKIt, 1.0 / 2.0));
    TF = GK(k15p * m + k15pp * SK, k16p + k16pp * m * CycB, J15, J16);
    Mad = 1.0;

    // REACTIONS
    // CycBt synthesis
    CycBt_synthesis = k1;

    // CycBt degradation
    CycBdegradation = k2p * CycBt;

    // CycBt degradation via Cdh1
    CycBdegradationviaCdh1 = k2pp * Cdh1 * CycBt;

    // CycBt degradation via Cdc20a
    CycBtdegradationviaCdc20a = k2ppp * Cdc20a * CycBt;

    // Cdh1 synthesis
    Cdh1synthesis = (k3p + k3pp * Cdc20a) * (1.0 - Cdh1) / (J3 + 1.0 - Cdh1);

    // Cdh1 degradation
    Cdh1degradation = (k4p * SK * Cdh1 + k4 * m * CycB * Cdh1) / (J4 + Cdh1);

    // Cdc20t synthesis
    Cdc20tsynthesis = k5p + k5pp * std::pow(CycB * m / J5, n) / (1.0 + std::pow(CycB * m / J5, n));

    // Cdc20t degradation
    Cdc20t_deg = k6 * Cdc20t;

    // Cdc20 activation
    Cdc20activation = k7 * IEP * (Cdc20t - Cdc20a) / (J7 + Cdc20t - Cdc20a);

    // Cdc20a inhibition
    Cdc20ainhibition = k8 * Mad * Cdc20a / (J8 + Cdc20a);

    // Cdc20a degradation
    Cdc20adegradation = k6 * Cdc20a;

    // IEP synthesis
    IEPsynthesis = k9 * m * CycB * (1.0 - IEP);

    // IEP degradation
    IEPdegradation = k10 * IEP;

    // growth
    growth = mu * m * (1.0 - m / mmax);

    // CKIt synthesis
    CKItsynthesis = k11;

    // CKIt degradation
    CKIdegradation = k12p * CKIt;

    // CKIt phosphorilation via SK
    CKItphosphorilationviaSK = k12pp * SK * CKIt;

    // CKIt Trimer sequestred
    eq_7 = k12ppp * m * CycB * CKIt;

    // SK synthesis
    SKsynthesis = k13 * TF;

    // SK degradation
    SKdegradation = k14 * SK;
}

// MODEL FUNCTIONS
inline double TysonNovak2001SbmlOdeSystem::GK(double A1, double A2, double A3, double A4)
{
    return 2.0 * A4 * A1 / (A2 - A1 + A3 * A2 + A4 * A1 + sm::root(2.0, std::pow(A2 - A1 + A3 * A2 + A4 * A1, 2.0) - 4.0 * (A2 - A1) * A4 * A1));
}

template <>
void CellwiseOdeSystemInformation<TysonNovak2001SbmlOdeSystem>::Initialise()
{
    // STATE VARIABLES
    this->mVariableNames.push_back("CycBt");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.001);

    this->mVariableNames.push_back("Cdc20a");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.001);

    this->mVariableNames.push_back("Cdh1");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.001);

    this->mVariableNames.push_back("m");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.5);

    this->mVariableNames.push_back("Cdc20t");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.001);

    this->mVariableNames.push_back("IEP");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.001);

    this->mVariableNames.push_back("CKIt");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.001);

    this->mVariableNames.push_back("SK");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.001);

    // DERIVED QUANTITIES
    this->mDerivedQuantityNames.push_back("CycB");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Trimer");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Mad");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("TF");
    this->mDerivedQuantityUnits.push_back("non-dim");

    // PARAMETERS
    this->mParameterNames.push_back("cell");
    this->mParameterUnits.push_back("non-dim");

    this->mInitialised = true;
}

// Register the ODE system with Boost serialization
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(TysonNovak2001SbmlOdeSystem)