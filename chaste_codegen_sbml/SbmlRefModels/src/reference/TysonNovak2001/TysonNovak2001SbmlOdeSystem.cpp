#include <cmath>
#include <limits>
#include <vector>

#include "CellwiseOdeSystemInformation.hpp"
#include "SbmlEventType.hpp"
#include "SbmlMath.hpp"

#include "TysonNovak2001SbmlOdeSystem.hpp"

namespace sm = sbmlmath;

TysonNovak2001SbmlOdeSystem::TysonNovak2001SbmlOdeSystem()
        : AbstractSbmlOdeSystem(8, 36, 1)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<TysonNovak2001SbmlOdeSystem>);

    Initialise();

    // EVENTS
    mEventType.resize(1, SbmlEventType::UNKNOWN);

    // Uncomment lines below for events that should trigger cell division

    mEventType[0] = SbmlEventType::CELL_DIVISION; // Cell division

    mEventSatisfied = { true }; // From SBML trigger initialValue
    mEventTriggered.resize(1, false);

    mEventAdjustedParameters.resize(36, false);
    mEventAdjustedParameterValues.resize(36, 0.0);

    mEventAdjustedStateVars.resize(8, false);
    mEventAdjustedStateValues.resize(8, 0.0);
}

TysonNovak2001SbmlOdeSystem::~TysonNovak2001SbmlOdeSystem()
{
}

std::vector<double> TysonNovak2001SbmlOdeSystem::ComputeDerivedQuantities(double time, const std::vector<double>& rY)
{
    std::vector<double> dqs;
    RunModelEquations(time, rY);

    // AMOUNT / CONCENTRATION CONVERSIONS
    double conc__CycBt = CycBt / cell;   //
    double conc__CycB = CycB / cell;     //
    double conc__Cdc20a = Cdc20a / cell; //
    double conc__Trimer = Trimer / cell; //
    double conc__Cdh1 = Cdh1 / cell;     //
    double conc__m = m / cell;           //
    double conc__Cdc20t = Cdc20t / cell; //
    double conc__IEP = IEP / cell;       //
    double conc__Mad = Mad / cell;       //
    double conc__CKIt = CKIt / cell;     //
    double conc__SK = SK / cell;         //

    dqs.push_back(cell);
    dqs.push_back(conc__CycBt);
    dqs.push_back(CycB);
    dqs.push_back(conc__CycB);
    dqs.push_back(conc__Cdc20a);
    dqs.push_back(Trimer);
    dqs.push_back(conc__Trimer);
    dqs.push_back(conc__Cdh1);
    dqs.push_back(conc__m);
    dqs.push_back(conc__Cdc20t);
    dqs.push_back(conc__IEP);
    dqs.push_back(Mad);
    dqs.push_back(conc__Mad);
    dqs.push_back(conc__CKIt);
    dqs.push_back(conc__SK);
    dqs.push_back(TF);

    return dqs;
}

void TysonNovak2001SbmlOdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    std::vector<double> derivatives = RunModelEquations(time, rY);
    for (unsigned i = 0; i < rDY.size(); ++i)
    {
        rDY[i] = derivatives[i];
    }

    // TODO: Scale time appropriately
}

void TysonNovak2001SbmlOdeSystem::Initialise(double time)
{
    cell = 1.0;                                                                                                                                               //
    CycBt = 0.001;                                                                                                                                            //
    Cdc20a = 0.001;                                                                                                                                           //
    Cdh1 = 0.001;                                                                                                                                             //
    m = 0.5;                                                                                                                                                  //
    Cdc20t = 0.001;                                                                                                                                           //
    IEP = 0.001;                                                                                                                                              //
    CKIt = 0.001;                                                                                                                                             //
    SK = 0.001;                                                                                                                                               //
    k1 = 0.04;                                                                                                                                                //
    k2p = 0.04;                                                                                                                                               //
    k2pp = 1.0;                                                                                                                                               //
    k2ppp = 1.0;                                                                                                                                              //
    k3p = 1.0;                                                                                                                                                //
    k3pp = 10.0;                                                                                                                                              //
    J3 = 0.04;                                                                                                                                                //
    k4 = 35.0;                                                                                                                                                //
    k5p = 0.005;                                                                                                                                              //
    k5pp = 0.2;                                                                                                                                               //
    J5 = 0.3;                                                                                                                                                 //
    k6 = 0.1;                                                                                                                                                 //
    n = 4.0;                                                                                                                                                  //
    k7 = 1.0;                                                                                                                                                 //
    J7 = 0.001;                                                                                                                                               //
    k8 = 0.5;                                                                                                                                                 //
    J8 = 0.001;                                                                                                                                               //
    k9 = 0.1;                                                                                                                                                 //
    k10 = 0.02;                                                                                                                                               //
    mu = 0.005;                                                                                                                                               //
    k11 = 1.0;                                                                                                                                                //
    k12p = 0.2;                                                                                                                                               //
    k12pp = 50.0;                                                                                                                                             //
    mmax = 10.0;                                                                                                                                              //
    k12ppp = 100.0;                                                                                                                                           //
    Keq = 1000.0;                                                                                                                                             //
    k13 = 1.0;                                                                                                                                                //
    k14 = 1.0;                                                                                                                                                //
    k15p = 1.5;                                                                                                                                               //
    k15pp = 0.05;                                                                                                                                             //
    k16p = 1.0;                                                                                                                                               //
    k16pp = 3.0;                                                                                                                                              //
    J15 = 0.01;                                                                                                                                               //
    J16 = 0.01;                                                                                                                                               //
    k4p = 2.0;                                                                                                                                                //
    J4 = 0.04;                                                                                                                                                //
    CycB = CycBt - 2.0 * CycBt * CKIt / (CycBt + CKIt + 1.0 / Keq + std::pow((std::pow((CycBt + CKIt + 1.0 / Keq), 2.0) - 4.0 * CycBt * CKIt), (1.0 / 2.0))); //
    Trimer = 2.0 * CycBt * CKIt / (CycBt + CKIt + 1.0 / Keq + std::pow((std::pow((CycBt + CKIt + 1.0 / Keq), 2.0) - 4.0 * CycBt * CKIt), (1.0 / 2.0)));       //
    Mad = 1.0;                                                                                                                                                //
    TF = GK(k15p * m + k15pp * SK, k16p + k16pp * m * CycB, J15, J16);                                                                                        //
    CycBt_synthesis = k1;                                                                                                                                     //
    CycBdegradation = k2p * CycBt;                                                                                                                            //
    CycBdegradationviaCdh1 = k2pp * Cdh1 * CycBt;                                                                                                             //
    CycBtdegradationviaCdc20a = k2ppp * Cdc20a * CycBt;                                                                                                       //
    Cdh1synthesis = (k3p + k3pp * Cdc20a) * (1.0 - Cdh1) / (J3 + 1.0 - Cdh1);                                                                                 //
    Cdh1degradation = (k4p * SK * Cdh1 + k4 * m * CycB * Cdh1) / (J4 + Cdh1);                                                                                 //
    Cdc20tsynthesis = k5p + k5pp * std::pow((CycB * m / J5), n) / (1.0 + std::pow((CycB * m / J5), n));                                                       //
    Cdc20t_deg = k6 * Cdc20t;                                                                                                                                 //
    Cdc20activation = k7 * IEP * (Cdc20t - Cdc20a) / (J7 + Cdc20t - Cdc20a);                                                                                  //
    Cdc20ainhibition = k8 * Mad * Cdc20a / (J8 + Cdc20a);                                                                                                     //
    Cdc20adegradation = k6 * Cdc20a;                                                                                                                          //
    IEPsynthesis = k9 * m * CycB * (1.0 - IEP);                                                                                                               //
    IEPdegradation = k10 * IEP;                                                                                                                               //
    growth = mu * m * (1.0 - m / mmax);                                                                                                                       //
    CKItsynthesis = k11;                                                                                                                                      //
    CKIdegradation = k12p * CKIt;                                                                                                                             //
    CKItphosphorilationviaSK = k12pp * SK * CKIt;                                                                                                             //
    eq_7 = k12ppp * m * CycB * CKIt;                                                                                                                          //
    SKsynthesis = k13 * TF;                                                                                                                                   //
    SKdegradation = k14 * SK;                                                                                                                                 //
    d_CycBt_dt = CycBt_synthesis - CycBdegradation - CycBdegradationviaCdh1 - CycBtdegradationviaCdc20a;                                                      //
    d_Cdc20a_dt = Cdc20activation - Cdc20ainhibition - Cdc20adegradation;                                                                                     //
    d_Cdh1_dt = Cdh1synthesis - Cdh1degradation;                                                                                                              //
    d_m_dt = growth;                                                                                                                                          //
    d_Cdc20t_dt = Cdc20tsynthesis - Cdc20t_deg;                                                                                                               //
    d_IEP_dt = IEPsynthesis - IEPdegradation;                                                                                                                 //
    d_CKIt_dt = CKItsynthesis - CKIdegradation - CKItphosphorilationviaSK - eq_7;                                                                             //
    d_SK_dt = SKsynthesis - SKdegradation;                                                                                                                    //

    mStateVariables.push_back(CycBt);
    mStateVariables.push_back(Cdc20a);
    mStateVariables.push_back(Cdh1);
    mStateVariables.push_back(m);
    mStateVariables.push_back(Cdc20t);
    mStateVariables.push_back(IEP);
    mStateVariables.push_back(CKIt);
    mStateVariables.push_back(SK);

    SetDefaultInitialCondition(0, CycBt);
    SetDefaultInitialCondition(1, Cdc20a);
    SetDefaultInitialCondition(2, Cdh1);
    SetDefaultInitialCondition(3, m);
    SetDefaultInitialCondition(4, Cdc20t);
    SetDefaultInitialCondition(5, IEP);
    SetDefaultInitialCondition(6, CKIt);
    SetDefaultInitialCondition(7, SK);

    mParameters.push_back(k1);
    mParameters.push_back(k2p);
    mParameters.push_back(k2pp);
    mParameters.push_back(k2ppp);
    mParameters.push_back(k3p);
    mParameters.push_back(k3pp);
    mParameters.push_back(J3);
    mParameters.push_back(k4);
    mParameters.push_back(k5p);
    mParameters.push_back(k5pp);
    mParameters.push_back(J5);
    mParameters.push_back(k6);
    mParameters.push_back(n);
    mParameters.push_back(k7);
    mParameters.push_back(J7);
    mParameters.push_back(k8);
    mParameters.push_back(J8);
    mParameters.push_back(k9);
    mParameters.push_back(k10);
    mParameters.push_back(mu);
    mParameters.push_back(k11);
    mParameters.push_back(k12p);
    mParameters.push_back(k12pp);
    mParameters.push_back(mmax);
    mParameters.push_back(k12ppp);
    mParameters.push_back(Keq);
    mParameters.push_back(k13);
    mParameters.push_back(k14);
    mParameters.push_back(k15p);
    mParameters.push_back(k15pp);
    mParameters.push_back(k16p);
    mParameters.push_back(k16pp);
    mParameters.push_back(J15);
    mParameters.push_back(J16);
    mParameters.push_back(k4p);
    mParameters.push_back(J4);
}

double TysonNovak2001SbmlOdeSystem::ProcessModelEvents(double time, const std::vector<double>& rY)
{
    // Ensure all member variables (state vars, parameters, derived quantities) reflect
    // the rY passed in. Without this, event triggers and assignments would use stale
    // values from the last EvaluateYDerivatives call, which may differ from rY when
    // called from CalculateRootFunction or CalculateStoppingEvent with a different state.
    RunModelEquations(time, rY);

    // Do NOT clear mEventAdjustedStateVars/Parameters here. Once set by an event fire,
    // they must persist across all CVODE bisection calls until AdjustParameters applies
    // them. Clearing here would erase the stored assignment when a later bisection call
    // lands in the clamped state (mEventSatisfied=true), causing the halving to be lost.
    // CalculateStoppingEvent (BackwardEuler path) clears these itself before calling.

    // Root function for CVODE: the maximum signed event distance, where each distance is
    // positive exactly when its event's trigger condition holds. Taking the MAXIMUM (not the
    // minimum absolute value) means the combined function crosses zero the moment ANY event
    // becomes triggered, and cannot be masked by another event that happens to sit just below
    // its own boundary (a small negative distance). A min-abs combination misses an event
    // whose rising edge coincides with another event re-arming near its threshold.
    double max_dist = -std::numeric_limits<double>::max();

    //========================================
    // EVENT: Cell division
    //========================================
    {
        double event_dist = (0.1) - (CycB)-std::numeric_limits<double>::epsilon();

        // Suppress an event whose trigger was already active when this Solve segment started
        // (a carried-over trigger) by forcing a large negative distance, so CVODE reports no
        // spurious root at the initial condition. mEventClampActive is frozen at segment start
        // (CalculateStoppingEvent) and cleared below the instant the trigger first goes false.
        // Using this monotonic per-segment flag rather than the live, in-step-mutated
        // mEventSatisfied keeps the root function stable across CVODE's root bracketing, so an
        // event localizes at its true crossing instead of the integration step endpoint.
        if (mEventClampActive[0] && (CycB < 0.1))
        {
            event_dist = -(std::abs(event_dist) + 1.0);
        }

        // Update max_dist (closest event to triggering)
        if (event_dist > max_dist)
        {
            max_dist = event_dist;
        }

        // Process the event
        if (CycB < 0.1)
        {
            if (!mEventSatisfied[0])
            {
                // The condition is transitioning from false -> true: trigger the event
                mEventTriggered[0] = true;

                // Adjust relevant state variables and parameters
                // m = m / 2.0
                mEventAdjustedStateVars[3] = true;
                mEventAdjustedStateValues[3] = m / 2.0;
            }
            mEventSatisfied[0] = true;
        }
        else if (!mEventTriggered[0])
        {
            // Trigger is false and the event has not fired in this segment, so it (re-)arms:
            // clear the satisfied latch and the clamp (the clamp permanently, monotonically,
            // so it stays stable across CVODE's in-step root bracketing and the next rising
            // edge is detected). Once the event HAS fired this segment we leave these sticky,
            // so a later root-bracketing evaluation that lands on the trigger-false side cannot
            // undo the fire and leave the event spuriously unsatisfied (which would re-fire it
            // at the next segment's initial condition).
            mEventSatisfied[0] = false;
            mEventClampActive[0] = false;
        }
    }

    return max_dist; // Signed distance of the event closest to triggering
}

// ASSIGNMENT RULES
void TysonNovak2001SbmlOdeSystem::RunAssignmentRules(double time)
{
}

// INITIAL ASSIGNMENTS
void TysonNovak2001SbmlOdeSystem::RunInitialAssignments(double time)
{
}

std::vector<double> TysonNovak2001SbmlOdeSystem::RunModelEquations(double time, const std::vector<double>& rStateVariables)
{
    CycBt = rStateVariables[0];
    Cdc20a = rStateVariables[1];
    Cdh1 = rStateVariables[2];
    m = rStateVariables[3];
    Cdc20t = rStateVariables[4];
    IEP = rStateVariables[5];
    CKIt = rStateVariables[6];
    SK = rStateVariables[7];

    k1 = GetParameter(0);
    k2p = GetParameter(1);
    k2pp = GetParameter(2);
    k2ppp = GetParameter(3);
    k3p = GetParameter(4);
    k3pp = GetParameter(5);
    J3 = GetParameter(6);
    k4 = GetParameter(7);
    k5p = GetParameter(8);
    k5pp = GetParameter(9);
    J5 = GetParameter(10);
    k6 = GetParameter(11);
    n = GetParameter(12);
    k7 = GetParameter(13);
    J7 = GetParameter(14);
    k8 = GetParameter(15);
    J8 = GetParameter(16);
    k9 = GetParameter(17);
    k10 = GetParameter(18);
    mu = GetParameter(19);
    k11 = GetParameter(20);
    k12p = GetParameter(21);
    k12pp = GetParameter(22);
    mmax = GetParameter(23);
    k12ppp = GetParameter(24);
    Keq = GetParameter(25);
    k13 = GetParameter(26);
    k14 = GetParameter(27);
    k15p = GetParameter(28);
    k15pp = GetParameter(29);
    k16p = GetParameter(30);
    k16pp = GetParameter(31);
    J15 = GetParameter(32);
    J16 = GetParameter(33);
    k4p = GetParameter(34);
    J4 = GetParameter(35);

    CycB = CycBt - 2.0 * CycBt * CKIt / (CycBt + CKIt + 1.0 / Keq + std::pow((std::pow((CycBt + CKIt + 1.0 / Keq), 2.0) - 4.0 * CycBt * CKIt), (1.0 / 2.0))); //
    Trimer = 2.0 * CycBt * CKIt / (CycBt + CKIt + 1.0 / Keq + std::pow((std::pow((CycBt + CKIt + 1.0 / Keq), 2.0) - 4.0 * CycBt * CKIt), (1.0 / 2.0)));       //
    Mad = 1.0;                                                                                                                                                //
    TF = GK(k15p * m + k15pp * SK, k16p + k16pp * m * CycB, J15, J16);                                                                                        //
    CycBt_synthesis = k1;                                                                                                                                     //
    CycBdegradation = k2p * CycBt;                                                                                                                            //
    CycBdegradationviaCdh1 = k2pp * Cdh1 * CycBt;                                                                                                             //
    CycBtdegradationviaCdc20a = k2ppp * Cdc20a * CycBt;                                                                                                       //
    Cdh1synthesis = (k3p + k3pp * Cdc20a) * (1.0 - Cdh1) / (J3 + 1.0 - Cdh1);                                                                                 //
    Cdh1degradation = (k4p * SK * Cdh1 + k4 * m * CycB * Cdh1) / (J4 + Cdh1);                                                                                 //
    Cdc20tsynthesis = k5p + k5pp * std::pow((CycB * m / J5), n) / (1.0 + std::pow((CycB * m / J5), n));                                                       //
    Cdc20t_deg = k6 * Cdc20t;                                                                                                                                 //
    Cdc20activation = k7 * IEP * (Cdc20t - Cdc20a) / (J7 + Cdc20t - Cdc20a);                                                                                  //
    Cdc20ainhibition = k8 * Mad * Cdc20a / (J8 + Cdc20a);                                                                                                     //
    Cdc20adegradation = k6 * Cdc20a;                                                                                                                          //
    IEPsynthesis = k9 * m * CycB * (1.0 - IEP);                                                                                                               //
    IEPdegradation = k10 * IEP;                                                                                                                               //
    growth = mu * m * (1.0 - m / mmax);                                                                                                                       //
    CKItsynthesis = k11;                                                                                                                                      //
    CKIdegradation = k12p * CKIt;                                                                                                                             //
    CKItphosphorilationviaSK = k12pp * SK * CKIt;                                                                                                             //
    eq_7 = k12ppp * m * CycB * CKIt;                                                                                                                          //
    SKsynthesis = k13 * TF;                                                                                                                                   //
    SKdegradation = k14 * SK;                                                                                                                                 //
    d_CycBt_dt = CycBt_synthesis - CycBdegradation - CycBdegradationviaCdh1 - CycBtdegradationviaCdc20a;                                                      //
    d_Cdc20a_dt = Cdc20activation - Cdc20ainhibition - Cdc20adegradation;                                                                                     //
    d_Cdh1_dt = Cdh1synthesis - Cdh1degradation;                                                                                                              //
    d_m_dt = growth;                                                                                                                                          //
    d_Cdc20t_dt = Cdc20tsynthesis - Cdc20t_deg;                                                                                                               //
    d_IEP_dt = IEPsynthesis - IEPdegradation;                                                                                                                 //
    d_CKIt_dt = CKItsynthesis - CKIdegradation - CKItphosphorilationviaSK - eq_7;                                                                             //
    d_SK_dt = SKsynthesis - SKdegradation;                                                                                                                    //

    std::vector<double> derivatives(8);
    derivatives[0] = d_CycBt_dt;
    derivatives[1] = d_Cdc20a_dt;
    derivatives[2] = d_Cdh1_dt;
    derivatives[3] = d_m_dt;
    derivatives[4] = d_Cdc20t_dt;
    derivatives[5] = d_IEP_dt;
    derivatives[6] = d_CKIt_dt;
    derivatives[7] = d_SK_dt;
    return derivatives;
}

// REACTIONS
void TysonNovak2001SbmlOdeSystem::RunReactions(double time)
{
}

// VARIABLE PARAMETERS
void TysonNovak2001SbmlOdeSystem::UpdateParameters(double time)
{
}

// STATE VARIABLES
void TysonNovak2001SbmlOdeSystem::UpdateStateVariables(double time, const std::vector<double>& rStateVariables)
{
    CycBt = rStateVariables[0];
    Cdc20a = rStateVariables[1];
    Cdh1 = rStateVariables[2];
    m = rStateVariables[3];
    Cdc20t = rStateVariables[4];
    IEP = rStateVariables[5];
    CKIt = rStateVariables[6];
    SK = rStateVariables[7];
}

// MODEL FUNCTIONS
inline double TysonNovak2001SbmlOdeSystem::GK(double A1, double A2, double A3, double A4)
{
    return 2.0 * A4 * A1 / ((A2 - A1) + A3 * A2 + A4 * A1 + sm::root(2.0, std::pow(((A2 - A1) + A3 * A2 + A4 * A1), 2.0) - 4.0 * (A2 - A1) * A4 * A1));
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
    this->mDerivedQuantityNames.push_back("cell");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__CycBt");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("CycB");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__CycB");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__Cdc20a");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Trimer");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__Trimer");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__Cdh1");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__m");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__Cdc20t");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__IEP");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Mad");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__Mad");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__CKIt");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__SK");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("TF");
    this->mDerivedQuantityUnits.push_back("non-dim");

    // PARAMETERS
    this->mParameterNames.push_back("k1");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("k2p");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("k2pp");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("k2ppp");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("k3p");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("k3pp");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("J3");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("k4");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("k5p");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("k5pp");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("J5");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("k6");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("n");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("k7");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("J7");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("k8");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("J8");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("k9");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("k10");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("mu");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("k11");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("k12p");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("k12pp");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("mmax");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("k12ppp");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Keq");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("k13");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("k14");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("k15p");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("k15pp");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("k16p");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("k16pp");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("J15");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("J16");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("k4p");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("J4");
    this->mParameterUnits.push_back("non-dim");

    this->mInitialised = true;
}

// Register the ODE system with Boost serialization
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(TysonNovak2001SbmlOdeSystem)