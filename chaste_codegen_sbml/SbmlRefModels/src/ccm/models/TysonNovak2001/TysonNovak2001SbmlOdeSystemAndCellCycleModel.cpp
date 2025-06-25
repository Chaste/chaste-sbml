#include <cmath>
#include <limits>
#include <vector>

#include "CellwiseOdeSystemInformation.hpp"
#include "SbmlMath.hpp"

#include "TysonNovak2001SbmlOdeSystemAndCellCycleModel.hpp"

namespace sm = sbmlmath;

TysonNovak2001SbmlOdeSystem::TysonNovak2001SbmlOdeSystem(std::vector<double> stateVariables)
    : AbstractOdeSystem(11)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<TysonNovak2001SbmlOdeSystem>);

    // COMPARTMENTS:
    cell = 1.0;

    // STATE VARIABLES:
    CycBt = 0.001; // CycBt
    Cdc20a = 0.001; // Cdc20a
    Cdh1 = 0.001; // Cdh1
    m = 0.5; // m
    Cdc20t = 0.001; // Cdc20t
    IEP = 0.001; // IEP
    CKIt = 0.001; // CKIt
    SK = 0.001; // SK

    CycB = CycBt - 2.0 * CycBt * CKIt / (CycBt + CKIt + 1.0 / Keq + std::pow(std::pow(CycBt + CKIt + 1.0 / Keq, 2.0) - 4.0 * CycBt * CKIt, 1.0 / 2.0));
    Trimer = 2.0 * CycBt * CKIt / (CycBt + CKIt + 1.0 / Keq + std::pow(std::pow(CycBt + CKIt + 1.0 / Keq, 2.0) - 4.0 * CycBt * CKIt, 1.0 / 2.0));
    Mad = 1.0;

    SetDefaultInitialCondition(0, CycBt);
    SetDefaultInitialCondition(1, CycB);
    SetDefaultInitialCondition(2, Cdc20a);
    SetDefaultInitialCondition(3, Trimer);
    SetDefaultInitialCondition(4, Cdh1);
    SetDefaultInitialCondition(5, m);
    SetDefaultInitialCondition(6, Cdc20t);
    SetDefaultInitialCondition(7, IEP);
    SetDefaultInitialCondition(8, Mad);
    SetDefaultInitialCondition(9, CKIt);
    SetDefaultInitialCondition(10, SK);

    if (stateVariables.size() == 11)
    {
        CycBt = stateVariables[0];
        CycB = stateVariables[1];
        Cdc20a = stateVariables[2];
        Trimer = stateVariables[3];
        Cdh1 = stateVariables[4];
        m = stateVariables[5];
        Cdc20t = stateVariables[6];
        IEP = stateVariables[7];
        Mad = stateVariables[8];
        CKIt = stateVariables[9];
        SK = stateVariables[10];
    }
    else if (stateVariables.size() != 0)
    {
        EXCEPTION("TysonNovak2001SbmlOdeSystem: Expected 11 state variables, got " + std::to_string(stateVariables.size()));
    }

    mStateVariables.push_back(CycBt);
    mStateVariables.push_back(CycB);
    mStateVariables.push_back(Cdc20a);
    mStateVariables.push_back(Trimer);
    mStateVariables.push_back(Cdh1);
    mStateVariables.push_back(m);
    mStateVariables.push_back(Cdc20t);
    mStateVariables.push_back(IEP);
    mStateVariables.push_back(Mad);
    mStateVariables.push_back(CKIt);
    mStateVariables.push_back(SK);

    // STATE PARAMETERS:

    TF = 0.0;

    TF = GK(k15p * m + k15pp * SK, k16p + k16pp * m * CycB, J15, J16);


    mParameters.push_back(TF);

    // EVENTS:
    eventsSatisfied.resize(1, false);
    eventsInitialised = false;
}

TysonNovak2001SbmlOdeSystem::~TysonNovak2001SbmlOdeSystem()
{
}

void TysonNovak2001SbmlOdeSystem::RefreshState(const std::vector<double> &rY)
{
    // STATE VARIABLES:
    CycBt = rY[0];
    Cdc20a = rY[2];
    Cdh1 = rY[4];
    m = rY[5];
    Cdc20t = rY[6];
    IEP = rY[7];
    CKIt = rY[9];
    SK = rY[10];

    // STATE PARAMETERS:

    TF = GetParameter("TF");
}

void TysonNovak2001SbmlOdeSystem::EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY)
{
    RefreshState(rY);

    // RULES:
    CycB = CycBt - 2.0 * CycBt * CKIt / (CycBt + CKIt + 1.0 / Keq + std::pow(std::pow(CycBt + CKIt + 1.0 / Keq, 2.0) - 4.0 * CycBt * CKIt, 1.0 / 2.0));
    Trimer = 2.0 * CycBt * CKIt / (CycBt + CKIt + 1.0 / Keq + std::pow(std::pow(CycBt + CKIt + 1.0 / Keq, 2.0) - 4.0 * CycBt * CKIt, 1.0 / 2.0));
    Mad = 1.0;

    TF = GK(k15p * m + k15pp * SK, k16p + k16pp * m * CycB, J15, J16);

    // UPDATE STATE PARAMETERS:

    SetParameter("TF", TF);

    // REACTIONS:

    // CycBt synthesis
    double CycBt_synthesis = k1;

    // CycBt degradation
    double CycBdegradation = k2p * CycBt;

    // CycBt degradation via Cdh1
    double CycBdegradationviaCdh1 = k2pp * Cdh1 * CycBt;

    // CycBt degradation via Cdc20a
    double CycBtdegradationviaCdc20a = k2ppp * Cdc20a * CycBt;

    // Cdh1 synthesis
    double Cdh1synthesis = (k3p + k3pp * Cdc20a) * (1.0 - Cdh1) / (J3 + 1.0 - Cdh1);

    // Cdh1 degradation
    double Cdh1degradation = (k4p * SK * Cdh1 + k4 * m * CycB * Cdh1) / (J4 + Cdh1);

    // Cdc20t synthesis
    double Cdc20tsynthesis = k5p + k5pp * std::pow(CycB * m / J5, n) / (1.0 + std::pow(CycB * m / J5, n));

    // Cdc20t degradation
    double Cdc20t_deg = k6 * Cdc20t;

    // Cdc20 activation
    double Cdc20activation = k7 * IEP * (Cdc20t - Cdc20a) / (J7 + Cdc20t - Cdc20a);

    // Cdc20a inhibition
    double Cdc20ainhibition = k8 * Mad * Cdc20a / (J8 + Cdc20a);

    // Cdc20a degradation
    double Cdc20adegradation = k6 * Cdc20a;

    // IEP synthesis
    double IEPsynthesis = k9 * m * CycB * (1.0 - IEP);

    // IEP degradation
    double IEPdegradation = k10 * IEP;

    // growth
    double growth = mu * m * (1.0 - m / mmax);

    // CKIt synthesis
    double CKItsynthesis = k11;

    // CKIt degradation
    double CKIdegradation = k12p * CKIt;

    // CKIt phosphorilation via SK
    double CKItphosphorilationviaSK = k12pp * SK * CKIt;

    // CKIt Trimer sequestred
    double eq_7 = k12ppp * m * CycB * CKIt;

    // SK synthesis
    double SKsynthesis = k13 * TF;

    // SK degradation
    double SKdegradation = k14 * SK;

    // ODES:
    rDY[0] = (CycBt_synthesis - CycBdegradation - CycBdegradationviaCdh1 - CycBtdegradationviaCdc20a) / cell; // d[CycBt]/dt
    rDY[1] = (CycB - rY[1]) / cell; // d[CycB]/dt
    rDY[2] = (Cdc20activation - Cdc20ainhibition - Cdc20adegradation) / cell; // d[Cdc20a]/dt
    rDY[3] = (Trimer - rY[3]) / cell; // d[Trimer]/dt
    rDY[4] = (Cdh1synthesis - Cdh1degradation) / cell; // d[Cdh1]/dt
    rDY[5] = (growth) / cell; // d[m]/dt
    rDY[6] = (Cdc20tsynthesis - Cdc20t_deg) / cell; // d[Cdc20t]/dt
    rDY[7] = (IEPsynthesis - IEPdegradation) / cell; // d[IEP]/dt
    rDY[8] = (Mad - rY[8]) / cell; // d[Mad]/dt
    rDY[9] = (CKItsynthesis - CKIdegradation - CKItphosphorilationviaSK - eq_7) / cell; // d[CKIt]/dt
    rDY[10] = (SKsynthesis - SKdegradation) / cell; // d[SK]/dt

    // Scale time appropriately
}

double TysonNovak2001SbmlOdeSystem::ProcessEvents(double time, const std::vector<double> &rY)
{
    RefreshState(rY);

    double min_dist = std::numeric_limits<double>::max();
    double event_dist = min_dist;

    // EVENT: sm::lt(CycB, 0.1)
    event_dist = (0.1) - (CycB) - std::numeric_limits<double>::epsilon();

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
        if (!eventsSatisfied[0] && eventsInitialised)
        {
            // The condition is transitioning from false to true,
            // and this is not the first time-step => trigger the event.
            event_dist = 0.0;
            min_dist = 0.0;

            UpdateDefaultInitialConditions(rY);
            SetStateVariable(5, m / 2.0);
            SetDefaultInitialCondition(5, m / 2.0);
        }
        eventsSatisfied[0] = true; // Flag the condition true
    }
    else
    {
        eventsSatisfied[0] = false; // Flag the condition false
    }

    eventsInitialised = true; // Flag that events have been processed at least once

    // Distance to closest event
    return min_dist;
}

double TysonNovak2001SbmlOdeSystem::CalculateRootFunction(double time, const std::vector<double> &rY)
{
    return ProcessEvents(time, rY);
}

bool TysonNovak2001SbmlOdeSystem::CalculateStoppingEvent(double time, const std::vector<double> &rY)
{
    return ProcessEvents(time, rY) == 0.0;
}

void TysonNovak2001SbmlOdeSystem::UpdateDefaultInitialConditions(const std::vector<double> &rY)
{
    SetDefaultInitialCondition(0, rY[0]); // CycBt
    SetDefaultInitialCondition(1, rY[1]); // CycB
    SetDefaultInitialCondition(2, rY[2]); // Cdc20a
    SetDefaultInitialCondition(3, rY[3]); // Trimer
    SetDefaultInitialCondition(4, rY[4]); // Cdh1
    SetDefaultInitialCondition(5, rY[5]); // m
    SetDefaultInitialCondition(6, rY[6]); // Cdc20t
    SetDefaultInitialCondition(7, rY[7]); // IEP
    SetDefaultInitialCondition(8, rY[8]); // Mad
    SetDefaultInitialCondition(9, rY[9]); // CKIt
    SetDefaultInitialCondition(10, rY[10]); // SK
}

// FUNCTION DEFINITIONS:
double TysonNovak2001SbmlOdeSystem::GK(double A1, double A2, double A3, double A4)
{
    return 2.0 * A4 * A1 / (A2 - A1 + A3 * A2 + A4 * A1 + sm::root(2.0, std::pow(A2 - A1 + A3 * A2 + A4 * A1, 2.0) - 4.0 * (A2 - A1) * A4 * A1));
}

template <>
void CellwiseOdeSystemInformation<TysonNovak2001SbmlOdeSystem>::Initialise()
{
    // STATE VARIABLES:
    this->mVariableNames.push_back("CycBt");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.001);

    this->mVariableNames.push_back("CycB");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("Cdc20a");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.001);

    this->mVariableNames.push_back("Trimer");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

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

    this->mVariableNames.push_back("Mad");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("CKIt");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.001);

    this->mVariableNames.push_back("SK");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.001);


    // STATE PARAMETERS:
    this->mParameterNames.push_back("TF");
    this->mParameterUnits.push_back("non-dim");

    this->mInitialised = true;
}

// Define SbmlCellCycleWrapperModel using wrappers
#include "SbmlCellCycleWrapperModel.hpp"
#include "SbmlCellCycleWrapperModel.cpp"

typedef SbmlCellCycleWrapperModel<TysonNovak2001SbmlOdeSystem, 11> TysonNovak2001SbmlCellCycleModel;

// Declare identifiers for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(TysonNovak2001SbmlOdeSystem)
EXPORT_TEMPLATE_CLASS2(SbmlCellCycleWrapperModel, TysonNovak2001SbmlOdeSystem, 11)

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(TysonNovak2001SbmlCellCycleModel)