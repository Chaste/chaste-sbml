#include <cmath>
#include <limits>

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
    CycBt = 0.001;  // CycBt
    CycB = 0.0;     // CycB
    Cdc20a = 0.001; // Cdc20a
    Trimer = 0.0;   // Trimer
    Cdh1 = 0.001;   // Cdh1
    m = 0.5;        // m
    Cdc20t = 0.001; // Cdc20t
    IEP = 0.001;    // IEP
    Mad = 0.0;      // Mad
    CKIt = 0.001;   // CKIt
    SK = 0.001;     // SK

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

    mParameters.push_back(TF);

    // EVENTS:
    satisfiedEvents.resize(1, false);
}

TysonNovak2001SbmlOdeSystem::~TysonNovak2001SbmlOdeSystem()
{
}

void TysonNovak2001SbmlOdeSystem::EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY)
{
    RefreshState(rY);

    // RULES:
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
    double Cdh1synthesis = (k3p + k3pp * Cdc20a) * (1 - Cdh1) / (J3 + 1 - Cdh1);

    // Cdh1 degradation
    double Cdh1degradation = (k4p * SK * Cdh1 + k4 * m * CycB * Cdh1) / (J4 + Cdh1);

    // Cdc20t synthesis
    double Cdc20tsynthesis = k5p + k5pp * std::pow(CycB * m / J5, n) / (1 + std::pow(CycB * m / J5, n));

    // Cdc20t degradation
    double Cdc20t_deg = k6 * Cdc20t;

    // Cdc20 activation
    double Cdc20activation = k7 * IEP * (Cdc20t - Cdc20a) / (J7 + Cdc20t - Cdc20a);

    // Cdc20a inhibition
    double Cdc20ainhibition = k8 * Mad * Cdc20a / (J8 + Cdc20a);

    // Cdc20a degradation
    double Cdc20adegradation = k6 * Cdc20a;

    // IEP synthesis
    double IEPsynthesis = k9 * m * CycB * (1 - IEP);

    // IEP degradation
    double IEPdegradation = k10 * IEP;

    // growth
    double growth = mu * m * (1 - m / mmax);

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
    rDY[0] = (CycBt_synthesis - CycBdegradation - CycBdegradationviaCdh1 - CycBtdegradationviaCdc20a) / cell;                                                  // d[CycBt]/dt
    rDY[1] = ((CycBt - 2 * CycBt * CKIt / (CycBt + CKIt + 1 / Keq + std::pow(std::pow(CycBt + CKIt + 1 / Keq, 2) - 4 * CycBt * CKIt, 1 / 2))) - rY[1]) / cell; // d[CycB]/dt
    rDY[2] = (Cdc20activation - Cdc20ainhibition - Cdc20adegradation) / cell;                                                                                  // d[Cdc20a]/dt
    rDY[3] = ((2 * CycBt * CKIt / (CycBt + CKIt + 1 / Keq + std::pow(std::pow(CycBt + CKIt + 1 / Keq, 2) - 4 * CycBt * CKIt, 1 / 2))) - rY[3]) / cell;         // d[Trimer]/dt
    rDY[4] = (Cdh1synthesis - Cdh1degradation) / cell;                                                                                                         // d[Cdh1]/dt
    rDY[5] = (growth) / cell;                                                                                                                                  // d[m]/dt
    rDY[6] = (Cdc20tsynthesis - Cdc20t_deg) / cell;                                                                                                            // d[Cdc20t]/dt
    rDY[7] = (IEPsynthesis - IEPdegradation) / cell;                                                                                                           // d[IEP]/dt
    rDY[8] = ((1) - rY[8]) / cell;                                                                                                                             // d[Mad]/dt
    rDY[9] = (CKItsynthesis - CKIdegradation - CKItphosphorilationviaSK - eq_7) / cell;                                                                        // d[CKIt]/dt
    rDY[10] = (SKsynthesis - SKdegradation) / cell;                                                                                                            // d[SK]/dt

    // Scale time appropriately
}

void TysonNovak2001SbmlOdeSystem::RefreshState(const std::vector<double> &rY)
{
    // STATE VARIABLES:
    CycBt = rY[0];
    CycB = rY[1];
    Cdc20a = rY[2];
    Trimer = rY[3];
    Cdh1 = rY[4];
    m = rY[5];
    Cdc20t = rY[6];
    IEP = rY[7];
    Mad = rY[8];
    CKIt = rY[9];
    SK = rY[10];

    // STATE PARAMETERS:

    TF = GetParameter("TF");
}

bool TysonNovak2001SbmlOdeSystem::CalculateStoppingEvent(double time, const std::vector<double> &rY)
{
    bool triggered = false;

    if (sm::lt(CycB, 0.1))
    {
        if (time <= 0.0)
        {
            // Event satisfied at first timestep: mark but don't trigger
            satisfiedEvents[0] = true;
        }
        else if (!satisfiedEvents[0])
        {
            // Event previously unsatisfied: mark and trigger
            satisfiedEvents[0] = true;
            triggered = true;
            SetDefaultInitialCondition(5, m / 2.0);
        }
        // else { Event previously satisfied: do nothing }
    }
    else
    {
        // Event not satisfied: unmark
        satisfiedEvents[0] = false;
    }

    return triggered;
}

double TysonNovak2001SbmlOdeSystem::CalculateRootFunction(double time, const std::vector<double> &rY)
{
    bool triggered = CalculateStoppingEvent(time, rY);
    if (triggered)
    {
        std::cout << "CycB root function triggered at time: " << time << std::endl;
    }
    if (sm::lt(CycB, 0.1))
    {
        return triggered ? 0.0 : 1.0;
    }
    return std::abs(CycB - 0.1);
}

// FUNCTION DEFINITIONS:
double TysonNovak2001SbmlOdeSystem::GK(double A1, double A2, double A3, double A4)
{
    return 2 * A4 * A1 / (A2 - A1 + A3 * A2 + A4 * A1 + sm::root(2, std::pow(A2 - A1 + A3 * A2 + A4 * A1, 2) - 4 * (A2 - A1) * A4 * A1));
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