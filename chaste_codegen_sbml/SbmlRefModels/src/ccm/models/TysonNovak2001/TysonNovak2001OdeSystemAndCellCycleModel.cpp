#include "CellwiseOdeSystemInformation.hpp"
#include "SbmlMath.hpp"

#include "TysonNovak2001OdeSystemAndCellCycleModel.hpp"

/* SBML ODE System */
TysonNovak2001OdeSystem::TysonNovak2001OdeSystem(std::vector<double> stateVariables)
    : AbstractOdeSystem(11)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<TysonNovak2001OdeSystem>);

    Init();

    SetDefaultInitialCondition(0, 0.001); // CycBt
    SetDefaultInitialCondition(1, 0.0); // CycB
    SetDefaultInitialCondition(2, 0.001); // Cdc20a
    SetDefaultInitialCondition(3, 0.0); // Trimer
    SetDefaultInitialCondition(4, 0.001); // Cdh1
    SetDefaultInitialCondition(5, 0.5); // m
    SetDefaultInitialCondition(6, 0.001); // Cdc20t
    SetDefaultInitialCondition(7, 0.001); // IEP
    SetDefaultInitialCondition(8, 0.0); // Mad
    SetDefaultInitialCondition(9, 0.001); // CKIt
    SetDefaultInitialCondition(10, 0.001); // SK


    if (stateVariables != std::vector<double>())
    {
        SetStateVariables(stateVariables);
    }
}

TysonNovak2001OdeSystem::~TysonNovak2001OdeSystem()
{
}

double TysonNovak2001OdeSystem::GK(double A1, double A2, double A3, double A4)
{
    return 2 * A4 * A1 / (A2 - A1 + A3 * A2 + A4 * A1 + sbmlmath::sm_root(2, std::pow(A2 - A1 + A3 * A2 + A4 * A1, 2) - 4 * (A2 - A1) * A4 * A1));
}

void TysonNovak2001OdeSystem::Init()
{
    /* Initialise model compartments. */
    cell = 1.0;

    /* Initialise model parameters. */
    k1 = 0.04;
    k2p = 0.04;
    k2pp = 1.0;
    k2ppp = 1.0;
    k3p = 1.0;
    k3pp = 10.0;
    J3 = 0.04;
    k4 = 35.0;
    k5p = 0.005;
    k5pp = 0.2;
    J5 = 0.3;
    k6 = 0.1;
    n = 4.0;
    k7 = 1.0;
    J7 = 0.001;
    k8 = 0.5;
    J8 = 0.001;
    k9 = 0.1;
    k10 = 0.02;
    mu = 0.005;
    k11 = 1.0;
    k12p = 0.2;
    k12pp = 50.0;
    mmax = 10.0;
    k12ppp = 100.0;
    Keq = 1000.0;
    k13 = 1.0;
    k14 = 1.0;
    k15p = 1.5;
    k15pp = 0.05;
    k16p = 1.0;
    k16pp = 3.0;
    J15 = 0.01;
    J16 = 0.01;
    k4p = 2.0;
    J4 = 0.04;
    TF = 0.0;

    /* Initialise vector to check if events have been triggered. */
    eventsSatisfied.resize(1, false);
}

void TysonNovak2001OdeSystem::EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY)
{
    /* Define algebraic rules. */
    double CycB = rY[0] - 2 * rY[0] * rY[9] / (rY[0] + rY[9] + 1 / Keq + std::pow(std::pow(rY[0] + rY[9] + 1 / Keq, 2) - 4 * rY[0] * rY[9], 1 / 2));
    double Trimer = 2 * rY[0] * rY[9] / (rY[0] + rY[9] + 1 / Keq + std::pow(std::pow(rY[0] + rY[9] + 1 / Keq, 2) - 4 * rY[0] * rY[9], 1 / 2));
    double TF = GK(k15p * rY[5] + k15pp * rY[10], k16p + k16pp * rY[5] * rY[1], J15, J16);
    double Mad = 1;

    /* Define the reactions in this model. */
    // CycBt synthesis
    double CycBt_synthesis = k1;

    // CycBt degradation
    double CycBdegradation = k2p * rY[0];

    // CycBt degradation via Cdh1
    double CycBdegradationviaCdh1 = k2pp * rY[4] * rY[0];

    // CycBt degradation via Cdc20a
    double CycBtdegradationviaCdc20a = k2ppp * rY[2] * rY[0];

    // Cdh1 synthesis
    double Cdh1synthesis = (k3p + k3pp * rY[2]) * (1 - rY[4]) / (J3 + 1 - rY[4]);

    // Cdh1 degradation
    double Cdh1degradation = (k4p * rY[10] * rY[4] + k4 * rY[5] * rY[1] * rY[4]) / (J4 + rY[4]);

    // Cdc20t synthesis
    double Cdc20tsynthesis = k5p + k5pp * std::pow(rY[1] * rY[5] / J5, n) / (1 + std::pow(rY[1] * rY[5] / J5, n));

    // Cdc20t degradation
    double Cdc20t_deg = k6 * rY[6];

    // Cdc20 activation
    double Cdc20activation = k7 * rY[7] * (rY[6] - rY[2]) / (J7 + rY[6] - rY[2]);

    // Cdc20a inhibition
    double Cdc20ainhibition = k8 * rY[8] * rY[2] / (J8 + rY[2]);

    // Cdc20a degradation
    double Cdc20adegradation = k6 * rY[2];

    // IEP synthesis
    double IEPsynthesis = k9 * rY[5] * rY[1] * (1 - rY[7]);

    // IEP degradation
    double IEPdegradation = k10 * rY[7];

    // growth
    double growth = mu * rY[5] * (1 - rY[5] / mmax);

    // CKIt synthesis
    double CKItsynthesis = k11;

    // CKIt degradation
    double CKIdegradation = k12p * rY[9];

    // CKIt phosphorilation via SK
    double CKItphosphorilationviaSK = k12pp * rY[10] * rY[9];

    // CKIt Trimer sequestred
    double eq_7 = k12ppp * rY[5] * rY[1] * rY[9];

    // SK synthesis
    double SKsynthesis = k13 * TF;

    // SK degradation
    double SKdegradation = k14 * rY[10];


    rDY[0] = (CycBt_synthesis - CycBdegradation - CycBdegradationviaCdh1 - CycBtdegradationviaCdc20a) / cell; // dCycBt/dt
    rDY[1] = (rDY[0] + rDY[0] + rDY[9] + rDY[9] + rDY[9] + rDY[0]) / cell; // dCycB/dt
    rDY[2] = (Cdc20activation - Cdc20ainhibition - Cdc20adegradation) / cell; // dCdc20a/dt
    rDY[3] = (rDY[0] + rDY[9] + rDY[9] + rDY[9] + rDY[0]) / cell; // dTrimer/dt
    rDY[4] = (Cdh1synthesis - Cdh1degradation) / cell; // dCdh1/dt
    rDY[5] = (growth) / cell; // dm/dt
    rDY[6] = (Cdc20tsynthesis - Cdc20t_deg) / cell; // dCdc20t/dt
    rDY[7] = (IEPsynthesis - IEPdegradation) / cell; // dIEP/dt
    rDY[8] = (CKItsynthesis - CKIdegradation - CKItphosphorilationviaSK - eq_7) / cell; // dCKIt/dt
    rDY[9] = (SKsynthesis - SKdegradation) / cell; // dSK/dt

    /* Account for the differences in timescales. */
}

bool TysonNovak2001OdeSystem::CalculateStoppingEvent(double time, const std::vector<double> & rY)
{
    // Return true if all events have been triggered.
    return AreAllEventsSatisfied(time, rY);
}

void TysonNovak2001OdeSystem::CheckAndUpdateEvents(double time, const std::vector<double> & rY)
{
    std::vector<double> dy(rY.size()); // Initialise derivatives vector
    EvaluateYDerivatives(time, rY, dy);

    if (sbmlmath::sm_lt(rY[1], 0.1))
    {
        this->SetStateVariable("m", static_cast<double>(rY[5] / 2));
        eventsSatisfied[0] = true;
    }
}

bool TysonNovak2001OdeSystem::AreAllEventsSatisfied(double time, const std::vector<double>& rY)
{
    CheckAndUpdateEvents(time, rY);
    bool events_satisfied = true;
    if (std::find(eventsSatisfied.begin(), eventsSatisfied.end(), false) != eventsSatisfied.end())
    {
        events_satisfied = false;
    }
    if (events_satisfied) // Reset events vector if cell division is triggered
    {
        std::fill(eventsSatisfied.begin(), eventsSatisfied.end(), false);
    }
    return events_satisfied;
}

template <>
void CellwiseOdeSystemInformation<TysonNovak2001OdeSystem>::Initialise()
{
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


    /* Define state parameters. */
    // Parameters without set values must be externally defined
    this->mInitialised = true;
}

/* Define SbmlCellCycleWrapperModel using wrappers. */
#include "SbmlCellCycleWrapperModel.hpp"
#include "SbmlCellCycleWrapperModel.cpp"

typedef SbmlCellCycleWrapperModel<TysonNovak2001OdeSystem, 11> TysonNovak2001CellCycleModel;

// Declare identifiers for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(TysonNovak2001OdeSystem)
EXPORT_TEMPLATE_CLASS2(SbmlCellCycleWrapperModel, TysonNovak2001OdeSystem, 11)

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(TysonNovak2001CellCycleModel)