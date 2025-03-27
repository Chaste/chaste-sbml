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
    return 2 * A4 * A1 / (A2 - A1 + A3 * A2 + A4 * A1 + sbmlmath::sm_root(2, pow(A2 - A1 + A3 * A2 + A4 * A1, 2) - 4 * (A2 - A1) * A4 * A1));
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
    /* Define state variables */
    double CycBt = rY[0]; // CycBt
    double CycB = rY[1]; // CycB
    double Cdc20a = rY[2]; // Cdc20a
    double Trimer = rY[3]; // Trimer
    double Cdh1 = rY[4]; // Cdh1
    double m = rY[5]; // m
    double Cdc20t = rY[6]; // Cdc20t
    double IEP = rY[7]; // IEP
    double Mad = rY[8]; // Mad
    double CKIt = rY[9]; // CKIt
    double SK = rY[10]; // SK

    /* Define state parameters. */


    /* Define algebraic rules. */
    CycB = this->GetStateVariable("CycBt") - 2 * this->GetStateVariable("CycBt") * this->GetStateVariable("CKIt") / (this->GetStateVariable("CycBt") + this->GetStateVariable("CKIt") + 1 / Keq + pow(pow(this->GetStateVariable("CycBt") + this->GetStateVariable("CKIt") + 1 / Keq, 2) - 4 * this->GetStateVariable("CycBt") * this->GetStateVariable("CKIt"), 1 / 2));
    Trimer = 2 * this->GetStateVariable("CycBt") * this->GetStateVariable("CKIt") / (this->GetStateVariable("CycBt") + this->GetStateVariable("CKIt") + 1 / Keq + pow(pow(this->GetStateVariable("CycBt") + this->GetStateVariable("CKIt") + 1 / Keq, 2) - 4 * this->GetStateVariable("CycBt") * this->GetStateVariable("CKIt"), 1 / 2));
    TF = GK(k15p * this->GetStateVariable("m") + k15pp * this->GetStateVariable("SK"), k16p + k16pp * this->GetStateVariable("m") * this->GetStateVariable("CycB"), J15, J16);
    Mad = 1;

    /* Define the reactions in this model. */
    // CycBt synthesis
    double CycBt_synthesis = k1;

    // CycBt degradation
    double CycBdegradation = k2p * this->GetStateVariable("CycBt");

    // CycBt degradation via Cdh1
    double CycBdegradationviaCdh1 = k2pp * this->GetStateVariable("Cdh1") * this->GetStateVariable("CycBt");

    // CycBt degradation via Cdc20a
    double CycBtdegradationviaCdc20a = k2ppp * this->GetStateVariable("Cdc20a") * this->GetStateVariable("CycBt");

    // Cdh1 synthesis
    double Cdh1synthesis = (k3p + k3pp * this->GetStateVariable("Cdc20a")) * (1 - this->GetStateVariable("Cdh1")) / (J3 + 1 - this->GetStateVariable("Cdh1"));

    // Cdh1 degradation
    double Cdh1degradation = (k4p * this->GetStateVariable("SK") * this->GetStateVariable("Cdh1") + k4 * this->GetStateVariable("m") * this->GetStateVariable("CycB") * this->GetStateVariable("Cdh1")) / (J4 + this->GetStateVariable("Cdh1"));

    // Cdc20t synthesis
    double Cdc20tsynthesis = k5p + k5pp * pow(this->GetStateVariable("CycB") * this->GetStateVariable("m") / J5, n) / (1 + pow(this->GetStateVariable("CycB") * this->GetStateVariable("m") / J5, n));

    // Cdc20t degradation
    double Cdc20t_deg = k6 * this->GetStateVariable("Cdc20t");

    // Cdc20 activation
    double Cdc20activation = k7 * this->GetStateVariable("IEP") * (this->GetStateVariable("Cdc20t") - this->GetStateVariable("Cdc20a")) / (J7 + this->GetStateVariable("Cdc20t") - this->GetStateVariable("Cdc20a"));

    // Cdc20a inhibition
    double Cdc20ainhibition = k8 * this->GetStateVariable("Mad") * this->GetStateVariable("Cdc20a") / (J8 + this->GetStateVariable("Cdc20a"));

    // Cdc20a degradation
    double Cdc20adegradation = k6 * this->GetStateVariable("Cdc20a");

    // IEP synthesis
    double IEPsynthesis = k9 * this->GetStateVariable("m") * this->GetStateVariable("CycB") * (1 - this->GetStateVariable("IEP"));

    // IEP degradation
    double IEPdegradation = k10 * this->GetStateVariable("IEP");

    // growth
    double growth = mu * this->GetStateVariable("m") * (1 - this->GetStateVariable("m") / mmax);

    // CKIt synthesis
    double CKItsynthesis = k11;

    // CKIt degradation
    double CKIdegradation = k12p * this->GetStateVariable("CKIt");

    // CKIt phosphorilation via SK
    double CKItphosphorilationviaSK = k12pp * this->GetStateVariable("SK") * this->GetStateVariable("CKIt");

    // CKIt Trimer sequestred
    double eq_7 = k12ppp * this->GetStateVariable("m") * this->GetStateVariable("CycB") * this->GetStateVariable("CKIt");

    // SK synthesis
    double SKsynthesis = k13 * TF;

    // SK degradation
    double SKdegradation = k14 * this->GetStateVariable("SK");


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

    if (sbmlmath::sm_lt(this->GetStateVariable("CycB"), 0.1))
    {
<<<<<<< Updated upstream
        this->rGetStateVariables()[5] = double(rY[5] / 2);
=======
        this->SetStateVariable("m", static_cast<double>(this->GetStateVariable("m") / 2));
>>>>>>> Stashed changes
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

/* Define SRN model using Wrappers. */
#include "SbmlCellCycleWrapperModel.hpp"
#include "SbmlCellCycleWrapperModel.cpp"

typedef SbmlCellCycleWrapperModel<TysonNovak2001OdeSystem, 11> TysonNovak2001CellCycleModel;

// Declare identifiers for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(TysonNovak2001OdeSystem)
EXPORT_TEMPLATE_CLASS2(SbmlCellCycleWrapperModel, TysonNovak2001OdeSystem, 11)

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(TysonNovak2001CellCycleModel)