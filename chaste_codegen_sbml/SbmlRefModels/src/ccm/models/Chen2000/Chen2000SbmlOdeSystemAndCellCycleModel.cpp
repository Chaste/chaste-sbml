#include <cmath>
#include <limits>
#include <vector>

#include "CellwiseOdeSystemInformation.hpp"
#include "SbmlEventType.hpp"
#include "SbmlMath.hpp"

#include "Chen2000SbmlOdeSystemAndCellCycleModel.hpp"

namespace sm = sbmlmath;

Chen2000SbmlOdeSystem::Chen2000SbmlOdeSystem(std::vector<double> stateVariables)
        : AbstractOdeSystem(13)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<Chen2000SbmlOdeSystem>);

    // STATE VARIABLES
    Cln2 = 0.0078;
    Clb2_T = 0.2342;
    Clb5_T = 0.0614;
    Sic1_T = 0.1231;
    Clb2_Sic1 = 0.079;
    Clb5_Sic1 = 0.0207;
    Cdc20_T = 0.8332;
    Cdc20 = 0.6848;
    Hct1 = 0.9946;
    mass = 0.6608;
    ORI = 0.0;
    BUD = 0.0;
    SPN = 0.0;

    SetDefaultInitialCondition(0, Cln2);
    SetDefaultInitialCondition(1, Clb2_T);
    SetDefaultInitialCondition(2, Clb5_T);
    SetDefaultInitialCondition(3, Sic1_T);
    SetDefaultInitialCondition(4, Clb2_Sic1);
    SetDefaultInitialCondition(5, Clb5_Sic1);
    SetDefaultInitialCondition(6, Cdc20_T);
    SetDefaultInitialCondition(7, Cdc20);
    SetDefaultInitialCondition(8, Hct1);
    SetDefaultInitialCondition(9, mass);
    SetDefaultInitialCondition(10, ORI);
    SetDefaultInitialCondition(11, BUD);
    SetDefaultInitialCondition(12, SPN);

    if (stateVariables.size() == 13)
    {
        Cln2 = stateVariables[0];
        Clb2_T = stateVariables[1];
        Clb5_T = stateVariables[2];
        Sic1_T = stateVariables[3];
        Clb2_Sic1 = stateVariables[4];
        Clb5_Sic1 = stateVariables[5];
        Cdc20_T = stateVariables[6];
        Cdc20 = stateVariables[7];
        Hct1 = stateVariables[8];
        mass = stateVariables[9];
        ORI = stateVariables[10];
        BUD = stateVariables[11];
        SPN = stateVariables[12];
    }
    else if (stateVariables.size() != 0)
    {
        EXCEPTION("Chen2000SbmlOdeSystem: Expected 13 state variables, got " + std::to_string(stateVariables.size()));
    }

    mStateVariables.push_back(Cln2);
    mStateVariables.push_back(Clb2_T);
    mStateVariables.push_back(Clb5_T);
    mStateVariables.push_back(Sic1_T);
    mStateVariables.push_back(Clb2_Sic1);
    mStateVariables.push_back(Clb5_Sic1);
    mStateVariables.push_back(Cdc20_T);
    mStateVariables.push_back(Cdc20);
    mStateVariables.push_back(Hct1);
    mStateVariables.push_back(mass);
    mStateVariables.push_back(ORI);
    mStateVariables.push_back(BUD);
    mStateVariables.push_back(SPN);

    // DERIVED QUANTITIES

    // VARIABLE PARAMETERS
    COMpartment = 1.0;

    mParameters.push_back(COMpartment);

    // RULE-BASED PARAMETERS
    Vd_b2 = 0.0;
    Clb2 = 0.0;
    Clb5 = 0.0;
    Sic1 = 0.0;
    Vd_b5 = 0.0;
    Bck2 = 0.0;
    Cln3 = 0.0;
    Vd2_c1 = 0.0;
    Vi_20 = 0.0;
    Vi_t1 = 0.0;
    SBF = 0.0;
    Va_sbf = 0.0;
    MBF = 0.0;
    Mcm1 = 0.0;
    Swi5 = 0.0;

    // REACTIONS

    // EVENTS
    mEventType.resize(0, SbmlEventType::UNKNOWN);

    // Uncomment lines below for events that should trigger cell division

    mEventSatisfied.resize(0, true); // Prevent events from triggering at the start
    mEventTriggered.resize(0, false);

    mEventAdjustedParameters.resize(1, false);
    mEventAdjustedParameterValues.resize(1, 0.0);

    mEventAdjustedStateVars.resize(13, false);
    mEventAdjustedStateValues.resize(13, 0.0);

    // Run model rules to complete state initialisation
    RunModelRules(0.0, mStateVariables);
}

Chen2000SbmlOdeSystem::Chen2000SbmlOdeSystem(const Chen2000SbmlOdeSystem& rOdeSystem)
        : Chen2000SbmlOdeSystem(rOdeSystem.mStateVariables)
{
    mEventSatisfied = rOdeSystem.mEventSatisfied;
    mEventTriggered = rOdeSystem.mEventTriggered;

    mEventAdjustedParameters = rOdeSystem.mEventAdjustedParameters;
    mEventAdjustedParameterValues = rOdeSystem.mEventAdjustedParameterValues;

    mEventAdjustedStateVars = rOdeSystem.mEventAdjustedStateVars;
    mEventAdjustedStateValues = rOdeSystem.mEventAdjustedStateValues;
}

Chen2000SbmlOdeSystem::~Chen2000SbmlOdeSystem()
{
}

void Chen2000SbmlOdeSystem::AdjustParameters(double time)
{
    for (unsigned i = 0; i < mEventAdjustedParameters.size(); ++i)
    {
        if (mEventAdjustedParameters[i])
        {
            SetParameter(i, mEventAdjustedParameterValues[i]);
        }
    }

    for (unsigned i = 0; i < mEventAdjustedStateVars.size(); ++i)
    {
        if (mEventAdjustedStateVars[i])
        {
            SetStateVariable(i, mEventAdjustedStateValues[i]);
            mEventAdjustedStateVars[i] = false;
        }
    }
}

double Chen2000SbmlOdeSystem::CalculateRootFunction(double time, const std::vector<double>& rY)
{
    return ProcessModelEvents(time, rY);
}

bool Chen2000SbmlOdeSystem::CalculateStoppingEvent(double time, const std::vector<double>& rY)
{
    return ProcessModelEvents(time, rY) == 0.0;
}

std::vector<double> Chen2000SbmlOdeSystem::ComputeDerivedQuantities(double time, const std::vector<double>& rY)
{
    RunModelRules(time, rY);

    std::vector<double> dqs;
    return dqs;
}

void Chen2000SbmlOdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    RunModelRules(time, rY);

    rDY[0] = mass * (ks_n2 + ks_n2_ * SBF) - kd_n2 * Cln2;                                                         // d[Cln2]/dt
    rDY[1] = mass * (ks_b2 + ks_b2_ * Mcm1) - Vd_b2 * Clb2_T;                                                      // d[Clb2_T]/dt
    rDY[2] = mass * (ks_b5 + ks_b5_ * MBF) - Vd_b5 * Clb5_T;                                                       // d[Clb5_T]/dt
    rDY[3] = ks_c1 + ks_c1_ * Swi5 - Sic1_T * (kd1_c1 + Vd2_c1 / (Jd2_c1 + Sic1_T));                               // d[Sic1_T]/dt
    rDY[4] = kas_b2 * Clb2 * Sic1 - Clb2_Sic1 * (kdi_b2 + Vd_b2 + kd1_c1 + Vd2_c1 / (Jd2_c1 + Sic1_T));            // d[Clb2_Sic1]/dt
    rDY[5] = kas_b5 * Clb5 * Sic1 - Clb5_Sic1 * (kdi_b5 + Vd_b5 + kd1_c1 + Vd2_c1 / (Jd2_c1 + Sic1_T));            // d[Clb5_Sic1]/dt
    rDY[6] = ks_20 + ks_20_ * Clb2 - kd_20 * Cdc20_T;                                                              // d[Cdc20_T]/dt
    rDY[7] = ka_20 * (Cdc20_T - Cdc20) - Cdc20 * (Vi_20 + kd_20);                                                  // d[Cdc20]/dt
    rDY[8] = (ka_t1 + ka_t1_ * Cdc20) * (Hct1_T - Hct1) / (Ja_t1 + Hct1_T - Hct1) - Vi_t1 * Hct1 / (Ji_t1 + Hct1); // d[Hct1]/dt
    rDY[9] = mu * mass;                                                                                            // d[mass]/dt
    rDY[10] = ks_ori * (Clb5 + epsilonori_b2 * Clb2) - kd_ori * ORI;                                               // d[ORI]/dt
    rDY[11] = ks_bud * (Cln2 + Cln3 + epsilonbud_b5 * Clb5) - kd_bud * BUD;                                        // d[BUD]/dt
    rDY[12] = ks_spn * Clb2 / (J_spn + Clb2) - kd_spn * SPN;                                                       // d[SPN]/dt

    // TODO: Scale time appropriately
}

bool Chen2000SbmlOdeSystem::HasEventOccurred(SbmlEventType eventType)
{
    for (unsigned i = 0; i < mEventTriggered.size(); ++i)
    {
        if (mEventTriggered[i] && mEventType[i] == eventType)
        {
            return true;
        }
    }
    return false;
}

double Chen2000SbmlOdeSystem::ProcessModelEvents(double time, const std::vector<double>& rY)
{
    std::fill(std::begin(mEventAdjustedParameters), std::end(mEventAdjustedParameters), false);
    std::fill(std::begin(mEventAdjustedStateVars), std::end(mEventAdjustedStateVars), false);

    double min_dist = std::numeric_limits<double>::max();

    return min_dist; // Distance to closest event
}

void Chen2000SbmlOdeSystem::ResetEventsOccurred()
{
    std::fill(mEventTriggered.begin(), mEventTriggered.end(), false);
}

void Chen2000SbmlOdeSystem::RunModelRules(double time, const std::vector<double>& rY)
{
    // STATE VARIABLES
    Cln2 = rY[0];
    Clb2_T = rY[1];
    Clb5_T = rY[2];
    Sic1_T = rY[3];
    Clb2_Sic1 = rY[4];
    Clb5_Sic1 = rY[5];
    Cdc20_T = rY[6];
    Cdc20 = rY[7];
    Hct1 = rY[8];
    mass = rY[9];
    ORI = rY[10];
    BUD = rY[11];
    SPN = rY[12];

    // VARIABLE PARAMETERS
    COMpartment = GetParameter(0);

    // ASSIGNMENT RULES
    Vd_b2 = kd_b2 * (Hct1_T - Hct1) + kd_b2_ * Hct1 + kd_b2__ * Cdc20;
    Clb2 = Clb2_T - Clb2_Sic1;
    Clb5 = Clb5_T - Clb5_Sic1;
    Sic1 = Sic1_T - (Clb2_Sic1 + Clb5_Sic1);
    Vd_b5 = kd_b5 + kd_b5_ * Cdc20;
    Bck2 = Bck2_0 * mass;
    Cln3 = Cln3_max * Dn3 * mass / (Jn3 + Dn3 * mass);
    Vd2_c1 = kd2_c1 * (epsilonc1_n3 * Cln3 + epsilonc1_k2 * Bck2 + Cln2 + epsilonc1_b5 * Clb5 + epsilonc1_b2 * Clb2);
    Vi_20 = sm::piecewise(ki_20_, sm::geq(ORI, 1.0), ki_20, sm::geq(SPN, 1.0), 0.1);
    Vi_t1 = ki_t1 + ki_t1_ * (Cln3 + epsiloni_t1_n2 * Cln2 + epsiloni_t1_b5 * Clb5 + epsiloni_t1_b2 * Clb2);
    Va_sbf = ka_sbf * (Cln2 + epsilonsbf_n3 * (Cln3 + Bck2) + epsilonsbf_b5 * Clb5);
    SBF = 2.0 * Va_sbf * Ji_sbf / (ki_sbf + ki_sbf_ * Clb2 + Va_sbf * Ji_sbf + (ki_sbf + ki_sbf_ * Clb2) * Ja_sbf - Va_sbf + std::sqrt(std::pow(ki_sbf + ki_sbf_ * Clb2 + Va_sbf * Ji_sbf + (ki_sbf + ki_sbf_ * Clb2) * Ja_sbf - Va_sbf, 2.0) - 4.0 * Va_sbf * Ji_sbf * (ki_sbf + ki_sbf_ * Clb2 - Va_sbf)));
    MBF = SBF;
    Mcm1 = 2.0 * ka_mcm * Clb2 * Ji_mcm / (ki_mcm + ka_mcm * Clb2 * Ji_mcm + ki_mcm * Ja_mcm - ka_mcm * Clb2 + std::sqrt(std::pow(ki_mcm + ka_mcm * Clb2 * Ji_mcm + ki_mcm * Ja_mcm - ka_mcm * Clb2, 2.0) - 4.0 * (ki_mcm - ka_mcm * Clb2) * ka_mcm * Clb2 * Ji_mcm));
    Swi5 = 2.0 * ka_swi * Cdc20 * Ji_swi / (ki_swi + ki_swi_ * Clb2 + ka_swi * Cdc20 * Ji_swi + (ki_swi + ki_swi_ * Clb2) * Ja_swi - ka_swi * Cdc20 + std::sqrt(std::pow(ki_swi + ki_swi_ * Clb2 + ka_swi * Cdc20 * Ji_swi + (ki_swi + ki_swi_ * Clb2) * Ja_swi - ka_swi * Cdc20, 2.0) - 4.0 * (ki_swi + ki_swi_ * Clb2 - ka_swi * Cdc20) * ka_swi * Cdc20 * Ji_swi));

    // REACTIONS
}

// MODEL FUNCTIONS

template <>
void CellwiseOdeSystemInformation<Chen2000SbmlOdeSystem>::Initialise()
{
    // STATE VARIABLES
    this->mVariableNames.push_back("Cln2");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0078);

    this->mVariableNames.push_back("Clb2_T");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.2342);

    this->mVariableNames.push_back("Clb5_T");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0614);

    this->mVariableNames.push_back("Sic1_T");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.1231);

    this->mVariableNames.push_back("Clb2_Sic1");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.079);

    this->mVariableNames.push_back("Clb5_Sic1");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0207);

    this->mVariableNames.push_back("Cdc20_T");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.8332);

    this->mVariableNames.push_back("Cdc20");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.6848);

    this->mVariableNames.push_back("Hct1");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.9946);

    this->mVariableNames.push_back("mass");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.6608);

    this->mVariableNames.push_back("ORI");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("BUD");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("SPN");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

    // DERIVED QUANTITIES

    // PARAMETERS
    this->mParameterNames.push_back("COMpartment");
    this->mParameterUnits.push_back("non-dim");

    this->mInitialised = true;
}

// Define SbmlCellCycleWrapperModel using wrappers
#include "SbmlCellCycleWrapperModel.cpp"
#include "SbmlCellCycleWrapperModel.hpp"

typedef SbmlCellCycleWrapperModel<Chen2000SbmlOdeSystem, 13> Chen2000SbmlCellCycleModel;

// Declare identifiers for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(Chen2000SbmlOdeSystem)
EXPORT_TEMPLATE_CLASS2(SbmlCellCycleWrapperModel, Chen2000SbmlOdeSystem, 13)

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(Chen2000SbmlCellCycleModel)