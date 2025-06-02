#include "CellwiseOdeSystemInformation.hpp"
#include "SbmlMath.hpp"

#include "Chen2000SbmlOdeSystemAndCellCycleModel.hpp"

namespace sm = sbmlmath;

Chen2000SbmlOdeSystem::Chen2000SbmlOdeSystem(std::vector<double> stateVariables)
    : AbstractOdeSystem(0)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<Chen2000SbmlOdeSystem>);

    // COMPARTMENTS:
    COMpartment = 1.0;

    // STATE VARIABLES:


    if (stateVariables.size() == 0)
    {
    }
    else if (stateVariables.size() != 0)
    {
        EXCEPTION("Chen2000SbmlOdeSystem: Expected 0 state variables, got " + std::to_string(stateVariables.size()));
    }


    // STATE PARAMETERS:

    Cln2 = 0.0078;
    Clb2_T = 0.2342;
    Vd_b2 = 0.0;
    Clb2 = 0.0;
    Clb5 = 0.0;
    Sic1 = 0.0;
    Clb5_T = 0.0614;
    Vd_b5 = 0.0;
    Bck2 = 0.0;
    Cln3 = 0.0;
    Sic1_T = 0.1231;
    Clb2_Sic1 = 0.079;
    Clb5_Sic1 = 0.0207;
    Vd2_c1 = 0.0;
    Cdc20_T = 0.8332;
    Cdc20 = 0.6848;
    Vi_20 = 0.0;
    Hct1 = 0.9946;
    Vi_t1 = 0.0;
    mass = 0.6608;
    ORI = 0.0;
    BUD = 0.0;
    SPN = 0.0;
    SBF = 0.0;
    Va_sbf = 0.0;
    MBF = 0.0;
    Mcm1 = 0.0;
    Swi5 = 0.0;


    mParameters.push_back(Cln2);
    mParameters.push_back(Clb2_T);
    mParameters.push_back(Vd_b2);
    mParameters.push_back(Clb2);
    mParameters.push_back(Clb5);
    mParameters.push_back(Sic1);
    mParameters.push_back(Clb5_T);
    mParameters.push_back(Vd_b5);
    mParameters.push_back(Bck2);
    mParameters.push_back(Cln3);
    mParameters.push_back(Sic1_T);
    mParameters.push_back(Clb2_Sic1);
    mParameters.push_back(Clb5_Sic1);
    mParameters.push_back(Vd2_c1);
    mParameters.push_back(Cdc20_T);
    mParameters.push_back(Cdc20);
    mParameters.push_back(Vi_20);
    mParameters.push_back(Hct1);
    mParameters.push_back(Vi_t1);
    mParameters.push_back(mass);
    mParameters.push_back(ORI);
    mParameters.push_back(BUD);
    mParameters.push_back(SPN);
    mParameters.push_back(SBF);
    mParameters.push_back(Va_sbf);
    mParameters.push_back(MBF);
    mParameters.push_back(Mcm1);
    mParameters.push_back(Swi5);

}

Chen2000SbmlOdeSystem::~Chen2000SbmlOdeSystem()
{
}


void Chen2000SbmlOdeSystem::EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY)
{
    // STATE VARIABLES:

    // STATE PARAMETERS:

    Cln2 = GetParameter("Cln2");
    Clb2_T = GetParameter("Clb2_T");
    Vd_b2 = GetParameter("Vd_b2");
    Clb2 = GetParameter("Clb2");
    Clb5 = GetParameter("Clb5");
    Sic1 = GetParameter("Sic1");
    Clb5_T = GetParameter("Clb5_T");
    Vd_b5 = GetParameter("Vd_b5");
    Bck2 = GetParameter("Bck2");
    Cln3 = GetParameter("Cln3");
    Sic1_T = GetParameter("Sic1_T");
    Clb2_Sic1 = GetParameter("Clb2_Sic1");
    Clb5_Sic1 = GetParameter("Clb5_Sic1");
    Vd2_c1 = GetParameter("Vd2_c1");
    Cdc20_T = GetParameter("Cdc20_T");
    Cdc20 = GetParameter("Cdc20");
    Vi_20 = GetParameter("Vi_20");
    Hct1 = GetParameter("Hct1");
    Vi_t1 = GetParameter("Vi_t1");
    mass = GetParameter("mass");
    ORI = GetParameter("ORI");
    BUD = GetParameter("BUD");
    SPN = GetParameter("SPN");
    SBF = GetParameter("SBF");
    Va_sbf = GetParameter("Va_sbf");
    MBF = GetParameter("MBF");
    Mcm1 = GetParameter("Mcm1");
    Swi5 = GetParameter("Swi5");

    // RULES:
    Cln2 = mass * (ks_n2 + ks_n2_ * SBF) - kd_n2 * Cln2;
    Clb2_T = mass * (ks_b2 + ks_b2_ * Mcm1) - Vd_b2 * Clb2_T;
    Clb5_T = mass * (ks_b5 + ks_b5_ * MBF) - Vd_b5 * Clb5_T;
    Sic1_T = ks_c1 + ks_c1_ * Swi5 - Sic1_T * (kd1_c1 + Vd2_c1 / (Jd2_c1 + Sic1_T));
    Clb2_Sic1 = kas_b2 * Clb2 * Sic1 - Clb2_Sic1 * (kdi_b2 + Vd_b2 + kd1_c1 + Vd2_c1 / (Jd2_c1 + Sic1_T));
    Clb5_Sic1 = kas_b5 * Clb5 * Sic1 - Clb5_Sic1 * (kdi_b5 + Vd_b5 + kd1_c1 + Vd2_c1 / (Jd2_c1 + Sic1_T));
    Cdc20_T = ks_20 + ks_20_ * Clb2 - kd_20 * Cdc20_T;
    Cdc20 = ka_20 * (Cdc20_T - Cdc20) - Cdc20 * (Vi_20 + kd_20);
    Hct1 = (ka_t1 + ka_t1_ * Cdc20) * (Hct1_T - Hct1) / (Ja_t1 + Hct1_T - Hct1) - Vi_t1 * Hct1 / (Ji_t1 + Hct1);
    mass = mu * mass;
    ORI = ks_ori * (Clb5 + epsilonori_b2 * Clb2) - kd_ori * ORI;
    BUD = ks_bud * (Cln2 + Cln3 + epsilonbud_b5 * Clb5) - kd_bud * BUD;
    SPN = ks_spn * Clb2 / (J_spn + Clb2) - kd_spn * SPN;
    Vd_b2 = kd_b2 * (Hct1_T - Hct1) + kd_b2_ * Hct1 + kd_b2__ * Cdc20;
    Clb2 = Clb2_T - Clb2_Sic1;
    Clb5 = Clb5_T - Clb5_Sic1;
    Sic1 = Sic1_T - (Clb2_Sic1 + Clb5_Sic1);
    Vd_b5 = kd_b5 + kd_b5_ * Cdc20;
    Bck2 = Bck2_0 * mass;
    Cln3 = Cln3_max * Dn3 * mass / (Jn3 + Dn3 * mass);
    Vd2_c1 = kd2_c1 * (epsilonc1_n3 * Cln3 + epsilonc1_k2 * Bck2 + Cln2 + epsilonc1_b5 * Clb5 + epsilonc1_b2 * Clb2);
    Vi_20 = sm::piecewise(ki_20_, sm::geq(ORI, 1), ki_20, sm::geq(SPN, 1), 0.1);
    Vi_t1 = ki_t1 + ki_t1_ * (Cln3 + epsiloni_t1_n2 * Cln2 + epsiloni_t1_b5 * Clb5 + epsiloni_t1_b2 * Clb2);
    SBF = 2 * Va_sbf * Ji_sbf / (ki_sbf + ki_sbf_ * Clb2 + Va_sbf * Ji_sbf + (ki_sbf + ki_sbf_ * Clb2) * Ja_sbf - Va_sbf + std::sqrt(std::pow(ki_sbf + ki_sbf_ * Clb2 + Va_sbf * Ji_sbf + (ki_sbf + ki_sbf_ * Clb2) * Ja_sbf - Va_sbf, 2) - 4 * Va_sbf * Ji_sbf * (ki_sbf + ki_sbf_ * Clb2 - Va_sbf)));
    Va_sbf = ka_sbf * (Cln2 + epsilonsbf_n3 * (Cln3 + Bck2) + epsilonsbf_b5 * Clb5);
    MBF = SBF;
    Mcm1 = 2 * ka_mcm * Clb2 * Ji_mcm / (ki_mcm + ka_mcm * Clb2 * Ji_mcm + ki_mcm * Ja_mcm - ka_mcm * Clb2 + std::sqrt(std::pow(ki_mcm + ka_mcm * Clb2 * Ji_mcm + ki_mcm * Ja_mcm - ka_mcm * Clb2, 2) - 4 * (ki_mcm - ka_mcm * Clb2) * ka_mcm * Clb2 * Ji_mcm));
    Swi5 = 2 * ka_swi * Cdc20 * Ji_swi / (ki_swi + ki_swi_ * Clb2 + ka_swi * Cdc20 * Ji_swi + (ki_swi + ki_swi_ * Clb2) * Ja_swi - ka_swi * Cdc20 + std::sqrt(std::pow(ki_swi + ki_swi_ * Clb2 + ka_swi * Cdc20 * Ji_swi + (ki_swi + ki_swi_ * Clb2) * Ja_swi - ka_swi * Cdc20, 2) - 4 * (ki_swi + ki_swi_ * Clb2 - ka_swi * Cdc20) * ka_swi * Cdc20 * Ji_swi));

    // UPDATE STATE PARAMETERS:

    SetParameter("Cln2", Cln2);
    SetParameter("Clb2_T", Clb2_T);
    SetParameter("Vd_b2", Vd_b2);
    SetParameter("Clb2", Clb2);
    SetParameter("Clb5", Clb5);
    SetParameter("Sic1", Sic1);
    SetParameter("Clb5_T", Clb5_T);
    SetParameter("Vd_b5", Vd_b5);
    SetParameter("Bck2", Bck2);
    SetParameter("Cln3", Cln3);
    SetParameter("Sic1_T", Sic1_T);
    SetParameter("Clb2_Sic1", Clb2_Sic1);
    SetParameter("Clb5_Sic1", Clb5_Sic1);
    SetParameter("Vd2_c1", Vd2_c1);
    SetParameter("Cdc20_T", Cdc20_T);
    SetParameter("Cdc20", Cdc20);
    SetParameter("Vi_20", Vi_20);
    SetParameter("Hct1", Hct1);
    SetParameter("Vi_t1", Vi_t1);
    SetParameter("mass", mass);
    SetParameter("ORI", ORI);
    SetParameter("BUD", BUD);
    SetParameter("SPN", SPN);
    SetParameter("SBF", SBF);
    SetParameter("Va_sbf", Va_sbf);
    SetParameter("MBF", MBF);
    SetParameter("Mcm1", Mcm1);
    SetParameter("Swi5", Swi5);

    // ODES:

    // Scale time appropriately
}


template <>
void CellwiseOdeSystemInformation<Chen2000SbmlOdeSystem>::Initialise()
{
    // STATE VARIABLES:

    // STATE PARAMETERS:
    this->mParameterNames.push_back("Cln2");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Clb2_T");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Vd_b2");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Clb2");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Clb5");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Sic1");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Clb5_T");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Vd_b5");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Bck2");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Cln3");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Sic1_T");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Clb2_Sic1");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Clb5_Sic1");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Vd2_c1");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Cdc20_T");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Cdc20");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Vi_20");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Hct1");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Vi_t1");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("mass");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ORI");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("BUD");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("SPN");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("SBF");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Va_sbf");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("MBF");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Mcm1");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Swi5");
    this->mParameterUnits.push_back("non-dim");

    this->mInitialised = true;
}

// Define SbmlCellCycleWrapperModel using wrappers
#include "SbmlCellCycleWrapperModel.hpp"
#include "SbmlCellCycleWrapperModel.cpp"

typedef SbmlCellCycleWrapperModel<Chen2000SbmlOdeSystem, 0> Chen2000SbmlCellCycleModel;

// Declare identifiers for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(Chen2000SbmlOdeSystem)
EXPORT_TEMPLATE_CLASS2(SbmlCellCycleWrapperModel, Chen2000SbmlOdeSystem, 0)

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(Chen2000SbmlCellCycleModel)