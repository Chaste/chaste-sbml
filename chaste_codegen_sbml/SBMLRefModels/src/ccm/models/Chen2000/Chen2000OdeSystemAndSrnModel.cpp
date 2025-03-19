#include "Chen2000OdeSystemAndSrnModel.hpp"
#include "CellwiseOdeSystemInformation.hpp"

/* SBML ODE System */
Chen2000OdeSystem::Chen2000OdeSystem(std::vector<double> stateVariables)
    : AbstractOdeSystem(0)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<Chen2000OdeSystem>);

    Init();



    if (stateVariables != std::vector<double>())
    {
        SetStateVariables(stateVariables);
    }
}

Chen2000OdeSystem::~Chen2000OdeSystem()
{
}


void Chen2000OdeSystem::Init()
{
    /* Initialise model compartments. */
    COMpartment = 1.0;

    /* Initialise model parameters. */
    Cln2 = 0.0078;
    ks_n2 = 0.0;
    ks_n2_ = 0.05;
    kd_n2 = 0.1;
    Clb2_T = 0.2342;
    Hct1_T = 1.0;
    ks_b2 = 0.002;
    ks_b2_ = 0.05;
    kd_b2 = 0.01;
    kd_b2_ = 2.0;
    kd_b2__ = 0.05;
    Vd_b2 = 0.0;
    Clb2 = 0.0;
    Clb5 = 0.0;
    Sic1 = 0.0;
    Clb5_T = 0.0614;
    ks_b5 = 0.006;
    ks_b5_ = 0.02;
    kd_b5 = 0.1;
    kd_b5_ = 0.25;
    Vd_b5 = 0.0;
    Bck2 = 0.0;
    Bck2_0 = 0.0027;
    Cln3 = 0.0;
    Jn3 = 6.0;
    Dn3 = 1.0;
    Cln3_max = 0.02;
    Sic1_T = 0.1231;
    ks_c1 = 0.02;
    ks_c1_ = 0.1;
    Clb2_Sic1 = 0.079;
    kas_b2 = 50.0;
    kdi_b2 = 0.05;
    Clb5_Sic1 = 0.0207;
    kas_b5 = 50.0;
    kdi_b5 = 0.05;
    Vd2_c1 = 0.0;
    kd2_c1 = 0.3;
    epsilonc1_n3 = 20.0;
    epsilonc1_k2 = 2.0;
    epsilonc1_b5 = 1.0;
    epsilonc1_b2 = 0.067;
    Cdc20_T = 0.8332;
    ks_20 = 0.005;
    ks_20_ = 0.06;
    Cdc20 = 0.6848;
    ka_20 = 1.0;
    ki_20 = 0.1;
    ki_20_ = 10.0;
    Vi_20 = 0.0;
    Hct1 = 0.9946;
    ka_t1 = 0.04;
    ka_t1_ = 2.0;
    ki_t1 = 0.0;
    ki_t1_ = 0.64;
    Vi_t1 = 0.0;
    Ji_t1 = 0.05;
    Ja_t1 = 0.05;
    epsiloni_t1_n2 = 1.0;
    epsiloni_t1_b5 = 0.5;
    epsiloni_t1_b2 = 1.0;
    mass = 0.6608;
    mu = 0.005776;
    ORI = 0.0;
    ks_ori = 2.0;
    kd_ori = 0.06;
    epsilonori_b2 = 0.4;
    BUD = 0.0;
    ks_bud = 0.3;
    kd_bud = 0.06;
    epsilonbud_b5 = 1.0;
    SPN = 0.0;
    ks_spn = 0.08;
    kd_spn = 0.06;
    J_spn = 0.2;
    SBF = 0.0;
    ka_sbf = 1.0;
    ki_sbf = 0.5;
    ki_sbf_ = 6.0;
    Va_sbf = 0.0;
    Ji_sbf = 0.01;
    Ja_sbf = 0.01;
    epsilonsbf_n3 = 75.0;
    epsilonsbf_b5 = 0.5;
    MBF = 0.0;
    Mcm1 = 0.0;
    ka_mcm = 1.0;
    ki_mcm = 0.15;
    Ji_mcm = 1.0;
    Ja_mcm = 1.0;
    Swi5 = 0.0;
    ka_swi = 1.0;
    ki_swi = 0.3;
    ki_swi_ = 0.2;
    Ji_swi = 0.1;
    Ja_swi = 0.1;
    kd1_c1 = 0.01;
    kd_20 = 0.08;
    Jd2_c1 = 0.05;

}

void Chen2000OdeSystem::EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY)
{
    /* Define state variables */

    /* Define state parameters. */


    /* Define algebraic rules. */
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
    Vi_20 = piecewise(ki_20_, >=(ORI, 1), ki_20, >=(SPN, 1), 0.1);
    Vi_t1 = ki_t1 + ki_t1_ * (Cln3 + epsiloni_t1_n2 * Cln2 + epsiloni_t1_b5 * Clb5 + epsiloni_t1_b2 * Clb2);
    SBF = 2 * Va_sbf * Ji_sbf / (ki_sbf + ki_sbf_ * Clb2 + Va_sbf * Ji_sbf + (ki_sbf + ki_sbf_ * Clb2) * Ja_sbf - Va_sbf + sqrt(pow(ki_sbf + ki_sbf_ * Clb2 + Va_sbf * Ji_sbf + (ki_sbf + ki_sbf_ * Clb2) * Ja_sbf - Va_sbf, 2) - 4 * Va_sbf * Ji_sbf * (ki_sbf + ki_sbf_ * Clb2 - Va_sbf)));
    Va_sbf = ka_sbf * (Cln2 + epsilonsbf_n3 * (Cln3 + Bck2) + epsilonsbf_b5 * Clb5);
    MBF = SBF;
    Mcm1 = 2 * ka_mcm * Clb2 * Ji_mcm / (ki_mcm + ka_mcm * Clb2 * Ji_mcm + ki_mcm * Ja_mcm - ka_mcm * Clb2 + sqrt(pow(ki_mcm + ka_mcm * Clb2 * Ji_mcm + ki_mcm * Ja_mcm - ka_mcm * Clb2, 2) - 4 * (ki_mcm - ka_mcm * Clb2) * ka_mcm * Clb2 * Ji_mcm));
    Swi5 = 2 * ka_swi * Cdc20 * Ji_swi / (ki_swi + ki_swi_ * Clb2 + ka_swi * Cdc20 * Ji_swi + (ki_swi + ki_swi_ * Clb2) * Ja_swi - ka_swi * Cdc20 + sqrt(pow(ki_swi + ki_swi_ * Clb2 + ka_swi * Cdc20 * Ji_swi + (ki_swi + ki_swi_ * Clb2) * Ja_swi - ka_swi * Cdc20, 2) - 4 * (ki_swi + ki_swi_ * Clb2 - ka_swi * Cdc20) * ka_swi * Cdc20 * Ji_swi));

    /* Define the reactions in this model. */


    /* Account for the differences in timescales. */
}


template <>
void CellwiseOdeSystemInformation<Chen2000OdeSystem>::Initialise()
{

    /* Define state parameters. */
    // Parameters without set values must be externally defined
    this->mInitialised = true;
}

/* Define SRN model using Wrappers. */
#include "SbmlSrnWrapperModel.hpp"
#include "SbmlSrnWrapperModel.cpp"

typedef SbmlSrnWrapperModel<Chen2000OdeSystem, 0> Chen2000SrnModel;

// Declare identifiers for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(Chen2000OdeSystem)
EXPORT_TEMPLATE_CLASS2(SbmlSrnWrapperModel, Chen2000OdeSystem, 0)

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(Chen2000SrnModel)