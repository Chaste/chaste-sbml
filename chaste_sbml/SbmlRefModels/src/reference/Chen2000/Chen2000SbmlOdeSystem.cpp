#include <cmath>
#include <limits>
#include <vector>

#include "CellwiseOdeSystemInformation.hpp"
#include "SbmlEventType.hpp"
#include "SbmlMath.hpp"

#include "Chen2000SbmlOdeSystem.hpp"

namespace sm = sbmlmath;

Chen2000SbmlOdeSystem::Chen2000SbmlOdeSystem()
        : AbstractSbmlOdeSystem(13, 71, 0)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<Chen2000SbmlOdeSystem>);

    Initialise();

    // EVENTS
}

Chen2000SbmlOdeSystem::~Chen2000SbmlOdeSystem()
{
}

std::vector<double> Chen2000SbmlOdeSystem::ComputeDerivedQuantities(double time, const std::vector<double>& rY)
{
    std::vector<double> dqs;
    RunModelEquations(time, rY);

    // AMOUNT / CONCENTRATION CONVERSIONS

    dqs.push_back(COMpartment);
    dqs.push_back(Vd_b2);
    dqs.push_back(Clb2);
    dqs.push_back(Clb5);
    dqs.push_back(Sic1);
    dqs.push_back(Vd_b5);
    dqs.push_back(Bck2);
    dqs.push_back(Cln3);
    dqs.push_back(Vd2_c1);
    dqs.push_back(Vi_20);
    dqs.push_back(Vi_t1);
    dqs.push_back(SBF);
    dqs.push_back(Va_sbf);
    dqs.push_back(MBF);
    dqs.push_back(Mcm1);
    dqs.push_back(Swi5);

    return dqs;
}

void Chen2000SbmlOdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    std::vector<double> derivatives = RunModelEquations(time, rY);
    for (unsigned i = 0; i < rDY.size(); ++i)
    {
        rDY[i] = derivatives[i];
    }

    // TODO: Scale time appropriately
}

void Chen2000SbmlOdeSystem::Initialise(double time)
{
    COMpartment = 1.0;                                                                                                                                                                                                                                                                                                                                                   //
    Cln2 = 0.0078;                                                                                                                                                                                                                                                                                                                                                       //
    ks_n2 = 0.0;                                                                                                                                                                                                                                                                                                                                                         //
    ks_n2_ = 0.05;                                                                                                                                                                                                                                                                                                                                                       //
    kd_n2 = 0.1;                                                                                                                                                                                                                                                                                                                                                         //
    Clb2_T = 0.2342;                                                                                                                                                                                                                                                                                                                                                     //
    Hct1_T = 1.0;                                                                                                                                                                                                                                                                                                                                                        //
    ks_b2 = 0.002;                                                                                                                                                                                                                                                                                                                                                       //
    ks_b2_ = 0.05;                                                                                                                                                                                                                                                                                                                                                       //
    kd_b2 = 0.01;                                                                                                                                                                                                                                                                                                                                                        //
    kd_b2_ = 2.0;                                                                                                                                                                                                                                                                                                                                                        //
    kd_b2__ = 0.05;                                                                                                                                                                                                                                                                                                                                                      //
    Clb5_T = 0.0614;                                                                                                                                                                                                                                                                                                                                                     //
    ks_b5 = 0.006;                                                                                                                                                                                                                                                                                                                                                       //
    ks_b5_ = 0.02;                                                                                                                                                                                                                                                                                                                                                       //
    kd_b5 = 0.1;                                                                                                                                                                                                                                                                                                                                                         //
    kd_b5_ = 0.25;                                                                                                                                                                                                                                                                                                                                                       //
    Bck2_0 = 0.0027;                                                                                                                                                                                                                                                                                                                                                     //
    Jn3 = 6.0;                                                                                                                                                                                                                                                                                                                                                           //
    Dn3 = 1.0;                                                                                                                                                                                                                                                                                                                                                           //
    Cln3_max = 0.02;                                                                                                                                                                                                                                                                                                                                                     //
    Sic1_T = 0.1231;                                                                                                                                                                                                                                                                                                                                                     //
    ks_c1 = 0.02;                                                                                                                                                                                                                                                                                                                                                        //
    ks_c1_ = 0.1;                                                                                                                                                                                                                                                                                                                                                        //
    Clb2_Sic1 = 0.079;                                                                                                                                                                                                                                                                                                                                                   //
    kas_b2 = 50.0;                                                                                                                                                                                                                                                                                                                                                       //
    kdi_b2 = 0.05;                                                                                                                                                                                                                                                                                                                                                       //
    Clb5_Sic1 = 0.0207;                                                                                                                                                                                                                                                                                                                                                  //
    kas_b5 = 50.0;                                                                                                                                                                                                                                                                                                                                                       //
    kdi_b5 = 0.05;                                                                                                                                                                                                                                                                                                                                                       //
    kd2_c1 = 0.3;                                                                                                                                                                                                                                                                                                                                                        //
    epsilonc1_n3 = 20.0;                                                                                                                                                                                                                                                                                                                                                 //
    epsilonc1_k2 = 2.0;                                                                                                                                                                                                                                                                                                                                                  //
    epsilonc1_b5 = 1.0;                                                                                                                                                                                                                                                                                                                                                  //
    epsilonc1_b2 = 0.067;                                                                                                                                                                                                                                                                                                                                                //
    Cdc20_T = 0.8332;                                                                                                                                                                                                                                                                                                                                                    //
    ks_20 = 0.005;                                                                                                                                                                                                                                                                                                                                                       //
    ks_20_ = 0.06;                                                                                                                                                                                                                                                                                                                                                       //
    Cdc20 = 0.6848;                                                                                                                                                                                                                                                                                                                                                      //
    ka_20 = 1.0;                                                                                                                                                                                                                                                                                                                                                         //
    ki_20 = 0.1;                                                                                                                                                                                                                                                                                                                                                         //
    ki_20_ = 10.0;                                                                                                                                                                                                                                                                                                                                                       //
    Hct1 = 0.9946;                                                                                                                                                                                                                                                                                                                                                       //
    ka_t1 = 0.04;                                                                                                                                                                                                                                                                                                                                                        //
    ka_t1_ = 2.0;                                                                                                                                                                                                                                                                                                                                                        //
    ki_t1 = 0.0;                                                                                                                                                                                                                                                                                                                                                         //
    ki_t1_ = 0.64;                                                                                                                                                                                                                                                                                                                                                       //
    Ji_t1 = 0.05;                                                                                                                                                                                                                                                                                                                                                        //
    Ja_t1 = 0.05;                                                                                                                                                                                                                                                                                                                                                        //
    epsiloni_t1_n2 = 1.0;                                                                                                                                                                                                                                                                                                                                                //
    epsiloni_t1_b5 = 0.5;                                                                                                                                                                                                                                                                                                                                                //
    epsiloni_t1_b2 = 1.0;                                                                                                                                                                                                                                                                                                                                                //
    mass = 0.6608;                                                                                                                                                                                                                                                                                                                                                       //
    mu = 0.005776;                                                                                                                                                                                                                                                                                                                                                       //
    ORI = 0.0;                                                                                                                                                                                                                                                                                                                                                           //
    ks_ori = 2.0;                                                                                                                                                                                                                                                                                                                                                        //
    kd_ori = 0.06;                                                                                                                                                                                                                                                                                                                                                       //
    epsilonori_b2 = 0.4;                                                                                                                                                                                                                                                                                                                                                 //
    BUD = 0.0;                                                                                                                                                                                                                                                                                                                                                           //
    ks_bud = 0.3;                                                                                                                                                                                                                                                                                                                                                        //
    kd_bud = 0.06;                                                                                                                                                                                                                                                                                                                                                       //
    epsilonbud_b5 = 1.0;                                                                                                                                                                                                                                                                                                                                                 //
    SPN = 0.0;                                                                                                                                                                                                                                                                                                                                                           //
    ks_spn = 0.08;                                                                                                                                                                                                                                                                                                                                                       //
    kd_spn = 0.06;                                                                                                                                                                                                                                                                                                                                                       //
    J_spn = 0.2;                                                                                                                                                                                                                                                                                                                                                         //
    ka_sbf = 1.0;                                                                                                                                                                                                                                                                                                                                                        //
    ki_sbf = 0.5;                                                                                                                                                                                                                                                                                                                                                        //
    ki_sbf_ = 6.0;                                                                                                                                                                                                                                                                                                                                                       //
    Ji_sbf = 0.01;                                                                                                                                                                                                                                                                                                                                                       //
    Ja_sbf = 0.01;                                                                                                                                                                                                                                                                                                                                                       //
    epsilonsbf_n3 = 75.0;                                                                                                                                                                                                                                                                                                                                                //
    epsilonsbf_b5 = 0.5;                                                                                                                                                                                                                                                                                                                                                 //
    ka_mcm = 1.0;                                                                                                                                                                                                                                                                                                                                                        //
    ki_mcm = 0.15;                                                                                                                                                                                                                                                                                                                                                       //
    Ji_mcm = 1.0;                                                                                                                                                                                                                                                                                                                                                        //
    Ja_mcm = 1.0;                                                                                                                                                                                                                                                                                                                                                        //
    ka_swi = 1.0;                                                                                                                                                                                                                                                                                                                                                        //
    ki_swi = 0.3;                                                                                                                                                                                                                                                                                                                                                        //
    ki_swi_ = 0.2;                                                                                                                                                                                                                                                                                                                                                       //
    Ji_swi = 0.1;                                                                                                                                                                                                                                                                                                                                                        //
    Ja_swi = 0.1;                                                                                                                                                                                                                                                                                                                                                        //
    kd1_c1 = 0.01;                                                                                                                                                                                                                                                                                                                                                       //
    kd_20 = 0.08;                                                                                                                                                                                                                                                                                                                                                        //
    Jd2_c1 = 0.05;                                                                                                                                                                                                                                                                                                                                                       //
    Vd_b2 = kd_b2 * (Hct1_T - Hct1) + kd_b2_ * Hct1 + kd_b2__ * Cdc20;                                                                                                                                                                                                                                                                                                   //
    Clb2 = Clb2_T - Clb2_Sic1;                                                                                                                                                                                                                                                                                                                                           //
    Clb5 = Clb5_T - Clb5_Sic1;                                                                                                                                                                                                                                                                                                                                           //
    Sic1 = Sic1_T - (Clb2_Sic1 + Clb5_Sic1);                                                                                                                                                                                                                                                                                                                             //
    Vd_b5 = kd_b5 + kd_b5_ * Cdc20;                                                                                                                                                                                                                                                                                                                                      //
    Bck2 = Bck2_0 * mass;                                                                                                                                                                                                                                                                                                                                                //
    Cln3 = Cln3_max * Dn3 * mass / (Jn3 + Dn3 * mass);                                                                                                                                                                                                                                                                                                                   //
    Vd2_c1 = kd2_c1 * (epsilonc1_n3 * Cln3 + epsilonc1_k2 * Bck2 + Cln2 + epsilonc1_b5 * Clb5 + epsilonc1_b2 * Clb2);                                                                                                                                                                                                                                                    //
    Vi_20 = sm::piecewise(ki_20_, ORI >= 1.0, ki_20, SPN >= 1.0, 0.1);                                                                                                                                                                                                                                                                                                   //
    Vi_t1 = ki_t1 + ki_t1_ * (Cln3 + epsiloni_t1_n2 * Cln2 + epsiloni_t1_b5 * Clb5 + epsiloni_t1_b2 * Clb2);                                                                                                                                                                                                                                                             //
    Va_sbf = ka_sbf * (Cln2 + epsilonsbf_n3 * (Cln3 + Bck2) + epsilonsbf_b5 * Clb5);                                                                                                                                                                                                                                                                                     //
    SBF = 2.0 * Va_sbf * Ji_sbf / ((ki_sbf + ki_sbf_ * Clb2 + Va_sbf * Ji_sbf + (ki_sbf + ki_sbf_ * Clb2) * Ja_sbf - Va_sbf) + std::sqrt(std::pow(ki_sbf + ki_sbf_ * Clb2 + Va_sbf * Ji_sbf + (ki_sbf + ki_sbf_ * Clb2) * Ja_sbf - Va_sbf, 2.0) - 4.0 * Va_sbf * Ji_sbf * (ki_sbf + ki_sbf_ * Clb2 - Va_sbf)));                                                          //
    MBF = SBF;                                                                                                                                                                                                                                                                                                                                                           //
    Mcm1 = 2.0 * ka_mcm * Clb2 * Ji_mcm / ((ki_mcm + ka_mcm * Clb2 * Ji_mcm + ki_mcm * Ja_mcm - ka_mcm * Clb2) + std::sqrt(std::pow(ki_mcm + ka_mcm * Clb2 * Ji_mcm + ki_mcm * Ja_mcm - ka_mcm * Clb2, 2.0) - 4.0 * (ki_mcm - ka_mcm * Clb2) * ka_mcm * Clb2 * Ji_mcm));                                                                                                 //
    Swi5 = 2.0 * ka_swi * Cdc20 * Ji_swi / ((ki_swi + ki_swi_ * Clb2 + ka_swi * Cdc20 * Ji_swi + (ki_swi + ki_swi_ * Clb2) * Ja_swi - ka_swi * Cdc20) + std::sqrt(std::pow(ki_swi + ki_swi_ * Clb2 + ka_swi * Cdc20 * Ji_swi + (ki_swi + ki_swi_ * Clb2) * Ja_swi - ka_swi * Cdc20, 2.0) - 4.0 * (ki_swi + ki_swi_ * Clb2 - ka_swi * Cdc20) * ka_swi * Cdc20 * Ji_swi)); //
    d_Cln2_dt = mass * (ks_n2 + ks_n2_ * SBF) - kd_n2 * Cln2;                                                                                                                                                                                                                                                                                                            //
    d_Clb2_T_dt = mass * (ks_b2 + ks_b2_ * Mcm1) - Vd_b2 * Clb2_T;                                                                                                                                                                                                                                                                                                       //
    d_Clb5_T_dt = mass * (ks_b5 + ks_b5_ * MBF) - Vd_b5 * Clb5_T;                                                                                                                                                                                                                                                                                                        //
    d_Sic1_T_dt = ks_c1 + ks_c1_ * Swi5 - Sic1_T * (kd1_c1 + Vd2_c1 / (Jd2_c1 + Sic1_T));                                                                                                                                                                                                                                                                                //
    d_Clb2_Sic1_dt = kas_b2 * Clb2 * Sic1 - Clb2_Sic1 * (kdi_b2 + Vd_b2 + kd1_c1 + Vd2_c1 / (Jd2_c1 + Sic1_T));                                                                                                                                                                                                                                                          //
    d_Clb5_Sic1_dt = kas_b5 * Clb5 * Sic1 - Clb5_Sic1 * (kdi_b5 + Vd_b5 + kd1_c1 + Vd2_c1 / (Jd2_c1 + Sic1_T));                                                                                                                                                                                                                                                          //
    d_Cdc20_T_dt = ks_20 + ks_20_ * Clb2 - kd_20 * Cdc20_T;                                                                                                                                                                                                                                                                                                              //
    d_Cdc20_dt = ka_20 * (Cdc20_T - Cdc20) - Cdc20 * (Vi_20 + kd_20);                                                                                                                                                                                                                                                                                                    //
    d_Hct1_dt = (ka_t1 + ka_t1_ * Cdc20) * (Hct1_T - Hct1) / (Ja_t1 + Hct1_T - Hct1) - Vi_t1 * Hct1 / (Ji_t1 + Hct1);                                                                                                                                                                                                                                                    //
    d_mass_dt = mu * mass;                                                                                                                                                                                                                                                                                                                                               //
    d_ORI_dt = ks_ori * (Clb5 + epsilonori_b2 * Clb2) - kd_ori * ORI;                                                                                                                                                                                                                                                                                                    //
    d_BUD_dt = ks_bud * (Cln2 + Cln3 + epsilonbud_b5 * Clb5) - kd_bud * BUD;                                                                                                                                                                                                                                                                                             //
    d_SPN_dt = ks_spn * Clb2 / (J_spn + Clb2) - kd_spn * SPN;                                                                                                                                                                                                                                                                                                            //

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

    mParameters.push_back(ks_n2);
    mParameters.push_back(ks_n2_);
    mParameters.push_back(kd_n2);
    mParameters.push_back(Hct1_T);
    mParameters.push_back(ks_b2);
    mParameters.push_back(ks_b2_);
    mParameters.push_back(kd_b2);
    mParameters.push_back(kd_b2_);
    mParameters.push_back(kd_b2__);
    mParameters.push_back(ks_b5);
    mParameters.push_back(ks_b5_);
    mParameters.push_back(kd_b5);
    mParameters.push_back(kd_b5_);
    mParameters.push_back(Bck2_0);
    mParameters.push_back(Jn3);
    mParameters.push_back(Dn3);
    mParameters.push_back(Cln3_max);
    mParameters.push_back(ks_c1);
    mParameters.push_back(ks_c1_);
    mParameters.push_back(kas_b2);
    mParameters.push_back(kdi_b2);
    mParameters.push_back(kas_b5);
    mParameters.push_back(kdi_b5);
    mParameters.push_back(kd2_c1);
    mParameters.push_back(epsilonc1_n3);
    mParameters.push_back(epsilonc1_k2);
    mParameters.push_back(epsilonc1_b5);
    mParameters.push_back(epsilonc1_b2);
    mParameters.push_back(ks_20);
    mParameters.push_back(ks_20_);
    mParameters.push_back(ka_20);
    mParameters.push_back(ki_20);
    mParameters.push_back(ki_20_);
    mParameters.push_back(ka_t1);
    mParameters.push_back(ka_t1_);
    mParameters.push_back(ki_t1);
    mParameters.push_back(ki_t1_);
    mParameters.push_back(Ji_t1);
    mParameters.push_back(Ja_t1);
    mParameters.push_back(epsiloni_t1_n2);
    mParameters.push_back(epsiloni_t1_b5);
    mParameters.push_back(epsiloni_t1_b2);
    mParameters.push_back(mu);
    mParameters.push_back(ks_ori);
    mParameters.push_back(kd_ori);
    mParameters.push_back(epsilonori_b2);
    mParameters.push_back(ks_bud);
    mParameters.push_back(kd_bud);
    mParameters.push_back(epsilonbud_b5);
    mParameters.push_back(ks_spn);
    mParameters.push_back(kd_spn);
    mParameters.push_back(J_spn);
    mParameters.push_back(ka_sbf);
    mParameters.push_back(ki_sbf);
    mParameters.push_back(ki_sbf_);
    mParameters.push_back(Ji_sbf);
    mParameters.push_back(Ja_sbf);
    mParameters.push_back(epsilonsbf_n3);
    mParameters.push_back(epsilonsbf_b5);
    mParameters.push_back(ka_mcm);
    mParameters.push_back(ki_mcm);
    mParameters.push_back(Ji_mcm);
    mParameters.push_back(Ja_mcm);
    mParameters.push_back(ka_swi);
    mParameters.push_back(ki_swi);
    mParameters.push_back(ki_swi_);
    mParameters.push_back(Ji_swi);
    mParameters.push_back(Ja_swi);
    mParameters.push_back(kd1_c1);
    mParameters.push_back(kd_20);
    mParameters.push_back(Jd2_c1);
}

double Chen2000SbmlOdeSystem::ProcessModelEvents(double time, const std::vector<double>& rY)
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

    return max_dist; // Signed distance of the event closest to triggering
}

std::vector<double> Chen2000SbmlOdeSystem::RunModelEquations(double time, const std::vector<double>& rStateVariables)
{
    Cln2 = rStateVariables[0];
    Clb2_T = rStateVariables[1];
    Clb5_T = rStateVariables[2];
    Sic1_T = rStateVariables[3];
    Clb2_Sic1 = rStateVariables[4];
    Clb5_Sic1 = rStateVariables[5];
    Cdc20_T = rStateVariables[6];
    Cdc20 = rStateVariables[7];
    Hct1 = rStateVariables[8];
    mass = rStateVariables[9];
    ORI = rStateVariables[10];
    BUD = rStateVariables[11];
    SPN = rStateVariables[12];

    ks_n2 = GetParameter(0);
    ks_n2_ = GetParameter(1);
    kd_n2 = GetParameter(2);
    Hct1_T = GetParameter(3);
    ks_b2 = GetParameter(4);
    ks_b2_ = GetParameter(5);
    kd_b2 = GetParameter(6);
    kd_b2_ = GetParameter(7);
    kd_b2__ = GetParameter(8);
    ks_b5 = GetParameter(9);
    ks_b5_ = GetParameter(10);
    kd_b5 = GetParameter(11);
    kd_b5_ = GetParameter(12);
    Bck2_0 = GetParameter(13);
    Jn3 = GetParameter(14);
    Dn3 = GetParameter(15);
    Cln3_max = GetParameter(16);
    ks_c1 = GetParameter(17);
    ks_c1_ = GetParameter(18);
    kas_b2 = GetParameter(19);
    kdi_b2 = GetParameter(20);
    kas_b5 = GetParameter(21);
    kdi_b5 = GetParameter(22);
    kd2_c1 = GetParameter(23);
    epsilonc1_n3 = GetParameter(24);
    epsilonc1_k2 = GetParameter(25);
    epsilonc1_b5 = GetParameter(26);
    epsilonc1_b2 = GetParameter(27);
    ks_20 = GetParameter(28);
    ks_20_ = GetParameter(29);
    ka_20 = GetParameter(30);
    ki_20 = GetParameter(31);
    ki_20_ = GetParameter(32);
    ka_t1 = GetParameter(33);
    ka_t1_ = GetParameter(34);
    ki_t1 = GetParameter(35);
    ki_t1_ = GetParameter(36);
    Ji_t1 = GetParameter(37);
    Ja_t1 = GetParameter(38);
    epsiloni_t1_n2 = GetParameter(39);
    epsiloni_t1_b5 = GetParameter(40);
    epsiloni_t1_b2 = GetParameter(41);
    mu = GetParameter(42);
    ks_ori = GetParameter(43);
    kd_ori = GetParameter(44);
    epsilonori_b2 = GetParameter(45);
    ks_bud = GetParameter(46);
    kd_bud = GetParameter(47);
    epsilonbud_b5 = GetParameter(48);
    ks_spn = GetParameter(49);
    kd_spn = GetParameter(50);
    J_spn = GetParameter(51);
    ka_sbf = GetParameter(52);
    ki_sbf = GetParameter(53);
    ki_sbf_ = GetParameter(54);
    Ji_sbf = GetParameter(55);
    Ja_sbf = GetParameter(56);
    epsilonsbf_n3 = GetParameter(57);
    epsilonsbf_b5 = GetParameter(58);
    ka_mcm = GetParameter(59);
    ki_mcm = GetParameter(60);
    Ji_mcm = GetParameter(61);
    Ja_mcm = GetParameter(62);
    ka_swi = GetParameter(63);
    ki_swi = GetParameter(64);
    ki_swi_ = GetParameter(65);
    Ji_swi = GetParameter(66);
    Ja_swi = GetParameter(67);
    kd1_c1 = GetParameter(68);
    kd_20 = GetParameter(69);
    Jd2_c1 = GetParameter(70);

    Vd_b2 = kd_b2 * (Hct1_T - Hct1) + kd_b2_ * Hct1 + kd_b2__ * Cdc20;                                                                                                                                                                                                                                                                                                   //
    Clb2 = Clb2_T - Clb2_Sic1;                                                                                                                                                                                                                                                                                                                                           //
    Clb5 = Clb5_T - Clb5_Sic1;                                                                                                                                                                                                                                                                                                                                           //
    Sic1 = Sic1_T - (Clb2_Sic1 + Clb5_Sic1);                                                                                                                                                                                                                                                                                                                             //
    Vd_b5 = kd_b5 + kd_b5_ * Cdc20;                                                                                                                                                                                                                                                                                                                                      //
    Bck2 = Bck2_0 * mass;                                                                                                                                                                                                                                                                                                                                                //
    Cln3 = Cln3_max * Dn3 * mass / (Jn3 + Dn3 * mass);                                                                                                                                                                                                                                                                                                                   //
    Vd2_c1 = kd2_c1 * (epsilonc1_n3 * Cln3 + epsilonc1_k2 * Bck2 + Cln2 + epsilonc1_b5 * Clb5 + epsilonc1_b2 * Clb2);                                                                                                                                                                                                                                                    //
    Vi_20 = sm::piecewise(ki_20_, ORI >= 1.0, ki_20, SPN >= 1.0, 0.1);                                                                                                                                                                                                                                                                                                   //
    Vi_t1 = ki_t1 + ki_t1_ * (Cln3 + epsiloni_t1_n2 * Cln2 + epsiloni_t1_b5 * Clb5 + epsiloni_t1_b2 * Clb2);                                                                                                                                                                                                                                                             //
    Va_sbf = ka_sbf * (Cln2 + epsilonsbf_n3 * (Cln3 + Bck2) + epsilonsbf_b5 * Clb5);                                                                                                                                                                                                                                                                                     //
    SBF = 2.0 * Va_sbf * Ji_sbf / ((ki_sbf + ki_sbf_ * Clb2 + Va_sbf * Ji_sbf + (ki_sbf + ki_sbf_ * Clb2) * Ja_sbf - Va_sbf) + std::sqrt(std::pow(ki_sbf + ki_sbf_ * Clb2 + Va_sbf * Ji_sbf + (ki_sbf + ki_sbf_ * Clb2) * Ja_sbf - Va_sbf, 2.0) - 4.0 * Va_sbf * Ji_sbf * (ki_sbf + ki_sbf_ * Clb2 - Va_sbf)));                                                          //
    MBF = SBF;                                                                                                                                                                                                                                                                                                                                                           //
    Mcm1 = 2.0 * ka_mcm * Clb2 * Ji_mcm / ((ki_mcm + ka_mcm * Clb2 * Ji_mcm + ki_mcm * Ja_mcm - ka_mcm * Clb2) + std::sqrt(std::pow(ki_mcm + ka_mcm * Clb2 * Ji_mcm + ki_mcm * Ja_mcm - ka_mcm * Clb2, 2.0) - 4.0 * (ki_mcm - ka_mcm * Clb2) * ka_mcm * Clb2 * Ji_mcm));                                                                                                 //
    Swi5 = 2.0 * ka_swi * Cdc20 * Ji_swi / ((ki_swi + ki_swi_ * Clb2 + ka_swi * Cdc20 * Ji_swi + (ki_swi + ki_swi_ * Clb2) * Ja_swi - ka_swi * Cdc20) + std::sqrt(std::pow(ki_swi + ki_swi_ * Clb2 + ka_swi * Cdc20 * Ji_swi + (ki_swi + ki_swi_ * Clb2) * Ja_swi - ka_swi * Cdc20, 2.0) - 4.0 * (ki_swi + ki_swi_ * Clb2 - ka_swi * Cdc20) * ka_swi * Cdc20 * Ji_swi)); //
    d_Cln2_dt = mass * (ks_n2 + ks_n2_ * SBF) - kd_n2 * Cln2;                                                                                                                                                                                                                                                                                                            //
    d_Clb2_T_dt = mass * (ks_b2 + ks_b2_ * Mcm1) - Vd_b2 * Clb2_T;                                                                                                                                                                                                                                                                                                       //
    d_Clb5_T_dt = mass * (ks_b5 + ks_b5_ * MBF) - Vd_b5 * Clb5_T;                                                                                                                                                                                                                                                                                                        //
    d_Sic1_T_dt = ks_c1 + ks_c1_ * Swi5 - Sic1_T * (kd1_c1 + Vd2_c1 / (Jd2_c1 + Sic1_T));                                                                                                                                                                                                                                                                                //
    d_Clb2_Sic1_dt = kas_b2 * Clb2 * Sic1 - Clb2_Sic1 * (kdi_b2 + Vd_b2 + kd1_c1 + Vd2_c1 / (Jd2_c1 + Sic1_T));                                                                                                                                                                                                                                                          //
    d_Clb5_Sic1_dt = kas_b5 * Clb5 * Sic1 - Clb5_Sic1 * (kdi_b5 + Vd_b5 + kd1_c1 + Vd2_c1 / (Jd2_c1 + Sic1_T));                                                                                                                                                                                                                                                          //
    d_Cdc20_T_dt = ks_20 + ks_20_ * Clb2 - kd_20 * Cdc20_T;                                                                                                                                                                                                                                                                                                              //
    d_Cdc20_dt = ka_20 * (Cdc20_T - Cdc20) - Cdc20 * (Vi_20 + kd_20);                                                                                                                                                                                                                                                                                                    //
    d_Hct1_dt = (ka_t1 + ka_t1_ * Cdc20) * (Hct1_T - Hct1) / (Ja_t1 + Hct1_T - Hct1) - Vi_t1 * Hct1 / (Ji_t1 + Hct1);                                                                                                                                                                                                                                                    //
    d_mass_dt = mu * mass;                                                                                                                                                                                                                                                                                                                                               //
    d_ORI_dt = ks_ori * (Clb5 + epsilonori_b2 * Clb2) - kd_ori * ORI;                                                                                                                                                                                                                                                                                                    //
    d_BUD_dt = ks_bud * (Cln2 + Cln3 + epsilonbud_b5 * Clb5) - kd_bud * BUD;                                                                                                                                                                                                                                                                                             //
    d_SPN_dt = ks_spn * Clb2 / (J_spn + Clb2) - kd_spn * SPN;                                                                                                                                                                                                                                                                                                            //

    std::vector<double> derivatives(13);
    derivatives[0] = d_Cln2_dt;
    derivatives[1] = d_Clb2_T_dt;
    derivatives[2] = d_Clb5_T_dt;
    derivatives[3] = d_Sic1_T_dt;
    derivatives[4] = d_Clb2_Sic1_dt;
    derivatives[5] = d_Clb5_Sic1_dt;
    derivatives[6] = d_Cdc20_T_dt;
    derivatives[7] = d_Cdc20_dt;
    derivatives[8] = d_Hct1_dt;
    derivatives[9] = d_mass_dt;
    derivatives[10] = d_ORI_dt;
    derivatives[11] = d_BUD_dt;
    derivatives[12] = d_SPN_dt;
    return derivatives;
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
    this->mDerivedQuantityNames.push_back("COMpartment");
    this->mDerivedQuantityUnits.push_back("litre");

    this->mDerivedQuantityNames.push_back("Vd_b2");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Clb2");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Clb5");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Sic1");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Vd_b5");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Bck2");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Cln3");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Vd2_c1");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Vi_20");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Vi_t1");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("SBF");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Va_sbf");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("MBF");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Mcm1");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Swi5");
    this->mDerivedQuantityUnits.push_back("non-dim");

    // PARAMETERS
    this->mParameterNames.push_back("ks_n2");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ks_n2_");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kd_n2");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Hct1_T");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ks_b2");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ks_b2_");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kd_b2");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kd_b2_");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kd_b2__");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ks_b5");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ks_b5_");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kd_b5");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kd_b5_");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Bck2_0");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Jn3");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Dn3");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Cln3_max");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ks_c1");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ks_c1_");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kas_b2");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kdi_b2");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kas_b5");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kdi_b5");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kd2_c1");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("epsilonc1_n3");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("epsilonc1_k2");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("epsilonc1_b5");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("epsilonc1_b2");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ks_20");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ks_20_");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ka_20");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ki_20");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ki_20_");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ka_t1");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ka_t1_");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ki_t1");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ki_t1_");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Ji_t1");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Ja_t1");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("epsiloni_t1_n2");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("epsiloni_t1_b5");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("epsiloni_t1_b2");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("mu");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ks_ori");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kd_ori");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("epsilonori_b2");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ks_bud");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kd_bud");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("epsilonbud_b5");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ks_spn");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kd_spn");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("J_spn");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ka_sbf");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ki_sbf");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ki_sbf_");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Ji_sbf");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Ja_sbf");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("epsilonsbf_n3");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("epsilonsbf_b5");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ka_mcm");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ki_mcm");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Ji_mcm");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Ja_mcm");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ka_swi");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ki_swi");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ki_swi_");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Ji_swi");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Ja_swi");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kd1_c1");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kd_20");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Jd2_c1");
    this->mParameterUnits.push_back("non-dim");

    this->mInitialised = true;
}

// Register the ODE system with Boost serialization
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(Chen2000SbmlOdeSystem)