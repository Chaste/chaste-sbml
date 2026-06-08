#include <cmath>
#include <limits>
#include <vector>

#include "CellwiseOdeSystemInformation.hpp"
#include "SbmlEventType.hpp"
#include "SbmlMath.hpp"

#include "VanLeeuwen2007SbmlOdeSystem.hpp"

namespace sm = sbmlmath;

VanLeeuwen2007SbmlOdeSystem::VanLeeuwen2007SbmlOdeSystem()
        : AbstractSbmlOdeSystem(11, 31, 0)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<VanLeeuwen2007SbmlOdeSystem>);

    Initialise();

    // EVENTS
}

VanLeeuwen2007SbmlOdeSystem::~VanLeeuwen2007SbmlOdeSystem()
{
}

std::vector<double> VanLeeuwen2007SbmlOdeSystem::ComputeDerivedQuantities(double time, const std::vector<double>& rY)
{
    std::vector<double> dqs;
    RunModelEquations(time, rY);

    // AMOUNTS
    double amt__X = X * cytosolmembraneandnucleus;       //
    double amt__D = D * cytosolmembraneandnucleus;       //
    double amt__C_o = C_o * cytosolmembraneandnucleus;   //
    double amt__C_u = C_u * cytosolmembraneandnucleus;   //
    double amt__C_c = C_c * cytosolmembraneandnucleus;   //
    double amt__A = A * cytosolmembraneandnucleus;       //
    double amt__C_A = C_A * cytosolmembraneandnucleus;   //
    double amt__T = T * cytosolmembraneandnucleus;       //
    double amt__C_oT = C_oT * cytosolmembraneandnucleus; //
    double amt__C_cT = C_cT * cytosolmembraneandnucleus; //
    double amt__Y = Y * cytosolmembraneandnucleus;       //
    double amt__C_F = C_F * cytosolmembraneandnucleus;   //
    double amt__C_T = C_T * cytosolmembraneandnucleus;   //
    double amt__drag = drag * cytosolmembraneandnucleus; //

    dqs.push_back(cytosolmembraneandnucleus);
    dqs.push_back(amt__X);
    dqs.push_back(amt__D);
    dqs.push_back(amt__C_o);
    dqs.push_back(amt__C_u);
    dqs.push_back(amt__C_c);
    dqs.push_back(amt__A);
    dqs.push_back(amt__C_A);
    dqs.push_back(amt__T);
    dqs.push_back(amt__C_oT);
    dqs.push_back(amt__C_cT);
    dqs.push_back(amt__Y);
    dqs.push_back(C_F);
    dqs.push_back(amt__C_F);
    dqs.push_back(C_T);
    dqs.push_back(amt__C_T);
    dqs.push_back(drag);
    dqs.push_back(amt__drag);

    return dqs;
}

void VanLeeuwen2007SbmlOdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    std::vector<double> derivatives = RunModelEquations(time, rY);
    for (unsigned i = 0; i < rDY.size(); ++i)
    {
        rDY[i] = derivatives[i];
    }

    // TODO: Scale time appropriately
}

void VanLeeuwen2007SbmlOdeSystem::Initialise(double time)
{
    cytosolmembraneandnucleus = 1.0;                                                                                                                                                                                                                                                                                                                                                     //
    X = 0.067;                                                                                                                                                                                                                                                                                                                                                                           //
    D = 0.67;                                                                                                                                                                                                                                                                                                                                                                            //
    C_o = 2.54;                                                                                                                                                                                                                                                                                                                                                                          //
    C_u = 0.45;                                                                                                                                                                                                                                                                                                                                                                          //
    C_c = 0.0;                                                                                                                                                                                                                                                                                                                                                                           //
    A = 10.0;                                                                                                                                                                                                                                                                                                                                                                            //
    C_A = 18.14;                                                                                                                                                                                                                                                                                                                                                                         //
    T = 25.0;                                                                                                                                                                                                                                                                                                                                                                            //
    C_oT = 2.54;                                                                                                                                                                                                                                                                                                                                                                         //
    C_cT = 0.0;                                                                                                                                                                                                                                                                                                                                                                          //
    Y = 0.48;                                                                                                                                                                                                                                                                                                                                                                            //
    C_F = 2.54;                                                                                                                                                                                                                                                                                                                                                                          //
    C_T = 2.54;                                                                                                                                                                                                                                                                                                                                                                          //
    drag = 1.0;                                                                                                                                                                                                                                                                                                                                                                          //
    K_T = 50.0;                                                                                                                                                                                                                                                                                                                                                                          //
    K_C = 200.0;                                                                                                                                                                                                                                                                                                                                                                         //
    K_D = 5.0;                                                                                                                                                                                                                                                                                                                                                                           //
    p_u = 100.0;                                                                                                                                                                                                                                                                                                                                                                         //
    wnt_level = 0.0;                                                                                                                                                                                                                                                                                                                                                                     //
    gamma1 = 1.0;                                                                                                                                                                                                                                                                                                                                                                        //
    gamma2 = 1.0;                                                                                                                                                                                                                                                                                                                                                                        //
    xi_D = 5.0;                                                                                                                                                                                                                                                                                                                                                                          //
    xi_Dx = 5.0;                                                                                                                                                                                                                                                                                                                                                                         //
    xi_X = 200.0;                                                                                                                                                                                                                                                                                                                                                                        //
    xi_C = 0.0;                                                                                                                                                                                                                                                                                                                                                                          //
    s_D = 100.0;                                                                                                                                                                                                                                                                                                                                                                         //
    d_Dx = 5.0;                                                                                                                                                                                                                                                                                                                                                                          //
    s_X = 10.0;                                                                                                                                                                                                                                                                                                                                                                          //
    d_X = 100.0;                                                                                                                                                                                                                                                                                                                                                                         //
    d_u = 50.0;                                                                                                                                                                                                                                                                                                                                                                          //
    s_c = 25.0;                                                                                                                                                                                                                                                                                                                                                                          //
    d_c = 1.0;                                                                                                                                                                                                                                                                                                                                                                           //
    s_CA = 250.0;                                                                                                                                                                                                                                                                                                                                                                        //
    d_CA = 350.0;                                                                                                                                                                                                                                                                                                                                                                        //
    s_CT = 30.0;                                                                                                                                                                                                                                                                                                                                                                         //
    d_CT = 750.0;                                                                                                                                                                                                                                                                                                                                                                        //
    p_c = 0.0;                                                                                                                                                                                                                                                                                                                                                                           //
    s_A = 20.0;                                                                                                                                                                                                                                                                                                                                                                          //
    d_A = 2.0;                                                                                                                                                                                                                                                                                                                                                                           //
    s_T = 10.0;                                                                                                                                                                                                                                                                                                                                                                          //
    d_T = 0.4;                                                                                                                                                                                                                                                                                                                                                                           //
    s_Y = 10.0;                                                                                                                                                                                                                                                                                                                                                                          //
    d_Y = 1.0;                                                                                                                                                                                                                                                                                                                                                                           //
    d_D = 5.0;                                                                                                                                                                                                                                                                                                                                                                           //
    ComplexTransitThreshold = 1.0;                                                                                                                                                                                                                                                                                                                                                       //
    C_F = C_o + C_c;                                                                                                                                                                                                                                                                                                                                                                     //
    C_T = C_oT + C_cT;                                                                                                                                                                                                                                                                                                                                                                   //
    drag = sm::piecewise((C_A - 100.0) / 3.0, ((C_A - 100.0) / 3.0) >= 1.0, 1.0);                                                                                                                                                                                                                                                                                                        //
    mwd6b35759_f098_484c_9c65_e84e7e4b61e4 = s_D * gamma1 * X;                                                                                                                                                                                                                                                                                                                           //
    mweddac6d0_231e_4c92_ba2a_c91edc682ff5 = (d_Dx + wnt_level * xi_Dx) * D;                                                                                                                                                                                                                                                                                                             //
    mwdf62dfed_ec88_4d81_bc9d_da0e10f41e4b = s_X * cytosolmembraneandnucleus;                                                                                                                                                                                                                                                                                                            //
    mwee9cc998_28e9_4173_a694_f3e278a639b7 = (d_X + wnt_level * xi_X) * X;                                                                                                                                                                                                                                                                                                               //
    mw661e341d_97d1_4e6f_8812_3be7ffc86d42 = p_u * gamma2 * C_o * D / (C_o + C_c + K_D);                                                                                                                                                                                                                                                                                                 //
    mw179aa33c_9a7e_43c0_9285_3d8f97719c60 = p_u * gamma2 * C_c * D / (C_c + C_o + K_D);                                                                                                                                                                                                                                                                                                 //
    mwff8d34f9_e036_49f1_b3b8_3706ecb98660 = d_u * C_u * cytosolmembraneandnucleus;                                                                                                                                                                                                                                                                                                      //
    mwd4fa317a_a484_4f38_b9b0_6ed404d9adcb = s_c * cytosolmembraneandnucleus;                                                                                                                                                                                                                                                                                                            //
    mwcb88a249_a200_4e95_9185_5654bf1ebfc0 = d_c * C_o * cytosolmembraneandnucleus;                                                                                                                                                                                                                                                                                                      //
    mw0aeac2fc_dd2f_4fe7_b60c_97f0b008ed89 = d_c * C_c * cytosolmembraneandnucleus;                                                                                                                                                                                                                                                                                                      //
    mw4b47c66d_37e6_4c33_b043_1f6b3b814449 = s_CA * C_o * A * cytosolmembraneandnucleus;                                                                                                                                                                                                                                                                                                 //
    mw69974db4_8ead_416c_a220_f6dc3be1f3b6 = d_CA * C_A * cytosolmembraneandnucleus;                                                                                                                                                                                                                                                                                                     //
    mwccb628c3_76c8_47be_9d6b_6d3f1ae30dcc = s_CT * C_o * T * cytosolmembraneandnucleus;                                                                                                                                                                                                                                                                                                 //
    mw581d69f1_60b3_4d21_9323_31b05ee89570 = s_CT * C_c * T * cytosolmembraneandnucleus;                                                                                                                                                                                                                                                                                                 //
    mwb17c2c57_279d_4e88_b9cf_896029135cc1 = d_CT * C_oT * cytosolmembraneandnucleus;                                                                                                                                                                                                                                                                                                    //
    mwe3236fc5_2118_40cb_8db3_ef9da29137cf = d_CT * C_cT * cytosolmembraneandnucleus;                                                                                                                                                                                                                                                                                                    //
    mw0be4a28b_e9c6_43da_8f95_d9c564a7caae = (p_c + wnt_level * xi_C) * C_o / (C_o + K_C);                                                                                                                                                                                                                                                                                               //
    mwc360befb_07da_4d19_bbec_523fbef47dc9 = s_A * cytosolmembraneandnucleus;                                                                                                                                                                                                                                                                                                            //
    mw985b39ad_50ff_4d55_95f7_42e7f1bf6a3d = d_A * A * cytosolmembraneandnucleus;                                                                                                                                                                                                                                                                                                        //
    mwf27b65bb_38d6_44cf_aa66_9b38ae5885eb = s_T * cytosolmembraneandnucleus;                                                                                                                                                                                                                                                                                                            //
    mwd7a7bdbb_8f2e_42ed_bf97_48e33e45dafb = d_T * T * cytosolmembraneandnucleus;                                                                                                                                                                                                                                                                                                        //
    mw988a8caf_bd68_462b_86d7_51844c1dcfd3 = s_Y * (C_oT + C_cT) / (C_oT + C_cT + K_T);                                                                                                                                                                                                                                                                                                  //
    mw9ab26a4c_bd70_45e0_bacc_f830ab28abca = d_Y * Y * cytosolmembraneandnucleus;                                                                                                                                                                                                                                                                                                        //
    mw931baf8f_6572_46f6_96eb_cae40ee267b7 = (d_D + wnt_level * xi_D) * D;                                                                                                                                                                                                                                                                                                               //
    d_X_dt = (-mwd6b35759_f098_484c_9c65_e84e7e4b61e4 + mweddac6d0_231e_4c92_ba2a_c91edc682ff5 + mwdf62dfed_ec88_4d81_bc9d_da0e10f41e4b - mwee9cc998_28e9_4173_a694_f3e278a639b7) / cytosolmembraneandnucleus;                                                                                                                                                                           //
    d_D_dt = (((mwd6b35759_f098_484c_9c65_e84e7e4b61e4 - mweddac6d0_231e_4c92_ba2a_c91edc682ff5) + mw661e341d_97d1_4e6f_8812_3be7ffc86d42 - mw661e341d_97d1_4e6f_8812_3be7ffc86d42) + mw179aa33c_9a7e_43c0_9285_3d8f97719c60 - mw179aa33c_9a7e_43c0_9285_3d8f97719c60 - mw931baf8f_6572_46f6_96eb_cae40ee267b7) / cytosolmembraneandnucleus;                                             //
    d_C_o_dt = (((-mw661e341d_97d1_4e6f_8812_3be7ffc86d42 + mwd4fa317a_a484_4f38_b9b0_6ed404d9adcb - mwcb88a249_a200_4e95_9185_5654bf1ebfc0 - mw4b47c66d_37e6_4c33_b043_1f6b3b814449) + mw69974db4_8ead_416c_a220_f6dc3be1f3b6 - mwccb628c3_76c8_47be_9d6b_6d3f1ae30dcc) + mwb17c2c57_279d_4e88_b9cf_896029135cc1 - mw0be4a28b_e9c6_43da_8f95_d9c564a7caae) / cytosolmembraneandnucleus; //
    d_C_u_dt = (mw661e341d_97d1_4e6f_8812_3be7ffc86d42 + mw179aa33c_9a7e_43c0_9285_3d8f97719c60 - mwff8d34f9_e036_49f1_b3b8_3706ecb98660) / cytosolmembraneandnucleus;                                                                                                                                                                                                                   //
    d_C_c_dt = ((-mw179aa33c_9a7e_43c0_9285_3d8f97719c60 - mw0aeac2fc_dd2f_4fe7_b60c_97f0b008ed89 - mw581d69f1_60b3_4d21_9323_31b05ee89570) + mwe3236fc5_2118_40cb_8db3_ef9da29137cf + mw0be4a28b_e9c6_43da_8f95_d9c564a7caae) / cytosolmembraneandnucleus;                                                                                                                              //
    d_A_dt = (-mw4b47c66d_37e6_4c33_b043_1f6b3b814449 + mw69974db4_8ead_416c_a220_f6dc3be1f3b6 + mwc360befb_07da_4d19_bbec_523fbef47dc9 - mw985b39ad_50ff_4d55_95f7_42e7f1bf6a3d) / cytosolmembraneandnucleus;                                                                                                                                                                           //
    d_C_A_dt = (mw4b47c66d_37e6_4c33_b043_1f6b3b814449 - mw69974db4_8ead_416c_a220_f6dc3be1f3b6) / cytosolmembraneandnucleus;                                                                                                                                                                                                                                                            //
    d_T_dt = ((-mwccb628c3_76c8_47be_9d6b_6d3f1ae30dcc - mw581d69f1_60b3_4d21_9323_31b05ee89570) + mwb17c2c57_279d_4e88_b9cf_896029135cc1 + mwe3236fc5_2118_40cb_8db3_ef9da29137cf + mwf27b65bb_38d6_44cf_aa66_9b38ae5885eb - mwd7a7bdbb_8f2e_42ed_bf97_48e33e45dafb) / cytosolmembraneandnucleus;                                                                                       //
    d_C_oT_dt = ((mwccb628c3_76c8_47be_9d6b_6d3f1ae30dcc - mwb17c2c57_279d_4e88_b9cf_896029135cc1) + mw988a8caf_bd68_462b_86d7_51844c1dcfd3 - mw988a8caf_bd68_462b_86d7_51844c1dcfd3) / cytosolmembraneandnucleus;                                                                                                                                                                       //
    d_C_cT_dt = (mw581d69f1_60b3_4d21_9323_31b05ee89570 - mwe3236fc5_2118_40cb_8db3_ef9da29137cf) / cytosolmembraneandnucleus;                                                                                                                                                                                                                                                           //
    d_Y_dt = (mw988a8caf_bd68_462b_86d7_51844c1dcfd3 - mw9ab26a4c_bd70_45e0_bacc_f830ab28abca) / cytosolmembraneandnucleus;                                                                                                                                                                                                                                                              //

    mStateVariables.push_back(X);
    mStateVariables.push_back(D);
    mStateVariables.push_back(C_o);
    mStateVariables.push_back(C_u);
    mStateVariables.push_back(C_c);
    mStateVariables.push_back(A);
    mStateVariables.push_back(C_A);
    mStateVariables.push_back(T);
    mStateVariables.push_back(C_oT);
    mStateVariables.push_back(C_cT);
    mStateVariables.push_back(Y);

    SetDefaultInitialCondition(0, X);
    SetDefaultInitialCondition(1, D);
    SetDefaultInitialCondition(2, C_o);
    SetDefaultInitialCondition(3, C_u);
    SetDefaultInitialCondition(4, C_c);
    SetDefaultInitialCondition(5, A);
    SetDefaultInitialCondition(6, C_A);
    SetDefaultInitialCondition(7, T);
    SetDefaultInitialCondition(8, C_oT);
    SetDefaultInitialCondition(9, C_cT);
    SetDefaultInitialCondition(10, Y);

    mParameters.push_back(K_T);
    mParameters.push_back(K_C);
    mParameters.push_back(K_D);
    mParameters.push_back(p_u);
    mParameters.push_back(wnt_level);
    mParameters.push_back(gamma1);
    mParameters.push_back(gamma2);
    mParameters.push_back(xi_D);
    mParameters.push_back(xi_Dx);
    mParameters.push_back(xi_X);
    mParameters.push_back(xi_C);
    mParameters.push_back(s_D);
    mParameters.push_back(d_Dx);
    mParameters.push_back(s_X);
    mParameters.push_back(d_X);
    mParameters.push_back(d_u);
    mParameters.push_back(s_c);
    mParameters.push_back(d_c);
    mParameters.push_back(s_CA);
    mParameters.push_back(d_CA);
    mParameters.push_back(s_CT);
    mParameters.push_back(d_CT);
    mParameters.push_back(p_c);
    mParameters.push_back(s_A);
    mParameters.push_back(d_A);
    mParameters.push_back(s_T);
    mParameters.push_back(d_T);
    mParameters.push_back(s_Y);
    mParameters.push_back(d_Y);
    mParameters.push_back(d_D);
    mParameters.push_back(ComplexTransitThreshold);
}

double VanLeeuwen2007SbmlOdeSystem::ProcessModelEvents(double time, const std::vector<double>& rY)
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

    double min_dist = std::numeric_limits<double>::max();

    return min_dist; // Distance to closest event
}

// ASSIGNMENT RULES
void VanLeeuwen2007SbmlOdeSystem::RunAssignmentRules(double time)
{
}

// INITIAL ASSIGNMENTS
void VanLeeuwen2007SbmlOdeSystem::RunInitialAssignments(double time)
{
}

std::vector<double> VanLeeuwen2007SbmlOdeSystem::RunModelEquations(double time, const std::vector<double>& rStateVariables)
{
    X = rStateVariables[0];
    D = rStateVariables[1];
    C_o = rStateVariables[2];
    C_u = rStateVariables[3];
    C_c = rStateVariables[4];
    A = rStateVariables[5];
    C_A = rStateVariables[6];
    T = rStateVariables[7];
    C_oT = rStateVariables[8];
    C_cT = rStateVariables[9];
    Y = rStateVariables[10];

    K_T = GetParameter(0);
    K_C = GetParameter(1);
    K_D = GetParameter(2);
    p_u = GetParameter(3);
    wnt_level = GetParameter(4);
    gamma1 = GetParameter(5);
    gamma2 = GetParameter(6);
    xi_D = GetParameter(7);
    xi_Dx = GetParameter(8);
    xi_X = GetParameter(9);
    xi_C = GetParameter(10);
    s_D = GetParameter(11);
    d_Dx = GetParameter(12);
    s_X = GetParameter(13);
    d_X = GetParameter(14);
    d_u = GetParameter(15);
    s_c = GetParameter(16);
    d_c = GetParameter(17);
    s_CA = GetParameter(18);
    d_CA = GetParameter(19);
    s_CT = GetParameter(20);
    d_CT = GetParameter(21);
    p_c = GetParameter(22);
    s_A = GetParameter(23);
    d_A = GetParameter(24);
    s_T = GetParameter(25);
    d_T = GetParameter(26);
    s_Y = GetParameter(27);
    d_Y = GetParameter(28);
    d_D = GetParameter(29);
    ComplexTransitThreshold = GetParameter(30);

    C_F = C_o + C_c;                                                                                                                                                                                                                                                                                                                                                                     //
    C_T = C_oT + C_cT;                                                                                                                                                                                                                                                                                                                                                                   //
    drag = sm::piecewise((C_A - 100.0) / 3.0, ((C_A - 100.0) / 3.0) >= 1.0, 1.0);                                                                                                                                                                                                                                                                                                        //
    mwd6b35759_f098_484c_9c65_e84e7e4b61e4 = s_D * gamma1 * X;                                                                                                                                                                                                                                                                                                                           //
    mweddac6d0_231e_4c92_ba2a_c91edc682ff5 = (d_Dx + wnt_level * xi_Dx) * D;                                                                                                                                                                                                                                                                                                             //
    mwdf62dfed_ec88_4d81_bc9d_da0e10f41e4b = s_X * cytosolmembraneandnucleus;                                                                                                                                                                                                                                                                                                            //
    mwee9cc998_28e9_4173_a694_f3e278a639b7 = (d_X + wnt_level * xi_X) * X;                                                                                                                                                                                                                                                                                                               //
    mw661e341d_97d1_4e6f_8812_3be7ffc86d42 = p_u * gamma2 * C_o * D / (C_o + C_c + K_D);                                                                                                                                                                                                                                                                                                 //
    mw179aa33c_9a7e_43c0_9285_3d8f97719c60 = p_u * gamma2 * C_c * D / (C_c + C_o + K_D);                                                                                                                                                                                                                                                                                                 //
    mwff8d34f9_e036_49f1_b3b8_3706ecb98660 = d_u * C_u * cytosolmembraneandnucleus;                                                                                                                                                                                                                                                                                                      //
    mwd4fa317a_a484_4f38_b9b0_6ed404d9adcb = s_c * cytosolmembraneandnucleus;                                                                                                                                                                                                                                                                                                            //
    mwcb88a249_a200_4e95_9185_5654bf1ebfc0 = d_c * C_o * cytosolmembraneandnucleus;                                                                                                                                                                                                                                                                                                      //
    mw0aeac2fc_dd2f_4fe7_b60c_97f0b008ed89 = d_c * C_c * cytosolmembraneandnucleus;                                                                                                                                                                                                                                                                                                      //
    mw4b47c66d_37e6_4c33_b043_1f6b3b814449 = s_CA * C_o * A * cytosolmembraneandnucleus;                                                                                                                                                                                                                                                                                                 //
    mw69974db4_8ead_416c_a220_f6dc3be1f3b6 = d_CA * C_A * cytosolmembraneandnucleus;                                                                                                                                                                                                                                                                                                     //
    mwccb628c3_76c8_47be_9d6b_6d3f1ae30dcc = s_CT * C_o * T * cytosolmembraneandnucleus;                                                                                                                                                                                                                                                                                                 //
    mw581d69f1_60b3_4d21_9323_31b05ee89570 = s_CT * C_c * T * cytosolmembraneandnucleus;                                                                                                                                                                                                                                                                                                 //
    mwb17c2c57_279d_4e88_b9cf_896029135cc1 = d_CT * C_oT * cytosolmembraneandnucleus;                                                                                                                                                                                                                                                                                                    //
    mwe3236fc5_2118_40cb_8db3_ef9da29137cf = d_CT * C_cT * cytosolmembraneandnucleus;                                                                                                                                                                                                                                                                                                    //
    mw0be4a28b_e9c6_43da_8f95_d9c564a7caae = (p_c + wnt_level * xi_C) * C_o / (C_o + K_C);                                                                                                                                                                                                                                                                                               //
    mwc360befb_07da_4d19_bbec_523fbef47dc9 = s_A * cytosolmembraneandnucleus;                                                                                                                                                                                                                                                                                                            //
    mw985b39ad_50ff_4d55_95f7_42e7f1bf6a3d = d_A * A * cytosolmembraneandnucleus;                                                                                                                                                                                                                                                                                                        //
    mwf27b65bb_38d6_44cf_aa66_9b38ae5885eb = s_T * cytosolmembraneandnucleus;                                                                                                                                                                                                                                                                                                            //
    mwd7a7bdbb_8f2e_42ed_bf97_48e33e45dafb = d_T * T * cytosolmembraneandnucleus;                                                                                                                                                                                                                                                                                                        //
    mw988a8caf_bd68_462b_86d7_51844c1dcfd3 = s_Y * (C_oT + C_cT) / (C_oT + C_cT + K_T);                                                                                                                                                                                                                                                                                                  //
    mw9ab26a4c_bd70_45e0_bacc_f830ab28abca = d_Y * Y * cytosolmembraneandnucleus;                                                                                                                                                                                                                                                                                                        //
    mw931baf8f_6572_46f6_96eb_cae40ee267b7 = (d_D + wnt_level * xi_D) * D;                                                                                                                                                                                                                                                                                                               //
    d_X_dt = (-mwd6b35759_f098_484c_9c65_e84e7e4b61e4 + mweddac6d0_231e_4c92_ba2a_c91edc682ff5 + mwdf62dfed_ec88_4d81_bc9d_da0e10f41e4b - mwee9cc998_28e9_4173_a694_f3e278a639b7) / cytosolmembraneandnucleus;                                                                                                                                                                           //
    d_D_dt = (((mwd6b35759_f098_484c_9c65_e84e7e4b61e4 - mweddac6d0_231e_4c92_ba2a_c91edc682ff5) + mw661e341d_97d1_4e6f_8812_3be7ffc86d42 - mw661e341d_97d1_4e6f_8812_3be7ffc86d42) + mw179aa33c_9a7e_43c0_9285_3d8f97719c60 - mw179aa33c_9a7e_43c0_9285_3d8f97719c60 - mw931baf8f_6572_46f6_96eb_cae40ee267b7) / cytosolmembraneandnucleus;                                             //
    d_C_o_dt = (((-mw661e341d_97d1_4e6f_8812_3be7ffc86d42 + mwd4fa317a_a484_4f38_b9b0_6ed404d9adcb - mwcb88a249_a200_4e95_9185_5654bf1ebfc0 - mw4b47c66d_37e6_4c33_b043_1f6b3b814449) + mw69974db4_8ead_416c_a220_f6dc3be1f3b6 - mwccb628c3_76c8_47be_9d6b_6d3f1ae30dcc) + mwb17c2c57_279d_4e88_b9cf_896029135cc1 - mw0be4a28b_e9c6_43da_8f95_d9c564a7caae) / cytosolmembraneandnucleus; //
    d_C_u_dt = (mw661e341d_97d1_4e6f_8812_3be7ffc86d42 + mw179aa33c_9a7e_43c0_9285_3d8f97719c60 - mwff8d34f9_e036_49f1_b3b8_3706ecb98660) / cytosolmembraneandnucleus;                                                                                                                                                                                                                   //
    d_C_c_dt = ((-mw179aa33c_9a7e_43c0_9285_3d8f97719c60 - mw0aeac2fc_dd2f_4fe7_b60c_97f0b008ed89 - mw581d69f1_60b3_4d21_9323_31b05ee89570) + mwe3236fc5_2118_40cb_8db3_ef9da29137cf + mw0be4a28b_e9c6_43da_8f95_d9c564a7caae) / cytosolmembraneandnucleus;                                                                                                                              //
    d_A_dt = (-mw4b47c66d_37e6_4c33_b043_1f6b3b814449 + mw69974db4_8ead_416c_a220_f6dc3be1f3b6 + mwc360befb_07da_4d19_bbec_523fbef47dc9 - mw985b39ad_50ff_4d55_95f7_42e7f1bf6a3d) / cytosolmembraneandnucleus;                                                                                                                                                                           //
    d_C_A_dt = (mw4b47c66d_37e6_4c33_b043_1f6b3b814449 - mw69974db4_8ead_416c_a220_f6dc3be1f3b6) / cytosolmembraneandnucleus;                                                                                                                                                                                                                                                            //
    d_T_dt = ((-mwccb628c3_76c8_47be_9d6b_6d3f1ae30dcc - mw581d69f1_60b3_4d21_9323_31b05ee89570) + mwb17c2c57_279d_4e88_b9cf_896029135cc1 + mwe3236fc5_2118_40cb_8db3_ef9da29137cf + mwf27b65bb_38d6_44cf_aa66_9b38ae5885eb - mwd7a7bdbb_8f2e_42ed_bf97_48e33e45dafb) / cytosolmembraneandnucleus;                                                                                       //
    d_C_oT_dt = ((mwccb628c3_76c8_47be_9d6b_6d3f1ae30dcc - mwb17c2c57_279d_4e88_b9cf_896029135cc1) + mw988a8caf_bd68_462b_86d7_51844c1dcfd3 - mw988a8caf_bd68_462b_86d7_51844c1dcfd3) / cytosolmembraneandnucleus;                                                                                                                                                                       //
    d_C_cT_dt = (mw581d69f1_60b3_4d21_9323_31b05ee89570 - mwe3236fc5_2118_40cb_8db3_ef9da29137cf) / cytosolmembraneandnucleus;                                                                                                                                                                                                                                                           //
    d_Y_dt = (mw988a8caf_bd68_462b_86d7_51844c1dcfd3 - mw9ab26a4c_bd70_45e0_bacc_f830ab28abca) / cytosolmembraneandnucleus;                                                                                                                                                                                                                                                              //

    std::vector<double> derivatives(11);
    derivatives[0] = d_X_dt;
    derivatives[1] = d_D_dt;
    derivatives[2] = d_C_o_dt;
    derivatives[3] = d_C_u_dt;
    derivatives[4] = d_C_c_dt;
    derivatives[5] = d_A_dt;
    derivatives[6] = d_C_A_dt;
    derivatives[7] = d_T_dt;
    derivatives[8] = d_C_oT_dt;
    derivatives[9] = d_C_cT_dt;
    derivatives[10] = d_Y_dt;
    return derivatives;
}

// REACTIONS
void VanLeeuwen2007SbmlOdeSystem::RunReactions(double time)
{
}

// VARIABLE PARAMETERS
void VanLeeuwen2007SbmlOdeSystem::UpdateParameters(double time)
{
}

// STATE VARIABLES
void VanLeeuwen2007SbmlOdeSystem::UpdateStateVariables(double time, const std::vector<double>& rStateVariables)
{
    X = rStateVariables[0];
    D = rStateVariables[1];
    C_o = rStateVariables[2];
    C_u = rStateVariables[3];
    C_c = rStateVariables[4];
    A = rStateVariables[5];
    C_A = rStateVariables[6];
    T = rStateVariables[7];
    C_oT = rStateVariables[8];
    C_cT = rStateVariables[9];
    Y = rStateVariables[10];
}

// MODEL FUNCTIONS

template <>
void CellwiseOdeSystemInformation<VanLeeuwen2007SbmlOdeSystem>::Initialise()
{
    // STATE VARIABLES
    this->mVariableNames.push_back("X");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.067);

    this->mVariableNames.push_back("D");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.67);

    this->mVariableNames.push_back("C_o");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(2.54);

    this->mVariableNames.push_back("C_u");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.45);

    this->mVariableNames.push_back("C_c");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("A");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(10.0);

    this->mVariableNames.push_back("C_A");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(18.14);

    this->mVariableNames.push_back("T");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(25.0);

    this->mVariableNames.push_back("C_oT");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(2.54);

    this->mVariableNames.push_back("C_cT");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("Y");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.48);

    // DERIVED QUANTITIES
    this->mDerivedQuantityNames.push_back("cytosolmembraneandnucleus");
    this->mDerivedQuantityUnits.push_back("MWDERIVEDUNIT_meter___3");

    this->mDerivedQuantityNames.push_back("amt__X");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("amt__D");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("amt__C_o");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("amt__C_u");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("amt__C_c");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("amt__A");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("amt__C_A");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("amt__T");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("amt__C_oT");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("amt__C_cT");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("amt__Y");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("C_F");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("amt__C_F");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("C_T");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("amt__C_T");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("drag");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("amt__drag");
    this->mDerivedQuantityUnits.push_back("non-dim");

    // PARAMETERS
    this->mParameterNames.push_back("K_T");
    this->mParameterUnits.push_back("MWDERIVEDUNIT_nanomole__meter___3");

    this->mParameterNames.push_back("K_C");
    this->mParameterUnits.push_back("MWDERIVEDUNIT_nanomole__meter___3");

    this->mParameterNames.push_back("K_D");
    this->mParameterUnits.push_back("MWDERIVEDUNIT_nanomole__meter___3");

    this->mParameterNames.push_back("p_u");
    this->mParameterUnits.push_back("MWDERIVEDUNIT_1__hour");

    this->mParameterNames.push_back("wnt_level");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("gamma1");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("gamma2");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("xi_D");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("xi_Dx");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("xi_X");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("xi_C");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("s_D");
    this->mParameterUnits.push_back("MWDERIVEDUNIT_1__hour");

    this->mParameterNames.push_back("d_Dx");
    this->mParameterUnits.push_back("MWDERIVEDUNIT_1__hour");

    this->mParameterNames.push_back("s_X");
    this->mParameterUnits.push_back("MWDERIVEDUNIT_nanomole__meter___3_hour");

    this->mParameterNames.push_back("d_X");
    this->mParameterUnits.push_back("MWDERIVEDUNIT_1__hour");

    this->mParameterNames.push_back("d_u");
    this->mParameterUnits.push_back("MWDERIVEDUNIT_1__hour");

    this->mParameterNames.push_back("s_c");
    this->mParameterUnits.push_back("MWDERIVEDUNIT_nanomole__meter___3_hour");

    this->mParameterNames.push_back("d_c");
    this->mParameterUnits.push_back("MWDERIVEDUNIT_1__hour");

    this->mParameterNames.push_back("s_CA");
    this->mParameterUnits.push_back("MWDERIVEDUNIT_1__nanomole__meter___3_hour");

    this->mParameterNames.push_back("d_CA");
    this->mParameterUnits.push_back("MWDERIVEDUNIT_1__hour");

    this->mParameterNames.push_back("s_CT");
    this->mParameterUnits.push_back("MWDERIVEDUNIT_1__nanomole__meter___3_hour");

    this->mParameterNames.push_back("d_CT");
    this->mParameterUnits.push_back("MWDERIVEDUNIT_1__hour");

    this->mParameterNames.push_back("p_c");
    this->mParameterUnits.push_back("MWDERIVEDUNIT_1__hour");

    this->mParameterNames.push_back("s_A");
    this->mParameterUnits.push_back("MWDERIVEDUNIT_nanomole__meter___3_hour");

    this->mParameterNames.push_back("d_A");
    this->mParameterUnits.push_back("MWDERIVEDUNIT_1__hour");

    this->mParameterNames.push_back("s_T");
    this->mParameterUnits.push_back("MWDERIVEDUNIT_nanomole__meter___3_hour");

    this->mParameterNames.push_back("d_T");
    this->mParameterUnits.push_back("MWDERIVEDUNIT_1__hour");

    this->mParameterNames.push_back("s_Y");
    this->mParameterUnits.push_back("MWDERIVEDUNIT_1__hour");

    this->mParameterNames.push_back("d_Y");
    this->mParameterUnits.push_back("MWDERIVEDUNIT_1__hour");

    this->mParameterNames.push_back("d_D");
    this->mParameterUnits.push_back("MWDERIVEDUNIT_1__hour");

    this->mParameterNames.push_back("ComplexTransitThreshold");
    this->mParameterUnits.push_back("non-dim");

    this->mInitialised = true;
}

// Register the ODE system with Boost serialization
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(VanLeeuwen2007SbmlOdeSystem)