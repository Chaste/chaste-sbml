#include <cmath>
#include <limits>
#include <vector>

#include "CellwiseOdeSystemInformation.hpp"
#include "SbmlEventType.hpp"
#include "SbmlMath.hpp"

#include "VanLeeuwen2007NonDimSbmlOdeSystem.hpp"

namespace sm = sbmlmath;

VanLeeuwen2007NonDimSbmlOdeSystem::VanLeeuwen2007NonDimSbmlOdeSystem()
        : AbstractSbmlOdeSystem(11, 31, 0)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<VanLeeuwen2007NonDimSbmlOdeSystem>);

    Initialise();

    // EVENTS
}

VanLeeuwen2007NonDimSbmlOdeSystem::~VanLeeuwen2007NonDimSbmlOdeSystem()
{
}

std::vector<double> VanLeeuwen2007NonDimSbmlOdeSystem::ComputeDerivedQuantities(double time, const std::vector<double>& rY)
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

void VanLeeuwen2007NonDimSbmlOdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    std::vector<double> derivatives = RunModelEquations(time, rY);
    for (unsigned i = 0; i < rDY.size(); ++i)
    {
        rDY[i] = derivatives[i];
    }

    // TODO: Scale time appropriately
}

void VanLeeuwen2007NonDimSbmlOdeSystem::Initialise(double time)
{
    cytosolmembraneandnucleus = 1.0;                                                                                                                                                                                                                                                                                                                                                     //
    X = 0.268;                                                                                                                                                                                                                                                                                                                                                                           //
    D = 2.68;                                                                                                                                                                                                                                                                                                                                                                            //
    C_o = 76.2;                                                                                                                                                                                                                                                                                                                                                                          //
    C_u = 13.5;                                                                                                                                                                                                                                                                                                                                                                          //
    C_c = 0.0;                                                                                                                                                                                                                                                                                                                                                                           //
    A = 300.0;                                                                                                                                                                                                                                                                                                                                                                           //
    C_A = 544.2;                                                                                                                                                                                                                                                                                                                                                                         //
    T = 750.0;                                                                                                                                                                                                                                                                                                                                                                           //
    C_oT = 76.2;                                                                                                                                                                                                                                                                                                                                                                         //
    C_cT = 0.0;                                                                                                                                                                                                                                                                                                                                                                          //
    Y = 36.0;                                                                                                                                                                                                                                                                                                                                                                            //
    C_F = 76.2;                                                                                                                                                                                                                                                                                                                                                                          //
    C_T = 76.2;                                                                                                                                                                                                                                                                                                                                                                          //
    drag = 1.0;                                                                                                                                                                                                                                                                                                                                                                          //
    K_T = 1500.0;                                                                                                                                                                                                                                                                                                                                                                        //
    K_C = 6000.0;                                                                                                                                                                                                                                                                                                                                                                        //
    K_D = 150.0;                                                                                                                                                                                                                                                                                                                                                                         //
    p_u = 1.0;                                                                                                                                                                                                                                                                                                                                                                           //
    wnt_level = 0.0;                                                                                                                                                                                                                                                                                                                                                                     //
    xi_D = 0.00666666666666667;                                                                                                                                                                                                                                                                                                                                                          //
    xi_Dx = 0.00666666666666667;                                                                                                                                                                                                                                                                                                                                                         //
    xi_X = 0.266666666666667;                                                                                                                                                                                                                                                                                                                                                            //
    xi_C = 0.0;                                                                                                                                                                                                                                                                                                                                                                          //
    d_Dx = 0.00666666666666667;                                                                                                                                                                                                                                                                                                                                                          //
    d_X = 0.133333333333333;                                                                                                                                                                                                                                                                                                                                                             //
    d_D = 0.00666666666666667;                                                                                                                                                                                                                                                                                                                                                           //
    p_c = 0.0;                                                                                                                                                                                                                                                                                                                                                                           //
    gamma1 = 1.0;                                                                                                                                                                                                                                                                                                                                                                        //
    gamma2 = 1.0;                                                                                                                                                                                                                                                                                                                                                                        //
    s_D = 0.133333333333333;                                                                                                                                                                                                                                                                                                                                                             //
    s_X = 0.0533333333333333;                                                                                                                                                                                                                                                                                                                                                            //
    d_u = 0.0666666666666667;                                                                                                                                                                                                                                                                                                                                                            //
    s_c = 1.0;                                                                                                                                                                                                                                                                                                                                                                           //
    d_c = 0.00133333333333333;                                                                                                                                                                                                                                                                                                                                                           //
    s_CA = 0.0111111111111111;                                                                                                                                                                                                                                                                                                                                                           //
    d_CA = 0.466666666666667;                                                                                                                                                                                                                                                                                                                                                            //
    s_CT = 0.00133333333333333;                                                                                                                                                                                                                                                                                                                                                          //
    d_CT = 1.0;                                                                                                                                                                                                                                                                                                                                                                          //
    s_A = 0.8;                                                                                                                                                                                                                                                                                                                                                                           //
    d_A = 0.00266666666666667;                                                                                                                                                                                                                                                                                                                                                           //
    s_T = 0.4;                                                                                                                                                                                                                                                                                                                                                                           //
    d_T = 0.000533333333333333;                                                                                                                                                                                                                                                                                                                                                          //
    s_Y = 1.0;                                                                                                                                                                                                                                                                                                                                                                           //
    d_Y = 0.00133333333333333;                                                                                                                                                                                                                                                                                                                                                           //
    ComplexTransitThreshold = 1.0;                                                                                                                                                                                                                                                                                                                                                       //
    X = X / cytosolmembraneandnucleus;                                                                                                                                                                                                                                                                                                                                                   //
    D = D / cytosolmembraneandnucleus;                                                                                                                                                                                                                                                                                                                                                   //
    C_o = C_o / cytosolmembraneandnucleus;                                                                                                                                                                                                                                                                                                                                               //
    C_u = C_u / cytosolmembraneandnucleus;                                                                                                                                                                                                                                                                                                                                               //
    C_c = C_c / cytosolmembraneandnucleus;                                                                                                                                                                                                                                                                                                                                               //
    A = A / cytosolmembraneandnucleus;                                                                                                                                                                                                                                                                                                                                                   //
    C_A = C_A / cytosolmembraneandnucleus;                                                                                                                                                                                                                                                                                                                                               //
    T = T / cytosolmembraneandnucleus;                                                                                                                                                                                                                                                                                                                                                   //
    C_oT = C_oT / cytosolmembraneandnucleus;                                                                                                                                                                                                                                                                                                                                             //
    C_cT = C_cT / cytosolmembraneandnucleus;                                                                                                                                                                                                                                                                                                                                             //
    Y = Y / cytosolmembraneandnucleus;                                                                                                                                                                                                                                                                                                                                                   //
    C_F = C_o + C_c;                                                                                                                                                                                                                                                                                                                                                                     //
    C_T = C_oT + C_cT;                                                                                                                                                                                                                                                                                                                                                                   //
    drag = sm::piecewise((C_A - 2300.0) / 36.0, ((C_A - 2300.0) / 36.0) > 1.0, 1.0);                                                                                                                                                                                                                                                                                                     //
    mwcfbf7716_cc13_473c_979a_033c57a28857 = s_D * gamma1 * X;                                                                                                                                                                                                                                                                                                                           //
    mwab347951_a286_432d_b03b_254dcdba4a25 = (d_Dx + wnt_level * xi_Dx) * D;                                                                                                                                                                                                                                                                                                             //
    mwc497befc_0edd_4b0d_8895_77dbfa05a4be = s_X;                                                                                                                                                                                                                                                                                                                                        //
    mwfa0f5940_d4e9_44e1_8a0c_379d6fdf9d0f = (d_X + wnt_level * xi_X) * X;                                                                                                                                                                                                                                                                                                               //
    mwbe108cd1_a5d9_4d55_97c9_ac584df6a005 = p_u * gamma2 * C_o * D / (C_o + C_c + K_D);                                                                                                                                                                                                                                                                                                 //
    mwad0ba91a_37fb_48f5_8392_6c043dfdd5dc = p_u * gamma2 * C_c * D / (C_c + C_o + K_D);                                                                                                                                                                                                                                                                                                 //
    mw1d3d836b_77bd_489b_9a37_acc874344652 = d_u * C_u;                                                                                                                                                                                                                                                                                                                                  //
    mw552840aa_6d3c_4188_850d_9a3caafcdaa5 = s_c;                                                                                                                                                                                                                                                                                                                                        //
    mw5a301d6f_ffa1_4362_98de_0919a68808bc = d_c * C_o;                                                                                                                                                                                                                                                                                                                                  //
    mw4b073704_f41a_43d1_8ff8_af03f1cb6e1d = d_c * C_c;                                                                                                                                                                                                                                                                                                                                  //
    mw1c7459d0_b5a9_4a89_8682_79fc8ca4ca6c = s_CA * C_o * A;                                                                                                                                                                                                                                                                                                                             //
    mw51e5ed29_82db_47fa_9eed_ab52b7786dd5 = d_CA * C_A;                                                                                                                                                                                                                                                                                                                                 //
    mw8704b9b0_f93d_405a_89a1_e000c5b66896 = s_CT * C_o * T;                                                                                                                                                                                                                                                                                                                             //
    mw322e8f78_68f6_4293_ae2a_ab238ec5cf9b = s_CT * C_c * T;                                                                                                                                                                                                                                                                                                                             //
    mwf5440ce3_586c_491f_ba71_bbf1c485027f = d_CT * C_oT;                                                                                                                                                                                                                                                                                                                                //
    mw08c332d6_1aaf_498f_8e1c_13f1b9da85bf = d_CT * C_cT;                                                                                                                                                                                                                                                                                                                                //
    mw3ece1442_a0b3_47eb_a8bf_cb317a46165f = (p_c + wnt_level * xi_C) * C_o / (C_o + K_C);                                                                                                                                                                                                                                                                                               //
    mwfb6f8dfa_2e07_4249_a8fe_af33ca27471e = s_A;                                                                                                                                                                                                                                                                                                                                        //
    mwaa4e7692_3290_45ab_8a38_b5ffd49ede87 = d_A * A;                                                                                                                                                                                                                                                                                                                                    //
    mw86ea3c28_4745_4a89_9cbd_0f23c603e01b = s_T;                                                                                                                                                                                                                                                                                                                                        //
    mwee62535f_931d_41e9_ad53_2cad54b94778 = d_T * T;                                                                                                                                                                                                                                                                                                                                    //
    mw1b250e04_290a_4689_9fc5_5f0bf6711e02 = s_Y * (C_oT + C_cT) / (C_oT + C_cT + K_T);                                                                                                                                                                                                                                                                                                  //
    mw925599eb_19a0_4434_8be3_67c40721b71d = d_Y * Y;                                                                                                                                                                                                                                                                                                                                    //
    mw321b3e5e_f6ed_4345_9346_55ffb1ff2137 = (d_D + wnt_level * xi_D) * D;                                                                                                                                                                                                                                                                                                               //
    d_X_dt = (-mwcfbf7716_cc13_473c_979a_033c57a28857 + mwab347951_a286_432d_b03b_254dcdba4a25 + mwc497befc_0edd_4b0d_8895_77dbfa05a4be - mwfa0f5940_d4e9_44e1_8a0c_379d6fdf9d0f) / cytosolmembraneandnucleus;                                                                                                                                                                           //
    d_D_dt = (((mwcfbf7716_cc13_473c_979a_033c57a28857 - mwab347951_a286_432d_b03b_254dcdba4a25) + mwbe108cd1_a5d9_4d55_97c9_ac584df6a005 - mwbe108cd1_a5d9_4d55_97c9_ac584df6a005) + mwad0ba91a_37fb_48f5_8392_6c043dfdd5dc - mwad0ba91a_37fb_48f5_8392_6c043dfdd5dc - mw321b3e5e_f6ed_4345_9346_55ffb1ff2137) / cytosolmembraneandnucleus;                                             //
    d_C_o_dt = (((-mwbe108cd1_a5d9_4d55_97c9_ac584df6a005 + mw552840aa_6d3c_4188_850d_9a3caafcdaa5 - mw5a301d6f_ffa1_4362_98de_0919a68808bc - mw1c7459d0_b5a9_4a89_8682_79fc8ca4ca6c) + mw51e5ed29_82db_47fa_9eed_ab52b7786dd5 - mw8704b9b0_f93d_405a_89a1_e000c5b66896) + mwf5440ce3_586c_491f_ba71_bbf1c485027f - mw3ece1442_a0b3_47eb_a8bf_cb317a46165f) / cytosolmembraneandnucleus; //
    d_C_u_dt = (mwbe108cd1_a5d9_4d55_97c9_ac584df6a005 + mwad0ba91a_37fb_48f5_8392_6c043dfdd5dc - mw1d3d836b_77bd_489b_9a37_acc874344652) / cytosolmembraneandnucleus;                                                                                                                                                                                                                   //
    d_C_c_dt = ((-mwad0ba91a_37fb_48f5_8392_6c043dfdd5dc - mw4b073704_f41a_43d1_8ff8_af03f1cb6e1d - mw322e8f78_68f6_4293_ae2a_ab238ec5cf9b) + mw08c332d6_1aaf_498f_8e1c_13f1b9da85bf + mw3ece1442_a0b3_47eb_a8bf_cb317a46165f) / cytosolmembraneandnucleus;                                                                                                                              //
    d_A_dt = (-mw1c7459d0_b5a9_4a89_8682_79fc8ca4ca6c + mw51e5ed29_82db_47fa_9eed_ab52b7786dd5 + mwfb6f8dfa_2e07_4249_a8fe_af33ca27471e - mwaa4e7692_3290_45ab_8a38_b5ffd49ede87) / cytosolmembraneandnucleus;                                                                                                                                                                           //
    d_C_A_dt = (mw1c7459d0_b5a9_4a89_8682_79fc8ca4ca6c - mw51e5ed29_82db_47fa_9eed_ab52b7786dd5) / cytosolmembraneandnucleus;                                                                                                                                                                                                                                                            //
    d_T_dt = ((-mw8704b9b0_f93d_405a_89a1_e000c5b66896 - mw322e8f78_68f6_4293_ae2a_ab238ec5cf9b) + mwf5440ce3_586c_491f_ba71_bbf1c485027f + mw08c332d6_1aaf_498f_8e1c_13f1b9da85bf + mw86ea3c28_4745_4a89_9cbd_0f23c603e01b - mwee62535f_931d_41e9_ad53_2cad54b94778) / cytosolmembraneandnucleus;                                                                                       //
    d_C_oT_dt = ((mw8704b9b0_f93d_405a_89a1_e000c5b66896 - mwf5440ce3_586c_491f_ba71_bbf1c485027f) + mw1b250e04_290a_4689_9fc5_5f0bf6711e02 - mw1b250e04_290a_4689_9fc5_5f0bf6711e02) / cytosolmembraneandnucleus;                                                                                                                                                                       //
    d_C_cT_dt = (mw322e8f78_68f6_4293_ae2a_ab238ec5cf9b - mw08c332d6_1aaf_498f_8e1c_13f1b9da85bf) / cytosolmembraneandnucleus;                                                                                                                                                                                                                                                           //
    d_Y_dt = (mw1b250e04_290a_4689_9fc5_5f0bf6711e02 - mw925599eb_19a0_4434_8be3_67c40721b71d) / cytosolmembraneandnucleus;                                                                                                                                                                                                                                                              //

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
    mParameters.push_back(xi_D);
    mParameters.push_back(xi_Dx);
    mParameters.push_back(xi_X);
    mParameters.push_back(xi_C);
    mParameters.push_back(d_Dx);
    mParameters.push_back(d_X);
    mParameters.push_back(d_D);
    mParameters.push_back(p_c);
    mParameters.push_back(gamma1);
    mParameters.push_back(gamma2);
    mParameters.push_back(s_D);
    mParameters.push_back(s_X);
    mParameters.push_back(d_u);
    mParameters.push_back(s_c);
    mParameters.push_back(d_c);
    mParameters.push_back(s_CA);
    mParameters.push_back(d_CA);
    mParameters.push_back(s_CT);
    mParameters.push_back(d_CT);
    mParameters.push_back(s_A);
    mParameters.push_back(d_A);
    mParameters.push_back(s_T);
    mParameters.push_back(d_T);
    mParameters.push_back(s_Y);
    mParameters.push_back(d_Y);
    mParameters.push_back(ComplexTransitThreshold);
}

double VanLeeuwen2007NonDimSbmlOdeSystem::ProcessModelEvents(double time, const std::vector<double>& rY)
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

// ASSIGNMENT RULES
void VanLeeuwen2007NonDimSbmlOdeSystem::RunAssignmentRules(double time)
{
}

// INITIAL ASSIGNMENTS
void VanLeeuwen2007NonDimSbmlOdeSystem::RunInitialAssignments(double time)
{
}

std::vector<double> VanLeeuwen2007NonDimSbmlOdeSystem::RunModelEquations(double time, const std::vector<double>& rStateVariables)
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
    xi_D = GetParameter(5);
    xi_Dx = GetParameter(6);
    xi_X = GetParameter(7);
    xi_C = GetParameter(8);
    d_Dx = GetParameter(9);
    d_X = GetParameter(10);
    d_D = GetParameter(11);
    p_c = GetParameter(12);
    gamma1 = GetParameter(13);
    gamma2 = GetParameter(14);
    s_D = GetParameter(15);
    s_X = GetParameter(16);
    d_u = GetParameter(17);
    s_c = GetParameter(18);
    d_c = GetParameter(19);
    s_CA = GetParameter(20);
    d_CA = GetParameter(21);
    s_CT = GetParameter(22);
    d_CT = GetParameter(23);
    s_A = GetParameter(24);
    d_A = GetParameter(25);
    s_T = GetParameter(26);
    d_T = GetParameter(27);
    s_Y = GetParameter(28);
    d_Y = GetParameter(29);
    ComplexTransitThreshold = GetParameter(30);

    C_F = C_o + C_c;                                                                                                                                                                                                                                                                                                                                                                     //
    C_T = C_oT + C_cT;                                                                                                                                                                                                                                                                                                                                                                   //
    drag = sm::piecewise((C_A - 2300.0) / 36.0, ((C_A - 2300.0) / 36.0) > 1.0, 1.0);                                                                                                                                                                                                                                                                                                     //
    mwcfbf7716_cc13_473c_979a_033c57a28857 = s_D * gamma1 * X;                                                                                                                                                                                                                                                                                                                           //
    mwab347951_a286_432d_b03b_254dcdba4a25 = (d_Dx + wnt_level * xi_Dx) * D;                                                                                                                                                                                                                                                                                                             //
    mwc497befc_0edd_4b0d_8895_77dbfa05a4be = s_X;                                                                                                                                                                                                                                                                                                                                        //
    mwfa0f5940_d4e9_44e1_8a0c_379d6fdf9d0f = (d_X + wnt_level * xi_X) * X;                                                                                                                                                                                                                                                                                                               //
    mwbe108cd1_a5d9_4d55_97c9_ac584df6a005 = p_u * gamma2 * C_o * D / (C_o + C_c + K_D);                                                                                                                                                                                                                                                                                                 //
    mwad0ba91a_37fb_48f5_8392_6c043dfdd5dc = p_u * gamma2 * C_c * D / (C_c + C_o + K_D);                                                                                                                                                                                                                                                                                                 //
    mw1d3d836b_77bd_489b_9a37_acc874344652 = d_u * C_u;                                                                                                                                                                                                                                                                                                                                  //
    mw552840aa_6d3c_4188_850d_9a3caafcdaa5 = s_c;                                                                                                                                                                                                                                                                                                                                        //
    mw5a301d6f_ffa1_4362_98de_0919a68808bc = d_c * C_o;                                                                                                                                                                                                                                                                                                                                  //
    mw4b073704_f41a_43d1_8ff8_af03f1cb6e1d = d_c * C_c;                                                                                                                                                                                                                                                                                                                                  //
    mw1c7459d0_b5a9_4a89_8682_79fc8ca4ca6c = s_CA * C_o * A;                                                                                                                                                                                                                                                                                                                             //
    mw51e5ed29_82db_47fa_9eed_ab52b7786dd5 = d_CA * C_A;                                                                                                                                                                                                                                                                                                                                 //
    mw8704b9b0_f93d_405a_89a1_e000c5b66896 = s_CT * C_o * T;                                                                                                                                                                                                                                                                                                                             //
    mw322e8f78_68f6_4293_ae2a_ab238ec5cf9b = s_CT * C_c * T;                                                                                                                                                                                                                                                                                                                             //
    mwf5440ce3_586c_491f_ba71_bbf1c485027f = d_CT * C_oT;                                                                                                                                                                                                                                                                                                                                //
    mw08c332d6_1aaf_498f_8e1c_13f1b9da85bf = d_CT * C_cT;                                                                                                                                                                                                                                                                                                                                //
    mw3ece1442_a0b3_47eb_a8bf_cb317a46165f = (p_c + wnt_level * xi_C) * C_o / (C_o + K_C);                                                                                                                                                                                                                                                                                               //
    mwfb6f8dfa_2e07_4249_a8fe_af33ca27471e = s_A;                                                                                                                                                                                                                                                                                                                                        //
    mwaa4e7692_3290_45ab_8a38_b5ffd49ede87 = d_A * A;                                                                                                                                                                                                                                                                                                                                    //
    mw86ea3c28_4745_4a89_9cbd_0f23c603e01b = s_T;                                                                                                                                                                                                                                                                                                                                        //
    mwee62535f_931d_41e9_ad53_2cad54b94778 = d_T * T;                                                                                                                                                                                                                                                                                                                                    //
    mw1b250e04_290a_4689_9fc5_5f0bf6711e02 = s_Y * (C_oT + C_cT) / (C_oT + C_cT + K_T);                                                                                                                                                                                                                                                                                                  //
    mw925599eb_19a0_4434_8be3_67c40721b71d = d_Y * Y;                                                                                                                                                                                                                                                                                                                                    //
    mw321b3e5e_f6ed_4345_9346_55ffb1ff2137 = (d_D + wnt_level * xi_D) * D;                                                                                                                                                                                                                                                                                                               //
    d_X_dt = (-mwcfbf7716_cc13_473c_979a_033c57a28857 + mwab347951_a286_432d_b03b_254dcdba4a25 + mwc497befc_0edd_4b0d_8895_77dbfa05a4be - mwfa0f5940_d4e9_44e1_8a0c_379d6fdf9d0f) / cytosolmembraneandnucleus;                                                                                                                                                                           //
    d_D_dt = (((mwcfbf7716_cc13_473c_979a_033c57a28857 - mwab347951_a286_432d_b03b_254dcdba4a25) + mwbe108cd1_a5d9_4d55_97c9_ac584df6a005 - mwbe108cd1_a5d9_4d55_97c9_ac584df6a005) + mwad0ba91a_37fb_48f5_8392_6c043dfdd5dc - mwad0ba91a_37fb_48f5_8392_6c043dfdd5dc - mw321b3e5e_f6ed_4345_9346_55ffb1ff2137) / cytosolmembraneandnucleus;                                             //
    d_C_o_dt = (((-mwbe108cd1_a5d9_4d55_97c9_ac584df6a005 + mw552840aa_6d3c_4188_850d_9a3caafcdaa5 - mw5a301d6f_ffa1_4362_98de_0919a68808bc - mw1c7459d0_b5a9_4a89_8682_79fc8ca4ca6c) + mw51e5ed29_82db_47fa_9eed_ab52b7786dd5 - mw8704b9b0_f93d_405a_89a1_e000c5b66896) + mwf5440ce3_586c_491f_ba71_bbf1c485027f - mw3ece1442_a0b3_47eb_a8bf_cb317a46165f) / cytosolmembraneandnucleus; //
    d_C_u_dt = (mwbe108cd1_a5d9_4d55_97c9_ac584df6a005 + mwad0ba91a_37fb_48f5_8392_6c043dfdd5dc - mw1d3d836b_77bd_489b_9a37_acc874344652) / cytosolmembraneandnucleus;                                                                                                                                                                                                                   //
    d_C_c_dt = ((-mwad0ba91a_37fb_48f5_8392_6c043dfdd5dc - mw4b073704_f41a_43d1_8ff8_af03f1cb6e1d - mw322e8f78_68f6_4293_ae2a_ab238ec5cf9b) + mw08c332d6_1aaf_498f_8e1c_13f1b9da85bf + mw3ece1442_a0b3_47eb_a8bf_cb317a46165f) / cytosolmembraneandnucleus;                                                                                                                              //
    d_A_dt = (-mw1c7459d0_b5a9_4a89_8682_79fc8ca4ca6c + mw51e5ed29_82db_47fa_9eed_ab52b7786dd5 + mwfb6f8dfa_2e07_4249_a8fe_af33ca27471e - mwaa4e7692_3290_45ab_8a38_b5ffd49ede87) / cytosolmembraneandnucleus;                                                                                                                                                                           //
    d_C_A_dt = (mw1c7459d0_b5a9_4a89_8682_79fc8ca4ca6c - mw51e5ed29_82db_47fa_9eed_ab52b7786dd5) / cytosolmembraneandnucleus;                                                                                                                                                                                                                                                            //
    d_T_dt = ((-mw8704b9b0_f93d_405a_89a1_e000c5b66896 - mw322e8f78_68f6_4293_ae2a_ab238ec5cf9b) + mwf5440ce3_586c_491f_ba71_bbf1c485027f + mw08c332d6_1aaf_498f_8e1c_13f1b9da85bf + mw86ea3c28_4745_4a89_9cbd_0f23c603e01b - mwee62535f_931d_41e9_ad53_2cad54b94778) / cytosolmembraneandnucleus;                                                                                       //
    d_C_oT_dt = ((mw8704b9b0_f93d_405a_89a1_e000c5b66896 - mwf5440ce3_586c_491f_ba71_bbf1c485027f) + mw1b250e04_290a_4689_9fc5_5f0bf6711e02 - mw1b250e04_290a_4689_9fc5_5f0bf6711e02) / cytosolmembraneandnucleus;                                                                                                                                                                       //
    d_C_cT_dt = (mw322e8f78_68f6_4293_ae2a_ab238ec5cf9b - mw08c332d6_1aaf_498f_8e1c_13f1b9da85bf) / cytosolmembraneandnucleus;                                                                                                                                                                                                                                                           //
    d_Y_dt = (mw1b250e04_290a_4689_9fc5_5f0bf6711e02 - mw925599eb_19a0_4434_8be3_67c40721b71d) / cytosolmembraneandnucleus;                                                                                                                                                                                                                                                              //

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
void VanLeeuwen2007NonDimSbmlOdeSystem::RunReactions(double time)
{
}

// VARIABLE PARAMETERS
void VanLeeuwen2007NonDimSbmlOdeSystem::UpdateParameters(double time)
{
}

// STATE VARIABLES
void VanLeeuwen2007NonDimSbmlOdeSystem::UpdateStateVariables(double time, const std::vector<double>& rStateVariables)
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
void CellwiseOdeSystemInformation<VanLeeuwen2007NonDimSbmlOdeSystem>::Initialise()
{
    // STATE VARIABLES
    this->mVariableNames.push_back("X");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.268);

    this->mVariableNames.push_back("D");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(2.68);

    this->mVariableNames.push_back("C_o");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(76.2);

    this->mVariableNames.push_back("C_u");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(13.5);

    this->mVariableNames.push_back("C_c");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("A");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(300.0);

    this->mVariableNames.push_back("C_A");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(544.2);

    this->mVariableNames.push_back("T");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(750.0);

    this->mVariableNames.push_back("C_oT");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(76.2);

    this->mVariableNames.push_back("C_cT");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("Y");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(36.0);

    // DERIVED QUANTITIES
    this->mDerivedQuantityNames.push_back("cytosolmembraneandnucleus");
    this->mDerivedQuantityUnits.push_back("non-dim");

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
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("K_C");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("K_D");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("p_u");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("wnt_level");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("xi_D");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("xi_Dx");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("xi_X");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("xi_C");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("d_Dx");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("d_X");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("d_D");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("p_c");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("gamma1");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("gamma2");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("s_D");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("s_X");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("d_u");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("s_c");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("d_c");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("s_CA");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("d_CA");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("s_CT");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("d_CT");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("s_A");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("d_A");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("s_T");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("d_T");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("s_Y");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("d_Y");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ComplexTransitThreshold");
    this->mParameterUnits.push_back("non-dim");

    this->mInitialised = true;
}

// Register the ODE system with Boost serialization
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(VanLeeuwen2007NonDimSbmlOdeSystem)