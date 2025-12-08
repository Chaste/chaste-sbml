#include <cmath>
#include <limits>
#include <vector>

#include "CellwiseOdeSystemInformation.hpp"
#include "SbmlEventType.hpp"
#include "SbmlMath.hpp"

#include "VanLeeuwen2007NonDimSbmlOdeSystem.hpp"

namespace sm = sbmlmath;

VanLeeuwen2007NonDimSbmlOdeSystem::VanLeeuwen2007NonDimSbmlOdeSystem()
        : AbstractSbmlOdeSystem(11, 4, 0)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<VanLeeuwen2007NonDimSbmlOdeSystem>);

    // VARIABLE PARAMETERS
    wnt_level = 0.0;
    gamma1 = 1.0;
    gamma2 = 1.0;
    ComplexTransitThreshold = 1.0;

    // STATE VARIABLES
    X = 0.268;
    D = 2.68;
    C_o = 76.2;
    C_u = 13.5;
    C_c = 0.0;
    A = 300.0;
    C_A = 544.2;
    T = 750.0;
    C_oT = 76.2;
    C_cT = 0.0;
    Y = 36.0;

    // DERIVED QUANTITIES
    cytosolmembraneandnucleus = 1.0;
    C_F = 76.2;
    C_T = 76.2;
    drag = 1.0;

    // INITIAL ASSIGNMENTS
    C_F = C_o + C_c;                            //
    C_T = C_oT + C_cT;                          //
    drag = sm::max((C_A - 2300.0) / 36.0, 1.0); //

    X = X / cytosolmembraneandnucleus;       // Convert X amount to concentration
    D = D / cytosolmembraneandnucleus;       // Convert D amount to concentration
    C_o = C_o / cytosolmembraneandnucleus;   // Convert C_o amount to concentration
    C_u = C_u / cytosolmembraneandnucleus;   // Convert C_u amount to concentration
    C_c = C_c / cytosolmembraneandnucleus;   // Convert C_c amount to concentration
    A = A / cytosolmembraneandnucleus;       // Convert A amount to concentration
    C_A = C_A / cytosolmembraneandnucleus;   // Convert C_A amount to concentration
    T = T / cytosolmembraneandnucleus;       // Convert T amount to concentration
    C_oT = C_oT / cytosolmembraneandnucleus; // Convert C_oT amount to concentration
    C_cT = C_cT / cytosolmembraneandnucleus; // Convert C_cT amount to concentration
    Y = Y / cytosolmembraneandnucleus;       // Convert Y amount to concentration

    // ODE SYSTEM INFORMATION
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

    mParameters.push_back(wnt_level);
    mParameters.push_back(gamma1);
    mParameters.push_back(gamma2);
    mParameters.push_back(ComplexTransitThreshold);

    // REACTIONS
    RunReactions(0.0);

    // EVENTS
}

VanLeeuwen2007NonDimSbmlOdeSystem::~VanLeeuwen2007NonDimSbmlOdeSystem()
{
}

std::vector<double> VanLeeuwen2007NonDimSbmlOdeSystem::ComputeDerivedQuantities(double time, const std::vector<double>& rY)
{
    std::vector<double> dqs;

    RunModelRules(time, rY);

    dqs.push_back(cytosolmembraneandnucleus);
    dqs.push_back(C_F);
    dqs.push_back(C_T);
    dqs.push_back(drag);

    // AMOUNTS
    double amt__X = X * cytosolmembraneandnucleus;
    double amt__D = D * cytosolmembraneandnucleus;
    double amt__C_o = C_o * cytosolmembraneandnucleus;
    double amt__C_u = C_u * cytosolmembraneandnucleus;
    double amt__C_c = C_c * cytosolmembraneandnucleus;
    double amt__A = A * cytosolmembraneandnucleus;
    double amt__C_A = C_A * cytosolmembraneandnucleus;
    double amt__T = T * cytosolmembraneandnucleus;
    double amt__C_oT = C_oT * cytosolmembraneandnucleus;
    double amt__C_cT = C_cT * cytosolmembraneandnucleus;
    double amt__Y = Y * cytosolmembraneandnucleus;
    double amt__C_F = C_F * cytosolmembraneandnucleus;
    double amt__C_T = C_T * cytosolmembraneandnucleus;
    double amt__drag = drag * cytosolmembraneandnucleus;

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
    dqs.push_back(amt__C_F);
    dqs.push_back(amt__C_T);
    dqs.push_back(amt__drag);
    return dqs;
}

void VanLeeuwen2007NonDimSbmlOdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    RunModelRules(time, rY);

    rDY[0] = (-mwcfbf7716_cc13_473c_979a_033c57a28857 + mwab347951_a286_432d_b03b_254dcdba4a25 + mwc497befc_0edd_4b0d_8895_77dbfa05a4be - mwfa0f5940_d4e9_44e1_8a0c_379d6fdf9d0f) / cytosolmembraneandnucleus;                                                                                                                                                                     // d[X]/dt
    rDY[1] = (mwcfbf7716_cc13_473c_979a_033c57a28857 - mwab347951_a286_432d_b03b_254dcdba4a25 + mwbe108cd1_a5d9_4d55_97c9_ac584df6a005 - mwbe108cd1_a5d9_4d55_97c9_ac584df6a005 + mwad0ba91a_37fb_48f5_8392_6c043dfdd5dc - mwad0ba91a_37fb_48f5_8392_6c043dfdd5dc - mw321b3e5e_f6ed_4345_9346_55ffb1ff2137) / cytosolmembraneandnucleus;                                           // d[D]/dt
    rDY[2] = (-mwbe108cd1_a5d9_4d55_97c9_ac584df6a005 + mw552840aa_6d3c_4188_850d_9a3caafcdaa5 - mw5a301d6f_ffa1_4362_98de_0919a68808bc - mw1c7459d0_b5a9_4a89_8682_79fc8ca4ca6c + mw51e5ed29_82db_47fa_9eed_ab52b7786dd5 - mw8704b9b0_f93d_405a_89a1_e000c5b66896 + mwf5440ce3_586c_491f_ba71_bbf1c485027f - mw3ece1442_a0b3_47eb_a8bf_cb317a46165f) / cytosolmembraneandnucleus; // d[C_o]/dt
    rDY[3] = (mwbe108cd1_a5d9_4d55_97c9_ac584df6a005 + mwad0ba91a_37fb_48f5_8392_6c043dfdd5dc - mw1d3d836b_77bd_489b_9a37_acc874344652) / cytosolmembraneandnucleus;                                                                                                                                                                                                               // d[C_u]/dt
    rDY[4] = (-mwad0ba91a_37fb_48f5_8392_6c043dfdd5dc - mw4b073704_f41a_43d1_8ff8_af03f1cb6e1d - mw322e8f78_68f6_4293_ae2a_ab238ec5cf9b + mw08c332d6_1aaf_498f_8e1c_13f1b9da85bf + mw3ece1442_a0b3_47eb_a8bf_cb317a46165f) / cytosolmembraneandnucleus;                                                                                                                            // d[C_c]/dt
    rDY[5] = (-mw1c7459d0_b5a9_4a89_8682_79fc8ca4ca6c + mw51e5ed29_82db_47fa_9eed_ab52b7786dd5 + mwfb6f8dfa_2e07_4249_a8fe_af33ca27471e - mwaa4e7692_3290_45ab_8a38_b5ffd49ede87) / cytosolmembraneandnucleus;                                                                                                                                                                     // d[A]/dt
    rDY[6] = (mw1c7459d0_b5a9_4a89_8682_79fc8ca4ca6c - mw51e5ed29_82db_47fa_9eed_ab52b7786dd5) / cytosolmembraneandnucleus;                                                                                                                                                                                                                                                        // d[C_A]/dt
    rDY[7] = (-mw8704b9b0_f93d_405a_89a1_e000c5b66896 - mw322e8f78_68f6_4293_ae2a_ab238ec5cf9b + mwf5440ce3_586c_491f_ba71_bbf1c485027f + mw08c332d6_1aaf_498f_8e1c_13f1b9da85bf + mw86ea3c28_4745_4a89_9cbd_0f23c603e01b - mwee62535f_931d_41e9_ad53_2cad54b94778) / cytosolmembraneandnucleus;                                                                                   // d[T]/dt
    rDY[8] = (mw8704b9b0_f93d_405a_89a1_e000c5b66896 - mwf5440ce3_586c_491f_ba71_bbf1c485027f + mw1b250e04_290a_4689_9fc5_5f0bf6711e02 - mw1b250e04_290a_4689_9fc5_5f0bf6711e02) / cytosolmembraneandnucleus;                                                                                                                                                                      // d[C_oT]/dt
    rDY[9] = (mw322e8f78_68f6_4293_ae2a_ab238ec5cf9b - mw08c332d6_1aaf_498f_8e1c_13f1b9da85bf) / cytosolmembraneandnucleus;                                                                                                                                                                                                                                                        // d[C_cT]/dt
    rDY[10] = (mw1b250e04_290a_4689_9fc5_5f0bf6711e02 - mw925599eb_19a0_4434_8be3_67c40721b71d) / cytosolmembraneandnucleus;                                                                                                                                                                                                                                                       // d[Y]/dt

    // TODO: Scale time appropriately
}

double VanLeeuwen2007NonDimSbmlOdeSystem::ProcessModelEvents(double time, const std::vector<double>& rY)
{
    std::fill(std::begin(mEventAdjustedParameters), std::end(mEventAdjustedParameters), false);
    std::fill(std::begin(mEventAdjustedStateVars), std::end(mEventAdjustedStateVars), false);

    double min_dist = std::numeric_limits<double>::max();

    return min_dist; // Distance to closest event
}

// ASSIGNMENT RULES
void VanLeeuwen2007NonDimSbmlOdeSystem::RunAssignmentRules(double time)
{
    C_F = C_o + C_c;
    C_T = C_oT + C_cT;
    drag = sm::max((C_A - 2300.0) / 36.0, 1.0);
}

// REACTIONS
void VanLeeuwen2007NonDimSbmlOdeSystem::RunReactions(double time)
{
    // r1
    mwcfbf7716_cc13_473c_979a_033c57a28857 = s_D * gamma1 * X;

    // r2
    mwab347951_a286_432d_b03b_254dcdba4a25 = (d_Dx + wnt_level * xi_Dx) * D;

    // r22
    mwc497befc_0edd_4b0d_8895_77dbfa05a4be = s_X;

    // r23
    mwfa0f5940_d4e9_44e1_8a0c_379d6fdf9d0f = (d_X + wnt_level * xi_X) * X;

    // r7
    mwbe108cd1_a5d9_4d55_97c9_ac584df6a005 = p_u * gamma2 * C_o * D / (C_o + C_c + K_D);

    // r16
    mwad0ba91a_37fb_48f5_8392_6c043dfdd5dc = p_u * gamma2 * C_c * D / (C_c + C_o + K_D);

    // r8
    mw1d3d836b_77bd_489b_9a37_acc874344652 = d_u * C_u;

    // r3
    mw552840aa_6d3c_4188_850d_9a3caafcdaa5 = s_c;

    // r4
    mw5a301d6f_ffa1_4362_98de_0919a68808bc = d_c * C_o;

    // r17
    mw4b073704_f41a_43d1_8ff8_af03f1cb6e1d = d_c * C_c;

    // r9
    mw1c7459d0_b5a9_4a89_8682_79fc8ca4ca6c = s_CA * C_o * A;

    // r10
    mw51e5ed29_82db_47fa_9eed_ab52b7786dd5 = d_CA * C_A;

    // r11
    mw8704b9b0_f93d_405a_89a1_e000c5b66896 = s_CT * C_o * T;

    // r18
    mw322e8f78_68f6_4293_ae2a_ab238ec5cf9b = s_CT * C_c * T;

    // r12
    mwf5440ce3_586c_491f_ba71_bbf1c485027f = d_CT * C_oT;

    // r19
    mw08c332d6_1aaf_498f_8e1c_13f1b9da85bf = d_CT * C_cT;

    // r15
    mw3ece1442_a0b3_47eb_a8bf_cb317a46165f = (p_c + wnt_level * xi_C) * C_o / (C_o + K_C);

    // r5
    mwfb6f8dfa_2e07_4249_a8fe_af33ca27471e = s_A;

    // r6
    mwaa4e7692_3290_45ab_8a38_b5ffd49ede87 = d_A * A;

    // r20
    mw86ea3c28_4745_4a89_9cbd_0f23c603e01b = s_T;

    // r21
    mwee62535f_931d_41e9_ad53_2cad54b94778 = d_T * T;

    // r13
    mw1b250e04_290a_4689_9fc5_5f0bf6711e02 = s_Y * (C_oT + C_cT) / (C_oT + C_cT + K_T);

    // r14
    mw925599eb_19a0_4434_8be3_67c40721b71d = d_Y * Y;

    // r24
    mw321b3e5e_f6ed_4345_9346_55ffb1ff2137 = (d_D + wnt_level * xi_D) * D;
}

// VARIABLE PARAMETERS
void VanLeeuwen2007NonDimSbmlOdeSystem::UpdateParameters(double time)
{
    wnt_level = GetParameter(0);
    gamma1 = GetParameter(1);
    gamma2 = GetParameter(2);
    ComplexTransitThreshold = GetParameter(3);
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

    this->mDerivedQuantityNames.push_back("C_F");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("C_T");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("drag");
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

    this->mDerivedQuantityNames.push_back("amt__C_F");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("amt__C_T");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("amt__drag");
    this->mDerivedQuantityUnits.push_back("non-dim");

    // PARAMETERS
    this->mParameterNames.push_back("wnt_level");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("gamma1");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("gamma2");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ComplexTransitThreshold");
    this->mParameterUnits.push_back("non-dim");

    this->mInitialised = true;
}

// Register the ODE system with Boost serialization
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(VanLeeuwen2007NonDimSbmlOdeSystem)