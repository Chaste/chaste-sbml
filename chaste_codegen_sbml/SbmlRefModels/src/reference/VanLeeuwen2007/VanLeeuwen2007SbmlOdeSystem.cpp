#include <cmath>
#include <limits>
#include <vector>

#include "CellwiseOdeSystemInformation.hpp"
#include "SbmlEventType.hpp"
#include "SbmlMath.hpp"

#include "VanLeeuwen2007SbmlOdeSystem.hpp"

namespace sm = sbmlmath;

VanLeeuwen2007SbmlOdeSystem::VanLeeuwen2007SbmlOdeSystem()
        : AbstractSbmlOdeSystem(11, 4, 0)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<VanLeeuwen2007SbmlOdeSystem>);

    // VARIABLE PARAMETERS
    wnt_level = 0.0;
    gamma1 = 1.0;
    gamma2 = 1.0;
    ComplexTransitThreshold = 1.0;

    // STATE VARIABLES
    X = 0.067;
    D = 0.67;
    C_o = 2.54;
    C_u = 0.45;
    C_c = 0.0;
    A = 10.0;
    C_A = 18.14;
    T = 25.0;
    C_oT = 2.54;
    C_cT = 0.0;
    Y = 0.48;

    // DERIVED QUANTITIES
    cytosolmembraneandnucleus = 1.0;
    C_F = 2.54;
    C_T = 2.54;
    drag = 1.0;

    // INITIAL ASSIGNMENTS
    C_F = C_o + C_c;                          //
    C_T = C_oT + C_cT;                        //
    drag = sm::max((C_A - 100.0) / 3.0, 1.0); //

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

VanLeeuwen2007SbmlOdeSystem::~VanLeeuwen2007SbmlOdeSystem()
{
}

std::vector<double> VanLeeuwen2007SbmlOdeSystem::ComputeDerivedQuantities(double time, const std::vector<double>& rY)
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

void VanLeeuwen2007SbmlOdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    RunModelRules(time, rY);

    rDY[0] = (-mwd6b35759_f098_484c_9c65_e84e7e4b61e4 + mweddac6d0_231e_4c92_ba2a_c91edc682ff5 + mwdf62dfed_ec88_4d81_bc9d_da0e10f41e4b - mwee9cc998_28e9_4173_a694_f3e278a639b7) / cytosolmembraneandnucleus;                                                                                                                                                                     // d[X]/dt
    rDY[1] = (mwd6b35759_f098_484c_9c65_e84e7e4b61e4 - mweddac6d0_231e_4c92_ba2a_c91edc682ff5 + mw661e341d_97d1_4e6f_8812_3be7ffc86d42 - mw661e341d_97d1_4e6f_8812_3be7ffc86d42 + mw179aa33c_9a7e_43c0_9285_3d8f97719c60 - mw179aa33c_9a7e_43c0_9285_3d8f97719c60 - mw931baf8f_6572_46f6_96eb_cae40ee267b7) / cytosolmembraneandnucleus;                                           // d[D]/dt
    rDY[2] = (-mw661e341d_97d1_4e6f_8812_3be7ffc86d42 + mwd4fa317a_a484_4f38_b9b0_6ed404d9adcb - mwcb88a249_a200_4e95_9185_5654bf1ebfc0 - mw4b47c66d_37e6_4c33_b043_1f6b3b814449 + mw69974db4_8ead_416c_a220_f6dc3be1f3b6 - mwccb628c3_76c8_47be_9d6b_6d3f1ae30dcc + mwb17c2c57_279d_4e88_b9cf_896029135cc1 - mw0be4a28b_e9c6_43da_8f95_d9c564a7caae) / cytosolmembraneandnucleus; // d[C_o]/dt
    rDY[3] = (mw661e341d_97d1_4e6f_8812_3be7ffc86d42 + mw179aa33c_9a7e_43c0_9285_3d8f97719c60 - mwff8d34f9_e036_49f1_b3b8_3706ecb98660) / cytosolmembraneandnucleus;                                                                                                                                                                                                               // d[C_u]/dt
    rDY[4] = (-mw179aa33c_9a7e_43c0_9285_3d8f97719c60 - mw0aeac2fc_dd2f_4fe7_b60c_97f0b008ed89 - mw581d69f1_60b3_4d21_9323_31b05ee89570 + mwe3236fc5_2118_40cb_8db3_ef9da29137cf + mw0be4a28b_e9c6_43da_8f95_d9c564a7caae) / cytosolmembraneandnucleus;                                                                                                                            // d[C_c]/dt
    rDY[5] = (-mw4b47c66d_37e6_4c33_b043_1f6b3b814449 + mw69974db4_8ead_416c_a220_f6dc3be1f3b6 + mwc360befb_07da_4d19_bbec_523fbef47dc9 - mw985b39ad_50ff_4d55_95f7_42e7f1bf6a3d) / cytosolmembraneandnucleus;                                                                                                                                                                     // d[A]/dt
    rDY[6] = (mw4b47c66d_37e6_4c33_b043_1f6b3b814449 - mw69974db4_8ead_416c_a220_f6dc3be1f3b6) / cytosolmembraneandnucleus;                                                                                                                                                                                                                                                        // d[C_A]/dt
    rDY[7] = (-mwccb628c3_76c8_47be_9d6b_6d3f1ae30dcc - mw581d69f1_60b3_4d21_9323_31b05ee89570 + mwb17c2c57_279d_4e88_b9cf_896029135cc1 + mwe3236fc5_2118_40cb_8db3_ef9da29137cf + mwf27b65bb_38d6_44cf_aa66_9b38ae5885eb - mwd7a7bdbb_8f2e_42ed_bf97_48e33e45dafb) / cytosolmembraneandnucleus;                                                                                   // d[T]/dt
    rDY[8] = (mwccb628c3_76c8_47be_9d6b_6d3f1ae30dcc - mwb17c2c57_279d_4e88_b9cf_896029135cc1 + mw988a8caf_bd68_462b_86d7_51844c1dcfd3 - mw988a8caf_bd68_462b_86d7_51844c1dcfd3) / cytosolmembraneandnucleus;                                                                                                                                                                      // d[C_oT]/dt
    rDY[9] = (mw581d69f1_60b3_4d21_9323_31b05ee89570 - mwe3236fc5_2118_40cb_8db3_ef9da29137cf) / cytosolmembraneandnucleus;                                                                                                                                                                                                                                                        // d[C_cT]/dt
    rDY[10] = (mw988a8caf_bd68_462b_86d7_51844c1dcfd3 - mw9ab26a4c_bd70_45e0_bacc_f830ab28abca) / cytosolmembraneandnucleus;                                                                                                                                                                                                                                                       // d[Y]/dt

    // TODO: Scale time appropriately
}

double VanLeeuwen2007SbmlOdeSystem::ProcessModelEvents(double time, const std::vector<double>& rY)
{
    std::fill(std::begin(mEventAdjustedParameters), std::end(mEventAdjustedParameters), false);
    std::fill(std::begin(mEventAdjustedStateVars), std::end(mEventAdjustedStateVars), false);

    double min_dist = std::numeric_limits<double>::max();

    return min_dist; // Distance to closest event
}

// ASSIGNMENT RULES
void VanLeeuwen2007SbmlOdeSystem::RunAssignmentRules(double time)
{
    C_F = C_o + C_c;
    C_T = C_oT + C_cT;
    drag = sm::max((C_A - 100.0) / 3.0, 1.0);
}

// REACTIONS
void VanLeeuwen2007SbmlOdeSystem::RunReactions(double time)
{
    // r1
    mwd6b35759_f098_484c_9c65_e84e7e4b61e4 = s_D * gamma1 * X;

    // r2
    mweddac6d0_231e_4c92_ba2a_c91edc682ff5 = (d_Dx + wnt_level * xi_Dx) * D;

    // r22
    mwdf62dfed_ec88_4d81_bc9d_da0e10f41e4b = s_X * cytosolmembraneandnucleus;

    // r23
    mwee9cc998_28e9_4173_a694_f3e278a639b7 = (d_X + wnt_level * xi_X) * X;

    // r7
    mw661e341d_97d1_4e6f_8812_3be7ffc86d42 = p_u * gamma2 * C_o * D / (C_o + C_c + K_D);

    // r16
    mw179aa33c_9a7e_43c0_9285_3d8f97719c60 = p_u * gamma2 * C_c * D / (C_c + C_o + K_D);

    // r8
    mwff8d34f9_e036_49f1_b3b8_3706ecb98660 = d_u * C_u * cytosolmembraneandnucleus;

    // r3
    mwd4fa317a_a484_4f38_b9b0_6ed404d9adcb = s_c * cytosolmembraneandnucleus;

    // r4
    mwcb88a249_a200_4e95_9185_5654bf1ebfc0 = d_c * C_o * cytosolmembraneandnucleus;

    // r17
    mw0aeac2fc_dd2f_4fe7_b60c_97f0b008ed89 = d_c * C_c * cytosolmembraneandnucleus;

    // r9
    mw4b47c66d_37e6_4c33_b043_1f6b3b814449 = s_CA * C_o * A * cytosolmembraneandnucleus;

    // r10
    mw69974db4_8ead_416c_a220_f6dc3be1f3b6 = d_CA * C_A * cytosolmembraneandnucleus;

    // r11
    mwccb628c3_76c8_47be_9d6b_6d3f1ae30dcc = s_CT * C_o * T * cytosolmembraneandnucleus;

    // r18
    mw581d69f1_60b3_4d21_9323_31b05ee89570 = s_CT * C_c * T * cytosolmembraneandnucleus;

    // r12
    mwb17c2c57_279d_4e88_b9cf_896029135cc1 = d_CT * C_oT * cytosolmembraneandnucleus;

    // r19
    mwe3236fc5_2118_40cb_8db3_ef9da29137cf = d_CT * C_cT * cytosolmembraneandnucleus;

    // r15
    mw0be4a28b_e9c6_43da_8f95_d9c564a7caae = (p_c + wnt_level * xi_C) * C_o / (C_o + K_C);

    // r5
    mwc360befb_07da_4d19_bbec_523fbef47dc9 = s_A * cytosolmembraneandnucleus;

    // r6
    mw985b39ad_50ff_4d55_95f7_42e7f1bf6a3d = d_A * A * cytosolmembraneandnucleus;

    // r20
    mwf27b65bb_38d6_44cf_aa66_9b38ae5885eb = s_T * cytosolmembraneandnucleus;

    // r21
    mwd7a7bdbb_8f2e_42ed_bf97_48e33e45dafb = d_T * T * cytosolmembraneandnucleus;

    // r13
    mw988a8caf_bd68_462b_86d7_51844c1dcfd3 = s_Y * (C_oT + C_cT) / (C_oT + C_cT + K_T);

    // r14
    mw9ab26a4c_bd70_45e0_bacc_f830ab28abca = d_Y * Y * cytosolmembraneandnucleus;

    // r24
    mw931baf8f_6572_46f6_96eb_cae40ee267b7 = (d_D + wnt_level * xi_D) * D;
}

// VARIABLE PARAMETERS
void VanLeeuwen2007SbmlOdeSystem::UpdateParameters(double time)
{
    wnt_level = GetParameter(0);
    gamma1 = GetParameter(1);
    gamma2 = GetParameter(2);
    ComplexTransitThreshold = GetParameter(3);
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
CHASTE_CLASS_EXPORT(VanLeeuwen2007SbmlOdeSystem)