#include <cmath>
#include <limits>

#include "CellwiseOdeSystemInformation.hpp"
#include "SbmlMath.hpp"

#include "VanLeeuwen2007NonDimSbmlOdeSystemAndSrnModel.hpp"

namespace sm = sbmlmath;

VanLeeuwen2007NonDimSbmlOdeSystem::VanLeeuwen2007NonDimSbmlOdeSystem(std::vector<double> stateVariables)
    : AbstractOdeSystem(14)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<VanLeeuwen2007NonDimSbmlOdeSystem>);

    // COMPARTMENTS:
    cytosolmembraneandnucleus = 1.0;

    // STATE VARIABLES:
    X = 0.268; // X
    D = 2.68; // D
    C_o = 76.2; // C_o
    C_u = 13.5; // C_u
    C_c = 0.0; // C_c
    A = 300.0; // A
    C_A = 544.2; // C_A
    T = 750.0; // T
    C_oT = 76.2; // C_oT
    C_cT = 0.0; // C_cT
    Y = 36.0; // Y
    C_F = 76.2; // C_F
    C_T = 76.2; // C_T
    drag = 1.0; // drag

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
    SetDefaultInitialCondition(11, C_F);
    SetDefaultInitialCondition(12, C_T);
    SetDefaultInitialCondition(13, drag);

    if (stateVariables.size() == 14)
    {
        X = stateVariables[0];
        D = stateVariables[1];
        C_o = stateVariables[2];
        C_u = stateVariables[3];
        C_c = stateVariables[4];
        A = stateVariables[5];
        C_A = stateVariables[6];
        T = stateVariables[7];
        C_oT = stateVariables[8];
        C_cT = stateVariables[9];
        Y = stateVariables[10];
        C_F = stateVariables[11];
        C_T = stateVariables[12];
        drag = stateVariables[13];
    }
    else if (stateVariables.size() != 0)
    {
        EXCEPTION("VanLeeuwen2007NonDimSbmlOdeSystem: Expected 14 state variables, got " + std::to_string(stateVariables.size()));
    }

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
    mStateVariables.push_back(C_F);
    mStateVariables.push_back(C_T);
    mStateVariables.push_back(drag);

    // STATE PARAMETERS:

    wnt_level = 0.0;
    gamma1 = 1.0;
    gamma2 = 1.0;
    ComplexTransitThreshold = 1.0;


    mParameters.push_back(wnt_level);
    mParameters.push_back(gamma1);
    mParameters.push_back(gamma2);
    mParameters.push_back(ComplexTransitThreshold);

}

VanLeeuwen2007NonDimSbmlOdeSystem::~VanLeeuwen2007NonDimSbmlOdeSystem()
{
}

void VanLeeuwen2007NonDimSbmlOdeSystem::RefreshState(const std::vector<double> &rY)
{
    // STATE VARIABLES:
    X = rY[0];
    D = rY[1];
    C_o = rY[2];
    C_u = rY[3];
    C_c = rY[4];
    A = rY[5];
    C_A = rY[6];
    T = rY[7];
    C_oT = rY[8];
    C_cT = rY[9];
    Y = rY[10];
    C_F = rY[11];
    C_T = rY[12];
    drag = rY[13];

    // STATE PARAMETERS:

    wnt_level = GetParameter("wnt_level");
    gamma1 = GetParameter("gamma1");
    gamma2 = GetParameter("gamma2");
    ComplexTransitThreshold = GetParameter("ComplexTransitThreshold");
}

void VanLeeuwen2007NonDimSbmlOdeSystem::EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY)
{
    RefreshState(rY);

    // RULES:

    // UPDATE STATE PARAMETERS:

    SetParameter("wnt_level", wnt_level);
    SetParameter("gamma1", gamma1);
    SetParameter("gamma2", gamma2);
    SetParameter("ComplexTransitThreshold", ComplexTransitThreshold);

    // REACTIONS:

    // r1
    double r1 = s_D * gamma1 * X;

    // r2
    double r2 = (d_Dx + wnt_level * xi_Dx) * D;

    // r22
    double r22 = s_X;

    // r23
    double r23 = (d_X + wnt_level * xi_X) * X;

    // r7
    double r7 = p_u * gamma2 * C_o * D / (C_o + C_c + K_D);

    // r16
    double r16 = p_u * gamma2 * C_c * D / (C_c + C_o + K_D);

    // r8
    double r8 = d_u * C_u;

    // r3
    double r3 = s_c;

    // r4
    double r4 = d_c * C_o;

    // r17
    double r17 = d_c * C_c;

    // r9
    double r9 = s_CA * C_o * A;

    // r10
    double r10 = d_CA * C_A;

    // r11
    double r11 = s_CT * C_o * T;

    // r18
    double r18 = s_CT * C_c * T;

    // r12
    double r12 = d_CT * C_oT;

    // r19
    double r19 = d_CT * C_cT;

    // r15
    double r15 = (p_c + wnt_level * xi_C) * C_o / (C_o + K_C);

    // r5
    double r5 = s_A;

    // r6
    double r6 = d_A * A;

    // r20
    double r20 = s_T;

    // r21
    double r21 = d_T * T;

    // r13
    double r13 = s_Y * (C_oT + C_cT) / (C_oT + C_cT + K_T);

    // r14
    double r14 = d_Y * Y;

    // r24
    double r24 = (d_D + wnt_level * xi_D) * D;

    // ODES:
    rDY[0] = (-r1 + r2 + r22 - r23) / cytosolmembraneandnucleus; // d[X]/dt
    rDY[1] = (r1 - r2 + r7 - r7 + r16 - r16 - r24) / cytosolmembraneandnucleus; // d[D]/dt
    rDY[2] = (-r7 + r3 - r4 - r9 + r10 - r11 + r12 - r15) / cytosolmembraneandnucleus; // d[C_o]/dt
    rDY[3] = (r7 + r16 - r8) / cytosolmembraneandnucleus; // d[C_u]/dt
    rDY[4] = (-r16 - r17 - r18 + r19 + r15) / cytosolmembraneandnucleus; // d[C_c]/dt
    rDY[5] = (-r9 + r10 + r5 - r6) / cytosolmembraneandnucleus; // d[A]/dt
    rDY[6] = (r9 - r10) / cytosolmembraneandnucleus; // d[C_A]/dt
    rDY[7] = (-r11 - r18 + r12 + r19 + r20 - r21) / cytosolmembraneandnucleus; // d[T]/dt
    rDY[8] = (r11 - r12 + r13 - r13) / cytosolmembraneandnucleus; // d[C_oT]/dt
    rDY[9] = (r18 - r19) / cytosolmembraneandnucleus; // d[C_cT]/dt
    rDY[10] = (r13 - r14) / cytosolmembraneandnucleus; // d[Y]/dt
    rDY[11] = ((C_o + C_c) - rY[11]) / cytosolmembraneandnucleus; // d[C_F]/dt
    rDY[12] = ((C_oT + C_cT) - rY[12]) / cytosolmembraneandnucleus; // d[C_T]/dt
    rDY[13] = ((sm::max((C_A - 2300.0) / 36.0, 1.0)) - rY[13]) / cytosolmembraneandnucleus; // d[drag]/dt

    // Scale time appropriately
}


// FUNCTION DEFINITIONS:

template <>
void CellwiseOdeSystemInformation<VanLeeuwen2007NonDimSbmlOdeSystem>::Initialise()
{
    // STATE VARIABLES:
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

    this->mVariableNames.push_back("C_F");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(76.2);

    this->mVariableNames.push_back("C_T");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(76.2);

    this->mVariableNames.push_back("drag");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(1.0);


    // STATE PARAMETERS:
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

// Define SbmlSrnWrapperModel using wrappers
#include "SbmlSrnWrapperModel.hpp"
#include "SbmlSrnWrapperModel.cpp"

typedef SbmlSrnWrapperModel<VanLeeuwen2007NonDimSbmlOdeSystem, 14> VanLeeuwen2007NonDimSbmlSrnModel;

// Declare identifiers for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(VanLeeuwen2007NonDimSbmlOdeSystem)
EXPORT_TEMPLATE_CLASS2(SbmlSrnWrapperModel, VanLeeuwen2007NonDimSbmlOdeSystem, 14)

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(VanLeeuwen2007NonDimSbmlSrnModel)