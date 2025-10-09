#include <cmath>
#include <limits>
#include <vector>

#include "CellwiseOdeSystemInformation.hpp"
#include "SbmlEventType.hpp"
#include "SbmlMath.hpp"

#include "Chen2004SbmlOdeSystem.hpp"

namespace sm = sbmlmath;

Chen2004SbmlOdeSystem::Chen2004SbmlOdeSystem()
        : AbstractSbmlOdeSystem(36, 14, 4)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<Chen2004SbmlOdeSystem>);

    // STATE VARIABLES
    BUD = 0.008473;
    C2 = 0.238404;
    C2P = 0.024034;
    C5 = 0.070081;
    C5P = 0.006878;
    CDC14 = 0.468344;
    CDC15 = 0.656533;
    CDC20 = 0.444296;
    CDC20i = 1.472044;
    CDC6 = 0.10758;
    CDC6P = 0.015486;
    CDH1 = 0.930499;
    CDH1i = 0.0695;
    CLB2 = 0.1469227;
    CLB5 = 0.0518014;
    CLN2 = 0.0652511;
    ESP1 = 0.301313;
    F2 = 0.236058;
    F2P = 0.0273938;
    F5 = 7.24e-05;
    F5P = 7.91e-05;
    IEP = 0.1015;
    MASS = 1.206019;
    NET1 = 0.018645;
    NET1P = 0.970271;
    ORI = 0.000909;
    PDS1 = 0.025612;
    PPX = 0.123179;
    RENT = 1.04954;
    RENTP = 0.6;
    SIC1 = 0.0228776;
    SIC1P = 0.00641;
    SPN = 0.03;
    SWI5 = 0.95;
    SWI5P = 0.02;
    TEM1GTP = 0.9;

    SetDefaultInitialCondition(0, BUD);
    SetDefaultInitialCondition(1, C2);
    SetDefaultInitialCondition(2, C2P);
    SetDefaultInitialCondition(3, C5);
    SetDefaultInitialCondition(4, C5P);
    SetDefaultInitialCondition(5, CDC14);
    SetDefaultInitialCondition(6, CDC15);
    SetDefaultInitialCondition(7, CDC20);
    SetDefaultInitialCondition(8, CDC20i);
    SetDefaultInitialCondition(9, CDC6);
    SetDefaultInitialCondition(10, CDC6P);
    SetDefaultInitialCondition(11, CDH1);
    SetDefaultInitialCondition(12, CDH1i);
    SetDefaultInitialCondition(13, CLB2);
    SetDefaultInitialCondition(14, CLB5);
    SetDefaultInitialCondition(15, CLN2);
    SetDefaultInitialCondition(16, ESP1);
    SetDefaultInitialCondition(17, F2);
    SetDefaultInitialCondition(18, F2P);
    SetDefaultInitialCondition(19, F5);
    SetDefaultInitialCondition(20, F5P);
    SetDefaultInitialCondition(21, IEP);
    SetDefaultInitialCondition(22, MASS);
    SetDefaultInitialCondition(23, NET1);
    SetDefaultInitialCondition(24, NET1P);
    SetDefaultInitialCondition(25, ORI);
    SetDefaultInitialCondition(26, PDS1);
    SetDefaultInitialCondition(27, PPX);
    SetDefaultInitialCondition(28, RENT);
    SetDefaultInitialCondition(29, RENTP);
    SetDefaultInitialCondition(30, SIC1);
    SetDefaultInitialCondition(31, SIC1P);
    SetDefaultInitialCondition(32, SPN);
    SetDefaultInitialCondition(33, SWI5);
    SetDefaultInitialCondition(34, SWI5P);
    SetDefaultInitialCondition(35, TEM1GTP);

    mStateVariables.push_back(BUD);
    mStateVariables.push_back(C2);
    mStateVariables.push_back(C2P);
    mStateVariables.push_back(C5);
    mStateVariables.push_back(C5P);
    mStateVariables.push_back(CDC14);
    mStateVariables.push_back(CDC15);
    mStateVariables.push_back(CDC20);
    mStateVariables.push_back(CDC20i);
    mStateVariables.push_back(CDC6);
    mStateVariables.push_back(CDC6P);
    mStateVariables.push_back(CDH1);
    mStateVariables.push_back(CDH1i);
    mStateVariables.push_back(CLB2);
    mStateVariables.push_back(CLB5);
    mStateVariables.push_back(CLN2);
    mStateVariables.push_back(ESP1);
    mStateVariables.push_back(F2);
    mStateVariables.push_back(F2P);
    mStateVariables.push_back(F5);
    mStateVariables.push_back(F5P);
    mStateVariables.push_back(IEP);
    mStateVariables.push_back(MASS);
    mStateVariables.push_back(NET1);
    mStateVariables.push_back(NET1P);
    mStateVariables.push_back(ORI);
    mStateVariables.push_back(PDS1);
    mStateVariables.push_back(PPX);
    mStateVariables.push_back(RENT);
    mStateVariables.push_back(RENTP);
    mStateVariables.push_back(SIC1);
    mStateVariables.push_back(SIC1P);
    mStateVariables.push_back(SPN);
    mStateVariables.push_back(SWI5);
    mStateVariables.push_back(SWI5P);
    mStateVariables.push_back(TEM1GTP);

    // DERIVED QUANTITIES
    BCK2 = 0.0;
    CDC14T = 2.0;
    CDC15i = 0.0;
    CDC6T = 0.0;
    CKIT = 0.0;
    CLB2T = 0.17;
    CLB5T = 0.12;
    CLN3 = 0.0;
    IE = 0.0;
    MCM1 = 0.0;
    NET1T = 2.8;
    PE = 0.0;
    SBF = 0.0;
    SIC1T = 0.0;
    TEM1GDP = 0.0;

    // VARIABLE PARAMETERS
    cell = 1.0;
    BUB2 = 0.2;
    LTE1 = 0.1;
    MAD2 = 0.01;
    bub2l = 0.2;
    CDC15T = 1.0;
    ESP1T = 1.0;
    IET = 1.0;
    KEZ = 0.3;
    KEZ2 = 0.2;
    lte1h = 1.0;
    lte1l = 0.1;
    mad2l = 0.01;
    TEM1T = 1.0;

    mParameters.push_back(cell);
    mParameters.push_back(BUB2);
    mParameters.push_back(LTE1);
    mParameters.push_back(MAD2);
    mParameters.push_back(bub2l);
    mParameters.push_back(CDC15T);
    mParameters.push_back(ESP1T);
    mParameters.push_back(IET);
    mParameters.push_back(KEZ);
    mParameters.push_back(KEZ2);
    mParameters.push_back(lte1h);
    mParameters.push_back(lte1l);
    mParameters.push_back(mad2l);
    mParameters.push_back(TEM1T);

    // RULE-BASED PARAMETERS
    D = 0.0;
    mu = 0.0;
    Vdb5 = 0.0;
    Vdb2 = 0.0;
    Vasbf = 0.0;
    Visbf = 0.0;
    Vkpc1 = 0.0;
    Vkpf6 = 0.0;
    Vacdh = 0.0;
    Vicdh = 0.0;
    Vppnet = 0.0;
    Vkpnet = 0.0;
    Vdppx = 0.0;
    Vdpds = 0.0;
    Vaiep = 0.0;
    Vd2c1 = 0.0;
    Vd2f6 = 0.0;
    Vppc1 = 0.0;
    Vppf6 = 0.0;
    F = 0.0;

    // REACTIONS
    Growth = 0.0;
    Synthesis_of_CLN2 = 0.0;
    Degradation_of_CLN2 = 0.0;
    Synthesis_of_CLB2 = 0.0;
    Degradation_of_CLB2 = 0.0;
    Synthesis_of_CLB5 = 0.0;
    Degradation_of_CLB5 = 0.0;
    Synthesis_of_SIC1 = 0.0;
    Phosphorylation_of_SIC1 = 0.0;
    Dephosphorylation_of_SIC1 = 0.0;
    Fast_Degradation_of_SIC1P = 0.0;
    Assoc_of_CLB2_and_SIC1 = 0.0;
    Dissoc_of_CLB2SIC1_complex = 0.0;
    Assoc_of_CLB5_and_SIC1 = 0.0;
    Dissoc_of_CLB5SIC1 = 0.0;
    Phosphorylation_of_C2 = 0.0;
    Dephosphorylation_of_C2P = 0.0;
    Phosphorylation_of_C5 = 0.0;
    Dephosphorylation_of_C5P = 0.0;
    Degradation_of_CLB2_in_C2 = 0.0;
    Degradation_of_CLB5_in_C5 = 0.0;
    Degradation_of_SIC1_in_C2P = 0.0;
    Degradation_of_SIC1P_in_C5P_ = 0.0;
    Degradation_of_CLB2_in_C2P = 0.0;
    Degradation_of_CLB5_in_C5P = 0.0;
    CDC6_synthesis = 0.0;
    Phosphorylation_of_CDC6 = 0.0;
    Dephosphorylation_of_CDC6 = 0.0;
    Degradation_of_CDC6P = 0.0;
    CLB2CDC6_complex_formation = 0.0;
    CLB2CDC6_dissociation = 0.0;
    CLB5CDC6_complex_formation = 0.0;
    CLB5CDC6_dissociation = 0.0;
    F2_phosphorylation = 0.0;
    F2P_dephosphorylation = 0.0;
    F5_phosphorylation = 0.0;
    F5P_dephosphorylation = 0.0;
    CLB2_degradation_in_F2 = 0.0;
    CLB5_degradation_in_F5 = 0.0;
    CDC6_degradation_in_F2P = 0.0;
    CDC6_degradation_in_F5P = 0.0;
    CLB2_degradation_in_F2P = 0.0;
    CLB5_degradation_in_F5P = 0.0;
    Synthesis_of_SWI5 = 0.0;
    Degradation_of_SWI5 = 0.0;
    Degradation_of_SWI5P = 0.0;
    Activation_of_SWI5 = 0.0;
    Inactivation_of_SWI5 = 0.0;
    Activation_of_IEP = 0.0;
    Inactivation_1 = 0.0;
    Synthesis_of_inactive_CDC20 = 0.0;
    Degradation_of_inactiveCDC20 = 0.0;
    Degradation_of_active_CDC20 = 0.0;
    Activation_of_CDC20 = 0.0;
    Inactivation_2 = 0.0;
    CDH1_synthesis = 0.0;
    CDH1_degradation = 0.0;
    CDH1i_degradation = 0.0;
    CDH1i_activation = 0.0;
    Inactivation_3 = 0.0;
    CDC14_synthesis = 0.0;
    CDC14_degradation = 0.0;
    Assoc_with_NET1_to_form_RENT = 0.0;
    Dissoc_from_RENT = 0.0;
    Assoc_with_NET1P_to_form_RENTP = 0.0;
    Dissoc_from_RENP = 0.0;
    Net1_synthesis = 0.0;
    Net1_degradation = 0.0;
    Net1P_degradation = 0.0;
    NET1_phosphorylation = 0.0;
    dephosphorylation_1 = 0.0;
    RENT_phosphorylation = 0.0;
    dephosphorylation_2 = 0.0;
    Degradation_of_NET1_in_RENT = 0.0;
    Degradation_of_NET1P_in_RENTP = 0.0;
    Degradation_of_CDC14_in_RENT = 0.0;
    Degradation_of_CDC14_in_RENTP = 0.0;
    TEM1_activation = 0.0;
    inactivation_1 = 0.0;
    CDC15_activation = 0.0;
    inactivation_2 = 0.0;
    PPX_synthesis = 0.0;
    degradation_1 = 0.0;
    PDS1_synthesis = 0.0;
    degradation_2 = 0.0;
    Degradation_of_PDS1_in_PE = 0.0;
    Assoc_with_ESP1_to_form_PE = 0.0;
    Disso_from_PE = 0.0;
    DNA_synthesis = 0.0;
    Negative_regulation_of_DNA_synthesis = 0.0;
    Budding = 0.0;
    Negative_regulation_of_Cell_budding = 0.0;
    Spindle_formation = 0.0;
    Spindle_disassembly = 0.0;

    // EVENTS
    mEventType.resize(4, SbmlEventType::UNKNOWN);

    // Uncomment lines below for events that should trigger cell division

    // mEventType[0] = SbmlEventType::CELL_DIVISION; // reset ORI
    // mEventType[1] = SbmlEventType::CELL_DIVISION; // start DNA synthesis
    // mEventType[2] = SbmlEventType::CELL_DIVISION; // spindle checkpoint
    mEventType[3] = SbmlEventType::CELL_DIVISION; // cell division

    mEventSatisfied.resize(4, true); // Prevent events from triggering at the start
    mEventTriggered.resize(4, false);

    mEventAdjustedParameters.resize(14, false);
    mEventAdjustedParameterValues.resize(14, 0.0);

    mEventAdjustedStateVars.resize(36, false);
    mEventAdjustedStateValues.resize(36, 0.0);

    // Run model rules to complete state initialisation
    RunModelRules(0.0, mStateVariables);
}

Chen2004SbmlOdeSystem::~Chen2004SbmlOdeSystem()
{
}

std::vector<double> Chen2004SbmlOdeSystem::ComputeDerivedQuantities(double time, const std::vector<double>& rY)
{
    RunModelRules(time, rY);

    std::vector<double> dqs;
    dqs.push_back(BCK2);
    dqs.push_back(CDC14T);
    dqs.push_back(CDC15i);
    dqs.push_back(CDC6T);
    dqs.push_back(CKIT);
    dqs.push_back(CLB2T);
    dqs.push_back(CLB5T);
    dqs.push_back(CLN3);
    dqs.push_back(IE);
    dqs.push_back(MCM1);
    dqs.push_back(NET1T);
    dqs.push_back(PE);
    dqs.push_back(SBF);
    dqs.push_back(SIC1T);
    dqs.push_back(TEM1GDP);
    return dqs;
}

void Chen2004SbmlOdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    RunModelRules(time, rY);

    rDY[0] = (Budding - Negative_regulation_of_Cell_budding) / cell;                                                                                                                                                                                        // d[BUD]/dt
    rDY[1] = (Assoc_of_CLB2_and_SIC1 - Dissoc_of_CLB2SIC1_complex - Phosphorylation_of_C2 + Dephosphorylation_of_C2P - Degradation_of_CLB2_in_C2) / cell;                                                                                                   // d[C2]/dt
    rDY[2] = (Phosphorylation_of_C2 - Dephosphorylation_of_C2P - Degradation_of_SIC1_in_C2P - Degradation_of_CLB2_in_C2P) / cell;                                                                                                                           // d[C2P]/dt
    rDY[3] = (Assoc_of_CLB5_and_SIC1 - Dissoc_of_CLB5SIC1 - Phosphorylation_of_C5 + Dephosphorylation_of_C5P - Degradation_of_CLB5_in_C5) / cell;                                                                                                           // d[C5]/dt
    rDY[4] = (Phosphorylation_of_C5 - Dephosphorylation_of_C5P - Degradation_of_SIC1P_in_C5P_ - Degradation_of_CLB5_in_C5P) / cell;                                                                                                                         // d[C5P]/dt
    rDY[5] = (CDC14_synthesis - CDC14_degradation - Assoc_with_NET1_to_form_RENT + Dissoc_from_RENT - Assoc_with_NET1P_to_form_RENTP + Dissoc_from_RENP + Degradation_of_NET1_in_RENT + Degradation_of_NET1P_in_RENTP) / cell;                              // d[CDC14]/dt
    rDY[6] = (CDC15_activation - inactivation_2) / cell;                                                                                                                                                                                                    // d[CDC15]/dt
    rDY[7] = (-Degradation_of_active_CDC20 + Activation_of_CDC20 - Inactivation_2) / cell;                                                                                                                                                                  // d[CDC20]/dt
    rDY[8] = (Synthesis_of_inactive_CDC20 - Degradation_of_inactiveCDC20 - Activation_of_CDC20 + Inactivation_2) / cell;                                                                                                                                    // d[CDC20i]/dt
    rDY[9] = (CDC6_synthesis - Phosphorylation_of_CDC6 + Dephosphorylation_of_CDC6 - CLB2CDC6_complex_formation + CLB2CDC6_dissociation - CLB5CDC6_complex_formation + CLB5CDC6_dissociation + CLB2_degradation_in_F2 + CLB5_degradation_in_F5) / cell;     // d[CDC6]/dt
    rDY[10] = (Phosphorylation_of_CDC6 - Dephosphorylation_of_CDC6 - Degradation_of_CDC6P + CLB2_degradation_in_F2P + CLB5_degradation_in_F5P) / cell;                                                                                                      // d[CDC6P]/dt
    rDY[11] = (CDH1_synthesis - CDH1_degradation + CDH1i_activation - Inactivation_3) / cell;                                                                                                                                                               // d[CDH1]/dt
    rDY[12] = (-CDH1i_degradation - CDH1i_activation + Inactivation_3) / cell;                                                                                                                                                                              // d[CDH1i]/dt
    rDY[13] = (Synthesis_of_CLB2 - Degradation_of_CLB2 - Assoc_of_CLB2_and_SIC1 + Dissoc_of_CLB2SIC1_complex + Degradation_of_SIC1_in_C2P - CLB2CDC6_complex_formation + CLB2CDC6_dissociation + CDC6_degradation_in_F2P) / cell;                           // d[CLB2]/dt
    rDY[14] = (Synthesis_of_CLB5 - Degradation_of_CLB5 - Assoc_of_CLB5_and_SIC1 + Dissoc_of_CLB5SIC1 + Degradation_of_SIC1P_in_C5P_ - CLB5CDC6_complex_formation + CLB5CDC6_dissociation + CDC6_degradation_in_F5P) / cell;                                 // d[CLB5]/dt
    rDY[15] = (Synthesis_of_CLN2 - Degradation_of_CLN2) / cell;                                                                                                                                                                                             // d[CLN2]/dt
    rDY[16] = (Degradation_of_PDS1_in_PE - Assoc_with_ESP1_to_form_PE + Disso_from_PE) / cell;                                                                                                                                                              // d[ESP1]/dt
    rDY[17] = (CLB2CDC6_complex_formation - CLB2CDC6_dissociation - F2_phosphorylation + F2P_dephosphorylation - CLB2_degradation_in_F2) / cell;                                                                                                            // d[F2]/dt
    rDY[18] = (F2_phosphorylation - F2P_dephosphorylation - CDC6_degradation_in_F2P - CLB2_degradation_in_F2P) / cell;                                                                                                                                      // d[F2P]/dt
    rDY[19] = (CLB5CDC6_complex_formation - CLB5CDC6_dissociation - F5_phosphorylation + F5P_dephosphorylation - CLB5_degradation_in_F5) / cell;                                                                                                            // d[F5]/dt
    rDY[20] = (F5_phosphorylation - F5P_dephosphorylation - CDC6_degradation_in_F5P - CLB5_degradation_in_F5P) / cell;                                                                                                                                      // d[F5P]/dt
    rDY[21] = (Activation_of_IEP - Inactivation_1) / cell;                                                                                                                                                                                                  // d[IEP]/dt
    rDY[22] = Growth / cell;                                                                                                                                                                                                                              // d[MASS]/dt
    rDY[23] = (-Assoc_with_NET1_to_form_RENT + Dissoc_from_RENT + Net1_synthesis - Net1_degradation - NET1_phosphorylation + dephosphorylation_1 + Degradation_of_CDC14_in_RENT) / cell;                                                                    // d[NET1]/dt
    rDY[24] = (-Assoc_with_NET1P_to_form_RENTP + Dissoc_from_RENP - Net1P_degradation + NET1_phosphorylation - dephosphorylation_1 + Degradation_of_CDC14_in_RENTP) / cell;                                                                                 // d[NET1P]/dt
    rDY[25] = (DNA_synthesis - Negative_regulation_of_DNA_synthesis) / cell;                                                                                                                                                                                // d[ORI]/dt
    rDY[26] = (PDS1_synthesis - degradation_2 - Assoc_with_ESP1_to_form_PE + Disso_from_PE) / cell;                                                                                                                                                         // d[PDS1]/dt
    rDY[27] = (PPX_synthesis - degradation_1) / cell;                                                                                                                                                                                                       // d[PPX]/dt
    rDY[28] = (Assoc_with_NET1_to_form_RENT - Dissoc_from_RENT - RENT_phosphorylation + dephosphorylation_2 - Degradation_of_NET1_in_RENT - Degradation_of_CDC14_in_RENT) / cell;                                                                           // d[RENT]/dt
    rDY[29] = (Assoc_with_NET1P_to_form_RENTP - Dissoc_from_RENP + RENT_phosphorylation - dephosphorylation_2 - Degradation_of_NET1P_in_RENTP - Degradation_of_CDC14_in_RENTP) / cell;                                                                      // d[RENTP]/dt
    rDY[30] = (Synthesis_of_SIC1 - Phosphorylation_of_SIC1 + Dephosphorylation_of_SIC1 - Assoc_of_CLB2_and_SIC1 + Dissoc_of_CLB2SIC1_complex - Assoc_of_CLB5_and_SIC1 + Dissoc_of_CLB5SIC1 + Degradation_of_CLB2_in_C2 + Degradation_of_CLB5_in_C5) / cell; // d[SIC1]/dt
    rDY[31] = (Phosphorylation_of_SIC1 - Dephosphorylation_of_SIC1 - Fast_Degradation_of_SIC1P + Degradation_of_CLB2_in_C2P + Degradation_of_CLB5_in_C5P) / cell;                                                                                           // d[SIC1P]/dt
    rDY[32] = (Spindle_formation - Spindle_disassembly) / cell;                                                                                                                                                                                             // d[SPN]/dt
    rDY[33] = (Synthesis_of_SWI5 - Degradation_of_SWI5 + Activation_of_SWI5 - Inactivation_of_SWI5) / cell;                                                                                                                                                 // d[SWI5]/dt
    rDY[34] = (-Degradation_of_SWI5P - Activation_of_SWI5 + Inactivation_of_SWI5) / cell;                                                                                                                                                                   // d[SWI5P]/dt
    rDY[35] = (TEM1_activation - inactivation_1) / cell;                                                                                                                                                                                                    // d[TEM1GTP]/dt

    // TODO: Scale time appropriately
}

double Chen2004SbmlOdeSystem::ProcessModelEvents(double time, const std::vector<double>& rY)
{
    std::fill(std::begin(mEventAdjustedParameters), std::end(mEventAdjustedParameters), false);
    std::fill(std::begin(mEventAdjustedStateVars), std::end(mEventAdjustedStateVars), false);

    double min_dist = std::numeric_limits<double>::max();

    //========================================
    // EVENT: reset ORI
    //========================================
    {
        double event_dist = (0.0) - (CLB2 + CLB5 - KEZ2) - std::numeric_limits<double>::epsilon();

        // Avoid oscillation by ensuring event_dist is not close to 0 unless triggered
        if (std::abs(event_dist) < 1.0)
        {
            event_dist = 1.0;
        }

        // Update min_dist
        if (std::abs(event_dist) < std::abs(min_dist))
        {
            min_dist = event_dist;
        }

        // Process the event
        if (sm::lt(CLB2 + CLB5 - KEZ2, 0.0))
        {
            if (!mEventSatisfied[0])
            {
                // The condition is transitioning from false -> true: trigger the event
                mEventTriggered[0] = true;
                event_dist = 0.0;
                min_dist = 0.0;

                // Adjust relevant state variables and parameters
                // ORI = 0.0
                mEventAdjustedStateVars[25] = true;
                mEventAdjustedStateValues[25] = 0.0;
            }
            mEventSatisfied[0] = true;
        }
        else
        {
            mEventSatisfied[0] = false;
            mEventTriggered[0] = false;
        }
    }

    //========================================
    // EVENT: start DNA synthesis
    //========================================
    {
        double event_dist = (ORI - 1.0) - (0.0) - std::numeric_limits<double>::epsilon();

        // Avoid oscillation by ensuring event_dist is not close to 0 unless triggered
        if (std::abs(event_dist) < 1.0)
        {
            event_dist = 1.0;
        }

        // Update min_dist
        if (std::abs(event_dist) < std::abs(min_dist))
        {
            min_dist = event_dist;
        }

        // Process the event
        if (sm::gt(ORI - 1.0, 0.0))
        {
            if (!mEventSatisfied[1])
            {
                // The condition is transitioning from false -> true: trigger the event
                mEventTriggered[1] = true;
                event_dist = 0.0;
                min_dist = 0.0;

                // Adjust relevant state variables and parameters
                // MAD2 = mad2h
                mEventAdjustedParameters[3] = true;
                mEventAdjustedParameterValues[3] = mad2h;

                // BUB2 = bub2h
                mEventAdjustedParameters[1] = true;
                mEventAdjustedParameterValues[1] = bub2h;
            }
            mEventSatisfied[1] = true;
        }
        else
        {
            mEventSatisfied[1] = false;
            mEventTriggered[1] = false;
        }
    }

    //========================================
    // EVENT: spindle checkpoint
    //========================================
    {
        double event_dist = (SPN - 1.0) - (0.0) - std::numeric_limits<double>::epsilon();

        // Avoid oscillation by ensuring event_dist is not close to 0 unless triggered
        if (std::abs(event_dist) < 1.0)
        {
            event_dist = 1.0;
        }

        // Update min_dist
        if (std::abs(event_dist) < std::abs(min_dist))
        {
            min_dist = event_dist;
        }

        // Process the event
        if (sm::gt(SPN - 1.0, 0.0))
        {
            if (!mEventSatisfied[2])
            {
                // The condition is transitioning from false -> true: trigger the event
                mEventTriggered[2] = true;
                event_dist = 0.0;
                min_dist = 0.0;

                // Adjust relevant state variables and parameters
                // MAD2 = mad2l
                mEventAdjustedParameters[3] = true;
                mEventAdjustedParameterValues[3] = mad2l;

                // LTE1 = lte1h
                mEventAdjustedParameters[2] = true;
                mEventAdjustedParameterValues[2] = lte1h;

                // BUB2 = bub2l
                mEventAdjustedParameters[1] = true;
                mEventAdjustedParameterValues[1] = bub2l;
            }
            mEventSatisfied[2] = true;
        }
        else
        {
            mEventSatisfied[2] = false;
            mEventTriggered[2] = false;
        }
    }

    //========================================
    // EVENT: cell division
    //========================================
    {
        double event_dist = (0.0) - (CLB2 - KEZ) - std::numeric_limits<double>::epsilon();

        // Avoid oscillation by ensuring event_dist is not close to 0 unless triggered
        if (std::abs(event_dist) < 1.0)
        {
            event_dist = 1.0;
        }

        // Update min_dist
        if (std::abs(event_dist) < std::abs(min_dist))
        {
            min_dist = event_dist;
        }

        // Process the event
        if (sm::lt(CLB2 - KEZ, 0.0))
        {
            if (!mEventSatisfied[3])
            {
                // The condition is transitioning from false -> true: trigger the event
                mEventTriggered[3] = true;
                event_dist = 0.0;
                min_dist = 0.0;

                // Adjust relevant state variables and parameters
                // MASS = F * MASS
                mEventAdjustedStateVars[22] = true;
                mEventAdjustedStateValues[22] = F * MASS;

                // LTE1 = lte1l
                mEventAdjustedParameters[2] = true;
                mEventAdjustedParameterValues[2] = lte1l;

                // BUD = 0.0
                mEventAdjustedStateVars[0] = true;
                mEventAdjustedStateValues[0] = 0.0;

                // SPN = 0.0
                mEventAdjustedStateVars[32] = true;
                mEventAdjustedStateValues[32] = 0.0;
            }
            mEventSatisfied[3] = true;
        }
        else
        {
            mEventSatisfied[3] = false;
            mEventTriggered[3] = false;
        }
    }

    return min_dist; // Distance to closest event
}

void Chen2004SbmlOdeSystem::RunModelRules(double time, const std::vector<double>& rY)
{
    // STATE VARIABLES
    BUD = rY[0];
    C2 = rY[1];
    C2P = rY[2];
    C5 = rY[3];
    C5P = rY[4];
    CDC14 = rY[5];
    CDC15 = rY[6];
    CDC20 = rY[7];
    CDC20i = rY[8];
    CDC6 = rY[9];
    CDC6P = rY[10];
    CDH1 = rY[11];
    CDH1i = rY[12];
    CLB2 = rY[13];
    CLB5 = rY[14];
    CLN2 = rY[15];
    ESP1 = rY[16];
    F2 = rY[17];
    F2P = rY[18];
    F5 = rY[19];
    F5P = rY[20];
    IEP = rY[21];
    MASS = rY[22];
    NET1 = rY[23];
    NET1P = rY[24];
    ORI = rY[25];
    PDS1 = rY[26];
    PPX = rY[27];
    RENT = rY[28];
    RENTP = rY[29];
    SIC1 = rY[30];
    SIC1P = rY[31];
    SPN = rY[32];
    SWI5 = rY[33];
    SWI5P = rY[34];
    TEM1GTP = rY[35];

    // VARIABLE PARAMETERS
    cell = GetParameter(0);
    BUB2 = GetParameter(1);
    LTE1 = GetParameter(2);
    MAD2 = GetParameter(3);
    bub2l = GetParameter(4);
    CDC15T = GetParameter(5);
    ESP1T = GetParameter(6);
    IET = GetParameter(7);
    KEZ = GetParameter(8);
    KEZ2 = GetParameter(9);
    lte1h = GetParameter(10);
    lte1l = GetParameter(11);
    mad2l = GetParameter(12);
    TEM1T = GetParameter(13);

    // ASSIGNMENT RULES
    BCK2 = b0 * MASS;
    Visbf = kisbf_p + kisbf_p_p * CLB2;
    CLN3 = C0 * Dn3 * MASS / (Jn3 + Dn3 * MASS);
    Vppc1 = kppc1 * CDC14;
    Vppf6 = kppf6 * CDC14;
    Vaiep = kaiep * CLB2;
    Vacdh = kacdh_p + kacdh_p_p * CDC14;
    Vicdh = kicdh_p + kicdh_p_p * (eicdhn3 * CLN3 + eicdhn2 * CLN2 + eicdhb5 * CLB5 + eicdhb2 * CLB2);
    Vkpnet = (kkpnet_p + kkpnet_p_p * CDC15) * MASS;
    Vppnet = kppnet_p + kppnet_p_p * PPX;
    Vasbf = kasbf * (esbfn2 * CLN2 + esbfn3 * (CLN3 + BCK2) + esbfb5 * CLB5);
    SBF = GK_219(Vasbf, Visbf, Jasbf, Jisbf);
    MCM1 = GK_219(kamcm * CLB2, kimcm, Jamcm, Jimcm);
    mu = sm::log(2.0) / mdt;
    D = 1.026 / mu - 32.0;
    F = std::exp(-mu * D);
    Vd2c1 = kd2c1 * (ec1n3 * CLN3 + ec1k2 * BCK2 + ec1n2 * CLN2 + ec1b5 * CLB5 + ec1b2 * CLB2);
    Vd2f6 = kd2f6 * (ef6n3 * CLN3 + ef6k2 * BCK2 + ef6n2 * CLN2 + ef6b5 * CLB5 + ef6b2 * CLB2);
    Vkpc1 = kd1c1 + Vd2c1 / (Jd2c1 + SIC1 + C2 + C5 + SIC1P + C2P + C5P);
    Vkpf6 = kd1f6 + Vd2f6 / (Jd2f6 + CDC6 + F2 + F5 + CDC6P + F2P + F5P);
    Vdb2 = kdb2_p + kdb2_p_p * CDH1 + kdb2p * CDC20;
    Vdb5 = kdb5_p + kdb5_p_p * CDC20;
    Vdpds = kd1pds_p + kd2pds_p_p * CDC20 + kd3pds_p_p * CDH1;
    Vdppx = kdppx_p + kdppx_p_p * (J20ppx + CDC20) * Jpds / (Jpds + PDS1);
    CLB2T = CLB2 + C2 + C2P + F2 + F2P;
    CLB5T = CLB5 + C5 + C5P + F5 + F5P;
    CDC14T = CDC14 + RENT + RENTP;
    NET1T = NET1 + NET1P + RENT + RENTP;
    SIC1T = SIC1 + C2 + C5 + SIC1P + C2P + C5P;
    CDC6T = CDC6 + F2 + F5 + CDC6P + F2P + F5P;
    CKIT = SIC1T + CDC6T;
    CDC15i = CDC15T - CDC15;
    IE = IET - IEP;
    PE = ESP1T - ESP1;
    TEM1GDP = TEM1T - TEM1GTP;

    // REACTIONS
    // Growth
    Growth = mu * MASS;

    // Synthesis of CLN2
    Synthesis_of_CLN2 = (ksn2_p + ksn2_p_p * SBF) * MASS;

    // Degradation of CLN2
    Degradation_of_CLN2 = Mass_Action_1_222(kdn2, CLN2);

    // Synthesis of CLB2
    Synthesis_of_CLB2 = (ksb2_p + ksb2_p_p * MCM1) * MASS;

    // Degradation of CLB2
    Degradation_of_CLB2 = Mass_Action_1_222(Vdb2, CLB2);

    // Synthesis of CLB5
    Synthesis_of_CLB5 = (ksb5_p + ksb5_p_p * SBF) * MASS;

    // Degradation of CLB5
    Degradation_of_CLB5 = Mass_Action_1_222(Vdb5, CLB5);

    // Synthesis of SIC1
    Synthesis_of_SIC1 = ksc1_p + ksc1_p_p * SWI5;

    // Phosphorylation of SIC1
    Phosphorylation_of_SIC1 = Mass_Action_1_222(Vkpc1, SIC1);

    // Dephosphorylation of SIC1
    Dephosphorylation_of_SIC1 = Mass_Action_1_222(Vppc1, SIC1P);

    // Fast Degradation of SIC1P
    Fast_Degradation_of_SIC1P = Mass_Action_1_222(kd3c1, SIC1P);

    // Assoc. of CLB2 and SIC1
    Assoc_of_CLB2_and_SIC1 = Mass_Action_2_221(kasb2, CLB2, SIC1);

    // Dissoc. of CLB2/SIC1 complex
    Dissoc_of_CLB2SIC1_complex = Mass_Action_1_222(kdib2, C2);

    // Assoc. of CLB5 and SIC1
    Assoc_of_CLB5_and_SIC1 = Mass_Action_2_221(kasb5, CLB5, SIC1);

    // Dissoc. of CLB5/SIC1
    Dissoc_of_CLB5SIC1 = Mass_Action_1_222(kdib5, C5);

    // Phosphorylation of C2
    Phosphorylation_of_C2 = Mass_Action_1_222(Vkpc1, C2);

    // Dephosphorylation of C2P
    Dephosphorylation_of_C2P = Mass_Action_1_222(Vppc1, C2P);

    // Phosphorylation of C5
    Phosphorylation_of_C5 = Mass_Action_1_222(Vkpc1, C5);

    // Dephosphorylation of C5P
    Dephosphorylation_of_C5P = Mass_Action_1_222(Vppc1, C5P);

    // Degradation of CLB2 in C2
    Degradation_of_CLB2_in_C2 = Mass_Action_1_222(Vdb2, C2);

    // Degradation of CLB5 in C5
    Degradation_of_CLB5_in_C5 = Mass_Action_1_222(Vdb5, C5);

    // Degradation of SIC1 in C2P
    Degradation_of_SIC1_in_C2P = Mass_Action_1_222(kd3c1, C2P);

    // Degradation of SIC1P in C5P
    Degradation_of_SIC1P_in_C5P_ = Mass_Action_1_222(kd3c1, C5P);

    // Degradation of CLB2 in C2P
    Degradation_of_CLB2_in_C2P = Mass_Action_1_222(Vdb2, C2P);

    // Degradation of CLB5 in C5P
    Degradation_of_CLB5_in_C5P = Mass_Action_1_222(Vdb5, C5P);

    // CDC6 synthesis
    CDC6_synthesis = ksf6_p + ksf6_p_p * SWI5 + ksf6_p_p_p * SBF;

    // Phosphorylation of CDC6
    Phosphorylation_of_CDC6 = Mass_Action_1_222(Vkpf6, CDC6);

    // Dephosphorylation of CDC6
    Dephosphorylation_of_CDC6 = Mass_Action_1_222(Vppf6, CDC6P);

    // Degradation of CDC6P
    Degradation_of_CDC6P = Mass_Action_1_222(kd3f6, CDC6P);

    // CLB2/CDC6 complex formation
    CLB2CDC6_complex_formation = Mass_Action_2_221(kasf2, CLB2, CDC6);

    // CLB2/CDC6 dissociation
    CLB2CDC6_dissociation = Mass_Action_1_222(kdif2, F2);

    // CLB5/CDC6 complex formation
    CLB5CDC6_complex_formation = Mass_Action_2_221(kasf5, CLB5, CDC6);

    // CLB5/CDC6 dissociation
    CLB5CDC6_dissociation = Mass_Action_1_222(kdif5, F5);

    // F2 phosphorylation
    F2_phosphorylation = Mass_Action_1_222(Vkpf6, F2);

    // F2P dephosphorylation
    F2P_dephosphorylation = Mass_Action_1_222(Vppf6, F2P);

    // F5 phosphorylation
    F5_phosphorylation = Mass_Action_1_222(Vkpf6, F5);

    // F5P dephosphorylation
    F5P_dephosphorylation = Mass_Action_1_222(Vppf6, F5P);

    // CLB2 degradation in F2
    CLB2_degradation_in_F2 = Mass_Action_1_222(Vdb2, F2);

    // CLB5 degradation in F5
    CLB5_degradation_in_F5 = Mass_Action_1_222(Vdb5, F5);

    // CDC6 degradation in F2P
    CDC6_degradation_in_F2P = Mass_Action_1_222(kd3f6, F2P);

    // CDC6 degradation in F5P
    CDC6_degradation_in_F5P = Mass_Action_1_222(kd3f6, F5P);

    // CLB2 degradation in F2P
    CLB2_degradation_in_F2P = Mass_Action_1_222(Vdb2, F2P);

    // CLB5 degradation in F5P
    CLB5_degradation_in_F5P = Mass_Action_1_222(Vdb5, F5P);

    // Synthesis of SWI5
    Synthesis_of_SWI5 = ksswi_p + ksswi_p_p * MCM1;

    // Degradation of SWI5
    Degradation_of_SWI5 = Mass_Action_1_222(kdswi, SWI5);

    // Degradation of SWI5P
    Degradation_of_SWI5P = Mass_Action_1_222(kdswi, SWI5P);

    // Activation of SWI5
    Activation_of_SWI5 = Mass_Action_1_222(kaswi * CDC14, SWI5P);

    // Inactivation of SWI5
    Inactivation_of_SWI5 = Mass_Action_1_222(kiswi * CLB2, SWI5);

    // Activation of IEP
    Activation_of_IEP = MichaelisMenten_220(Vaiep, Jaiep, 1.0, IE);

    // Inactivation
    Inactivation_1 = MichaelisMenten_220(1.0, Jiiep, kiiep, IEP);

    // Synthesis of inactive CDC20
    Synthesis_of_inactive_CDC20 = ks20_p + ks20_p_p * MCM1;

    // Degradation of inactiveCDC20
    Degradation_of_inactiveCDC20 = Mass_Action_1_222(kd20, CDC20i);

    // Degradation of active CDC20
    Degradation_of_active_CDC20 = Mass_Action_1_222(kd20, CDC20);

    // Activation of CDC20
    Activation_of_CDC20 = Mass_Action_1_222(ka20_p + ka20_p_p * IEP, CDC20i);

    // Inactivation
    Inactivation_2 = 0.0;
    {
        double k = 1.0;
        Inactivation_2 = k * Mass_Action_1_222(MAD2, CDC20);
    }

    // CDH1 synthesis
    CDH1_synthesis = kscdh;

    // CDH1 degradation
    CDH1_degradation = Mass_Action_1_222(kdcdh, CDH1);

    // CDH1i degradation
    CDH1i_degradation = Mass_Action_1_222(kdcdh, CDH1i);

    // CDH1i activation
    CDH1i_activation = MichaelisMenten_220(Vacdh, Jacdh, 1.0, CDH1i);

    // Inactivation
    Inactivation_3 = MichaelisMenten_220(Vicdh, Jicdh, 1.0, CDH1);

    // CDC14 synthesis
    CDC14_synthesis = ks14;

    // CDC14 degradation
    CDC14_degradation = Mass_Action_1_222(kd14, CDC14);

    // Assoc. with NET1 to form RENT
    Assoc_with_NET1_to_form_RENT = Mass_Action_2_221(kasrent, CDC14, NET1);

    // Dissoc. from RENT
    Dissoc_from_RENT = Mass_Action_1_222(kdirent, RENT);

    // Assoc with NET1P to form RENTP
    Assoc_with_NET1P_to_form_RENTP = Mass_Action_2_221(kasrentp, CDC14, NET1P);

    // Dissoc. from RENP
    Dissoc_from_RENP = Mass_Action_1_222(kdirentp, RENTP);

    // Net1 synthesis
    Net1_synthesis = ksnet;

    // Net1 degradation
    Net1_degradation = Mass_Action_1_222(kdnet, NET1);

    // Net1P degradation
    Net1P_degradation = Mass_Action_1_222(kdnet, NET1P);

    // NET1 phosphorylation
    NET1_phosphorylation = Mass_Action_1_222(Vkpnet, NET1);

    // dephosphorylation
    dephosphorylation_1 = Mass_Action_1_222(Vppnet, NET1P);

    // RENT phosphorylation
    RENT_phosphorylation = Mass_Action_1_222(Vkpnet, RENT);

    // dephosphorylation
    dephosphorylation_2 = Mass_Action_1_222(Vppnet, RENTP);

    // Degradation of NET1 in RENT
    Degradation_of_NET1_in_RENT = Mass_Action_1_222(kdnet, RENT);

    // Degradation of NET1P in RENTP
    Degradation_of_NET1P_in_RENTP = Mass_Action_1_222(kdnet, RENTP);

    // Degradation of CDC14 in RENT
    Degradation_of_CDC14_in_RENT = Mass_Action_1_222(kd14, RENT);

    // Degradation of CDC14 in RENTP
    Degradation_of_CDC14_in_RENTP = Mass_Action_1_222(kd14, RENTP);

    // TEM1 activation
    TEM1_activation = MichaelisMenten_220(LTE1, Jatem, 1.0, TEM1GDP);

    // inactivation
    inactivation_1 = MichaelisMenten_220(BUB2, Jitem, 1.0, TEM1GTP);

    // CDC15 activation
    CDC15_activation = Mass_Action_1_222(ka15_p * TEM1GDP + ka15_p_p * TEM1GTP + ka15p * CDC14, CDC15i);

    // inactivation
    inactivation_2 = Mass_Action_1_222(ki15, CDC15);

    // PPX synthesis
    PPX_synthesis = ksppx;

    // degradation
    degradation_1 = Mass_Action_1_222(Vdppx, PPX);

    // PDS1 synthesis
    PDS1_synthesis = kspds_p + ks1pds_p_p * SBF + ks2pds_p_p * MCM1;

    // degradation
    degradation_2 = Mass_Action_1_222(Vdpds, PDS1);

    // Degradation of PDS1 in PE
    Degradation_of_PDS1_in_PE = Mass_Action_1_222(Vdpds, PE);

    // Assoc. with ESP1 to form PE
    Assoc_with_ESP1_to_form_PE = Mass_Action_2_221(kasesp, PDS1, ESP1);

    // Disso. from PE
    Disso_from_PE = Mass_Action_1_222(kdiesp, PE);

    // DNA synthesis
    DNA_synthesis = ksori * (eorib5 * CLB5 + eorib2 * CLB2);

    // Negative regulation of DNA synthesis
    Negative_regulation_of_DNA_synthesis = Mass_Action_1_222(kdori, ORI);

    // Budding
    Budding = ksbud * (ebudn2 * CLN2 + ebudn3 * CLN3 + ebudb5 * CLB5);

    // Negative regulation of Cell budding
    Negative_regulation_of_Cell_budding = Mass_Action_1_222(kdbud, BUD);

    // Spindle formation
    Spindle_formation = ksspn * CLB2 / (Jspn + CLB2);

    // Spindle disassembly
    Spindle_disassembly = Mass_Action_1_222(kdspn, SPN);
}

// MODEL FUNCTIONS
inline double Chen2004SbmlOdeSystem::BB_218(double A1, double A2, double A3, double A4)
{
    return A2 - A1 + A3 * A2 + A4 * A1;
}

inline double Chen2004SbmlOdeSystem::GK_219(double A1, double A2, double A3, double A4)
{
    return 2.0 * A4 * A1 / (A2 - A1 + A3 * A2 + A4 * A1 + sm::root(2.0, std::pow(A2 - A1 + A3 * A2 + A4 * A1, 2.0) - 4.0 * (A2 - A1) * A4 * A1));
}

inline double Chen2004SbmlOdeSystem::MichaelisMenten_220(double M1, double J1, double k1, double S1)
{
    return k1 * S1 * M1 / (J1 + S1);
}

inline double Chen2004SbmlOdeSystem::Mass_Action_2_221(double k1, double S1, double S2)
{
    return k1 * S1 * S2;
}

inline double Chen2004SbmlOdeSystem::Mass_Action_1_222(double k1, double S1)
{
    return k1 * S1;
}

template <>
void CellwiseOdeSystemInformation<Chen2004SbmlOdeSystem>::Initialise()
{
    // STATE VARIABLES
    this->mVariableNames.push_back("BUD");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.008473);

    this->mVariableNames.push_back("C2");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.238404);

    this->mVariableNames.push_back("C2P");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.024034);

    this->mVariableNames.push_back("C5");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.070081);

    this->mVariableNames.push_back("C5P");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.006878);

    this->mVariableNames.push_back("CDC14");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.468344);

    this->mVariableNames.push_back("CDC15");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.656533);

    this->mVariableNames.push_back("CDC20");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.444296);

    this->mVariableNames.push_back("CDC20i");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(1.472044);

    this->mVariableNames.push_back("CDC6");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.10758);

    this->mVariableNames.push_back("CDC6P");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.015486);

    this->mVariableNames.push_back("CDH1");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.930499);

    this->mVariableNames.push_back("CDH1i");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0695);

    this->mVariableNames.push_back("CLB2");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.1469227);

    this->mVariableNames.push_back("CLB5");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0518014);

    this->mVariableNames.push_back("CLN2");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0652511);

    this->mVariableNames.push_back("ESP1");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.301313);

    this->mVariableNames.push_back("F2");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.236058);

    this->mVariableNames.push_back("F2P");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0273938);

    this->mVariableNames.push_back("F5");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(7.24e-05);

    this->mVariableNames.push_back("F5P");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(7.91e-05);

    this->mVariableNames.push_back("IEP");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.1015);

    this->mVariableNames.push_back("MASS");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(1.206019);

    this->mVariableNames.push_back("NET1");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.018645);

    this->mVariableNames.push_back("NET1P");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.970271);

    this->mVariableNames.push_back("ORI");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.000909);

    this->mVariableNames.push_back("PDS1");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.025612);

    this->mVariableNames.push_back("PPX");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.123179);

    this->mVariableNames.push_back("RENT");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(1.04954);

    this->mVariableNames.push_back("RENTP");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.6);

    this->mVariableNames.push_back("SIC1");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0228776);

    this->mVariableNames.push_back("SIC1P");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.00641);

    this->mVariableNames.push_back("SPN");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.03);

    this->mVariableNames.push_back("SWI5");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.95);

    this->mVariableNames.push_back("SWI5P");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.02);

    this->mVariableNames.push_back("TEM1GTP");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.9);

    // DERIVED QUANTITIES
    this->mDerivedQuantityNames.push_back("BCK2");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("CDC14T");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("CDC15i");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("CDC6T");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("CKIT");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("CLB2T");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("CLB5T");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("CLN3");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("IE");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("MCM1");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("NET1T");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("PE");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("SBF");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("SIC1T");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("TEM1GDP");
    this->mDerivedQuantityUnits.push_back("non-dim");

    // PARAMETERS
    this->mParameterNames.push_back("cell");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("BUB2");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("LTE1");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("MAD2");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("bub2l");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("CDC15T");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ESP1T");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("IET");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("KEZ");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("KEZ2");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("lte1h");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("lte1l");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("mad2l");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("TEM1T");
    this->mParameterUnits.push_back("non-dim");

    this->mInitialised = true;
}

// Register the ODE system with Boost serialization
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(Chen2004SbmlOdeSystem)