#include <cmath>
#include <limits>

#include "CellwiseOdeSystemInformation.hpp"
#include "SbmlMath.hpp"

#include "Chen2004SbmlOdeSystemAndCellCycleModel.hpp"

namespace sm = sbmlmath;

Chen2004SbmlOdeSystem::Chen2004SbmlOdeSystem(std::vector<double> stateVariables)
    : AbstractOdeSystem(51)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<Chen2004SbmlOdeSystem>);

    // COMPARTMENTS:
    cell = 1.0;

    // STATE VARIABLES:
    BCK2 = 0.0; // BCK2
    BUD = 0.008473; // BUD
    C2 = 0.238404; // C2
    C2P = 0.024034; // C2P
    C5 = 0.070081; // C5
    C5P = 0.006878; // C5P
    CDC14 = 0.468344; // CDC14
    CDC14T = 2.0; // CDC14T
    CDC15 = 0.656533; // CDC15
    CDC15i = 0.0; // CDC15i
    CDC20 = 0.444296; // CDC20
    CDC20i = 1.472044; // CDC20i
    CDC6 = 0.10758; // CDC6
    CDC6P = 0.015486; // CDC6P
    CDC6T = 0.0; // CDC6T
    CDH1 = 0.930499; // CDH1
    CDH1i = 0.0695; // CDH1i
    CKIT = 0.0; // CKIT
    CLB2 = 0.1469227; // CLB2
    CLB2T = 0.17; // CLB2T
    CLB5 = 0.0518014; // CLB5
    CLB5T = 0.12; // CLB5T
    CLN2 = 0.0652511; // CLN2
    CLN3 = 0.0; // CLN3
    ESP1 = 0.301313; // ESP1
    F2 = 0.236058; // F2
    F2P = 0.0273938; // F2P
    F5 = 7.24e-05; // F5
    F5P = 7.91e-05; // F5P
    IE = 0.0; // IE
    IEP = 0.1015; // IEP
    MASS = 1.206019; // MASS
    MCM1 = 0.0; // MCM1
    NET1 = 0.018645; // NET1
    NET1P = 0.970271; // NET1P
    NET1T = 2.8; // NET1T
    ORI = 0.000909; // ORI
    PDS1 = 0.025612; // PDS1
    PE = 0.0; // PE
    PPX = 0.123179; // PPX
    RENT = 1.04954; // RENT
    RENTP = 0.6; // RENTP
    SBF = 0.0; // SBF
    SIC1 = 0.0228776; // SIC1
    SIC1P = 0.00641; // SIC1P
    SIC1T = 0.0; // SIC1T
    SPN = 0.03; // SPN
    SWI5 = 0.95; // SWI5
    SWI5P = 0.02; // SWI5P
    TEM1GDP = 0.0; // TEM1GDP
    TEM1GTP = 0.9; // TEM1GTP

    SetDefaultInitialCondition(0, BCK2);
    SetDefaultInitialCondition(1, BUD);
    SetDefaultInitialCondition(2, C2);
    SetDefaultInitialCondition(3, C2P);
    SetDefaultInitialCondition(4, C5);
    SetDefaultInitialCondition(5, C5P);
    SetDefaultInitialCondition(6, CDC14);
    SetDefaultInitialCondition(7, CDC14T);
    SetDefaultInitialCondition(8, CDC15);
    SetDefaultInitialCondition(9, CDC15i);
    SetDefaultInitialCondition(10, CDC20);
    SetDefaultInitialCondition(11, CDC20i);
    SetDefaultInitialCondition(12, CDC6);
    SetDefaultInitialCondition(13, CDC6P);
    SetDefaultInitialCondition(14, CDC6T);
    SetDefaultInitialCondition(15, CDH1);
    SetDefaultInitialCondition(16, CDH1i);
    SetDefaultInitialCondition(17, CKIT);
    SetDefaultInitialCondition(18, CLB2);
    SetDefaultInitialCondition(19, CLB2T);
    SetDefaultInitialCondition(20, CLB5);
    SetDefaultInitialCondition(21, CLB5T);
    SetDefaultInitialCondition(22, CLN2);
    SetDefaultInitialCondition(23, CLN3);
    SetDefaultInitialCondition(24, ESP1);
    SetDefaultInitialCondition(25, F2);
    SetDefaultInitialCondition(26, F2P);
    SetDefaultInitialCondition(27, F5);
    SetDefaultInitialCondition(28, F5P);
    SetDefaultInitialCondition(29, IE);
    SetDefaultInitialCondition(30, IEP);
    SetDefaultInitialCondition(31, MASS);
    SetDefaultInitialCondition(32, MCM1);
    SetDefaultInitialCondition(33, NET1);
    SetDefaultInitialCondition(34, NET1P);
    SetDefaultInitialCondition(35, NET1T);
    SetDefaultInitialCondition(36, ORI);
    SetDefaultInitialCondition(37, PDS1);
    SetDefaultInitialCondition(38, PE);
    SetDefaultInitialCondition(39, PPX);
    SetDefaultInitialCondition(40, RENT);
    SetDefaultInitialCondition(41, RENTP);
    SetDefaultInitialCondition(42, SBF);
    SetDefaultInitialCondition(43, SIC1);
    SetDefaultInitialCondition(44, SIC1P);
    SetDefaultInitialCondition(45, SIC1T);
    SetDefaultInitialCondition(46, SPN);
    SetDefaultInitialCondition(47, SWI5);
    SetDefaultInitialCondition(48, SWI5P);
    SetDefaultInitialCondition(49, TEM1GDP);
    SetDefaultInitialCondition(50, TEM1GTP);

    if (stateVariables.size() == 51)
    {
        BCK2 = stateVariables[0];
        BUD = stateVariables[1];
        C2 = stateVariables[2];
        C2P = stateVariables[3];
        C5 = stateVariables[4];
        C5P = stateVariables[5];
        CDC14 = stateVariables[6];
        CDC14T = stateVariables[7];
        CDC15 = stateVariables[8];
        CDC15i = stateVariables[9];
        CDC20 = stateVariables[10];
        CDC20i = stateVariables[11];
        CDC6 = stateVariables[12];
        CDC6P = stateVariables[13];
        CDC6T = stateVariables[14];
        CDH1 = stateVariables[15];
        CDH1i = stateVariables[16];
        CKIT = stateVariables[17];
        CLB2 = stateVariables[18];
        CLB2T = stateVariables[19];
        CLB5 = stateVariables[20];
        CLB5T = stateVariables[21];
        CLN2 = stateVariables[22];
        CLN3 = stateVariables[23];
        ESP1 = stateVariables[24];
        F2 = stateVariables[25];
        F2P = stateVariables[26];
        F5 = stateVariables[27];
        F5P = stateVariables[28];
        IE = stateVariables[29];
        IEP = stateVariables[30];
        MASS = stateVariables[31];
        MCM1 = stateVariables[32];
        NET1 = stateVariables[33];
        NET1P = stateVariables[34];
        NET1T = stateVariables[35];
        ORI = stateVariables[36];
        PDS1 = stateVariables[37];
        PE = stateVariables[38];
        PPX = stateVariables[39];
        RENT = stateVariables[40];
        RENTP = stateVariables[41];
        SBF = stateVariables[42];
        SIC1 = stateVariables[43];
        SIC1P = stateVariables[44];
        SIC1T = stateVariables[45];
        SPN = stateVariables[46];
        SWI5 = stateVariables[47];
        SWI5P = stateVariables[48];
        TEM1GDP = stateVariables[49];
        TEM1GTP = stateVariables[50];
    }
    else if (stateVariables.size() != 0)
    {
        EXCEPTION("Chen2004SbmlOdeSystem: Expected 51 state variables, got " + std::to_string(stateVariables.size()));
    }

    mStateVariables.push_back(BCK2);
    mStateVariables.push_back(BUD);
    mStateVariables.push_back(C2);
    mStateVariables.push_back(C2P);
    mStateVariables.push_back(C5);
    mStateVariables.push_back(C5P);
    mStateVariables.push_back(CDC14);
    mStateVariables.push_back(CDC14T);
    mStateVariables.push_back(CDC15);
    mStateVariables.push_back(CDC15i);
    mStateVariables.push_back(CDC20);
    mStateVariables.push_back(CDC20i);
    mStateVariables.push_back(CDC6);
    mStateVariables.push_back(CDC6P);
    mStateVariables.push_back(CDC6T);
    mStateVariables.push_back(CDH1);
    mStateVariables.push_back(CDH1i);
    mStateVariables.push_back(CKIT);
    mStateVariables.push_back(CLB2);
    mStateVariables.push_back(CLB2T);
    mStateVariables.push_back(CLB5);
    mStateVariables.push_back(CLB5T);
    mStateVariables.push_back(CLN2);
    mStateVariables.push_back(CLN3);
    mStateVariables.push_back(ESP1);
    mStateVariables.push_back(F2);
    mStateVariables.push_back(F2P);
    mStateVariables.push_back(F5);
    mStateVariables.push_back(F5P);
    mStateVariables.push_back(IE);
    mStateVariables.push_back(IEP);
    mStateVariables.push_back(MASS);
    mStateVariables.push_back(MCM1);
    mStateVariables.push_back(NET1);
    mStateVariables.push_back(NET1P);
    mStateVariables.push_back(NET1T);
    mStateVariables.push_back(ORI);
    mStateVariables.push_back(PDS1);
    mStateVariables.push_back(PE);
    mStateVariables.push_back(PPX);
    mStateVariables.push_back(RENT);
    mStateVariables.push_back(RENTP);
    mStateVariables.push_back(SBF);
    mStateVariables.push_back(SIC1);
    mStateVariables.push_back(SIC1P);
    mStateVariables.push_back(SIC1T);
    mStateVariables.push_back(SPN);
    mStateVariables.push_back(SWI5);
    mStateVariables.push_back(SWI5P);
    mStateVariables.push_back(TEM1GDP);
    mStateVariables.push_back(TEM1GTP);

    // STATE PARAMETERS:
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
    mParameters.push_back(D);
    mParameters.push_back(mu);
    mParameters.push_back(Vdb5);
    mParameters.push_back(Vdb2);
    mParameters.push_back(Vasbf);
    mParameters.push_back(Visbf);
    mParameters.push_back(Vkpc1);
    mParameters.push_back(Vkpf6);
    mParameters.push_back(Vacdh);
    mParameters.push_back(Vicdh);
    mParameters.push_back(Vppnet);
    mParameters.push_back(Vkpnet);
    mParameters.push_back(Vdppx);
    mParameters.push_back(Vdpds);
    mParameters.push_back(Vaiep);
    mParameters.push_back(Vd2c1);
    mParameters.push_back(Vd2f6);
    mParameters.push_back(Vppc1);
    mParameters.push_back(Vppf6);
    mParameters.push_back(F);

    // EVENTS:
    eventsSatisfied.resize(4, false);
    eventsInitialised = false;
}

Chen2004SbmlOdeSystem::~Chen2004SbmlOdeSystem()
{
}

void Chen2004SbmlOdeSystem::RefreshState(const std::vector<double> &rY)
{
    // STATE VARIABLES:
    BCK2 = rY[0];
    BUD = rY[1];
    C2 = rY[2];
    C2P = rY[3];
    C5 = rY[4];
    C5P = rY[5];
    CDC14 = rY[6];
    CDC14T = rY[7];
    CDC15 = rY[8];
    CDC15i = rY[9];
    CDC20 = rY[10];
    CDC20i = rY[11];
    CDC6 = rY[12];
    CDC6P = rY[13];
    CDC6T = rY[14];
    CDH1 = rY[15];
    CDH1i = rY[16];
    CKIT = rY[17];
    CLB2 = rY[18];
    CLB2T = rY[19];
    CLB5 = rY[20];
    CLB5T = rY[21];
    CLN2 = rY[22];
    CLN3 = rY[23];
    ESP1 = rY[24];
    F2 = rY[25];
    F2P = rY[26];
    F5 = rY[27];
    F5P = rY[28];
    IE = rY[29];
    IEP = rY[30];
    MASS = rY[31];
    MCM1 = rY[32];
    NET1 = rY[33];
    NET1P = rY[34];
    NET1T = rY[35];
    ORI = rY[36];
    PDS1 = rY[37];
    PE = rY[38];
    PPX = rY[39];
    RENT = rY[40];
    RENTP = rY[41];
    SBF = rY[42];
    SIC1 = rY[43];
    SIC1P = rY[44];
    SIC1T = rY[45];
    SPN = rY[46];
    SWI5 = rY[47];
    SWI5P = rY[48];
    TEM1GDP = rY[49];
    TEM1GTP = rY[50];

    // STATE PARAMETERS:
    BUB2 = GetParameter("BUB2");
    LTE1 = GetParameter("LTE1");
    MAD2 = GetParameter("MAD2");

    bub2l = GetParameter("bub2l");
    CDC15T = GetParameter("CDC15T");
    ESP1T = GetParameter("ESP1T");
    IET = GetParameter("IET");
    KEZ = GetParameter("KEZ");
    KEZ2 = GetParameter("KEZ2");
    lte1h = GetParameter("lte1h");
    lte1l = GetParameter("lte1l");
    mad2l = GetParameter("mad2l");
    TEM1T = GetParameter("TEM1T");
    D = GetParameter("D");
    mu = GetParameter("mu");
    Vdb5 = GetParameter("Vdb5");
    Vdb2 = GetParameter("Vdb2");
    Vasbf = GetParameter("Vasbf");
    Visbf = GetParameter("Visbf");
    Vkpc1 = GetParameter("Vkpc1");
    Vkpf6 = GetParameter("Vkpf6");
    Vacdh = GetParameter("Vacdh");
    Vicdh = GetParameter("Vicdh");
    Vppnet = GetParameter("Vppnet");
    Vkpnet = GetParameter("Vkpnet");
    Vdppx = GetParameter("Vdppx");
    Vdpds = GetParameter("Vdpds");
    Vaiep = GetParameter("Vaiep");
    Vd2c1 = GetParameter("Vd2c1");
    Vd2f6 = GetParameter("Vd2f6");
    Vppc1 = GetParameter("Vppc1");
    Vppf6 = GetParameter("Vppf6");
    F = GetParameter("F");
}

void Chen2004SbmlOdeSystem::EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY)
{
    RefreshState(rY);

    // RULES:
    Visbf = kisbf_p + kisbf_p_p * CLB2;
    Vppc1 = kppc1 * CDC14;
    Vppf6 = kppf6 * CDC14;
    Vaiep = kaiep * CLB2;
    Vacdh = kacdh_p + kacdh_p_p * CDC14;
    Vicdh = kicdh_p + kicdh_p_p * (eicdhn3 * CLN3 + eicdhn2 * CLN2 + eicdhb5 * CLB5 + eicdhb2 * CLB2);
    Vkpnet = (kkpnet_p + kkpnet_p_p * CDC15) * MASS;
    Vppnet = kppnet_p + kppnet_p_p * PPX;
    Vasbf = kasbf * (esbfn2 * CLN2 + esbfn3 * (CLN3 + BCK2) + esbfb5 * CLB5);
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

    // UPDATE STATE PARAMETERS:
    SetParameter("BUB2", BUB2);
    SetParameter("LTE1", LTE1);
    SetParameter("MAD2", MAD2);

    SetParameter("bub2l", bub2l);
    SetParameter("CDC15T", CDC15T);
    SetParameter("ESP1T", ESP1T);
    SetParameter("IET", IET);
    SetParameter("KEZ", KEZ);
    SetParameter("KEZ2", KEZ2);
    SetParameter("lte1h", lte1h);
    SetParameter("lte1l", lte1l);
    SetParameter("mad2l", mad2l);
    SetParameter("TEM1T", TEM1T);
    SetParameter("D", D);
    SetParameter("mu", mu);
    SetParameter("Vdb5", Vdb5);
    SetParameter("Vdb2", Vdb2);
    SetParameter("Vasbf", Vasbf);
    SetParameter("Visbf", Visbf);
    SetParameter("Vkpc1", Vkpc1);
    SetParameter("Vkpf6", Vkpf6);
    SetParameter("Vacdh", Vacdh);
    SetParameter("Vicdh", Vicdh);
    SetParameter("Vppnet", Vppnet);
    SetParameter("Vkpnet", Vkpnet);
    SetParameter("Vdppx", Vdppx);
    SetParameter("Vdpds", Vdpds);
    SetParameter("Vaiep", Vaiep);
    SetParameter("Vd2c1", Vd2c1);
    SetParameter("Vd2f6", Vd2f6);
    SetParameter("Vppc1", Vppc1);
    SetParameter("Vppf6", Vppf6);
    SetParameter("F", F);

    // REACTIONS:

    // Growth
    double Growth = mu * MASS;

    // Synthesis of CLN2
    double Synthesis_of_CLN2 = (ksn2_p + ksn2_p_p * SBF) * MASS;

    // Degradation of CLN2
    double Degradation_of_CLN2 = Mass_Action_1_222(kdn2, CLN2);

    // Synthesis of CLB2
    double Synthesis_of_CLB2 = (ksb2_p + ksb2_p_p * MCM1) * MASS;

    // Degradation of CLB2
    double Degradation_of_CLB2 = Mass_Action_1_222(Vdb2, CLB2);

    // Synthesis of CLB5
    double Synthesis_of_CLB5 = (ksb5_p + ksb5_p_p * SBF) * MASS;

    // Degradation of CLB5
    double Degradation_of_CLB5 = Mass_Action_1_222(Vdb5, CLB5);

    // Synthesis of SIC1
    double Synthesis_of_SIC1 = ksc1_p + ksc1_p_p * SWI5;

    // Phosphorylation of SIC1
    double Phosphorylation_of_SIC1 = Mass_Action_1_222(Vkpc1, SIC1);

    // Dephosphorylation of SIC1
    double Dephosphorylation_of_SIC1 = Mass_Action_1_222(Vppc1, SIC1P);

    // Fast Degradation of SIC1P
    double Fast_Degradation_of_SIC1P = Mass_Action_1_222(kd3c1, SIC1P);

    // Assoc. of CLB2 and SIC1
    double Assoc_of_CLB2_and_SIC1 = Mass_Action_2_221(kasb2, CLB2, SIC1);

    // Dissoc. of CLB2/SIC1 complex
    double Dissoc_of_CLB2SIC1_complex = Mass_Action_1_222(kdib2, C2);

    // Assoc. of CLB5 and SIC1
    double Assoc_of_CLB5_and_SIC1 = Mass_Action_2_221(kasb5, CLB5, SIC1);

    // Dissoc. of CLB5/SIC1
    double Dissoc_of_CLB5SIC1 = Mass_Action_1_222(kdib5, C5);

    // Phosphorylation of C2
    double Phosphorylation_of_C2 = Mass_Action_1_222(Vkpc1, C2);

    // Dephosphorylation of C2P
    double Dephosphorylation_of_C2P = Mass_Action_1_222(Vppc1, C2P);

    // Phosphorylation of C5
    double Phosphorylation_of_C5 = Mass_Action_1_222(Vkpc1, C5);

    // Dephosphorylation of C5P
    double Dephosphorylation_of_C5P = Mass_Action_1_222(Vppc1, C5P);

    // Degradation of CLB2 in C2
    double Degradation_of_CLB2_in_C2 = Mass_Action_1_222(Vdb2, C2);

    // Degradation of CLB5 in C5
    double Degradation_of_CLB5_in_C5 = Mass_Action_1_222(Vdb5, C5);

    // Degradation of SIC1 in C2P
    double Degradation_of_SIC1_in_C2P = Mass_Action_1_222(kd3c1, C2P);

    // Degradation of SIC1P in C5P
    double Degradation_of_SIC1P_in_C5P = Mass_Action_1_222(kd3c1, C5P);

    // Degradation of CLB2 in C2P
    double Degradation_of_CLB2_in_C2P = Mass_Action_1_222(Vdb2, C2P);

    // Degradation of CLB5 in C5P
    double Degradation_of_CLB5_in_C5P = Mass_Action_1_222(Vdb5, C5P);

    // CDC6 synthesis
    double CDC6_synthesis = ksf6_p + ksf6_p_p * SWI5 + ksf6_p_p_p * SBF;

    // Phosphorylation of CDC6
    double Phosphorylation_of_CDC6 = Mass_Action_1_222(Vkpf6, CDC6);

    // Dephosphorylation of CDC6
    double Dephosphorylation_of_CDC6 = Mass_Action_1_222(Vppf6, CDC6P);

    // Degradation of CDC6P
    double Degradation_of_CDC6P = Mass_Action_1_222(kd3f6, CDC6P);

    // CLB2/CDC6 complex formation
    double CLB2CDC6_complex_formation = Mass_Action_2_221(kasf2, CLB2, CDC6);

    // CLB2/CDC6 dissociation
    double CLB2CDC6_dissociation = Mass_Action_1_222(kdif2, F2);

    // CLB5/CDC6 complex formation
    double CLB5CDC6_complex_formation = Mass_Action_2_221(kasf5, CLB5, CDC6);

    // CLB5/CDC6 dissociation
    double CLB5CDC6_dissociation = Mass_Action_1_222(kdif5, F5);

    // F2 phosphorylation
    double F2_phosphorylation = Mass_Action_1_222(Vkpf6, F2);

    // F2P dephosphorylation
    double F2P_dephosphorylation = Mass_Action_1_222(Vppf6, F2P);

    // F5 phosphorylation
    double F5_phosphorylation = Mass_Action_1_222(Vkpf6, F5);

    // F5P dephosphorylation
    double F5P_dephosphorylation = Mass_Action_1_222(Vppf6, F5P);

    // CLB2 degradation in F2
    double CLB2_degradation_in_F2 = Mass_Action_1_222(Vdb2, F2);

    // CLB5 degradation in F5
    double CLB5_degradation_in_F5 = Mass_Action_1_222(Vdb5, F5);

    // CDC6 degradation in F2P
    double CDC6_degradation_in_F2P = Mass_Action_1_222(kd3f6, F2P);

    // CDC6 degradation in F5P
    double CDC6_degradation_in_F5P = Mass_Action_1_222(kd3f6, F5P);

    // CLB2 degradation in F2P
    double CLB2_degradation_in_F2P = Mass_Action_1_222(Vdb2, F2P);

    // CLB5 degradation in F5P
    double CLB5_degradation_in_F5P = Mass_Action_1_222(Vdb5, F5P);

    // Synthesis of SWI5
    double Synthesis_of_SWI5 = ksswi_p + ksswi_p_p * MCM1;

    // Degradation of SWI5
    double Degradation_of_SWI5 = Mass_Action_1_222(kdswi, SWI5);

    // Degradation of SWI5P
    double Degradation_of_SWI5P = Mass_Action_1_222(kdswi, SWI5P);

    // Activation of SWI5
    double Activation_of_SWI5 = Mass_Action_1_222(kaswi * CDC14, SWI5P);

    // Inactivation of SWI5
    double Inactivation_of_SWI5 = Mass_Action_1_222(kiswi * CLB2, SWI5);

    // Activation of IEP
    double Activation_of_IEP = MichaelisMenten_220(Vaiep, Jaiep, 1.0, IE);

    // Inactivation
    double Inactivation = MichaelisMenten_220(1.0, Jiiep, kiiep, IEP);

    // Synthesis of inactive CDC20
    double Synthesis_of_inactive_CDC20 = ks20_p + ks20_p_p * MCM1;

    // Degradation of inactiveCDC20
    double Degradation_of_inactiveCDC20 = Mass_Action_1_222(kd20, CDC20i);

    // Degradation of active CDC20
    double Degradation_of_active_CDC20 = Mass_Action_1_222(kd20, CDC20);

    // Activation of CDC20
    double Activation_of_CDC20 = Mass_Action_1_222(ka20_p + ka20_p_p * IEP, CDC20i);

    // Inactivation
    double Inactivation_0 = 0.0;
    {
        double k = 1.0;
        Inactivation_0 = k * Mass_Action_1_222(MAD2, CDC20);
    }

    // CDH1 synthesis
    double CDH1_synthesis = kscdh;

    // CDH1 degradation
    double CDH1_degradation = Mass_Action_1_222(kdcdh, CDH1);

    // CDH1i degradation
    double CDH1i_degradation = Mass_Action_1_222(kdcdh, CDH1i);

    // CDH1i activation
    double CDH1i_activation = MichaelisMenten_220(Vacdh, Jacdh, 1.0, CDH1i);

    // Inactivation
    double Inactivation_1 = MichaelisMenten_220(Vicdh, Jicdh, 1.0, CDH1);

    // CDC14 synthesis
    double CDC14_synthesis = ks14;

    // CDC14 degradation
    double CDC14_degradation = Mass_Action_1_222(kd14, CDC14);

    // Assoc. with NET1 to form RENT
    double Assoc_with_NET1_to_form_RENT = Mass_Action_2_221(kasrent, CDC14, NET1);

    // Dissoc. from RENT
    double Dissoc_from_RENT = Mass_Action_1_222(kdirent, RENT);

    // Assoc with NET1P to form RENTP
    double Assoc_with_NET1P_to_form_RENTP = Mass_Action_2_221(kasrentp, CDC14, NET1P);

    // Dissoc. from RENP
    double Dissoc_from_RENP = Mass_Action_1_222(kdirentp, RENTP);

    // Net1 synthesis
    double Net1_synthesis = ksnet;

    // Net1 degradation
    double Net1_degradation = Mass_Action_1_222(kdnet, NET1);

    // Net1P degradation
    double Net1P_degradation = Mass_Action_1_222(kdnet, NET1P);

    // NET1 phosphorylation
    double NET1_phosphorylation = Mass_Action_1_222(Vkpnet, NET1);

    // dephosphorylation
    double dephosphorylation = Mass_Action_1_222(Vppnet, NET1P);

    // RENT phosphorylation
    double RENT_phosphorylation = Mass_Action_1_222(Vkpnet, RENT);

    // dephosphorylation
    double dephosphorylation_0 = Mass_Action_1_222(Vppnet, RENTP);

    // Degradation of NET1 in RENT
    double Degradation_of_NET1_in_RENT = Mass_Action_1_222(kdnet, RENT);

    // Degradation of NET1P in RENTP
    double Degradation_of_NET1P_in_RENTP = Mass_Action_1_222(kdnet, RENTP);

    // Degradation of CDC14 in RENT
    double Degradation_of_CDC14_in_RENT = Mass_Action_1_222(kd14, RENT);

    // Degradation of CDC14 in RENTP
    double Degradation_of_CDC14_in_RENTP = Mass_Action_1_222(kd14, RENTP);

    // TEM1 activation
    double TEM1_activation = MichaelisMenten_220(LTE1, Jatem, 1.0, TEM1GDP);

    // inactivation
    double inactivation = MichaelisMenten_220(BUB2, Jitem, 1.0, TEM1GTP);

    // CDC15 activation
    double CDC15_activation = Mass_Action_1_222(ka15_p * TEM1GDP + ka15_p_p * TEM1GTP + ka15p * CDC14, CDC15i);

    // inactivation
    double inactivation_0 = Mass_Action_1_222(ki15, CDC15);

    // PPX synthesis
    double PPX_synthesis = ksppx;

    // degradation
    double degradation = Mass_Action_1_222(Vdppx, PPX);

    // PDS1 synthesis
    double PDS1_synthesis = kspds_p + ks1pds_p_p * SBF + ks2pds_p_p * MCM1;

    // degradation
    double degradation_0 = Mass_Action_1_222(Vdpds, PDS1);

    // Degradation of PDS1 in PE
    double Degradation_of_PDS1_in_PE = Mass_Action_1_222(Vdpds, PE);

    // Assoc. with ESP1 to form PE
    double Assoc_with_ESP1_to_form_PE = Mass_Action_2_221(kasesp, PDS1, ESP1);

    // Disso. from PE
    double Disso_from_PE = Mass_Action_1_222(kdiesp, PE);

    // DNA synthesis
    double DNA_synthesis = ksori * (eorib5 * CLB5 + eorib2 * CLB2);

    // Negative regulation of DNA synthesis
    double Negative_regulation_of_DNA_synthesis = Mass_Action_1_222(kdori, ORI);

    // Budding
    double Budding = ksbud * (ebudn2 * CLN2 + ebudn3 * CLN3 + ebudb5 * CLB5);

    // Negative regulation of Cell budding
    double Negative_regulation_of_Cell_budding = Mass_Action_1_222(kdbud, BUD);

    // Spindle formation
    double Spindle_formation = ksspn * CLB2 / (Jspn + CLB2);

    // Spindle disassembly
    double Spindle_disassembly = Mass_Action_1_222(kdspn, SPN);

    // ODES:
    rDY[0] = ((b0 * MASS) - rY[0]) / cell; // d[BCK2]/dt
    rDY[1] = (Budding - Negative_regulation_of_Cell_budding) / cell; // d[BUD]/dt
    rDY[2] = (Assoc_of_CLB2_and_SIC1 - Dissoc_of_CLB2SIC1_complex - Phosphorylation_of_C2 + Dephosphorylation_of_C2P - Degradation_of_CLB2_in_C2) / cell; // d[C2]/dt
    rDY[3] = (Phosphorylation_of_C2 - Dephosphorylation_of_C2P - Degradation_of_SIC1_in_C2P - Degradation_of_CLB2_in_C2P) / cell; // d[C2P]/dt
    rDY[4] = (Assoc_of_CLB5_and_SIC1 - Dissoc_of_CLB5SIC1 - Phosphorylation_of_C5 + Dephosphorylation_of_C5P - Degradation_of_CLB5_in_C5) / cell; // d[C5]/dt
    rDY[5] = (Phosphorylation_of_C5 - Dephosphorylation_of_C5P - Degradation_of_SIC1P_in_C5P - Degradation_of_CLB5_in_C5P) / cell; // d[C5P]/dt
    rDY[6] = (CDC14_synthesis - CDC14_degradation - Assoc_with_NET1_to_form_RENT + Dissoc_from_RENT - Assoc_with_NET1P_to_form_RENTP + Dissoc_from_RENP + Degradation_of_NET1_in_RENT + Degradation_of_NET1P_in_RENTP) / cell; // d[CDC14]/dt
    rDY[7] = ((CDC14 + RENT + RENTP) - rY[7]) / cell; // d[CDC14T]/dt
    rDY[8] = (CDC15_activation - inactivation_0) / cell; // d[CDC15]/dt
    rDY[9] = (-CDC15_activation + inactivation_0) / cell; // d[CDC15i]/dt
    rDY[10] = (-Degradation_of_active_CDC20 + Activation_of_CDC20 - Inactivation_0) / cell; // d[CDC20]/dt
    rDY[11] = (Synthesis_of_inactive_CDC20 - Degradation_of_inactiveCDC20 - Activation_of_CDC20 + Inactivation_0) / cell; // d[CDC20i]/dt
    rDY[12] = (CDC6_synthesis - Phosphorylation_of_CDC6 + Dephosphorylation_of_CDC6 - CLB2CDC6_complex_formation + CLB2CDC6_dissociation - CLB5CDC6_complex_formation + CLB5CDC6_dissociation + CLB2_degradation_in_F2 + CLB5_degradation_in_F5) / cell; // d[CDC6]/dt
    rDY[13] = (Phosphorylation_of_CDC6 - Dephosphorylation_of_CDC6 - Degradation_of_CDC6P + CLB2_degradation_in_F2P + CLB5_degradation_in_F5P) / cell; // d[CDC6P]/dt
    rDY[14] = ((CDC6 + F2 + F5 + CDC6P + F2P + F5P) - rY[14]) / cell; // d[CDC6T]/dt
    rDY[15] = (CDH1_synthesis - CDH1_degradation + CDH1i_activation - Inactivation_1) / cell; // d[CDH1]/dt
    rDY[16] = (-CDH1i_degradation - CDH1i_activation + Inactivation_1) / cell; // d[CDH1i]/dt
    rDY[17] = ((SIC1T + CDC6T) - rY[17]) / cell; // d[CKIT]/dt
    rDY[18] = (Synthesis_of_CLB2 - Degradation_of_CLB2 - Assoc_of_CLB2_and_SIC1 + Dissoc_of_CLB2SIC1_complex + Degradation_of_SIC1_in_C2P - CLB2CDC6_complex_formation + CLB2CDC6_dissociation + CDC6_degradation_in_F2P) / cell; // d[CLB2]/dt
    rDY[19] = ((CLB2 + C2 + C2P + F2 + F2P) - rY[19]) / cell; // d[CLB2T]/dt
    rDY[20] = (Synthesis_of_CLB5 - Degradation_of_CLB5 - Assoc_of_CLB5_and_SIC1 + Dissoc_of_CLB5SIC1 + Degradation_of_SIC1P_in_C5P - CLB5CDC6_complex_formation + CLB5CDC6_dissociation + CDC6_degradation_in_F5P) / cell; // d[CLB5]/dt
    rDY[21] = ((CLB5 + C5 + C5P + F5 + F5P) - rY[21]) / cell; // d[CLB5T]/dt
    rDY[22] = (Synthesis_of_CLN2 - Degradation_of_CLN2) / cell; // d[CLN2]/dt
    rDY[23] = ((C0 * Dn3 * MASS / (Jn3 + Dn3 * MASS)) - rY[23]) / cell; // d[CLN3]/dt
    rDY[24] = (Degradation_of_PDS1_in_PE - Assoc_with_ESP1_to_form_PE + Disso_from_PE) / cell; // d[ESP1]/dt
    rDY[25] = (CLB2CDC6_complex_formation - CLB2CDC6_dissociation - F2_phosphorylation + F2P_dephosphorylation - CLB2_degradation_in_F2) / cell; // d[F2]/dt
    rDY[26] = (F2_phosphorylation - F2P_dephosphorylation - CDC6_degradation_in_F2P - CLB2_degradation_in_F2P) / cell; // d[F2P]/dt
    rDY[27] = (CLB5CDC6_complex_formation - CLB5CDC6_dissociation - F5_phosphorylation + F5P_dephosphorylation - CLB5_degradation_in_F5) / cell; // d[F5]/dt
    rDY[28] = (F5_phosphorylation - F5P_dephosphorylation - CDC6_degradation_in_F5P - CLB5_degradation_in_F5P) / cell; // d[F5P]/dt
    rDY[29] = (-Activation_of_IEP + Inactivation) / cell; // d[IE]/dt
    rDY[30] = (Activation_of_IEP - Inactivation) / cell; // d[IEP]/dt
    rDY[31] = (Growth) / cell; // d[MASS]/dt
    rDY[32] = ((GK_219(kamcm * CLB2, kimcm, Jamcm, Jimcm)) - rY[32]) / cell; // d[MCM1]/dt
    rDY[33] = (-Assoc_with_NET1_to_form_RENT + Dissoc_from_RENT + Net1_synthesis - Net1_degradation - NET1_phosphorylation + dephosphorylation + Degradation_of_CDC14_in_RENT) / cell; // d[NET1]/dt
    rDY[34] = (-Assoc_with_NET1P_to_form_RENTP + Dissoc_from_RENP - Net1P_degradation + NET1_phosphorylation - dephosphorylation + Degradation_of_CDC14_in_RENTP) / cell; // d[NET1P]/dt
    rDY[35] = ((NET1 + NET1P + RENT + RENTP) - rY[35]) / cell; // d[NET1T]/dt
    rDY[36] = (DNA_synthesis - Negative_regulation_of_DNA_synthesis) / cell; // d[ORI]/dt
    rDY[37] = (PDS1_synthesis - degradation_0 - Assoc_with_ESP1_to_form_PE + Disso_from_PE) / cell; // d[PDS1]/dt
    rDY[38] = (-Degradation_of_PDS1_in_PE + Assoc_with_ESP1_to_form_PE - Disso_from_PE) / cell; // d[PE]/dt
    rDY[39] = (PPX_synthesis - degradation) / cell; // d[PPX]/dt
    rDY[40] = (Assoc_with_NET1_to_form_RENT - Dissoc_from_RENT - RENT_phosphorylation + dephosphorylation_0 - Degradation_of_NET1_in_RENT - Degradation_of_CDC14_in_RENT) / cell; // d[RENT]/dt
    rDY[41] = (Assoc_with_NET1P_to_form_RENTP - Dissoc_from_RENP + RENT_phosphorylation - dephosphorylation_0 - Degradation_of_NET1P_in_RENTP - Degradation_of_CDC14_in_RENTP) / cell; // d[RENTP]/dt
    rDY[42] = ((GK_219(Vasbf, Visbf, Jasbf, Jisbf)) - rY[42]) / cell; // d[SBF]/dt
    rDY[43] = (Synthesis_of_SIC1 - Phosphorylation_of_SIC1 + Dephosphorylation_of_SIC1 - Assoc_of_CLB2_and_SIC1 + Dissoc_of_CLB2SIC1_complex - Assoc_of_CLB5_and_SIC1 + Dissoc_of_CLB5SIC1 + Degradation_of_CLB2_in_C2 + Degradation_of_CLB5_in_C5) / cell; // d[SIC1]/dt
    rDY[44] = (Phosphorylation_of_SIC1 - Dephosphorylation_of_SIC1 - Fast_Degradation_of_SIC1P + Degradation_of_CLB2_in_C2P + Degradation_of_CLB5_in_C5P) / cell; // d[SIC1P]/dt
    rDY[45] = ((SIC1 + C2 + C5 + SIC1P + C2P + C5P) - rY[45]) / cell; // d[SIC1T]/dt
    rDY[46] = (Spindle_formation - Spindle_disassembly) / cell; // d[SPN]/dt
    rDY[47] = (Synthesis_of_SWI5 - Degradation_of_SWI5 + Activation_of_SWI5 - Inactivation_of_SWI5) / cell; // d[SWI5]/dt
    rDY[48] = (-Degradation_of_SWI5P - Activation_of_SWI5 + Inactivation_of_SWI5) / cell; // d[SWI5P]/dt
    rDY[49] = (-TEM1_activation + inactivation) / cell; // d[TEM1GDP]/dt
    rDY[50] = (TEM1_activation - inactivation) / cell; // d[TEM1GTP]/dt

    // Scale time appropriately
}

double Chen2004SbmlOdeSystem::CalculateRootFunction(double time, const std::vector<double> &rY)
{
    RefreshState(rY);

    double dist = std::numeric_limits<double>::max();

    if (sm::lt(CLB2 + CLB5 - KEZ2, 0.0))
    {
        if (!eventsInitialised)
        {
            // Condition true at first timestep: don't trigger
            dist = std::abs(dist) < 1.0 ? dist : 1.0;
        }
        else if (eventsSatisfied[0])
        {
            // Condition already true: don't trigger (again)
            dist = std::abs(dist) < 1.0 ? dist : 1.0;
        }
        else // if (!eventsSatisfied[0])
        {
            // Condition transitioning from false to true: trigger
            dist = 0.0;
            
            UpdateDefaultInitialConditions(rY);
            SetStateVariable(36, 0.0);
            SetDefaultInitialCondition(36, 0.0);
        }
        // Mark condition true
        eventsSatisfied[0] = true;
    }
    else
    {
        double event_dist = CLB2 + CLB5 - KEZ2 - 0.0 + std::numeric_limits<double>::epsilon();
        dist = std::abs(dist) < std::abs(event_dist) ? dist : event_dist;

        // Mark condition false
        eventsSatisfied[0] = false;
    }
    if (sm::gt(ORI - 1.0, 0.0))
    {
        if (!eventsInitialised)
        {
            // Condition true at first timestep: don't trigger
            dist = std::abs(dist) < 1.0 ? dist : 1.0;
        }
        else if (eventsSatisfied[1])
        {
            // Condition already true: don't trigger (again)
            dist = std::abs(dist) < 1.0 ? dist : 1.0;
        }
        else // if (!eventsSatisfied[1])
        {
            // Condition transitioning from false to true: trigger
            dist = 0.0;
            
            UpdateDefaultInitialConditions(rY);
            SetParameter("MAD2", mad2h);
            SetParameter("BUB2", bub2h);
        }
        // Mark condition true
        eventsSatisfied[1] = true;
    }
    else
    {
        double event_dist = 0.0 - ORI - 1.0 + std::numeric_limits<double>::epsilon();
        dist = std::abs(dist) < std::abs(event_dist) ? dist : event_dist;

        // Mark condition false
        eventsSatisfied[1] = false;
    }
    if (sm::gt(SPN - 1.0, 0.0))
    {
        if (!eventsInitialised)
        {
            // Condition true at first timestep: don't trigger
            dist = std::abs(dist) < 1.0 ? dist : 1.0;
        }
        else if (eventsSatisfied[2])
        {
            // Condition already true: don't trigger (again)
            dist = std::abs(dist) < 1.0 ? dist : 1.0;
        }
        else // if (!eventsSatisfied[2])
        {
            // Condition transitioning from false to true: trigger
            dist = 0.0;
            
            UpdateDefaultInitialConditions(rY);
            SetParameter("MAD2", mad2l);
            SetParameter("LTE1", lte1h);
            SetParameter("BUB2", bub2l);
        }
        // Mark condition true
        eventsSatisfied[2] = true;
    }
    else
    {
        double event_dist = 0.0 - SPN - 1.0 + std::numeric_limits<double>::epsilon();
        dist = std::abs(dist) < std::abs(event_dist) ? dist : event_dist;

        // Mark condition false
        eventsSatisfied[2] = false;
    }
    if (sm::lt(CLB2 - KEZ, 0.0))
    {
        if (!eventsInitialised)
        {
            // Condition true at first timestep: don't trigger
            dist = std::abs(dist) < 1.0 ? dist : 1.0;
        }
        else if (eventsSatisfied[3])
        {
            // Condition already true: don't trigger (again)
            dist = std::abs(dist) < 1.0 ? dist : 1.0;
        }
        else // if (!eventsSatisfied[3])
        {
            // Condition transitioning from false to true: trigger
            dist = 0.0;
            
            UpdateDefaultInitialConditions(rY);
            SetStateVariable(31, F * MASS);
            SetDefaultInitialCondition(31, F * MASS);
            SetParameter("LTE1", lte1l);
            SetStateVariable(1, 0.0);
            SetDefaultInitialCondition(1, 0.0);
            SetStateVariable(46, 0.0);
            SetDefaultInitialCondition(46, 0.0);
        }
        // Mark condition true
        eventsSatisfied[3] = true;
    }
    else
    {
        double event_dist = CLB2 - KEZ - 0.0 + std::numeric_limits<double>::epsilon();
        dist = std::abs(dist) < std::abs(event_dist) ? dist : event_dist;

        // Mark condition false
        eventsSatisfied[3] = false;
    }

    eventsInitialised = true;

    return dist;
}

bool Chen2004SbmlOdeSystem::CalculateStoppingEvent(double time, const std::vector<double> &rY)
{
    return CalculateRootFunction(time, rY) == 0.0;
}

void Chen2004SbmlOdeSystem::UpdateDefaultInitialConditions(const std::vector<double> &rY)
{
    SetDefaultInitialCondition(0, rY[0]); // BCK2
    SetDefaultInitialCondition(1, rY[1]); // BUD
    SetDefaultInitialCondition(2, rY[2]); // C2
    SetDefaultInitialCondition(3, rY[3]); // C2P
    SetDefaultInitialCondition(4, rY[4]); // C5
    SetDefaultInitialCondition(5, rY[5]); // C5P
    SetDefaultInitialCondition(6, rY[6]); // CDC14
    SetDefaultInitialCondition(7, rY[7]); // CDC14T
    SetDefaultInitialCondition(8, rY[8]); // CDC15
    SetDefaultInitialCondition(9, rY[9]); // CDC15i
    SetDefaultInitialCondition(10, rY[10]); // CDC20
    SetDefaultInitialCondition(11, rY[11]); // CDC20i
    SetDefaultInitialCondition(12, rY[12]); // CDC6
    SetDefaultInitialCondition(13, rY[13]); // CDC6P
    SetDefaultInitialCondition(14, rY[14]); // CDC6T
    SetDefaultInitialCondition(15, rY[15]); // CDH1
    SetDefaultInitialCondition(16, rY[16]); // CDH1i
    SetDefaultInitialCondition(17, rY[17]); // CKIT
    SetDefaultInitialCondition(18, rY[18]); // CLB2
    SetDefaultInitialCondition(19, rY[19]); // CLB2T
    SetDefaultInitialCondition(20, rY[20]); // CLB5
    SetDefaultInitialCondition(21, rY[21]); // CLB5T
    SetDefaultInitialCondition(22, rY[22]); // CLN2
    SetDefaultInitialCondition(23, rY[23]); // CLN3
    SetDefaultInitialCondition(24, rY[24]); // ESP1
    SetDefaultInitialCondition(25, rY[25]); // F2
    SetDefaultInitialCondition(26, rY[26]); // F2P
    SetDefaultInitialCondition(27, rY[27]); // F5
    SetDefaultInitialCondition(28, rY[28]); // F5P
    SetDefaultInitialCondition(29, rY[29]); // IE
    SetDefaultInitialCondition(30, rY[30]); // IEP
    SetDefaultInitialCondition(31, rY[31]); // MASS
    SetDefaultInitialCondition(32, rY[32]); // MCM1
    SetDefaultInitialCondition(33, rY[33]); // NET1
    SetDefaultInitialCondition(34, rY[34]); // NET1P
    SetDefaultInitialCondition(35, rY[35]); // NET1T
    SetDefaultInitialCondition(36, rY[36]); // ORI
    SetDefaultInitialCondition(37, rY[37]); // PDS1
    SetDefaultInitialCondition(38, rY[38]); // PE
    SetDefaultInitialCondition(39, rY[39]); // PPX
    SetDefaultInitialCondition(40, rY[40]); // RENT
    SetDefaultInitialCondition(41, rY[41]); // RENTP
    SetDefaultInitialCondition(42, rY[42]); // SBF
    SetDefaultInitialCondition(43, rY[43]); // SIC1
    SetDefaultInitialCondition(44, rY[44]); // SIC1P
    SetDefaultInitialCondition(45, rY[45]); // SIC1T
    SetDefaultInitialCondition(46, rY[46]); // SPN
    SetDefaultInitialCondition(47, rY[47]); // SWI5
    SetDefaultInitialCondition(48, rY[48]); // SWI5P
    SetDefaultInitialCondition(49, rY[49]); // TEM1GDP
    SetDefaultInitialCondition(50, rY[50]); // TEM1GTP
}

// FUNCTION DEFINITIONS:
double Chen2004SbmlOdeSystem::BB_218(double A1, double A2, double A3, double A4)
{
    return A2 - A1 + A3 * A2 + A4 * A1;
}
double Chen2004SbmlOdeSystem::GK_219(double A1, double A2, double A3, double A4)
{
    return 2.0 * A4 * A1 / (A2 - A1 + A3 * A2 + A4 * A1 + sm::root(2.0, std::pow(A2 - A1 + A3 * A2 + A4 * A1, 2.0) - 4.0 * (A2 - A1) * A4 * A1));
}
double Chen2004SbmlOdeSystem::MichaelisMenten_220(double M1, double J1, double k1, double S1)
{
    return k1 * S1 * M1 / (J1 + S1);
}
double Chen2004SbmlOdeSystem::Mass_Action_2_221(double k1, double S1, double S2)
{
    return k1 * S1 * S2;
}
double Chen2004SbmlOdeSystem::Mass_Action_1_222(double k1, double S1)
{
    return k1 * S1;
}

template <>
void CellwiseOdeSystemInformation<Chen2004SbmlOdeSystem>::Initialise()
{
    // STATE VARIABLES:
    this->mVariableNames.push_back("BCK2");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

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

    this->mVariableNames.push_back("CDC14T");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(2.0);

    this->mVariableNames.push_back("CDC15");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.656533);

    this->mVariableNames.push_back("CDC15i");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

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

    this->mVariableNames.push_back("CDC6T");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("CDH1");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.930499);

    this->mVariableNames.push_back("CDH1i");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0695);

    this->mVariableNames.push_back("CKIT");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("CLB2");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.1469227);

    this->mVariableNames.push_back("CLB2T");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.17);

    this->mVariableNames.push_back("CLB5");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0518014);

    this->mVariableNames.push_back("CLB5T");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.12);

    this->mVariableNames.push_back("CLN2");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0652511);

    this->mVariableNames.push_back("CLN3");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

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

    this->mVariableNames.push_back("IE");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("IEP");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.1015);

    this->mVariableNames.push_back("MASS");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(1.206019);

    this->mVariableNames.push_back("MCM1");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("NET1");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.018645);

    this->mVariableNames.push_back("NET1P");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.970271);

    this->mVariableNames.push_back("NET1T");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(2.8);

    this->mVariableNames.push_back("ORI");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.000909);

    this->mVariableNames.push_back("PDS1");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.025612);

    this->mVariableNames.push_back("PE");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("PPX");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.123179);

    this->mVariableNames.push_back("RENT");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(1.04954);

    this->mVariableNames.push_back("RENTP");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.6);

    this->mVariableNames.push_back("SBF");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("SIC1");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0228776);

    this->mVariableNames.push_back("SIC1P");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.00641);

    this->mVariableNames.push_back("SIC1T");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("SPN");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.03);

    this->mVariableNames.push_back("SWI5");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.95);

    this->mVariableNames.push_back("SWI5P");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.02);

    this->mVariableNames.push_back("TEM1GDP");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("TEM1GTP");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.9);


    // STATE PARAMETERS:
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

    this->mParameterNames.push_back("D");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("mu");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Vdb5");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Vdb2");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Vasbf");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Visbf");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Vkpc1");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Vkpf6");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Vacdh");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Vicdh");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Vppnet");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Vkpnet");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Vdppx");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Vdpds");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Vaiep");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Vd2c1");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Vd2f6");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Vppc1");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Vppf6");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("F");
    this->mParameterUnits.push_back("non-dim");

    this->mInitialised = true;
}

// Define SbmlCellCycleWrapperModel using wrappers
#include "SbmlCellCycleWrapperModel.hpp"
#include "SbmlCellCycleWrapperModel.cpp"

typedef SbmlCellCycleWrapperModel<Chen2004SbmlOdeSystem, 51> Chen2004SbmlCellCycleModel;

// Declare identifiers for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(Chen2004SbmlOdeSystem)
EXPORT_TEMPLATE_CLASS2(SbmlCellCycleWrapperModel, Chen2004SbmlOdeSystem, 51)

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(Chen2004SbmlCellCycleModel)