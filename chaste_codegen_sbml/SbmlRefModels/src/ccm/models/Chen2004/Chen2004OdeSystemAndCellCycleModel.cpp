#include "CellwiseOdeSystemInformation.hpp"
#include "SbmlMath.hpp"

#include "Chen2004OdeSystemAndCellCycleModel.hpp"

/* SBML ODE System */
Chen2004OdeSystem::Chen2004OdeSystem(std::vector<double> stateVariables)
    : AbstractOdeSystem(51)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<Chen2004OdeSystem>);

    Init();

    SetDefaultInitialCondition(0, 0.0); // BCK2
    SetDefaultInitialCondition(1, 0.008473); // BUD
    SetDefaultInitialCondition(2, 0.238404); // C2
    SetDefaultInitialCondition(3, 0.024034); // C2P
    SetDefaultInitialCondition(4, 0.070081); // C5
    SetDefaultInitialCondition(5, 0.006878); // C5P
    SetDefaultInitialCondition(6, 0.468344); // CDC14
    SetDefaultInitialCondition(7, 2.0); // CDC14T
    SetDefaultInitialCondition(8, 0.656533); // CDC15
    SetDefaultInitialCondition(9, 0.0); // CDC15i
    SetDefaultInitialCondition(10, 0.444296); // CDC20
    SetDefaultInitialCondition(11, 1.472044); // CDC20i
    SetDefaultInitialCondition(12, 0.10758); // CDC6
    SetDefaultInitialCondition(13, 0.015486); // CDC6P
    SetDefaultInitialCondition(14, 0.0); // CDC6T
    SetDefaultInitialCondition(15, 0.930499); // CDH1
    SetDefaultInitialCondition(16, 0.0695); // CDH1i
    SetDefaultInitialCondition(17, 0.0); // CKIT
    SetDefaultInitialCondition(18, 0.1469227); // CLB2
    SetDefaultInitialCondition(19, 0.17); // CLB2T
    SetDefaultInitialCondition(20, 0.0518014); // CLB5
    SetDefaultInitialCondition(21, 0.12); // CLB5T
    SetDefaultInitialCondition(22, 0.0652511); // CLN2
    SetDefaultInitialCondition(23, 0.0); // CLN3
    SetDefaultInitialCondition(24, 0.301313); // ESP1
    SetDefaultInitialCondition(25, 0.236058); // F2
    SetDefaultInitialCondition(26, 0.0273938); // F2P
    SetDefaultInitialCondition(27, 7.24e-05); // F5
    SetDefaultInitialCondition(28, 7.91e-05); // F5P
    SetDefaultInitialCondition(29, 0.0); // IE
    SetDefaultInitialCondition(30, 0.1015); // IEP
    SetDefaultInitialCondition(31, 1.206019); // MASS
    SetDefaultInitialCondition(32, 0.0); // MCM1
    SetDefaultInitialCondition(33, 0.018645); // NET1
    SetDefaultInitialCondition(34, 0.970271); // NET1P
    SetDefaultInitialCondition(35, 2.8); // NET1T
    SetDefaultInitialCondition(36, 0.000909); // ORI
    SetDefaultInitialCondition(37, 0.025612); // PDS1
    SetDefaultInitialCondition(38, 0.0); // PE
    SetDefaultInitialCondition(39, 0.123179); // PPX
    SetDefaultInitialCondition(40, 1.04954); // RENT
    SetDefaultInitialCondition(41, 0.6); // RENTP
    SetDefaultInitialCondition(42, 0.0); // SBF
    SetDefaultInitialCondition(43, 0.0228776); // SIC1
    SetDefaultInitialCondition(44, 0.00641); // SIC1P
    SetDefaultInitialCondition(45, 0.0); // SIC1T
    SetDefaultInitialCondition(46, 0.03); // SPN
    SetDefaultInitialCondition(47, 0.95); // SWI5
    SetDefaultInitialCondition(48, 0.02); // SWI5P
    SetDefaultInitialCondition(49, 0.0); // TEM1GDP
    SetDefaultInitialCondition(50, 0.9); // TEM1GTP


    if (stateVariables != std::vector<double>())
    {
        SetStateVariables(stateVariables);
    }
}

Chen2004OdeSystem::~Chen2004OdeSystem()
{
}

double Chen2004OdeSystem::BB_218(double A1, double A2, double A3, double A4)
{
    return A2 - A1 + A3 * A2 + A4 * A1;
}
double Chen2004OdeSystem::GK_219(double A1, double A2, double A3, double A4)
{
    return 2 * A4 * A1 / (A2 - A1 + A3 * A2 + A4 * A1 + sbmlmath::sm_root(2, std::pow(A2 - A1 + A3 * A2 + A4 * A1, 2) - 4 * (A2 - A1) * A4 * A1));
}
double Chen2004OdeSystem::MichaelisMenten_220(double M1, double J1, double k1, double S1)
{
    return k1 * S1 * M1 / (J1 + S1);
}
double Chen2004OdeSystem::Mass_Action_2_221(double k1, double S1, double S2)
{
    return k1 * S1 * S2;
}
double Chen2004OdeSystem::Mass_Action_1_222(double k1, double S1)
{
    return k1 * S1;
}

void Chen2004OdeSystem::Init()
{
    /* Initialise model compartments. */
    cell = 1.0;

    /* Initialise model parameters. */
    b0 = 0.054;
    bub2h = 1.0;
    bub2l = 0.2;
    C0 = 0.4;
    CDC15T = 1.0;
    Dn3 = 1.0;
    ebudb5 = 1.0;
    ebudn2 = 0.25;
    ebudn3 = 0.05;
    ec1b2 = 0.45;
    ec1b5 = 0.1;
    ec1k2 = 0.03;
    ec1n2 = 0.06;
    ec1n3 = 0.3;
    ef6b2 = 0.55;
    ef6b5 = 0.1;
    ef6k2 = 0.03;
    ef6n2 = 0.06;
    ef6n3 = 0.3;
    eicdhb2 = 1.2;
    eicdhb5 = 8.0;
    eicdhn2 = 0.4;
    eicdhn3 = 0.25;
    eorib2 = 0.45;
    eorib5 = 0.9;
    esbfb5 = 2.0;
    esbfn2 = 2.0;
    esbfn3 = 10.0;
    ESP1T = 1.0;
    IET = 1.0;
    J20ppx = 0.15;
    Jacdh = 0.03;
    Jaiep = 0.1;
    Jamcm = 0.1;
    Jasbf = 0.01;
    Jatem = 0.1;
    Jd2c1 = 0.05;
    Jd2f6 = 0.05;
    Jicdh = 0.03;
    Jiiep = 0.1;
    Jimcm = 0.1;
    Jisbf = 0.01;
    Jitem = 0.1;
    Jn3 = 6.0;
    Jpds = 0.04;
    Jspn = 0.14;
    ka15_p = 0.002;
    ka15_p_p = 1.0;
    ka15p = 0.001;
    ka20_p = 0.05;
    ka20_p_p = 0.2;
    kacdh_p = 0.01;
    kacdh_p_p = 0.8;
    kaiep = 0.1;
    kamcm = 1.0;
    kasb2 = 50.0;
    kasb5 = 50.0;
    kasbf = 0.38;
    kasesp = 50.0;
    kasf2 = 15.0;
    kasf5 = 0.01;
    kasrent = 200.0;
    kasrentp = 1.0;
    kaswi = 2.0;
    kd14 = 0.1;
    kd1c1 = 0.01;
    kd1f6 = 0.01;
    kd1pds_p = 0.01;
    kd20 = 0.3;
    kd2c1 = 1.0;
    kd2f6 = 1.0;
    kd2pds_p_p = 0.2;
    kd3c1 = 1.0;
    kd3f6 = 1.0;
    kd3pds_p_p = 0.04;
    kdb2_p = 0.003;
    kdb2_p_p = 0.4;
    kdb2p = 0.15;
    kdb5_p = 0.01;
    kdb5_p_p = 0.16;
    kdbud = 0.06;
    kdcdh = 0.01;
    kdib2 = 0.05;
    kdib5 = 0.06;
    kdiesp = 0.5;
    kdif2 = 0.5;
    kdif5 = 0.01;
    kdirent = 1.0;
    kdirentp = 2.0;
    kdn2 = 0.12;
    kdnet = 0.03;
    kdori = 0.06;
    kdppx_p = 0.17;
    kdppx_p_p = 2.0;
    kdspn = 0.06;
    kdswi = 0.08;
    KEZ = 0.3;
    KEZ2 = 0.2;
    ki15 = 0.5;
    kicdh_p = 0.001;
    kicdh_p_p = 0.08;
    kiiep = 0.15;
    kimcm = 0.15;
    kisbf_p = 0.6;
    kisbf_p_p = 8.0;
    kiswi = 0.05;
    kkpnet_p = 0.01;
    kkpnet_p_p = 0.6;
    kppc1 = 4.0;
    kppf6 = 4.0;
    kppnet_p = 0.05;
    kppnet_p_p = 3.0;
    ks14 = 0.2;
    ks1pds_p_p = 0.03;
    ks20_p = 0.006;
    ks20_p_p = 0.6;
    ks2pds_p_p = 0.055;
    ksb2_p = 0.001;
    ksb2_p_p = 0.04;
    ksb5_p = 0.0008;
    ksb5_p_p = 0.005;
    ksbud = 0.2;
    ksc1_p = 0.012;
    ksc1_p_p = 0.12;
    kscdh = 0.01;
    ksf6_p = 0.024;
    ksf6_p_p = 0.12;
    ksf6_p_p_p = 0.004;
    ksn2_p = 0.0;
    ksn2_p_p = 0.15;
    ksnet = 0.084;
    ksori = 2.0;
    kspds_p = 0.0;
    ksppx = 0.1;
    ksspn = 0.1;
    ksswi_p = 0.005;
    ksswi_p_p = 0.08;
    lte1h = 1.0;
    lte1l = 0.1;
    mad2h = 8.0;
    mad2l = 0.01;
    mdt = 90.0;
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

    /* Initialise vector to check if events have been triggered. */
    eventsSatisfied.resize(4, false);
}

void Chen2004OdeSystem::EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY)
{
    /* Define algebraic rules. */
    double BCK2 = b0 * rY[31];
    double Visbf = kisbf_p + kisbf_p_p * rY[18];
    double CLN3 = C0 * Dn3 * rY[31] / (Jn3 + Dn3 * rY[31]);
    double Vppc1 = kppc1 * rY[6];
    double Vppf6 = kppf6 * rY[6];
    double Vaiep = kaiep * rY[18];
    double Vacdh = kacdh_p + kacdh_p_p * rY[6];
    double Vicdh = kicdh_p + kicdh_p_p * (eicdhn3 * rY[23] + eicdhn2 * rY[22] + eicdhb5 * rY[20] + eicdhb2 * rY[18]);
    double Vkpnet = (kkpnet_p + kkpnet_p_p * rY[8]) * rY[31];
    double Vppnet = kppnet_p + kppnet_p_p * rY[39];
    double Vasbf = kasbf * (esbfn2 * rY[22] + esbfn3 * (rY[23] + rY[0]) + esbfb5 * rY[20]);
    double SBF = GK_219(Vasbf, Visbf, Jasbf, Jisbf);
    double MCM1 = GK_219(kamcm * rY[18], kimcm, Jamcm, Jimcm);
    double mu = sbmlmath::sm_log(2) / mdt;
    double D = 1.026 / mu - 32;
    double F = std::exp(-mu * D);
    double Vd2c1 = kd2c1 * (ec1n3 * rY[23] + ec1k2 * rY[0] + ec1n2 * rY[22] + ec1b5 * rY[20] + ec1b2 * rY[18]);
    double Vd2f6 = kd2f6 * (ef6n3 * rY[23] + ef6k2 * rY[0] + ef6n2 * rY[22] + ef6b5 * rY[20] + ef6b2 * rY[18]);
    double Vkpc1 = kd1c1 + Vd2c1 / (Jd2c1 + rY[43] + rY[2] + rY[4] + rY[44] + rY[3] + rY[5]);
    double Vkpf6 = kd1f6 + Vd2f6 / (Jd2f6 + rY[12] + rY[25] + rY[27] + rY[13] + rY[26] + rY[28]);
    double Vdb2 = kdb2_p + kdb2_p_p * rY[15] + kdb2p * rY[10];
    double Vdb5 = kdb5_p + kdb5_p_p * rY[10];
    double Vdpds = kd1pds_p + kd2pds_p_p * rY[10] + kd3pds_p_p * rY[15];
    double Vdppx = kdppx_p + kdppx_p_p * (J20ppx + rY[10]) * Jpds / (Jpds + rY[37]);
    double CLB2T = rY[18] + rY[2] + rY[3] + rY[25] + rY[26];
    double CLB5T = rY[20] + rY[4] + rY[5] + rY[27] + rY[28];
    double CDC14T = rY[6] + rY[40] + rY[41];
    double NET1T = rY[33] + rY[34] + rY[40] + rY[41];
    double SIC1T = rY[43] + rY[2] + rY[4] + rY[44] + rY[3] + rY[5];
    double CDC6T = rY[12] + rY[25] + rY[27] + rY[13] + rY[26] + rY[28];
    double CKIT = rY[45] + rY[14];
    double CDC15i = CDC15T - rY[8];
    double IE = IET - rY[30];
    double PE = ESP1T - rY[24];
    double TEM1GDP = TEM1T - rY[50];

    /* Define the reactions in this model. */
    // Growth
    double Growth = mu * rY[31];

    // Synthesis of CLN2
    double Synthesis_of_CLN2 = (ksn2_p + ksn2_p_p * rY[42]) * rY[31];

    // Degradation of CLN2
    double Degradation_of_CLN2 = Mass_Action_1_222(kdn2, rY[22]);

    // Synthesis of CLB2
    double Synthesis_of_CLB2 = (ksb2_p + ksb2_p_p * rY[32]) * rY[31];

    // Degradation of CLB2
    double Degradation_of_CLB2 = Mass_Action_1_222(Vdb2, rY[18]);

    // Synthesis of CLB5
    double Synthesis_of_CLB5 = (ksb5_p + ksb5_p_p * rY[42]) * rY[31];

    // Degradation of CLB5
    double Degradation_of_CLB5 = Mass_Action_1_222(Vdb5, rY[20]);

    // Synthesis of SIC1
    double Synthesis_of_SIC1 = ksc1_p + ksc1_p_p * rY[47];

    // Phosphorylation of SIC1
    double Phosphorylation_of_SIC1 = Mass_Action_1_222(Vkpc1, rY[43]);

    // Dephosphorylation of SIC1
    double Dephosphorylation_of_SIC1 = Mass_Action_1_222(Vppc1, rY[44]);

    // Fast Degradation of SIC1P
    double Fast_Degradation_of_SIC1P = Mass_Action_1_222(kd3c1, rY[44]);

    // Assoc. of CLB2 and SIC1
    double Assoc_of_CLB2_and_SIC1 = Mass_Action_2_221(kasb2, rY[18], rY[43]);

    // Dissoc. of CLB2/SIC1 complex
    double Dissoc_of_CLB2SIC1_complex = Mass_Action_1_222(kdib2, rY[2]);

    // Assoc. of CLB5 and SIC1
    double Assoc_of_CLB5_and_SIC1 = Mass_Action_2_221(kasb5, rY[20], rY[43]);

    // Dissoc. of CLB5/SIC1
    double Dissoc_of_CLB5SIC1 = Mass_Action_1_222(kdib5, rY[4]);

    // Phosphorylation of C2
    double Phosphorylation_of_C2 = Mass_Action_1_222(Vkpc1, rY[2]);

    // Dephosphorylation of C2P
    double Dephosphorylation_of_C2P = Mass_Action_1_222(Vppc1, rY[3]);

    // Phosphorylation of C5
    double Phosphorylation_of_C5 = Mass_Action_1_222(Vkpc1, rY[4]);

    // Dephosphorylation of C5P
    double Dephosphorylation_of_C5P = Mass_Action_1_222(Vppc1, rY[5]);

    // Degradation of CLB2 in C2
    double Degradation_of_CLB2_in_C2 = Mass_Action_1_222(Vdb2, rY[2]);

    // Degradation of CLB5 in C5
    double Degradation_of_CLB5_in_C5 = Mass_Action_1_222(Vdb5, rY[4]);

    // Degradation of SIC1 in C2P
    double Degradation_of_SIC1_in_C2P = Mass_Action_1_222(kd3c1, rY[3]);

    // Degradation of SIC1P in C5P
    double Degradation_of_SIC1P_in_C5P = Mass_Action_1_222(kd3c1, rY[5]);

    // Degradation of CLB2 in C2P
    double Degradation_of_CLB2_in_C2P = Mass_Action_1_222(Vdb2, rY[3]);

    // Degradation of CLB5 in C5P
    double Degradation_of_CLB5_in_C5P = Mass_Action_1_222(Vdb5, rY[5]);

    // CDC6 synthesis
    double CDC6_synthesis = ksf6_p + ksf6_p_p * rY[47] + ksf6_p_p_p * rY[42];

    // Phosphorylation of CDC6
    double Phosphorylation_of_CDC6 = Mass_Action_1_222(Vkpf6, rY[12]);

    // Dephosphorylation of CDC6
    double Dephosphorylation_of_CDC6 = Mass_Action_1_222(Vppf6, rY[13]);

    // Degradation of CDC6P
    double Degradation_of_CDC6P = Mass_Action_1_222(kd3f6, rY[13]);

    // CLB2/CDC6 complex formation
    double CLB2CDC6_complex_formation = Mass_Action_2_221(kasf2, rY[18], rY[12]);

    // CLB2/CDC6 dissociation
    double CLB2CDC6_dissociation = Mass_Action_1_222(kdif2, rY[25]);

    // CLB5/CDC6 complex formation
    double CLB5CDC6_complex_formation = Mass_Action_2_221(kasf5, rY[20], rY[12]);

    // CLB5/CDC6 dissociation
    double CLB5CDC6_dissociation = Mass_Action_1_222(kdif5, rY[27]);

    // F2 phosphorylation
    double F2_phosphorylation = Mass_Action_1_222(Vkpf6, rY[25]);

    // F2P dephosphorylation
    double F2P_dephosphorylation = Mass_Action_1_222(Vppf6, rY[26]);

    // F5 phosphorylation
    double F5_phosphorylation = Mass_Action_1_222(Vkpf6, rY[27]);

    // F5P dephosphorylation
    double F5P_dephosphorylation = Mass_Action_1_222(Vppf6, rY[28]);

    // CLB2 degradation in F2
    double CLB2_degradation_in_F2 = Mass_Action_1_222(Vdb2, rY[25]);

    // CLB5 degradation in F5
    double CLB5_degradation_in_F5 = Mass_Action_1_222(Vdb5, rY[27]);

    // CDC6 degradation in F2P
    double CDC6_degradation_in_F2P = Mass_Action_1_222(kd3f6, rY[26]);

    // CDC6 degradation in F5P
    double CDC6_degradation_in_F5P = Mass_Action_1_222(kd3f6, rY[28]);

    // CLB2 degradation in F2P
    double CLB2_degradation_in_F2P = Mass_Action_1_222(Vdb2, rY[26]);

    // CLB5 degradation in F5P
    double CLB5_degradation_in_F5P = Mass_Action_1_222(Vdb5, rY[28]);

    // Synthesis of SWI5
    double Synthesis_of_SWI5 = ksswi_p + ksswi_p_p * rY[32];

    // Degradation of SWI5
    double Degradation_of_SWI5 = Mass_Action_1_222(kdswi, rY[47]);

    // Degradation of SWI5P
    double Degradation_of_SWI5P = Mass_Action_1_222(kdswi, rY[48]);

    // Activation of SWI5
    double Activation_of_SWI5 = Mass_Action_1_222(kaswi * rY[6], rY[48]);

    // Inactivation of SWI5
    double Inactivation_of_SWI5 = Mass_Action_1_222(kiswi * rY[18], rY[47]);

    // Activation of IEP
    double Activation_of_IEP = MichaelisMenten_220(Vaiep, Jaiep, 1, rY[29]);

    // Inactivation
    double Inactivation = MichaelisMenten_220(1, Jiiep, kiiep, rY[30]);

    // Synthesis of inactive CDC20
    double Synthesis_of_inactive_CDC20 = ks20_p + ks20_p_p * rY[32];

    // Degradation of inactiveCDC20
    double Degradation_of_inactiveCDC20 = Mass_Action_1_222(kd20, rY[11]);

    // Degradation of active CDC20
    double Degradation_of_active_CDC20 = Mass_Action_1_222(kd20, rY[10]);

    // Activation of CDC20
    double Activation_of_CDC20 = Mass_Action_1_222(ka20_p + ka20_p_p * rY[30], rY[11]);

    // Inactivation
    double Inactivation_0 = 0.0;
    {
        double k = 1.0;
        Inactivation_0 = k * Mass_Action_1_222(this->mParameters[2], rY[10]);
    }

    // CDH1 synthesis
    double CDH1_synthesis = kscdh;

    // CDH1 degradation
    double CDH1_degradation = Mass_Action_1_222(kdcdh, rY[15]);

    // CDH1i degradation
    double CDH1i_degradation = Mass_Action_1_222(kdcdh, rY[16]);

    // CDH1i activation
    double CDH1i_activation = MichaelisMenten_220(Vacdh, Jacdh, 1, rY[16]);

    // Inactivation
    double Inactivation_1 = MichaelisMenten_220(Vicdh, Jicdh, 1, rY[15]);

    // CDC14 synthesis
    double CDC14_synthesis = ks14;

    // CDC14 degradation
    double CDC14_degradation = Mass_Action_1_222(kd14, rY[6]);

    // Assoc. with NET1 to form RENT
    double Assoc_with_NET1_to_form_RENT = Mass_Action_2_221(kasrent, rY[6], rY[33]);

    // Dissoc. from RENT
    double Dissoc_from_RENT = Mass_Action_1_222(kdirent, rY[40]);

    // Assoc with NET1P to form RENTP
    double Assoc_with_NET1P_to_form_RENTP = Mass_Action_2_221(kasrentp, rY[6], rY[34]);

    // Dissoc. from RENP
    double Dissoc_from_RENP = Mass_Action_1_222(kdirentp, rY[41]);

    // Net1 synthesis
    double Net1_synthesis = ksnet;

    // Net1 degradation
    double Net1_degradation = Mass_Action_1_222(kdnet, rY[33]);

    // Net1P degradation
    double Net1P_degradation = Mass_Action_1_222(kdnet, rY[34]);

    // NET1 phosphorylation
    double NET1_phosphorylation = Mass_Action_1_222(Vkpnet, rY[33]);

    // dephosphorylation
    double dephosphorylation = Mass_Action_1_222(Vppnet, rY[34]);

    // RENT phosphorylation
    double RENT_phosphorylation = Mass_Action_1_222(Vkpnet, rY[40]);

    // dephosphorylation
    double dephosphorylation_0 = Mass_Action_1_222(Vppnet, rY[41]);

    // Degradation of NET1 in RENT
    double Degradation_of_NET1_in_RENT = Mass_Action_1_222(kdnet, rY[40]);

    // Degradation of NET1P in RENTP
    double Degradation_of_NET1P_in_RENTP = Mass_Action_1_222(kdnet, rY[41]);

    // Degradation of CDC14 in RENT
    double Degradation_of_CDC14_in_RENT = Mass_Action_1_222(kd14, rY[40]);

    // Degradation of CDC14 in RENTP
    double Degradation_of_CDC14_in_RENTP = Mass_Action_1_222(kd14, rY[41]);

    // TEM1 activation
    double TEM1_activation = MichaelisMenten_220(this->mParameters[1], Jatem, 1, rY[49]);

    // inactivation
    double inactivation = MichaelisMenten_220(this->mParameters[0], Jitem, 1, rY[50]);

    // CDC15 activation
    double CDC15_activation = Mass_Action_1_222(ka15_p * rY[49] + ka15_p_p * rY[50] + ka15p * rY[6], rY[9]);

    // inactivation
    double inactivation_0 = Mass_Action_1_222(ki15, rY[8]);

    // PPX synthesis
    double PPX_synthesis = ksppx;

    // degradation
    double degradation = Mass_Action_1_222(Vdppx, rY[39]);

    // PDS1 synthesis
    double PDS1_synthesis = kspds_p + ks1pds_p_p * rY[42] + ks2pds_p_p * rY[32];

    // degradation
    double degradation_0 = Mass_Action_1_222(Vdpds, rY[37]);

    // Degradation of PDS1 in PE
    double Degradation_of_PDS1_in_PE = Mass_Action_1_222(Vdpds, rY[38]);

    // Assoc. with ESP1 to form PE
    double Assoc_with_ESP1_to_form_PE = Mass_Action_2_221(kasesp, rY[37], rY[24]);

    // Disso. from PE
    double Disso_from_PE = Mass_Action_1_222(kdiesp, rY[38]);

    // DNA synthesis
    double DNA_synthesis = ksori * (eorib5 * rY[20] + eorib2 * rY[18]);

    // Negative regulation of DNA synthesis
    double Negative_regulation_of_DNA_synthesis = Mass_Action_1_222(kdori, rY[36]);

    // Budding
    double Budding = ksbud * (ebudn2 * rY[22] + ebudn3 * rY[23] + ebudb5 * rY[20]);

    // Negative regulation of Cell budding
    double Negative_regulation_of_Cell_budding = Mass_Action_1_222(kdbud, rY[1]);

    // Spindle formation
    double Spindle_formation = ksspn * rY[18] / (Jspn + rY[18]);

    // Spindle disassembly
    double Spindle_disassembly = Mass_Action_1_222(kdspn, rY[46]);


    rDY[0] = (rDY[31]) / cell; // dBCK2/dt
    rDY[1] = (Budding - Negative_regulation_of_Cell_budding) / cell; // dBUD/dt
    rDY[2] = (Assoc_of_CLB2_and_SIC1 - Dissoc_of_CLB2SIC1_complex - Phosphorylation_of_C2 + Dephosphorylation_of_C2P - Degradation_of_CLB2_in_C2) / cell; // dC2/dt
    rDY[3] = (Phosphorylation_of_C2 - Dephosphorylation_of_C2P - Degradation_of_SIC1_in_C2P - Degradation_of_CLB2_in_C2P) / cell; // dC2P/dt
    rDY[4] = (Assoc_of_CLB5_and_SIC1 - Dissoc_of_CLB5SIC1 - Phosphorylation_of_C5 + Dephosphorylation_of_C5P - Degradation_of_CLB5_in_C5) / cell; // dC5/dt
    rDY[5] = (Phosphorylation_of_C5 - Dephosphorylation_of_C5P - Degradation_of_SIC1P_in_C5P - Degradation_of_CLB5_in_C5P) / cell; // dC5P/dt
    rDY[6] = (CDC14_synthesis - CDC14_degradation - Assoc_with_NET1_to_form_RENT + Dissoc_from_RENT - Assoc_with_NET1P_to_form_RENTP + Dissoc_from_RENP + Degradation_of_NET1_in_RENT + Degradation_of_NET1P_in_RENTP) / cell; // dCDC14/dt
    rDY[7] = (rDY[6] + rDY[40] + rDY[41]) / cell; // dCDC14T/dt
    rDY[8] = (CDC15_activation - inactivation_0) / cell; // dCDC15/dt
    rDY[9] = (-Degradation_of_active_CDC20 + Activation_of_CDC20 - Inactivation_0) / cell; // dCDC20/dt
    rDY[10] = (Synthesis_of_inactive_CDC20 - Degradation_of_inactiveCDC20 - Activation_of_CDC20 + Inactivation_0) / cell; // dCDC20i/dt
    rDY[11] = (CDC6_synthesis - Phosphorylation_of_CDC6 + Dephosphorylation_of_CDC6 - CLB2CDC6_complex_formation + CLB2CDC6_dissociation - CLB5CDC6_complex_formation + CLB5CDC6_dissociation + CLB2_degradation_in_F2 + CLB5_degradation_in_F5) / cell; // dCDC6/dt
    rDY[12] = (Phosphorylation_of_CDC6 - Dephosphorylation_of_CDC6 - Degradation_of_CDC6P + CLB2_degradation_in_F2P + CLB5_degradation_in_F5P) / cell; // dCDC6P/dt
    rDY[13] = (rDY[12] + rDY[25] + rDY[27] + rDY[13] + rDY[26] + rDY[28]) / cell; // dCDC6T/dt
    rDY[14] = (CDH1_synthesis - CDH1_degradation + CDH1i_activation - Inactivation_1) / cell; // dCDH1/dt
    rDY[15] = (-CDH1i_degradation - CDH1i_activation + Inactivation_1) / cell; // dCDH1i/dt
    rDY[16] = (rDY[45] + rDY[14]) / cell; // dCKIT/dt
    rDY[17] = (Synthesis_of_CLB2 - Degradation_of_CLB2 - Assoc_of_CLB2_and_SIC1 + Dissoc_of_CLB2SIC1_complex + Degradation_of_SIC1_in_C2P - CLB2CDC6_complex_formation + CLB2CDC6_dissociation + CDC6_degradation_in_F2P) / cell; // dCLB2/dt
    rDY[18] = (rDY[18] + rDY[2] + rDY[3] + rDY[25] + rDY[26]) / cell; // dCLB2T/dt
    rDY[19] = (Synthesis_of_CLB5 - Degradation_of_CLB5 - Assoc_of_CLB5_and_SIC1 + Dissoc_of_CLB5SIC1 + Degradation_of_SIC1P_in_C5P - CLB5CDC6_complex_formation + CLB5CDC6_dissociation + CDC6_degradation_in_F5P) / cell; // dCLB5/dt
    rDY[20] = (rDY[20] + rDY[4] + rDY[5] + rDY[27] + rDY[28]) / cell; // dCLB5T/dt
    rDY[21] = (Synthesis_of_CLN2 - Degradation_of_CLN2) / cell; // dCLN2/dt
    rDY[22] = (rDY[31]) / cell; // dCLN3/dt
    rDY[23] = (Degradation_of_PDS1_in_PE - Assoc_with_ESP1_to_form_PE + Disso_from_PE) / cell; // dESP1/dt
    rDY[24] = (CLB2CDC6_complex_formation - CLB2CDC6_dissociation - F2_phosphorylation + F2P_dephosphorylation - CLB2_degradation_in_F2) / cell; // dF2/dt
    rDY[25] = (F2_phosphorylation - F2P_dephosphorylation - CDC6_degradation_in_F2P - CLB2_degradation_in_F2P) / cell; // dF2P/dt
    rDY[26] = (CLB5CDC6_complex_formation - CLB5CDC6_dissociation - F5_phosphorylation + F5P_dephosphorylation - CLB5_degradation_in_F5) / cell; // dF5/dt
    rDY[27] = (F5_phosphorylation - F5P_dephosphorylation - CDC6_degradation_in_F5P - CLB5_degradation_in_F5P) / cell; // dF5P/dt
    rDY[28] = (Activation_of_IEP - Inactivation) / cell; // dIEP/dt
    rDY[29] = (Growth) / cell; // dMASS/dt
    rDY[30] = (-Assoc_with_NET1_to_form_RENT + Dissoc_from_RENT + Net1_synthesis - Net1_degradation - NET1_phosphorylation + dephosphorylation + Degradation_of_CDC14_in_RENT) / cell; // dNET1/dt
    rDY[31] = (-Assoc_with_NET1P_to_form_RENTP + Dissoc_from_RENP - Net1P_degradation + NET1_phosphorylation - dephosphorylation + Degradation_of_CDC14_in_RENTP) / cell; // dNET1P/dt
    rDY[32] = (rDY[33] + rDY[34] + rDY[40] + rDY[41]) / cell; // dNET1T/dt
    rDY[33] = (DNA_synthesis - Negative_regulation_of_DNA_synthesis) / cell; // dORI/dt
    rDY[34] = (PDS1_synthesis - degradation_0 - Assoc_with_ESP1_to_form_PE + Disso_from_PE) / cell; // dPDS1/dt
    rDY[35] = (PPX_synthesis - degradation) / cell; // dPPX/dt
    rDY[36] = (Assoc_with_NET1_to_form_RENT - Dissoc_from_RENT - RENT_phosphorylation + dephosphorylation_0 - Degradation_of_NET1_in_RENT - Degradation_of_CDC14_in_RENT) / cell; // dRENT/dt
    rDY[37] = (Assoc_with_NET1P_to_form_RENTP - Dissoc_from_RENP + RENT_phosphorylation - dephosphorylation_0 - Degradation_of_NET1P_in_RENTP - Degradation_of_CDC14_in_RENTP) / cell; // dRENTP/dt
    rDY[38] = (Synthesis_of_SIC1 - Phosphorylation_of_SIC1 + Dephosphorylation_of_SIC1 - Assoc_of_CLB2_and_SIC1 + Dissoc_of_CLB2SIC1_complex - Assoc_of_CLB5_and_SIC1 + Dissoc_of_CLB5SIC1 + Degradation_of_CLB2_in_C2 + Degradation_of_CLB5_in_C5) / cell; // dSIC1/dt
    rDY[39] = (Phosphorylation_of_SIC1 - Dephosphorylation_of_SIC1 - Fast_Degradation_of_SIC1P + Degradation_of_CLB2_in_C2P + Degradation_of_CLB5_in_C5P) / cell; // dSIC1P/dt
    rDY[40] = (rDY[43] + rDY[2] + rDY[4] + rDY[44] + rDY[3] + rDY[5]) / cell; // dSIC1T/dt
    rDY[41] = (Spindle_formation - Spindle_disassembly) / cell; // dSPN/dt
    rDY[42] = (Synthesis_of_SWI5 - Degradation_of_SWI5 + Activation_of_SWI5 - Inactivation_of_SWI5) / cell; // dSWI5/dt
    rDY[43] = (-Degradation_of_SWI5P - Activation_of_SWI5 + Inactivation_of_SWI5) / cell; // dSWI5P/dt
    rDY[44] = (TEM1_activation - inactivation) / cell; // dTEM1GTP/dt

    /* Account for the differences in timescales. */
}

bool Chen2004OdeSystem::CalculateStoppingEvent(double time, const std::vector<double> & rY)
{
    // Return true if all events have been triggered.
    return AreAllEventsSatisfied(time, rY);
}

void Chen2004OdeSystem::CheckAndUpdateEvents(double time, const std::vector<double> & rY)
{
    std::vector<double> dy(rY.size()); // Initialise derivatives vector
    EvaluateYDerivatives(time, rY, dy);

    if (sbmlmath::sm_lt(rY[18] + rY[20] - KEZ2, 0))
    {
        this->SetStateVariable("ORI", static_cast<double>(0));
        eventsSatisfied[0] = true;
    }
    if (sbmlmath::sm_gt(rY[36] - 1, 0))
    {
        this->SetParameter("MAD2", static_cast<double>(mad2h));
        this->SetParameter("BUB2", static_cast<double>(bub2h));
        eventsSatisfied[1] = true;
    }
    if (sbmlmath::sm_gt(rY[46] - 1, 0))
    {
        this->SetParameter("MAD2", static_cast<double>(mad2l));
        this->SetParameter("LTE1", static_cast<double>(lte1h));
        this->SetParameter("BUB2", static_cast<double>(bub2l));
        eventsSatisfied[2] = true;
    }
    if (sbmlmath::sm_lt(rY[18] - KEZ, 0))
    {
        this->SetStateVariable("MASS", static_cast<double>(F * rY[31]));
        this->SetParameter("LTE1", static_cast<double>(lte1l));
        this->SetStateVariable("BUD", static_cast<double>(0));
        this->SetStateVariable("SPN", static_cast<double>(0));
        eventsSatisfied[3] = true;
    }
}

bool Chen2004OdeSystem::AreAllEventsSatisfied(double time, const std::vector<double>& rY)
{
    CheckAndUpdateEvents(time, rY);
    bool events_satisfied = true;
    if (std::find(eventsSatisfied.begin(), eventsSatisfied.end(), false) != eventsSatisfied.end())
    {
        events_satisfied = false;
    }
    if (events_satisfied) // Reset events vector if cell division is triggered
    {
        std::fill(eventsSatisfied.begin(), eventsSatisfied.end(), false);
    }
    return events_satisfied;
}

template <>
void CellwiseOdeSystemInformation<Chen2004OdeSystem>::Initialise()
{
    this->mVariableNames.push_back("BCK2");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

    this->mParameterNames.push_back("BUB2");
    this->mParameterUnits.push_back("non-dim");

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

    this->mParameterNames.push_back("LTE1");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("MAD2");
    this->mParameterUnits.push_back("non-dim");

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


    /* Define state parameters. */
    // Parameters without set values must be externally defined
    this->mInitialised = true;
}

/* Define SbmlCellCycleWrapperModel using wrappers. */
#include "SbmlCellCycleWrapperModel.hpp"
#include "SbmlCellCycleWrapperModel.cpp"

typedef SbmlCellCycleWrapperModel<Chen2004OdeSystem, 51> Chen2004CellCycleModel;

// Declare identifiers for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(Chen2004OdeSystem)
EXPORT_TEMPLATE_CLASS2(SbmlCellCycleWrapperModel, Chen2004OdeSystem, 51)

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(Chen2004CellCycleModel)