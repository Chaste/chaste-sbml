#include <cmath>
#include <limits>
#include <vector>

#include "CellwiseOdeSystemInformation.hpp"
#include "SbmlEventType.hpp"
#include "SbmlMath.hpp"

#include "Chen2004SbmlOdeSystem.hpp"

namespace sm = sbmlmath;

namespace
{
// Convert the model's native time units (minutes) to Chaste default (hours) and
// scale the derivatives by this factor (minutes per hour).
constexpr double TIMESCALE_MULTIPLIER = 60.0;
} // namespace

Chen2004SbmlOdeSystem::Chen2004SbmlOdeSystem()
        : AbstractSbmlOdeSystem(36, 146, 4)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<Chen2004SbmlOdeSystem>);

    Initialise();

    // EVENTS
    mEventType.resize(4, SbmlEventType::UNKNOWN);

    // Uncomment lines below for events that should trigger cell division

    // mEventType[0] = SbmlEventType::CELL_DIVISION; // reset ORI
    // mEventType[1] = SbmlEventType::CELL_DIVISION; // start DNA synthesis
    // mEventType[2] = SbmlEventType::CELL_DIVISION; // spindle checkpoint
    mEventType[3] = SbmlEventType::CELL_DIVISION; // cell division

    mEventSatisfied = { true, true, true, true }; // From SBML trigger initialValue
    mEventTriggered.resize(4, false);

    mEventAdjustedParameters.resize(146, false);
    mEventAdjustedParameterValues.resize(146, 0.0);

    mEventAdjustedStateVars.resize(36, false);
    mEventAdjustedStateValues.resize(36, 0.0);
}

Chen2004SbmlOdeSystem::~Chen2004SbmlOdeSystem()
{
}

std::vector<double> Chen2004SbmlOdeSystem::ComputeDerivedQuantities(double time, const std::vector<double>& rY)
{
    std::vector<double> dqs;
    dqs.reserve(184);
    time *= TIMESCALE_MULTIPLIER; // Chaste integrates in hours; use the model's native time units
    RunModelEquations(time, rY);

    // AMOUNT / CONCENTRATION CONVERSIONS
    double conc__BCK2 = BCK2 / cell;       //
    double conc__BUB2 = BUB2 / cell;       //
    double conc__BUD = BUD / cell;         //
    double conc__C2 = C2 / cell;           //
    double conc__C2P = C2P / cell;         //
    double conc__C5 = C5 / cell;           //
    double conc__C5P = C5P / cell;         //
    double conc__CDC14 = CDC14 / cell;     //
    double conc__CDC14T = CDC14T / cell;   //
    double conc__CDC15 = CDC15 / cell;     //
    double conc__CDC15i = CDC15i / cell;   //
    double conc__CDC20 = CDC20 / cell;     //
    double conc__CDC20i = CDC20i / cell;   //
    double conc__CDC6 = CDC6 / cell;       //
    double conc__CDC6P = CDC6P / cell;     //
    double conc__CDC6T = CDC6T / cell;     //
    double conc__CDH1 = CDH1 / cell;       //
    double conc__CDH1i = CDH1i / cell;     //
    double conc__CKIT = CKIT / cell;       //
    double conc__CLB2 = CLB2 / cell;       //
    double conc__CLB2T = CLB2T / cell;     //
    double conc__CLB5 = CLB5 / cell;       //
    double conc__CLB5T = CLB5T / cell;     //
    double conc__CLN2 = CLN2 / cell;       //
    double conc__CLN3 = CLN3 / cell;       //
    double conc__ESP1 = ESP1 / cell;       //
    double conc__F2 = F2 / cell;           //
    double conc__F2P = F2P / cell;         //
    double conc__F5 = F5 / cell;           //
    double conc__F5P = F5P / cell;         //
    double conc__IE = IE / cell;           //
    double conc__IEP = IEP / cell;         //
    double conc__LTE1 = LTE1 / cell;       //
    double conc__MAD2 = MAD2 / cell;       //
    double conc__MASS = MASS / cell;       //
    double conc__MCM1 = MCM1 / cell;       //
    double conc__NET1 = NET1 / cell;       //
    double conc__NET1P = NET1P / cell;     //
    double conc__NET1T = NET1T / cell;     //
    double conc__ORI = ORI / cell;         //
    double conc__PDS1 = PDS1 / cell;       //
    double conc__PE = PE / cell;           //
    double conc__PPX = PPX / cell;         //
    double conc__RENT = RENT / cell;       //
    double conc__RENTP = RENTP / cell;     //
    double conc__SBF = SBF / cell;         //
    double conc__SIC1 = SIC1 / cell;       //
    double conc__SIC1P = SIC1P / cell;     //
    double conc__SIC1T = SIC1T / cell;     //
    double conc__SPN = SPN / cell;         //
    double conc__SWI5 = SWI5 / cell;       //
    double conc__SWI5P = SWI5P / cell;     //
    double conc__TEM1GDP = TEM1GDP / cell; //
    double conc__TEM1GTP = TEM1GTP / cell; //

    dqs.push_back(cell);
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
    dqs.push_back(D);
    dqs.push_back(mu);
    dqs.push_back(Vdb5);
    dqs.push_back(Vdb2);
    dqs.push_back(Vasbf);
    dqs.push_back(Visbf);
    dqs.push_back(Vkpc1);
    dqs.push_back(Vkpf6);
    dqs.push_back(Vacdh);
    dqs.push_back(Vicdh);
    dqs.push_back(Vppnet);
    dqs.push_back(Vkpnet);
    dqs.push_back(Vdppx);
    dqs.push_back(Vdpds);
    dqs.push_back(Vaiep);
    dqs.push_back(Vd2c1);
    dqs.push_back(Vd2f6);
    dqs.push_back(Vppc1);
    dqs.push_back(Vppf6);
    dqs.push_back(F);
    dqs.push_back(Growth);
    dqs.push_back(Synthesis_of_CLN2);
    dqs.push_back(Degradation_of_CLN2);
    dqs.push_back(Synthesis_of_CLB2);
    dqs.push_back(Degradation_of_CLB2);
    dqs.push_back(Synthesis_of_CLB5);
    dqs.push_back(Degradation_of_CLB5);
    dqs.push_back(Synthesis_of_SIC1);
    dqs.push_back(Phosphorylation_of_SIC1);
    dqs.push_back(Dephosphorylation_of_SIC1);
    dqs.push_back(Fast_Degradation_of_SIC1P);
    dqs.push_back(Assoc_of_CLB2_and_SIC1);
    dqs.push_back(Dissoc_of_CLB2SIC1_complex);
    dqs.push_back(Assoc_of_CLB5_and_SIC1);
    dqs.push_back(Dissoc_of_CLB5SIC1);
    dqs.push_back(Phosphorylation_of_C2);
    dqs.push_back(Dephosphorylation_of_C2P);
    dqs.push_back(Phosphorylation_of_C5);
    dqs.push_back(Dephosphorylation_of_C5P);
    dqs.push_back(Degradation_of_CLB2_in_C2);
    dqs.push_back(Degradation_of_CLB5_in_C5);
    dqs.push_back(Degradation_of_SIC1_in_C2P);
    dqs.push_back(Degradation_of_SIC1P_in_C5P_);
    dqs.push_back(Degradation_of_CLB2_in_C2P);
    dqs.push_back(Degradation_of_CLB5_in_C5P);
    dqs.push_back(CDC6_synthesis);
    dqs.push_back(Phosphorylation_of_CDC6);
    dqs.push_back(Dephosphorylation_of_CDC6);
    dqs.push_back(Degradation_of_CDC6P);
    dqs.push_back(CLB2CDC6_complex_formation);
    dqs.push_back(CLB2CDC6_dissociation);
    dqs.push_back(CLB5CDC6_complex_formation);
    dqs.push_back(CLB5CDC6_dissociation);
    dqs.push_back(F2_phosphorylation);
    dqs.push_back(F2P_dephosphorylation);
    dqs.push_back(F5_phosphorylation);
    dqs.push_back(F5P_dephosphorylation);
    dqs.push_back(CLB2_degradation_in_F2);
    dqs.push_back(CLB5_degradation_in_F5);
    dqs.push_back(CDC6_degradation_in_F2P);
    dqs.push_back(CDC6_degradation_in_F5P);
    dqs.push_back(CLB2_degradation_in_F2P);
    dqs.push_back(CLB5_degradation_in_F5P);
    dqs.push_back(Synthesis_of_SWI5);
    dqs.push_back(Degradation_of_SWI5);
    dqs.push_back(Degradation_of_SWI5P);
    dqs.push_back(Activation_of_SWI5);
    dqs.push_back(Inactivation_of_SWI5);
    dqs.push_back(Activation_of_IEP);
    dqs.push_back(Inactivation_1);
    dqs.push_back(Synthesis_of_inactive_CDC20);
    dqs.push_back(Degradation_of_inactiveCDC20);
    dqs.push_back(Degradation_of_active_CDC20);
    dqs.push_back(Activation_of_CDC20);
    dqs.push_back(Inactivation_2);
    dqs.push_back(CDH1_synthesis);
    dqs.push_back(CDH1_degradation);
    dqs.push_back(CDH1i_degradation);
    dqs.push_back(CDH1i_activation);
    dqs.push_back(Inactivation_3);
    dqs.push_back(CDC14_synthesis);
    dqs.push_back(CDC14_degradation);
    dqs.push_back(Assoc_with_NET1_to_form_RENT);
    dqs.push_back(Dissoc_from_RENT);
    dqs.push_back(Assoc_with_NET1P_to_form_RENTP);
    dqs.push_back(Dissoc_from_RENP);
    dqs.push_back(Net1_synthesis);
    dqs.push_back(Net1_degradation);
    dqs.push_back(Net1P_degradation);
    dqs.push_back(NET1_phosphorylation);
    dqs.push_back(dephosphorylation_1);
    dqs.push_back(RENT_phosphorylation);
    dqs.push_back(dephosphorylation_2);
    dqs.push_back(Degradation_of_NET1_in_RENT);
    dqs.push_back(Degradation_of_NET1P_in_RENTP);
    dqs.push_back(Degradation_of_CDC14_in_RENT);
    dqs.push_back(Degradation_of_CDC14_in_RENTP);
    dqs.push_back(TEM1_activation);
    dqs.push_back(inactivation_1);
    dqs.push_back(CDC15_activation);
    dqs.push_back(inactivation_2);
    dqs.push_back(PPX_synthesis);
    dqs.push_back(degradation_1);
    dqs.push_back(PDS1_synthesis);
    dqs.push_back(degradation_2);
    dqs.push_back(Degradation_of_PDS1_in_PE);
    dqs.push_back(Assoc_with_ESP1_to_form_PE);
    dqs.push_back(Disso_from_PE);
    dqs.push_back(DNA_synthesis);
    dqs.push_back(Negative_regulation_of_DNA_synthesis);
    dqs.push_back(Budding);
    dqs.push_back(Negative_regulation_of_Cell_budding);
    dqs.push_back(Spindle_formation);
    dqs.push_back(Spindle_disassembly);
    dqs.push_back(conc__BCK2);
    dqs.push_back(conc__BUB2);
    dqs.push_back(conc__BUD);
    dqs.push_back(conc__C2);
    dqs.push_back(conc__C2P);
    dqs.push_back(conc__C5);
    dqs.push_back(conc__C5P);
    dqs.push_back(conc__CDC14);
    dqs.push_back(conc__CDC14T);
    dqs.push_back(conc__CDC15);
    dqs.push_back(conc__CDC15i);
    dqs.push_back(conc__CDC20);
    dqs.push_back(conc__CDC20i);
    dqs.push_back(conc__CDC6);
    dqs.push_back(conc__CDC6P);
    dqs.push_back(conc__CDC6T);
    dqs.push_back(conc__CDH1);
    dqs.push_back(conc__CDH1i);
    dqs.push_back(conc__CKIT);
    dqs.push_back(conc__CLB2);
    dqs.push_back(conc__CLB2T);
    dqs.push_back(conc__CLB5);
    dqs.push_back(conc__CLB5T);
    dqs.push_back(conc__CLN2);
    dqs.push_back(conc__CLN3);
    dqs.push_back(conc__ESP1);
    dqs.push_back(conc__F2);
    dqs.push_back(conc__F2P);
    dqs.push_back(conc__F5);
    dqs.push_back(conc__F5P);
    dqs.push_back(conc__IE);
    dqs.push_back(conc__IEP);
    dqs.push_back(conc__LTE1);
    dqs.push_back(conc__MAD2);
    dqs.push_back(conc__MASS);
    dqs.push_back(conc__MCM1);
    dqs.push_back(conc__NET1);
    dqs.push_back(conc__NET1P);
    dqs.push_back(conc__NET1T);
    dqs.push_back(conc__ORI);
    dqs.push_back(conc__PDS1);
    dqs.push_back(conc__PE);
    dqs.push_back(conc__PPX);
    dqs.push_back(conc__RENT);
    dqs.push_back(conc__RENTP);
    dqs.push_back(conc__SBF);
    dqs.push_back(conc__SIC1);
    dqs.push_back(conc__SIC1P);
    dqs.push_back(conc__SIC1T);
    dqs.push_back(conc__SPN);
    dqs.push_back(conc__SWI5);
    dqs.push_back(conc__SWI5P);
    dqs.push_back(conc__TEM1GDP);
    dqs.push_back(conc__TEM1GTP);

    return dqs;
} // LCOV_EXCL_LINE (gcov marks this closing brace of a std::vector-returning function as uncovered)

void Chen2004SbmlOdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    // Convert the model's native time units to Chaste default (hours) and scale the derivatives.
    time *= TIMESCALE_MULTIPLIER;
    std::vector<double> derivatives = RunModelEquations(time, rY);
    for (unsigned i = 0; i < rDY.size(); ++i)
    {
        rDY[i] = TIMESCALE_MULTIPLIER * derivatives[i];
    }
}

void Chen2004SbmlOdeSystem::Initialise(double time)
{
    // This does NOT scale time as Initialise only runs at time=0 from the constructor.
    cell = 1.0;                                                                                        //
    BUB2 = 0.2;                                                                                        //
    BUD = 0.008473;                                                                                    //
    C2 = 0.238404;                                                                                     //
    C2P = 0.024034;                                                                                    //
    C5 = 0.070081;                                                                                     //
    C5P = 0.006878;                                                                                    //
    CDC14 = 0.468344;                                                                                  //
    CDC14T = 2.0;                                                                                      //
    CDC15 = 0.656533;                                                                                  //
    CDC20 = 0.444296;                                                                                  //
    CDC20i = 1.472044;                                                                                 //
    CDC6 = 0.10758;                                                                                    //
    CDC6P = 0.015486;                                                                                  //
    CDH1 = 0.930499;                                                                                   //
    CDH1i = 0.0695;                                                                                    //
    CLB2 = 0.1469227;                                                                                  //
    CLB2T = 0.17;                                                                                      //
    CLB5 = 0.0518014;                                                                                  //
    CLB5T = 0.12;                                                                                      //
    CLN2 = 0.0652511;                                                                                  //
    ESP1 = 0.301313;                                                                                   //
    F2 = 0.236058;                                                                                     //
    F2P = 0.0273938;                                                                                   //
    F5 = 7.24e-5;                                                                                      //
    F5P = 7.91e-5;                                                                                     //
    IEP = 0.1015;                                                                                      //
    LTE1 = 0.1;                                                                                        //
    MAD2 = 0.01;                                                                                       //
    MASS = 1.206019;                                                                                   //
    NET1 = 0.018645;                                                                                   //
    NET1P = 0.970271;                                                                                  //
    NET1T = 2.8;                                                                                       //
    ORI = 0.000909;                                                                                    //
    PDS1 = 0.025612;                                                                                   //
    PPX = 0.123179;                                                                                    //
    RENT = 1.04954;                                                                                    //
    RENTP = 0.6;                                                                                       //
    SIC1 = 0.0228776;                                                                                  //
    SIC1P = 0.00641;                                                                                   //
    SPN = 0.03;                                                                                        //
    SWI5 = 0.95;                                                                                       //
    SWI5P = 0.02;                                                                                      //
    TEM1GTP = 0.9;                                                                                     //
    b0 = 0.054;                                                                                        //
    bub2h = 1.0;                                                                                       //
    bub2l = 0.2;                                                                                       //
    C0 = 0.4;                                                                                          //
    CDC15T = 1.0;                                                                                      //
    Dn3 = 1.0;                                                                                         //
    ebudb5 = 1.0;                                                                                      //
    ebudn2 = 0.25;                                                                                     //
    ebudn3 = 0.05;                                                                                     //
    ec1b2 = 0.45;                                                                                      //
    ec1b5 = 0.1;                                                                                       //
    ec1k2 = 0.03;                                                                                      //
    ec1n2 = 0.06;                                                                                      //
    ec1n3 = 0.3;                                                                                       //
    ef6b2 = 0.55;                                                                                      //
    ef6b5 = 0.1;                                                                                       //
    ef6k2 = 0.03;                                                                                      //
    ef6n2 = 0.06;                                                                                      //
    ef6n3 = 0.3;                                                                                       //
    eicdhb2 = 1.2;                                                                                     //
    eicdhb5 = 8.0;                                                                                     //
    eicdhn2 = 0.4;                                                                                     //
    eicdhn3 = 0.25;                                                                                    //
    eorib2 = 0.45;                                                                                     //
    eorib5 = 0.9;                                                                                      //
    esbfb5 = 2.0;                                                                                      //
    esbfn2 = 2.0;                                                                                      //
    esbfn3 = 10.0;                                                                                     //
    ESP1T = 1.0;                                                                                       //
    IET = 1.0;                                                                                         //
    J20ppx = 0.15;                                                                                     //
    Jacdh = 0.03;                                                                                      //
    Jaiep = 0.1;                                                                                       //
    Jamcm = 0.1;                                                                                       //
    Jasbf = 0.01;                                                                                      //
    Jatem = 0.1;                                                                                       //
    Jd2c1 = 0.05;                                                                                      //
    Jd2f6 = 0.05;                                                                                      //
    Jicdh = 0.03;                                                                                      //
    Jiiep = 0.1;                                                                                       //
    Jimcm = 0.1;                                                                                       //
    Jisbf = 0.01;                                                                                      //
    Jitem = 0.1;                                                                                       //
    Jn3 = 6.0;                                                                                         //
    Jpds = 0.04;                                                                                       //
    Jspn = 0.14;                                                                                       //
    ka15_p = 0.002;                                                                                    //
    ka15_p_p = 1.0;                                                                                    //
    ka15p = 0.001;                                                                                     //
    ka20_p = 0.05;                                                                                     //
    ka20_p_p = 0.2;                                                                                    //
    kacdh_p = 0.01;                                                                                    //
    kacdh_p_p = 0.8;                                                                                   //
    kaiep = 0.1;                                                                                       //
    kamcm = 1.0;                                                                                       //
    kasb2 = 50.0;                                                                                      //
    kasb5 = 50.0;                                                                                      //
    kasbf = 0.38;                                                                                      //
    kasesp = 50.0;                                                                                     //
    kasf2 = 15.0;                                                                                      //
    kasf5 = 0.01;                                                                                      //
    kasrent = 200.0;                                                                                   //
    kasrentp = 1.0;                                                                                    //
    kaswi = 2.0;                                                                                       //
    kd14 = 0.1;                                                                                        //
    kd1c1 = 0.01;                                                                                      //
    kd1f6 = 0.01;                                                                                      //
    kd1pds_p = 0.01;                                                                                   //
    kd20 = 0.3;                                                                                        //
    kd2c1 = 1.0;                                                                                       //
    kd2f6 = 1.0;                                                                                       //
    kd2pds_p_p = 0.2;                                                                                  //
    kd3c1 = 1.0;                                                                                       //
    kd3f6 = 1.0;                                                                                       //
    kd3pds_p_p = 0.04;                                                                                 //
    kdb2_p = 0.003;                                                                                    //
    kdb2_p_p = 0.4;                                                                                    //
    kdb2p = 0.15;                                                                                      //
    kdb5_p = 0.01;                                                                                     //
    kdb5_p_p = 0.16;                                                                                   //
    kdbud = 0.06;                                                                                      //
    kdcdh = 0.01;                                                                                      //
    kdib2 = 0.05;                                                                                      //
    kdib5 = 0.06;                                                                                      //
    kdiesp = 0.5;                                                                                      //
    kdif2 = 0.5;                                                                                       //
    kdif5 = 0.01;                                                                                      //
    kdirent = 1.0;                                                                                     //
    kdirentp = 2.0;                                                                                    //
    kdn2 = 0.12;                                                                                       //
    kdnet = 0.03;                                                                                      //
    kdori = 0.06;                                                                                      //
    kdppx_p = 0.17;                                                                                    //
    kdppx_p_p = 2.0;                                                                                   //
    kdspn = 0.06;                                                                                      //
    kdswi = 0.08;                                                                                      //
    KEZ = 0.3;                                                                                         //
    KEZ2 = 0.2;                                                                                        //
    ki15 = 0.5;                                                                                        //
    kicdh_p = 0.001;                                                                                   //
    kicdh_p_p = 0.08;                                                                                  //
    kiiep = 0.15;                                                                                      //
    kimcm = 0.15;                                                                                      //
    kisbf_p = 0.6;                                                                                     //
    kisbf_p_p = 8.0;                                                                                   //
    kiswi = 0.05;                                                                                      //
    kkpnet_p = 0.01;                                                                                   //
    kkpnet_p_p = 0.6;                                                                                  //
    kppc1 = 4.0;                                                                                       //
    kppf6 = 4.0;                                                                                       //
    kppnet_p = 0.05;                                                                                   //
    kppnet_p_p = 3.0;                                                                                  //
    ks14 = 0.2;                                                                                        //
    ks1pds_p_p = 0.03;                                                                                 //
    ks20_p = 0.006;                                                                                    //
    ks20_p_p = 0.6;                                                                                    //
    ks2pds_p_p = 0.055;                                                                                //
    ksb2_p = 0.001;                                                                                    //
    ksb2_p_p = 0.04;                                                                                   //
    ksb5_p = 0.0008;                                                                                   //
    ksb5_p_p = 0.005;                                                                                  //
    ksbud = 0.2;                                                                                       //
    ksc1_p = 0.012;                                                                                    //
    ksc1_p_p = 0.12;                                                                                   //
    kscdh = 0.01;                                                                                      //
    ksf6_p = 0.024;                                                                                    //
    ksf6_p_p = 0.12;                                                                                   //
    ksf6_p_p_p = 0.004;                                                                                //
    ksn2_p = 0.0;                                                                                      //
    ksn2_p_p = 0.15;                                                                                   //
    ksnet = 0.084;                                                                                     //
    ksori = 2.0;                                                                                       //
    kspds_p = 0.0;                                                                                     //
    ksppx = 0.1;                                                                                       //
    ksspn = 0.1;                                                                                       //
    ksswi_p = 0.005;                                                                                   //
    ksswi_p_p = 0.08;                                                                                  //
    lte1h = 1.0;                                                                                       //
    lte1l = 0.1;                                                                                       //
    mad2h = 8.0;                                                                                       //
    mad2l = 0.01;                                                                                      //
    mdt = 90.0;                                                                                        //
    TEM1T = 1.0;                                                                                       //
    BCK2 = b0 * MASS;                                                                                  //
    CDC14T = CDC14 + RENT + RENTP;                                                                     //
    CDC15i = CDC15T - CDC15;                                                                           //
    CDC6T = CDC6 + F2 + F5 + CDC6P + F2P + F5P;                                                        //
    SIC1T = SIC1 + C2 + C5 + SIC1P + C2P + C5P;                                                        //
    CKIT = SIC1T + CDC6T;                                                                              //
    CLB2T = CLB2 + C2 + C2P + F2 + F2P;                                                                //
    CLB5T = CLB5 + C5 + C5P + F5 + F5P;                                                                //
    CLN3 = C0 * Dn3 * MASS / (Jn3 + Dn3 * MASS);                                                       //
    IE = IET - IEP;                                                                                    //
    MCM1 = GK_219(kamcm * CLB2, kimcm, Jamcm, Jimcm);                                                  //
    NET1T = NET1 + NET1P + RENT + RENTP;                                                               //
    PE = ESP1T - ESP1;                                                                                 //
    Vasbf = kasbf * (esbfn2 * CLN2 + esbfn3 * (CLN3 + BCK2) + esbfb5 * CLB5);                          //
    Visbf = kisbf_p + kisbf_p_p * CLB2;                                                                //
    SBF = GK_219(Vasbf, Visbf, Jasbf, Jisbf);                                                          //
    TEM1GDP = TEM1T - TEM1GTP;                                                                         //
    mu = std::log(2.0) / mdt;                                                                          //
    D = 1.026 / mu - 32.0;                                                                             //
    Vdb5 = kdb5_p + kdb5_p_p * CDC20;                                                                  //
    Vdb2 = kdb2_p + kdb2_p_p * CDH1 + kdb2p * CDC20;                                                   //
    Vd2c1 = kd2c1 * (ec1n3 * CLN3 + ec1k2 * BCK2 + ec1n2 * CLN2 + ec1b5 * CLB5 + ec1b2 * CLB2);        //
    Vkpc1 = kd1c1 + Vd2c1 / (Jd2c1 + SIC1 + C2 + C5 + SIC1P + C2P + C5P);                              //
    Vd2f6 = kd2f6 * (ef6n3 * CLN3 + ef6k2 * BCK2 + ef6n2 * CLN2 + ef6b5 * CLB5 + ef6b2 * CLB2);        //
    Vkpf6 = kd1f6 + Vd2f6 / (Jd2f6 + CDC6 + F2 + F5 + CDC6P + F2P + F5P);                              //
    Vacdh = kacdh_p + kacdh_p_p * CDC14;                                                               //
    Vicdh = kicdh_p + kicdh_p_p * (eicdhn3 * CLN3 + eicdhn2 * CLN2 + eicdhb5 * CLB5 + eicdhb2 * CLB2); //
    Vppnet = kppnet_p + kppnet_p_p * PPX;                                                              //
    Vkpnet = (kkpnet_p + kkpnet_p_p * CDC15) * MASS;                                                   //
    Vdppx = kdppx_p + kdppx_p_p * (J20ppx + CDC20) * Jpds / (Jpds + PDS1);                             //
    Vdpds = kd1pds_p + kd2pds_p_p * CDC20 + kd3pds_p_p * CDH1;                                         //
    Vaiep = kaiep * CLB2;                                                                              //
    Vppc1 = kppc1 * CDC14;                                                                             //
    Vppf6 = kppf6 * CDC14;                                                                             //
    F = std::exp(-mu * D);                                                                             //
    Growth = mu * MASS;                                                                                //
    Synthesis_of_CLN2 = (ksn2_p + ksn2_p_p * SBF) * MASS;                                              //
    Degradation_of_CLN2 = Mass_Action_1_222(kdn2, CLN2);                                               //
    Synthesis_of_CLB2 = (ksb2_p + ksb2_p_p * MCM1) * MASS;                                             //
    Degradation_of_CLB2 = Mass_Action_1_222(Vdb2, CLB2);                                               //
    Synthesis_of_CLB5 = (ksb5_p + ksb5_p_p * SBF) * MASS;                                              //
    Degradation_of_CLB5 = Mass_Action_1_222(Vdb5, CLB5);                                               //
    Synthesis_of_SIC1 = ksc1_p + ksc1_p_p * SWI5;                                                      //
    Phosphorylation_of_SIC1 = Mass_Action_1_222(Vkpc1, SIC1);                                          //
    Dephosphorylation_of_SIC1 = Mass_Action_1_222(Vppc1, SIC1P);                                       //
    Fast_Degradation_of_SIC1P = Mass_Action_1_222(kd3c1, SIC1P);                                       //
    Assoc_of_CLB2_and_SIC1 = Mass_Action_2_221(kasb2, CLB2, SIC1);                                     //
    Dissoc_of_CLB2SIC1_complex = Mass_Action_1_222(kdib2, C2);                                         //
    Assoc_of_CLB5_and_SIC1 = Mass_Action_2_221(kasb5, CLB5, SIC1);                                     //
    Dissoc_of_CLB5SIC1 = Mass_Action_1_222(kdib5, C5);                                                 //
    Phosphorylation_of_C2 = Mass_Action_1_222(Vkpc1, C2);                                              //
    Dephosphorylation_of_C2P = Mass_Action_1_222(Vppc1, C2P);                                          //
    Phosphorylation_of_C5 = Mass_Action_1_222(Vkpc1, C5);                                              //
    Dephosphorylation_of_C5P = Mass_Action_1_222(Vppc1, C5P);                                          //
    Degradation_of_CLB2_in_C2 = Mass_Action_1_222(Vdb2, C2);                                           //
    Degradation_of_CLB5_in_C5 = Mass_Action_1_222(Vdb5, C5);                                           //
    Degradation_of_SIC1_in_C2P = Mass_Action_1_222(kd3c1, C2P);                                        //
    Degradation_of_SIC1P_in_C5P_ = Mass_Action_1_222(kd3c1, C5P);                                      //
    Degradation_of_CLB2_in_C2P = Mass_Action_1_222(Vdb2, C2P);                                         //
    Degradation_of_CLB5_in_C5P = Mass_Action_1_222(Vdb5, C5P);                                         //
    CDC6_synthesis = ksf6_p + ksf6_p_p * SWI5 + ksf6_p_p_p * SBF;                                      //
    Phosphorylation_of_CDC6 = Mass_Action_1_222(Vkpf6, CDC6);                                          //
    Dephosphorylation_of_CDC6 = Mass_Action_1_222(Vppf6, CDC6P);                                       //
    Degradation_of_CDC6P = Mass_Action_1_222(kd3f6, CDC6P);                                            //
    CLB2CDC6_complex_formation = Mass_Action_2_221(kasf2, CLB2, CDC6);                                 //
    CLB2CDC6_dissociation = Mass_Action_1_222(kdif2, F2);                                              //
    CLB5CDC6_complex_formation = Mass_Action_2_221(kasf5, CLB5, CDC6);                                 //
    CLB5CDC6_dissociation = Mass_Action_1_222(kdif5, F5);                                              //
    F2_phosphorylation = Mass_Action_1_222(Vkpf6, F2);                                                 //
    F2P_dephosphorylation = Mass_Action_1_222(Vppf6, F2P);                                             //
    F5_phosphorylation = Mass_Action_1_222(Vkpf6, F5);                                                 //
    F5P_dephosphorylation = Mass_Action_1_222(Vppf6, F5P);                                             //
    CLB2_degradation_in_F2 = Mass_Action_1_222(Vdb2, F2);                                              //
    CLB5_degradation_in_F5 = Mass_Action_1_222(Vdb5, F5);                                              //
    CDC6_degradation_in_F2P = Mass_Action_1_222(kd3f6, F2P);                                           //
    CDC6_degradation_in_F5P = Mass_Action_1_222(kd3f6, F5P);                                           //
    CLB2_degradation_in_F2P = Mass_Action_1_222(Vdb2, F2P);                                            //
    CLB5_degradation_in_F5P = Mass_Action_1_222(Vdb5, F5P);                                            //
    Synthesis_of_SWI5 = ksswi_p + ksswi_p_p * MCM1;                                                    //
    Degradation_of_SWI5 = Mass_Action_1_222(kdswi, SWI5);                                              //
    Degradation_of_SWI5P = Mass_Action_1_222(kdswi, SWI5P);                                            //
    Activation_of_SWI5 = Mass_Action_1_222(kaswi * CDC14, SWI5P);                                      //
    Inactivation_of_SWI5 = Mass_Action_1_222(kiswi * CLB2, SWI5);                                      //
    Activation_of_IEP = MichaelisMenten_220(Vaiep, Jaiep, 1.0, IE);                                    //
    Inactivation_1 = MichaelisMenten_220(1.0, Jiiep, kiiep, IEP);                                      //
    Synthesis_of_inactive_CDC20 = ks20_p + ks20_p_p * MCM1;                                            //
    Degradation_of_inactiveCDC20 = Mass_Action_1_222(kd20, CDC20i);                                    //
    Degradation_of_active_CDC20 = Mass_Action_1_222(kd20, CDC20);                                      //
    Activation_of_CDC20 = Mass_Action_1_222(ka20_p + ka20_p_p * IEP, CDC20i);                          //
    // Inactivation_2:
    {
        [[maybe_unused]] double k = 1.0;
        // Qualify with this-> so the assignment targets the reaction member even when a local
        // parameter shadows its name (an SBML local parameter may shadow the reaction ID, e.g. a
        // reaction J1 with a local parameter also named J1).
        this->Inactivation_2 = k * Mass_Action_1_222(MAD2, CDC20);
    }
    CDH1_synthesis = kscdh;                                                                                                                                                                                                                                //
    CDH1_degradation = Mass_Action_1_222(kdcdh, CDH1);                                                                                                                                                                                                     //
    CDH1i_degradation = Mass_Action_1_222(kdcdh, CDH1i);                                                                                                                                                                                                   //
    CDH1i_activation = MichaelisMenten_220(Vacdh, Jacdh, 1.0, CDH1i);                                                                                                                                                                                      //
    Inactivation_3 = MichaelisMenten_220(Vicdh, Jicdh, 1.0, CDH1);                                                                                                                                                                                         //
    CDC14_synthesis = ks14;                                                                                                                                                                                                                                //
    CDC14_degradation = Mass_Action_1_222(kd14, CDC14);                                                                                                                                                                                                    //
    Assoc_with_NET1_to_form_RENT = Mass_Action_2_221(kasrent, CDC14, NET1);                                                                                                                                                                                //
    Dissoc_from_RENT = Mass_Action_1_222(kdirent, RENT);                                                                                                                                                                                                   //
    Assoc_with_NET1P_to_form_RENTP = Mass_Action_2_221(kasrentp, CDC14, NET1P);                                                                                                                                                                            //
    Dissoc_from_RENP = Mass_Action_1_222(kdirentp, RENTP);                                                                                                                                                                                                 //
    Net1_synthesis = ksnet;                                                                                                                                                                                                                                //
    Net1_degradation = Mass_Action_1_222(kdnet, NET1);                                                                                                                                                                                                     //
    Net1P_degradation = Mass_Action_1_222(kdnet, NET1P);                                                                                                                                                                                                   //
    NET1_phosphorylation = Mass_Action_1_222(Vkpnet, NET1);                                                                                                                                                                                                //
    dephosphorylation_1 = Mass_Action_1_222(Vppnet, NET1P);                                                                                                                                                                                                //
    RENT_phosphorylation = Mass_Action_1_222(Vkpnet, RENT);                                                                                                                                                                                                //
    dephosphorylation_2 = Mass_Action_1_222(Vppnet, RENTP);                                                                                                                                                                                                //
    Degradation_of_NET1_in_RENT = Mass_Action_1_222(kdnet, RENT);                                                                                                                                                                                          //
    Degradation_of_NET1P_in_RENTP = Mass_Action_1_222(kdnet, RENTP);                                                                                                                                                                                       //
    Degradation_of_CDC14_in_RENT = Mass_Action_1_222(kd14, RENT);                                                                                                                                                                                          //
    Degradation_of_CDC14_in_RENTP = Mass_Action_1_222(kd14, RENTP);                                                                                                                                                                                        //
    TEM1_activation = MichaelisMenten_220(LTE1, Jatem, 1.0, TEM1GDP);                                                                                                                                                                                      //
    inactivation_1 = MichaelisMenten_220(BUB2, Jitem, 1.0, TEM1GTP);                                                                                                                                                                                       //
    CDC15_activation = Mass_Action_1_222(ka15_p * TEM1GDP + ka15_p_p * TEM1GTP + ka15p * CDC14, CDC15i);                                                                                                                                                   //
    inactivation_2 = Mass_Action_1_222(ki15, CDC15);                                                                                                                                                                                                       //
    PPX_synthesis = ksppx;                                                                                                                                                                                                                                 //
    degradation_1 = Mass_Action_1_222(Vdppx, PPX);                                                                                                                                                                                                         //
    PDS1_synthesis = kspds_p + ks1pds_p_p * SBF + ks2pds_p_p * MCM1;                                                                                                                                                                                       //
    degradation_2 = Mass_Action_1_222(Vdpds, PDS1);                                                                                                                                                                                                        //
    Degradation_of_PDS1_in_PE = Mass_Action_1_222(Vdpds, PE);                                                                                                                                                                                              //
    Assoc_with_ESP1_to_form_PE = Mass_Action_2_221(kasesp, PDS1, ESP1);                                                                                                                                                                                    //
    Disso_from_PE = Mass_Action_1_222(kdiesp, PE);                                                                                                                                                                                                         //
    DNA_synthesis = ksori * (eorib5 * CLB5 + eorib2 * CLB2);                                                                                                                                                                                               //
    Negative_regulation_of_DNA_synthesis = Mass_Action_1_222(kdori, ORI);                                                                                                                                                                                  //
    Budding = ksbud * (ebudn2 * CLN2 + ebudn3 * CLN3 + ebudb5 * CLB5);                                                                                                                                                                                     //
    Negative_regulation_of_Cell_budding = Mass_Action_1_222(kdbud, BUD);                                                                                                                                                                                   //
    Spindle_formation = ksspn * CLB2 / (Jspn + CLB2);                                                                                                                                                                                                      //
    Spindle_disassembly = Mass_Action_1_222(kdspn, SPN);                                                                                                                                                                                                   //
    d_BUD_dt = Budding - Negative_regulation_of_Cell_budding;                                                                                                                                                                                              //
    d_C2_dt = (Assoc_of_CLB2_and_SIC1 - Dissoc_of_CLB2SIC1_complex - Phosphorylation_of_C2) + Dephosphorylation_of_C2P - Degradation_of_CLB2_in_C2;                                                                                                        //
    d_C2P_dt = Phosphorylation_of_C2 - Dephosphorylation_of_C2P - Degradation_of_SIC1_in_C2P - Degradation_of_CLB2_in_C2P;                                                                                                                                 //
    d_C5_dt = (Assoc_of_CLB5_and_SIC1 - Dissoc_of_CLB5SIC1 - Phosphorylation_of_C5) + Dephosphorylation_of_C5P - Degradation_of_CLB5_in_C5;                                                                                                                //
    d_C5P_dt = Phosphorylation_of_C5 - Dephosphorylation_of_C5P - Degradation_of_SIC1P_in_C5P_ - Degradation_of_CLB5_in_C5P;                                                                                                                               //
    d_CDC14_dt = ((CDC14_synthesis - CDC14_degradation - Assoc_with_NET1_to_form_RENT) + Dissoc_from_RENT - Assoc_with_NET1P_to_form_RENTP) + Dissoc_from_RENP + Degradation_of_NET1_in_RENT + Degradation_of_NET1P_in_RENTP;                              //
    d_CDC15_dt = CDC15_activation - inactivation_2;                                                                                                                                                                                                        //
    d_CDC20_dt = -Degradation_of_active_CDC20 + Activation_of_CDC20 - Inactivation_2;                                                                                                                                                                      //
    d_CDC20i_dt = (Synthesis_of_inactive_CDC20 - Degradation_of_inactiveCDC20 - Activation_of_CDC20) + Inactivation_2;                                                                                                                                     //
    d_CDC6_dt = (((CDC6_synthesis - Phosphorylation_of_CDC6) + Dephosphorylation_of_CDC6 - CLB2CDC6_complex_formation) + CLB2CDC6_dissociation - CLB5CDC6_complex_formation) + CLB5CDC6_dissociation + CLB2_degradation_in_F2 + CLB5_degradation_in_F5;    //
    d_CDC6P_dt = (Phosphorylation_of_CDC6 - Dephosphorylation_of_CDC6 - Degradation_of_CDC6P) + CLB2_degradation_in_F2P + CLB5_degradation_in_F5P;                                                                                                         //
    d_CDH1_dt = (CDH1_synthesis - CDH1_degradation) + CDH1i_activation - Inactivation_3;                                                                                                                                                                   //
    d_CDH1i_dt = (-CDH1i_degradation - CDH1i_activation) + Inactivation_3;                                                                                                                                                                                 //
    d_CLB2_dt = ((Synthesis_of_CLB2 - Degradation_of_CLB2 - Assoc_of_CLB2_and_SIC1) + Dissoc_of_CLB2SIC1_complex + Degradation_of_SIC1_in_C2P - CLB2CDC6_complex_formation) + CLB2CDC6_dissociation + CDC6_degradation_in_F2P;                             //
    d_CLB5_dt = ((Synthesis_of_CLB5 - Degradation_of_CLB5 - Assoc_of_CLB5_and_SIC1) + Dissoc_of_CLB5SIC1 + Degradation_of_SIC1P_in_C5P_ - CLB5CDC6_complex_formation) + CLB5CDC6_dissociation + CDC6_degradation_in_F5P;                                   //
    d_CLN2_dt = Synthesis_of_CLN2 - Degradation_of_CLN2;                                                                                                                                                                                                   //
    d_ESP1_dt = (Degradation_of_PDS1_in_PE - Assoc_with_ESP1_to_form_PE) + Disso_from_PE;                                                                                                                                                                  //
    d_F2_dt = (CLB2CDC6_complex_formation - CLB2CDC6_dissociation - F2_phosphorylation) + F2P_dephosphorylation - CLB2_degradation_in_F2;                                                                                                                  //
    d_F2P_dt = F2_phosphorylation - F2P_dephosphorylation - CDC6_degradation_in_F2P - CLB2_degradation_in_F2P;                                                                                                                                             //
    d_F5_dt = (CLB5CDC6_complex_formation - CLB5CDC6_dissociation - F5_phosphorylation) + F5P_dephosphorylation - CLB5_degradation_in_F5;                                                                                                                  //
    d_F5P_dt = F5_phosphorylation - F5P_dephosphorylation - CDC6_degradation_in_F5P - CLB5_degradation_in_F5P;                                                                                                                                             //
    d_IEP_dt = Activation_of_IEP - Inactivation_1;                                                                                                                                                                                                         //
    d_MASS_dt = Growth;                                                                                                                                                                                                                                    //
    d_NET1_dt = (-Assoc_with_NET1_to_form_RENT + Dissoc_from_RENT + Net1_synthesis - Net1_degradation - NET1_phosphorylation) + dephosphorylation_1 + Degradation_of_CDC14_in_RENT;                                                                        //
    d_NET1P_dt = ((-Assoc_with_NET1P_to_form_RENTP + Dissoc_from_RENP - Net1P_degradation) + NET1_phosphorylation - dephosphorylation_1) + Degradation_of_CDC14_in_RENTP;                                                                                  //
    d_ORI_dt = DNA_synthesis - Negative_regulation_of_DNA_synthesis;                                                                                                                                                                                       //
    d_PDS1_dt = (PDS1_synthesis - degradation_2 - Assoc_with_ESP1_to_form_PE) + Disso_from_PE;                                                                                                                                                             //
    d_PPX_dt = PPX_synthesis - degradation_1;                                                                                                                                                                                                              //
    d_RENT_dt = (Assoc_with_NET1_to_form_RENT - Dissoc_from_RENT - RENT_phosphorylation) + dephosphorylation_2 - Degradation_of_NET1_in_RENT - Degradation_of_CDC14_in_RENT;                                                                               //
    d_RENTP_dt = (Assoc_with_NET1P_to_form_RENTP - Dissoc_from_RENP) + RENT_phosphorylation - dephosphorylation_2 - Degradation_of_NET1P_in_RENTP - Degradation_of_CDC14_in_RENTP;                                                                         //
    d_SIC1_dt = (((Synthesis_of_SIC1 - Phosphorylation_of_SIC1) + Dephosphorylation_of_SIC1 - Assoc_of_CLB2_and_SIC1) + Dissoc_of_CLB2SIC1_complex - Assoc_of_CLB5_and_SIC1) + Dissoc_of_CLB5SIC1 + Degradation_of_CLB2_in_C2 + Degradation_of_CLB5_in_C5; //
    d_SIC1P_dt = (Phosphorylation_of_SIC1 - Dephosphorylation_of_SIC1 - Fast_Degradation_of_SIC1P) + Degradation_of_CLB2_in_C2P + Degradation_of_CLB5_in_C5P;                                                                                              //
    d_SPN_dt = Spindle_formation - Spindle_disassembly;                                                                                                                                                                                                    //
    d_SWI5_dt = (Synthesis_of_SWI5 - Degradation_of_SWI5) + Activation_of_SWI5 - Inactivation_of_SWI5;                                                                                                                                                     //
    d_SWI5P_dt = (-Degradation_of_SWI5P - Activation_of_SWI5) + Inactivation_of_SWI5;                                                                                                                                                                      //
    d_TEM1GTP_dt = TEM1_activation - inactivation_1;                                                                                                                                                                                                       //

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

    mParameters.push_back(BUB2);
    mParameters.push_back(LTE1);
    mParameters.push_back(MAD2);
    mParameters.push_back(b0);
    mParameters.push_back(bub2h);
    mParameters.push_back(bub2l);
    mParameters.push_back(C0);
    mParameters.push_back(CDC15T);
    mParameters.push_back(Dn3);
    mParameters.push_back(ebudb5);
    mParameters.push_back(ebudn2);
    mParameters.push_back(ebudn3);
    mParameters.push_back(ec1b2);
    mParameters.push_back(ec1b5);
    mParameters.push_back(ec1k2);
    mParameters.push_back(ec1n2);
    mParameters.push_back(ec1n3);
    mParameters.push_back(ef6b2);
    mParameters.push_back(ef6b5);
    mParameters.push_back(ef6k2);
    mParameters.push_back(ef6n2);
    mParameters.push_back(ef6n3);
    mParameters.push_back(eicdhb2);
    mParameters.push_back(eicdhb5);
    mParameters.push_back(eicdhn2);
    mParameters.push_back(eicdhn3);
    mParameters.push_back(eorib2);
    mParameters.push_back(eorib5);
    mParameters.push_back(esbfb5);
    mParameters.push_back(esbfn2);
    mParameters.push_back(esbfn3);
    mParameters.push_back(ESP1T);
    mParameters.push_back(IET);
    mParameters.push_back(J20ppx);
    mParameters.push_back(Jacdh);
    mParameters.push_back(Jaiep);
    mParameters.push_back(Jamcm);
    mParameters.push_back(Jasbf);
    mParameters.push_back(Jatem);
    mParameters.push_back(Jd2c1);
    mParameters.push_back(Jd2f6);
    mParameters.push_back(Jicdh);
    mParameters.push_back(Jiiep);
    mParameters.push_back(Jimcm);
    mParameters.push_back(Jisbf);
    mParameters.push_back(Jitem);
    mParameters.push_back(Jn3);
    mParameters.push_back(Jpds);
    mParameters.push_back(Jspn);
    mParameters.push_back(ka15_p);
    mParameters.push_back(ka15_p_p);
    mParameters.push_back(ka15p);
    mParameters.push_back(ka20_p);
    mParameters.push_back(ka20_p_p);
    mParameters.push_back(kacdh_p);
    mParameters.push_back(kacdh_p_p);
    mParameters.push_back(kaiep);
    mParameters.push_back(kamcm);
    mParameters.push_back(kasb2);
    mParameters.push_back(kasb5);
    mParameters.push_back(kasbf);
    mParameters.push_back(kasesp);
    mParameters.push_back(kasf2);
    mParameters.push_back(kasf5);
    mParameters.push_back(kasrent);
    mParameters.push_back(kasrentp);
    mParameters.push_back(kaswi);
    mParameters.push_back(kd14);
    mParameters.push_back(kd1c1);
    mParameters.push_back(kd1f6);
    mParameters.push_back(kd1pds_p);
    mParameters.push_back(kd20);
    mParameters.push_back(kd2c1);
    mParameters.push_back(kd2f6);
    mParameters.push_back(kd2pds_p_p);
    mParameters.push_back(kd3c1);
    mParameters.push_back(kd3f6);
    mParameters.push_back(kd3pds_p_p);
    mParameters.push_back(kdb2_p);
    mParameters.push_back(kdb2_p_p);
    mParameters.push_back(kdb2p);
    mParameters.push_back(kdb5_p);
    mParameters.push_back(kdb5_p_p);
    mParameters.push_back(kdbud);
    mParameters.push_back(kdcdh);
    mParameters.push_back(kdib2);
    mParameters.push_back(kdib5);
    mParameters.push_back(kdiesp);
    mParameters.push_back(kdif2);
    mParameters.push_back(kdif5);
    mParameters.push_back(kdirent);
    mParameters.push_back(kdirentp);
    mParameters.push_back(kdn2);
    mParameters.push_back(kdnet);
    mParameters.push_back(kdori);
    mParameters.push_back(kdppx_p);
    mParameters.push_back(kdppx_p_p);
    mParameters.push_back(kdspn);
    mParameters.push_back(kdswi);
    mParameters.push_back(KEZ);
    mParameters.push_back(KEZ2);
    mParameters.push_back(ki15);
    mParameters.push_back(kicdh_p);
    mParameters.push_back(kicdh_p_p);
    mParameters.push_back(kiiep);
    mParameters.push_back(kimcm);
    mParameters.push_back(kisbf_p);
    mParameters.push_back(kisbf_p_p);
    mParameters.push_back(kiswi);
    mParameters.push_back(kkpnet_p);
    mParameters.push_back(kkpnet_p_p);
    mParameters.push_back(kppc1);
    mParameters.push_back(kppf6);
    mParameters.push_back(kppnet_p);
    mParameters.push_back(kppnet_p_p);
    mParameters.push_back(ks14);
    mParameters.push_back(ks1pds_p_p);
    mParameters.push_back(ks20_p);
    mParameters.push_back(ks20_p_p);
    mParameters.push_back(ks2pds_p_p);
    mParameters.push_back(ksb2_p);
    mParameters.push_back(ksb2_p_p);
    mParameters.push_back(ksb5_p);
    mParameters.push_back(ksb5_p_p);
    mParameters.push_back(ksbud);
    mParameters.push_back(ksc1_p);
    mParameters.push_back(ksc1_p_p);
    mParameters.push_back(kscdh);
    mParameters.push_back(ksf6_p);
    mParameters.push_back(ksf6_p_p);
    mParameters.push_back(ksf6_p_p_p);
    mParameters.push_back(ksn2_p);
    mParameters.push_back(ksn2_p_p);
    mParameters.push_back(ksnet);
    mParameters.push_back(ksori);
    mParameters.push_back(kspds_p);
    mParameters.push_back(ksppx);
    mParameters.push_back(ksspn);
    mParameters.push_back(ksswi_p);
    mParameters.push_back(ksswi_p_p);
    mParameters.push_back(lte1h);
    mParameters.push_back(lte1l);
    mParameters.push_back(mad2h);
    mParameters.push_back(mad2l);
    mParameters.push_back(mdt);
    mParameters.push_back(TEM1T);
}

double Chen2004SbmlOdeSystem::ProcessModelEvents(double time, const std::vector<double>& rY)
{
    time *= TIMESCALE_MULTIPLIER; // Chaste integrates in hours; use the model's native time units
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

    //========================================
    // EVENT: reset ORI
    //========================================
    {
        double event_dist = (0.0) - (CLB2 + CLB5 - KEZ2) - std::numeric_limits<double>::epsilon();

        // active: the raw SBML trigger condition. detected: the signed distance has reached zero,
        // the point CVODE roots on. For a >=/<= trigger detected lags active by an epsilon at the
        // boundary; for a non-relational trigger the distance is a constant >= 0 so detected is
        // always true and the logic below reduces to the raw condition.
        // The fire is gated on detected so an event whose crossing lands on a sample grid point is
        // not latched by an uncommitted evaluation at the grid point before it can be applied. The
        // satisfied/re-arm state tracks active, so a trigger that is already true at the initial
        // condition (initialValue=true) stays satisfied and does not spuriously fire.
        bool active = (CLB2 + CLB5 - KEZ2) < 0.0;
        bool detected = event_dist >= 0.0;

        // Suppress an event whose trigger was already active when this Solve segment started
        // (a carried-over trigger) by forcing a large negative distance, so CVODE reports no
        // spurious root at the initial condition. mEventClampActive is frozen at segment start
        // (CalculateStoppingEvent) and cleared below the instant the trigger first goes false.
        // Using this monotonic per-segment flag rather than the live, in-step-mutated
        // mEventSatisfied keeps the root function stable across CVODE's root bracketing, so an
        // event localizes at its true crossing instead of the integration step endpoint.
        if (mEventClampActive[0] && active)
        {
            event_dist = -(std::abs(event_dist) + 1.0);
        }

        // Update max_dist (closest event to triggering)
        if (event_dist > max_dist)
        {
            max_dist = event_dist;
        }

        // Process the event
        if (active)
        {
            if (!mEventSatisfied[0] && detected)
            {
                // The condition is transitioning from false -> true: trigger the event. The
                // assignment values are recorded below (not here) so they can be refreshed at the
                // committed event point.
                mEventTriggered[0] = true;
            }
            // Latch only once the distance has crossed zero. Until then (active but not yet
            // detected, i.e. within the epsilon boundary) leave the satisfied state untouched, so
            // a crossing exactly on a grid point is neither prematurely latched nor, when the
            // trigger is already true at t=0, re-armed.
            if (detected)
            {
                mEventSatisfied[0] = true;
            }
        }
        else if (!mEventTriggered[0])
        {
            // Trigger is false and the event has not fired in this segment, so it (re-)arms:
            // clear the satisfied latch and the clamp (the clamp permanently, monotonically,
            // so it stays stable across CVODE's in-step root bracketing and the next rising
            // edge is detected). Once the event HAS fired this segment we leave these sticky,
            // so a later root-bracketing evaluation that lands on the trigger-false side cannot
            // undo the fire and leave the event spuriously unsatisfied (which would re-fire it
            // at the next segment's initial condition).
            mEventSatisfied[0] = false;
            mEventClampActive[0] = false;
        }

        // Record this event's assignments while it has fired this segment - re-evaluated on every
        // call rather than only at first detection. The assignment is still deferred (applied by
        // AdjustParameters at the committed point), but recording it here lets the harness re-run
        // ProcessModelEvents at the localized root before applying, so a state-dependent value
        // (e.g. a compartment-resize rescale S * C_old / C_new) uses the root state rather than the
        // integration step where the event was first detected. event_priority orders simultaneous
        // events: an assignment only overwrites one already recorded this firing if its event has
        // lower-or-equal priority, so the lowest-priority event - which SBML executes last - wins a
        // conflict. Events with no priority use +inf, reducing to last-writer-wins.
        if (mEventTriggered[0])
        {
            [[maybe_unused]] double event_priority = std::numeric_limits<double>::max();
            // ORI = 0.0
            if (!mEventAdjustedStateVars[25]
                || event_priority <= mEventAdjustedStatePriority[25])
            {
                mEventAdjustedStateVars[25] = true;
                mEventAdjustedStateValues[25] = 0.0;
                mEventAdjustedStatePriority[25] = event_priority;
            }
        }
    }

    //========================================
    // EVENT: start DNA synthesis
    //========================================
    {
        double event_dist = (ORI - 1.0) - (0.0) - std::numeric_limits<double>::epsilon();

        // active: the raw SBML trigger condition. detected: the signed distance has reached zero,
        // the point CVODE roots on. For a >=/<= trigger detected lags active by an epsilon at the
        // boundary; for a non-relational trigger the distance is a constant >= 0 so detected is
        // always true and the logic below reduces to the raw condition.
        // The fire is gated on detected so an event whose crossing lands on a sample grid point is
        // not latched by an uncommitted evaluation at the grid point before it can be applied. The
        // satisfied/re-arm state tracks active, so a trigger that is already true at the initial
        // condition (initialValue=true) stays satisfied and does not spuriously fire.
        bool active = (ORI - 1.0) > 0.0;
        bool detected = event_dist >= 0.0;

        // Suppress an event whose trigger was already active when this Solve segment started
        // (a carried-over trigger) by forcing a large negative distance, so CVODE reports no
        // spurious root at the initial condition. mEventClampActive is frozen at segment start
        // (CalculateStoppingEvent) and cleared below the instant the trigger first goes false.
        // Using this monotonic per-segment flag rather than the live, in-step-mutated
        // mEventSatisfied keeps the root function stable across CVODE's root bracketing, so an
        // event localizes at its true crossing instead of the integration step endpoint.
        if (mEventClampActive[1] && active)
        {
            event_dist = -(std::abs(event_dist) + 1.0);
        }

        // Update max_dist (closest event to triggering)
        if (event_dist > max_dist)
        {
            max_dist = event_dist;
        }

        // Process the event
        if (active)
        {
            if (!mEventSatisfied[1] && detected)
            {
                // The condition is transitioning from false -> true: trigger the event. The
                // assignment values are recorded below (not here) so they can be refreshed at the
                // committed event point.
                mEventTriggered[1] = true;
            }
            // Latch only once the distance has crossed zero. Until then (active but not yet
            // detected, i.e. within the epsilon boundary) leave the satisfied state untouched, so
            // a crossing exactly on a grid point is neither prematurely latched nor, when the
            // trigger is already true at t=0, re-armed.
            if (detected)
            {
                mEventSatisfied[1] = true;
            }
        }
        else if (!mEventTriggered[1])
        {
            // Trigger is false and the event has not fired in this segment, so it (re-)arms:
            // clear the satisfied latch and the clamp (the clamp permanently, monotonically,
            // so it stays stable across CVODE's in-step root bracketing and the next rising
            // edge is detected). Once the event HAS fired this segment we leave these sticky,
            // so a later root-bracketing evaluation that lands on the trigger-false side cannot
            // undo the fire and leave the event spuriously unsatisfied (which would re-fire it
            // at the next segment's initial condition).
            mEventSatisfied[1] = false;
            mEventClampActive[1] = false;
        }

        // Record this event's assignments while it has fired this segment - re-evaluated on every
        // call rather than only at first detection. The assignment is still deferred (applied by
        // AdjustParameters at the committed point), but recording it here lets the harness re-run
        // ProcessModelEvents at the localized root before applying, so a state-dependent value
        // (e.g. a compartment-resize rescale S * C_old / C_new) uses the root state rather than the
        // integration step where the event was first detected. event_priority orders simultaneous
        // events: an assignment only overwrites one already recorded this firing if its event has
        // lower-or-equal priority, so the lowest-priority event - which SBML executes last - wins a
        // conflict. Events with no priority use +inf, reducing to last-writer-wins.
        if (mEventTriggered[1])
        {
            [[maybe_unused]] double event_priority = std::numeric_limits<double>::max();
            // MAD2 = mad2h
            if (!mEventAdjustedParameters[2]
                || event_priority <= mEventAdjustedParameterPriority[2])
            {
                mEventAdjustedParameters[2] = true;
                mEventAdjustedParameterValues[2] = mad2h;
                mEventAdjustedParameterPriority[2] = event_priority;
            }

            // BUB2 = bub2h
            if (!mEventAdjustedParameters[0]
                || event_priority <= mEventAdjustedParameterPriority[0])
            {
                mEventAdjustedParameters[0] = true;
                mEventAdjustedParameterValues[0] = bub2h;
                mEventAdjustedParameterPriority[0] = event_priority;
            }
        }
    }

    //========================================
    // EVENT: spindle checkpoint
    //========================================
    {
        double event_dist = (SPN - 1.0) - (0.0) - std::numeric_limits<double>::epsilon();

        // active: the raw SBML trigger condition. detected: the signed distance has reached zero,
        // the point CVODE roots on. For a >=/<= trigger detected lags active by an epsilon at the
        // boundary; for a non-relational trigger the distance is a constant >= 0 so detected is
        // always true and the logic below reduces to the raw condition.
        // The fire is gated on detected so an event whose crossing lands on a sample grid point is
        // not latched by an uncommitted evaluation at the grid point before it can be applied. The
        // satisfied/re-arm state tracks active, so a trigger that is already true at the initial
        // condition (initialValue=true) stays satisfied and does not spuriously fire.
        bool active = (SPN - 1.0) > 0.0;
        bool detected = event_dist >= 0.0;

        // Suppress an event whose trigger was already active when this Solve segment started
        // (a carried-over trigger) by forcing a large negative distance, so CVODE reports no
        // spurious root at the initial condition. mEventClampActive is frozen at segment start
        // (CalculateStoppingEvent) and cleared below the instant the trigger first goes false.
        // Using this monotonic per-segment flag rather than the live, in-step-mutated
        // mEventSatisfied keeps the root function stable across CVODE's root bracketing, so an
        // event localizes at its true crossing instead of the integration step endpoint.
        if (mEventClampActive[2] && active)
        {
            event_dist = -(std::abs(event_dist) + 1.0);
        }

        // Update max_dist (closest event to triggering)
        if (event_dist > max_dist)
        {
            max_dist = event_dist;
        }

        // Process the event
        if (active)
        {
            if (!mEventSatisfied[2] && detected)
            {
                // The condition is transitioning from false -> true: trigger the event. The
                // assignment values are recorded below (not here) so they can be refreshed at the
                // committed event point.
                mEventTriggered[2] = true;
            }
            // Latch only once the distance has crossed zero. Until then (active but not yet
            // detected, i.e. within the epsilon boundary) leave the satisfied state untouched, so
            // a crossing exactly on a grid point is neither prematurely latched nor, when the
            // trigger is already true at t=0, re-armed.
            if (detected)
            {
                mEventSatisfied[2] = true;
            }
        }
        else if (!mEventTriggered[2])
        {
            // Trigger is false and the event has not fired in this segment, so it (re-)arms:
            // clear the satisfied latch and the clamp (the clamp permanently, monotonically,
            // so it stays stable across CVODE's in-step root bracketing and the next rising
            // edge is detected). Once the event HAS fired this segment we leave these sticky,
            // so a later root-bracketing evaluation that lands on the trigger-false side cannot
            // undo the fire and leave the event spuriously unsatisfied (which would re-fire it
            // at the next segment's initial condition).
            mEventSatisfied[2] = false;
            mEventClampActive[2] = false;
        }

        // Record this event's assignments while it has fired this segment - re-evaluated on every
        // call rather than only at first detection. The assignment is still deferred (applied by
        // AdjustParameters at the committed point), but recording it here lets the harness re-run
        // ProcessModelEvents at the localized root before applying, so a state-dependent value
        // (e.g. a compartment-resize rescale S * C_old / C_new) uses the root state rather than the
        // integration step where the event was first detected. event_priority orders simultaneous
        // events: an assignment only overwrites one already recorded this firing if its event has
        // lower-or-equal priority, so the lowest-priority event - which SBML executes last - wins a
        // conflict. Events with no priority use +inf, reducing to last-writer-wins.
        if (mEventTriggered[2])
        {
            [[maybe_unused]] double event_priority = std::numeric_limits<double>::max();
            // MAD2 = mad2l
            if (!mEventAdjustedParameters[2]
                || event_priority <= mEventAdjustedParameterPriority[2])
            {
                mEventAdjustedParameters[2] = true;
                mEventAdjustedParameterValues[2] = mad2l;
                mEventAdjustedParameterPriority[2] = event_priority;
            }

            // LTE1 = lte1h
            if (!mEventAdjustedParameters[1]
                || event_priority <= mEventAdjustedParameterPriority[1])
            {
                mEventAdjustedParameters[1] = true;
                mEventAdjustedParameterValues[1] = lte1h;
                mEventAdjustedParameterPriority[1] = event_priority;
            }

            // BUB2 = bub2l
            if (!mEventAdjustedParameters[0]
                || event_priority <= mEventAdjustedParameterPriority[0])
            {
                mEventAdjustedParameters[0] = true;
                mEventAdjustedParameterValues[0] = bub2l;
                mEventAdjustedParameterPriority[0] = event_priority;
            }
        }
    }

    //========================================
    // EVENT: cell division
    //========================================
    {
        double event_dist = (0.0) - (CLB2 - KEZ) - std::numeric_limits<double>::epsilon();

        // active: the raw SBML trigger condition. detected: the signed distance has reached zero,
        // the point CVODE roots on. For a >=/<= trigger detected lags active by an epsilon at the
        // boundary; for a non-relational trigger the distance is a constant >= 0 so detected is
        // always true and the logic below reduces to the raw condition.
        // The fire is gated on detected so an event whose crossing lands on a sample grid point is
        // not latched by an uncommitted evaluation at the grid point before it can be applied. The
        // satisfied/re-arm state tracks active, so a trigger that is already true at the initial
        // condition (initialValue=true) stays satisfied and does not spuriously fire.
        bool active = (CLB2 - KEZ) < 0.0;
        bool detected = event_dist >= 0.0;

        // Suppress an event whose trigger was already active when this Solve segment started
        // (a carried-over trigger) by forcing a large negative distance, so CVODE reports no
        // spurious root at the initial condition. mEventClampActive is frozen at segment start
        // (CalculateStoppingEvent) and cleared below the instant the trigger first goes false.
        // Using this monotonic per-segment flag rather than the live, in-step-mutated
        // mEventSatisfied keeps the root function stable across CVODE's root bracketing, so an
        // event localizes at its true crossing instead of the integration step endpoint.
        if (mEventClampActive[3] && active)
        {
            event_dist = -(std::abs(event_dist) + 1.0);
        }

        // Update max_dist (closest event to triggering)
        if (event_dist > max_dist)
        {
            max_dist = event_dist;
        }

        // Process the event
        if (active)
        {
            if (!mEventSatisfied[3] && detected)
            {
                // The condition is transitioning from false -> true: trigger the event. The
                // assignment values are recorded below (not here) so they can be refreshed at the
                // committed event point.
                mEventTriggered[3] = true;
            }
            // Latch only once the distance has crossed zero. Until then (active but not yet
            // detected, i.e. within the epsilon boundary) leave the satisfied state untouched, so
            // a crossing exactly on a grid point is neither prematurely latched nor, when the
            // trigger is already true at t=0, re-armed.
            if (detected)
            {
                mEventSatisfied[3] = true;
            }
        }
        else if (!mEventTriggered[3])
        {
            // Trigger is false and the event has not fired in this segment, so it (re-)arms:
            // clear the satisfied latch and the clamp (the clamp permanently, monotonically,
            // so it stays stable across CVODE's in-step root bracketing and the next rising
            // edge is detected). Once the event HAS fired this segment we leave these sticky,
            // so a later root-bracketing evaluation that lands on the trigger-false side cannot
            // undo the fire and leave the event spuriously unsatisfied (which would re-fire it
            // at the next segment's initial condition).
            mEventSatisfied[3] = false;
            mEventClampActive[3] = false;
        }

        // Record this event's assignments while it has fired this segment - re-evaluated on every
        // call rather than only at first detection. The assignment is still deferred (applied by
        // AdjustParameters at the committed point), but recording it here lets the harness re-run
        // ProcessModelEvents at the localized root before applying, so a state-dependent value
        // (e.g. a compartment-resize rescale S * C_old / C_new) uses the root state rather than the
        // integration step where the event was first detected. event_priority orders simultaneous
        // events: an assignment only overwrites one already recorded this firing if its event has
        // lower-or-equal priority, so the lowest-priority event - which SBML executes last - wins a
        // conflict. Events with no priority use +inf, reducing to last-writer-wins.
        if (mEventTriggered[3])
        {
            [[maybe_unused]] double event_priority = std::numeric_limits<double>::max();
            // MASS = F * MASS
            if (!mEventAdjustedStateVars[22]
                || event_priority <= mEventAdjustedStatePriority[22])
            {
                mEventAdjustedStateVars[22] = true;
                mEventAdjustedStateValues[22] = F * MASS;
                mEventAdjustedStatePriority[22] = event_priority;
            }

            // LTE1 = lte1l
            if (!mEventAdjustedParameters[1]
                || event_priority <= mEventAdjustedParameterPriority[1])
            {
                mEventAdjustedParameters[1] = true;
                mEventAdjustedParameterValues[1] = lte1l;
                mEventAdjustedParameterPriority[1] = event_priority;
            }

            // BUD = 0.0
            if (!mEventAdjustedStateVars[0]
                || event_priority <= mEventAdjustedStatePriority[0])
            {
                mEventAdjustedStateVars[0] = true;
                mEventAdjustedStateValues[0] = 0.0;
                mEventAdjustedStatePriority[0] = event_priority;
            }

            // SPN = 0.0
            if (!mEventAdjustedStateVars[32]
                || event_priority <= mEventAdjustedStatePriority[32])
            {
                mEventAdjustedStateVars[32] = true;
                mEventAdjustedStateValues[32] = 0.0;
                mEventAdjustedStatePriority[32] = event_priority;
            }
        }
    }

    return max_dist; // Signed distance of the event closest to triggering
}

std::vector<double> Chen2004SbmlOdeSystem::RunModelEquations(double time, const std::vector<double>& rStateVariables)
{
    BUD = rStateVariables[0];
    C2 = rStateVariables[1];
    C2P = rStateVariables[2];
    C5 = rStateVariables[3];
    C5P = rStateVariables[4];
    CDC14 = rStateVariables[5];
    CDC15 = rStateVariables[6];
    CDC20 = rStateVariables[7];
    CDC20i = rStateVariables[8];
    CDC6 = rStateVariables[9];
    CDC6P = rStateVariables[10];
    CDH1 = rStateVariables[11];
    CDH1i = rStateVariables[12];
    CLB2 = rStateVariables[13];
    CLB5 = rStateVariables[14];
    CLN2 = rStateVariables[15];
    ESP1 = rStateVariables[16];
    F2 = rStateVariables[17];
    F2P = rStateVariables[18];
    F5 = rStateVariables[19];
    F5P = rStateVariables[20];
    IEP = rStateVariables[21];
    MASS = rStateVariables[22];
    NET1 = rStateVariables[23];
    NET1P = rStateVariables[24];
    ORI = rStateVariables[25];
    PDS1 = rStateVariables[26];
    PPX = rStateVariables[27];
    RENT = rStateVariables[28];
    RENTP = rStateVariables[29];
    SIC1 = rStateVariables[30];
    SIC1P = rStateVariables[31];
    SPN = rStateVariables[32];
    SWI5 = rStateVariables[33];
    SWI5P = rStateVariables[34];
    TEM1GTP = rStateVariables[35];

    BUB2 = GetParameter(0);
    LTE1 = GetParameter(1);
    MAD2 = GetParameter(2);
    b0 = GetParameter(3);
    bub2h = GetParameter(4);
    bub2l = GetParameter(5);
    C0 = GetParameter(6);
    CDC15T = GetParameter(7);
    Dn3 = GetParameter(8);
    ebudb5 = GetParameter(9);
    ebudn2 = GetParameter(10);
    ebudn3 = GetParameter(11);
    ec1b2 = GetParameter(12);
    ec1b5 = GetParameter(13);
    ec1k2 = GetParameter(14);
    ec1n2 = GetParameter(15);
    ec1n3 = GetParameter(16);
    ef6b2 = GetParameter(17);
    ef6b5 = GetParameter(18);
    ef6k2 = GetParameter(19);
    ef6n2 = GetParameter(20);
    ef6n3 = GetParameter(21);
    eicdhb2 = GetParameter(22);
    eicdhb5 = GetParameter(23);
    eicdhn2 = GetParameter(24);
    eicdhn3 = GetParameter(25);
    eorib2 = GetParameter(26);
    eorib5 = GetParameter(27);
    esbfb5 = GetParameter(28);
    esbfn2 = GetParameter(29);
    esbfn3 = GetParameter(30);
    ESP1T = GetParameter(31);
    IET = GetParameter(32);
    J20ppx = GetParameter(33);
    Jacdh = GetParameter(34);
    Jaiep = GetParameter(35);
    Jamcm = GetParameter(36);
    Jasbf = GetParameter(37);
    Jatem = GetParameter(38);
    Jd2c1 = GetParameter(39);
    Jd2f6 = GetParameter(40);
    Jicdh = GetParameter(41);
    Jiiep = GetParameter(42);
    Jimcm = GetParameter(43);
    Jisbf = GetParameter(44);
    Jitem = GetParameter(45);
    Jn3 = GetParameter(46);
    Jpds = GetParameter(47);
    Jspn = GetParameter(48);
    ka15_p = GetParameter(49);
    ka15_p_p = GetParameter(50);
    ka15p = GetParameter(51);
    ka20_p = GetParameter(52);
    ka20_p_p = GetParameter(53);
    kacdh_p = GetParameter(54);
    kacdh_p_p = GetParameter(55);
    kaiep = GetParameter(56);
    kamcm = GetParameter(57);
    kasb2 = GetParameter(58);
    kasb5 = GetParameter(59);
    kasbf = GetParameter(60);
    kasesp = GetParameter(61);
    kasf2 = GetParameter(62);
    kasf5 = GetParameter(63);
    kasrent = GetParameter(64);
    kasrentp = GetParameter(65);
    kaswi = GetParameter(66);
    kd14 = GetParameter(67);
    kd1c1 = GetParameter(68);
    kd1f6 = GetParameter(69);
    kd1pds_p = GetParameter(70);
    kd20 = GetParameter(71);
    kd2c1 = GetParameter(72);
    kd2f6 = GetParameter(73);
    kd2pds_p_p = GetParameter(74);
    kd3c1 = GetParameter(75);
    kd3f6 = GetParameter(76);
    kd3pds_p_p = GetParameter(77);
    kdb2_p = GetParameter(78);
    kdb2_p_p = GetParameter(79);
    kdb2p = GetParameter(80);
    kdb5_p = GetParameter(81);
    kdb5_p_p = GetParameter(82);
    kdbud = GetParameter(83);
    kdcdh = GetParameter(84);
    kdib2 = GetParameter(85);
    kdib5 = GetParameter(86);
    kdiesp = GetParameter(87);
    kdif2 = GetParameter(88);
    kdif5 = GetParameter(89);
    kdirent = GetParameter(90);
    kdirentp = GetParameter(91);
    kdn2 = GetParameter(92);
    kdnet = GetParameter(93);
    kdori = GetParameter(94);
    kdppx_p = GetParameter(95);
    kdppx_p_p = GetParameter(96);
    kdspn = GetParameter(97);
    kdswi = GetParameter(98);
    KEZ = GetParameter(99);
    KEZ2 = GetParameter(100);
    ki15 = GetParameter(101);
    kicdh_p = GetParameter(102);
    kicdh_p_p = GetParameter(103);
    kiiep = GetParameter(104);
    kimcm = GetParameter(105);
    kisbf_p = GetParameter(106);
    kisbf_p_p = GetParameter(107);
    kiswi = GetParameter(108);
    kkpnet_p = GetParameter(109);
    kkpnet_p_p = GetParameter(110);
    kppc1 = GetParameter(111);
    kppf6 = GetParameter(112);
    kppnet_p = GetParameter(113);
    kppnet_p_p = GetParameter(114);
    ks14 = GetParameter(115);
    ks1pds_p_p = GetParameter(116);
    ks20_p = GetParameter(117);
    ks20_p_p = GetParameter(118);
    ks2pds_p_p = GetParameter(119);
    ksb2_p = GetParameter(120);
    ksb2_p_p = GetParameter(121);
    ksb5_p = GetParameter(122);
    ksb5_p_p = GetParameter(123);
    ksbud = GetParameter(124);
    ksc1_p = GetParameter(125);
    ksc1_p_p = GetParameter(126);
    kscdh = GetParameter(127);
    ksf6_p = GetParameter(128);
    ksf6_p_p = GetParameter(129);
    ksf6_p_p_p = GetParameter(130);
    ksn2_p = GetParameter(131);
    ksn2_p_p = GetParameter(132);
    ksnet = GetParameter(133);
    ksori = GetParameter(134);
    kspds_p = GetParameter(135);
    ksppx = GetParameter(136);
    ksspn = GetParameter(137);
    ksswi_p = GetParameter(138);
    ksswi_p_p = GetParameter(139);
    lte1h = GetParameter(140);
    lte1l = GetParameter(141);
    mad2h = GetParameter(142);
    mad2l = GetParameter(143);
    mdt = GetParameter(144);
    TEM1T = GetParameter(145);

    BCK2 = b0 * MASS;                                                                                  //
    CDC14T = CDC14 + RENT + RENTP;                                                                     //
    CDC15i = CDC15T - CDC15;                                                                           //
    CDC6T = CDC6 + F2 + F5 + CDC6P + F2P + F5P;                                                        //
    SIC1T = SIC1 + C2 + C5 + SIC1P + C2P + C5P;                                                        //
    CKIT = SIC1T + CDC6T;                                                                              //
    CLB2T = CLB2 + C2 + C2P + F2 + F2P;                                                                //
    CLB5T = CLB5 + C5 + C5P + F5 + F5P;                                                                //
    CLN3 = C0 * Dn3 * MASS / (Jn3 + Dn3 * MASS);                                                       //
    IE = IET - IEP;                                                                                    //
    MCM1 = GK_219(kamcm * CLB2, kimcm, Jamcm, Jimcm);                                                  //
    NET1T = NET1 + NET1P + RENT + RENTP;                                                               //
    PE = ESP1T - ESP1;                                                                                 //
    Vasbf = kasbf * (esbfn2 * CLN2 + esbfn3 * (CLN3 + BCK2) + esbfb5 * CLB5);                          //
    Visbf = kisbf_p + kisbf_p_p * CLB2;                                                                //
    SBF = GK_219(Vasbf, Visbf, Jasbf, Jisbf);                                                          //
    TEM1GDP = TEM1T - TEM1GTP;                                                                         //
    mu = std::log(2.0) / mdt;                                                                          //
    D = 1.026 / mu - 32.0;                                                                             //
    Vdb5 = kdb5_p + kdb5_p_p * CDC20;                                                                  //
    Vdb2 = kdb2_p + kdb2_p_p * CDH1 + kdb2p * CDC20;                                                   //
    Vd2c1 = kd2c1 * (ec1n3 * CLN3 + ec1k2 * BCK2 + ec1n2 * CLN2 + ec1b5 * CLB5 + ec1b2 * CLB2);        //
    Vkpc1 = kd1c1 + Vd2c1 / (Jd2c1 + SIC1 + C2 + C5 + SIC1P + C2P + C5P);                              //
    Vd2f6 = kd2f6 * (ef6n3 * CLN3 + ef6k2 * BCK2 + ef6n2 * CLN2 + ef6b5 * CLB5 + ef6b2 * CLB2);        //
    Vkpf6 = kd1f6 + Vd2f6 / (Jd2f6 + CDC6 + F2 + F5 + CDC6P + F2P + F5P);                              //
    Vacdh = kacdh_p + kacdh_p_p * CDC14;                                                               //
    Vicdh = kicdh_p + kicdh_p_p * (eicdhn3 * CLN3 + eicdhn2 * CLN2 + eicdhb5 * CLB5 + eicdhb2 * CLB2); //
    Vppnet = kppnet_p + kppnet_p_p * PPX;                                                              //
    Vkpnet = (kkpnet_p + kkpnet_p_p * CDC15) * MASS;                                                   //
    Vdppx = kdppx_p + kdppx_p_p * (J20ppx + CDC20) * Jpds / (Jpds + PDS1);                             //
    Vdpds = kd1pds_p + kd2pds_p_p * CDC20 + kd3pds_p_p * CDH1;                                         //
    Vaiep = kaiep * CLB2;                                                                              //
    Vppc1 = kppc1 * CDC14;                                                                             //
    Vppf6 = kppf6 * CDC14;                                                                             //
    F = std::exp(-mu * D);                                                                             //
    Growth = mu * MASS;                                                                                //
    Synthesis_of_CLN2 = (ksn2_p + ksn2_p_p * SBF) * MASS;                                              //
    Degradation_of_CLN2 = Mass_Action_1_222(kdn2, CLN2);                                               //
    Synthesis_of_CLB2 = (ksb2_p + ksb2_p_p * MCM1) * MASS;                                             //
    Degradation_of_CLB2 = Mass_Action_1_222(Vdb2, CLB2);                                               //
    Synthesis_of_CLB5 = (ksb5_p + ksb5_p_p * SBF) * MASS;                                              //
    Degradation_of_CLB5 = Mass_Action_1_222(Vdb5, CLB5);                                               //
    Synthesis_of_SIC1 = ksc1_p + ksc1_p_p * SWI5;                                                      //
    Phosphorylation_of_SIC1 = Mass_Action_1_222(Vkpc1, SIC1);                                          //
    Dephosphorylation_of_SIC1 = Mass_Action_1_222(Vppc1, SIC1P);                                       //
    Fast_Degradation_of_SIC1P = Mass_Action_1_222(kd3c1, SIC1P);                                       //
    Assoc_of_CLB2_and_SIC1 = Mass_Action_2_221(kasb2, CLB2, SIC1);                                     //
    Dissoc_of_CLB2SIC1_complex = Mass_Action_1_222(kdib2, C2);                                         //
    Assoc_of_CLB5_and_SIC1 = Mass_Action_2_221(kasb5, CLB5, SIC1);                                     //
    Dissoc_of_CLB5SIC1 = Mass_Action_1_222(kdib5, C5);                                                 //
    Phosphorylation_of_C2 = Mass_Action_1_222(Vkpc1, C2);                                              //
    Dephosphorylation_of_C2P = Mass_Action_1_222(Vppc1, C2P);                                          //
    Phosphorylation_of_C5 = Mass_Action_1_222(Vkpc1, C5);                                              //
    Dephosphorylation_of_C5P = Mass_Action_1_222(Vppc1, C5P);                                          //
    Degradation_of_CLB2_in_C2 = Mass_Action_1_222(Vdb2, C2);                                           //
    Degradation_of_CLB5_in_C5 = Mass_Action_1_222(Vdb5, C5);                                           //
    Degradation_of_SIC1_in_C2P = Mass_Action_1_222(kd3c1, C2P);                                        //
    Degradation_of_SIC1P_in_C5P_ = Mass_Action_1_222(kd3c1, C5P);                                      //
    Degradation_of_CLB2_in_C2P = Mass_Action_1_222(Vdb2, C2P);                                         //
    Degradation_of_CLB5_in_C5P = Mass_Action_1_222(Vdb5, C5P);                                         //
    CDC6_synthesis = ksf6_p + ksf6_p_p * SWI5 + ksf6_p_p_p * SBF;                                      //
    Phosphorylation_of_CDC6 = Mass_Action_1_222(Vkpf6, CDC6);                                          //
    Dephosphorylation_of_CDC6 = Mass_Action_1_222(Vppf6, CDC6P);                                       //
    Degradation_of_CDC6P = Mass_Action_1_222(kd3f6, CDC6P);                                            //
    CLB2CDC6_complex_formation = Mass_Action_2_221(kasf2, CLB2, CDC6);                                 //
    CLB2CDC6_dissociation = Mass_Action_1_222(kdif2, F2);                                              //
    CLB5CDC6_complex_formation = Mass_Action_2_221(kasf5, CLB5, CDC6);                                 //
    CLB5CDC6_dissociation = Mass_Action_1_222(kdif5, F5);                                              //
    F2_phosphorylation = Mass_Action_1_222(Vkpf6, F2);                                                 //
    F2P_dephosphorylation = Mass_Action_1_222(Vppf6, F2P);                                             //
    F5_phosphorylation = Mass_Action_1_222(Vkpf6, F5);                                                 //
    F5P_dephosphorylation = Mass_Action_1_222(Vppf6, F5P);                                             //
    CLB2_degradation_in_F2 = Mass_Action_1_222(Vdb2, F2);                                              //
    CLB5_degradation_in_F5 = Mass_Action_1_222(Vdb5, F5);                                              //
    CDC6_degradation_in_F2P = Mass_Action_1_222(kd3f6, F2P);                                           //
    CDC6_degradation_in_F5P = Mass_Action_1_222(kd3f6, F5P);                                           //
    CLB2_degradation_in_F2P = Mass_Action_1_222(Vdb2, F2P);                                            //
    CLB5_degradation_in_F5P = Mass_Action_1_222(Vdb5, F5P);                                            //
    Synthesis_of_SWI5 = ksswi_p + ksswi_p_p * MCM1;                                                    //
    Degradation_of_SWI5 = Mass_Action_1_222(kdswi, SWI5);                                              //
    Degradation_of_SWI5P = Mass_Action_1_222(kdswi, SWI5P);                                            //
    Activation_of_SWI5 = Mass_Action_1_222(kaswi * CDC14, SWI5P);                                      //
    Inactivation_of_SWI5 = Mass_Action_1_222(kiswi * CLB2, SWI5);                                      //
    Activation_of_IEP = MichaelisMenten_220(Vaiep, Jaiep, 1.0, IE);                                    //
    Inactivation_1 = MichaelisMenten_220(1.0, Jiiep, kiiep, IEP);                                      //
    Synthesis_of_inactive_CDC20 = ks20_p + ks20_p_p * MCM1;                                            //
    Degradation_of_inactiveCDC20 = Mass_Action_1_222(kd20, CDC20i);                                    //
    Degradation_of_active_CDC20 = Mass_Action_1_222(kd20, CDC20);                                      //
    Activation_of_CDC20 = Mass_Action_1_222(ka20_p + ka20_p_p * IEP, CDC20i);                          //
    // Inactivation_2:
    {
        [[maybe_unused]] double k = 1.0;
        // Qualify with this-> so the assignment targets the reaction member even when a local
        // parameter shadows its name (an SBML local parameter may shadow the reaction ID, e.g. a
        // reaction J1 with a local parameter also named J1).
        this->Inactivation_2 = k * Mass_Action_1_222(MAD2, CDC20);
    }
    CDH1_synthesis = kscdh;                                                                                                                                                                                                                                //
    CDH1_degradation = Mass_Action_1_222(kdcdh, CDH1);                                                                                                                                                                                                     //
    CDH1i_degradation = Mass_Action_1_222(kdcdh, CDH1i);                                                                                                                                                                                                   //
    CDH1i_activation = MichaelisMenten_220(Vacdh, Jacdh, 1.0, CDH1i);                                                                                                                                                                                      //
    Inactivation_3 = MichaelisMenten_220(Vicdh, Jicdh, 1.0, CDH1);                                                                                                                                                                                         //
    CDC14_synthesis = ks14;                                                                                                                                                                                                                                //
    CDC14_degradation = Mass_Action_1_222(kd14, CDC14);                                                                                                                                                                                                    //
    Assoc_with_NET1_to_form_RENT = Mass_Action_2_221(kasrent, CDC14, NET1);                                                                                                                                                                                //
    Dissoc_from_RENT = Mass_Action_1_222(kdirent, RENT);                                                                                                                                                                                                   //
    Assoc_with_NET1P_to_form_RENTP = Mass_Action_2_221(kasrentp, CDC14, NET1P);                                                                                                                                                                            //
    Dissoc_from_RENP = Mass_Action_1_222(kdirentp, RENTP);                                                                                                                                                                                                 //
    Net1_synthesis = ksnet;                                                                                                                                                                                                                                //
    Net1_degradation = Mass_Action_1_222(kdnet, NET1);                                                                                                                                                                                                     //
    Net1P_degradation = Mass_Action_1_222(kdnet, NET1P);                                                                                                                                                                                                   //
    NET1_phosphorylation = Mass_Action_1_222(Vkpnet, NET1);                                                                                                                                                                                                //
    dephosphorylation_1 = Mass_Action_1_222(Vppnet, NET1P);                                                                                                                                                                                                //
    RENT_phosphorylation = Mass_Action_1_222(Vkpnet, RENT);                                                                                                                                                                                                //
    dephosphorylation_2 = Mass_Action_1_222(Vppnet, RENTP);                                                                                                                                                                                                //
    Degradation_of_NET1_in_RENT = Mass_Action_1_222(kdnet, RENT);                                                                                                                                                                                          //
    Degradation_of_NET1P_in_RENTP = Mass_Action_1_222(kdnet, RENTP);                                                                                                                                                                                       //
    Degradation_of_CDC14_in_RENT = Mass_Action_1_222(kd14, RENT);                                                                                                                                                                                          //
    Degradation_of_CDC14_in_RENTP = Mass_Action_1_222(kd14, RENTP);                                                                                                                                                                                        //
    TEM1_activation = MichaelisMenten_220(LTE1, Jatem, 1.0, TEM1GDP);                                                                                                                                                                                      //
    inactivation_1 = MichaelisMenten_220(BUB2, Jitem, 1.0, TEM1GTP);                                                                                                                                                                                       //
    CDC15_activation = Mass_Action_1_222(ka15_p * TEM1GDP + ka15_p_p * TEM1GTP + ka15p * CDC14, CDC15i);                                                                                                                                                   //
    inactivation_2 = Mass_Action_1_222(ki15, CDC15);                                                                                                                                                                                                       //
    PPX_synthesis = ksppx;                                                                                                                                                                                                                                 //
    degradation_1 = Mass_Action_1_222(Vdppx, PPX);                                                                                                                                                                                                         //
    PDS1_synthesis = kspds_p + ks1pds_p_p * SBF + ks2pds_p_p * MCM1;                                                                                                                                                                                       //
    degradation_2 = Mass_Action_1_222(Vdpds, PDS1);                                                                                                                                                                                                        //
    Degradation_of_PDS1_in_PE = Mass_Action_1_222(Vdpds, PE);                                                                                                                                                                                              //
    Assoc_with_ESP1_to_form_PE = Mass_Action_2_221(kasesp, PDS1, ESP1);                                                                                                                                                                                    //
    Disso_from_PE = Mass_Action_1_222(kdiesp, PE);                                                                                                                                                                                                         //
    DNA_synthesis = ksori * (eorib5 * CLB5 + eorib2 * CLB2);                                                                                                                                                                                               //
    Negative_regulation_of_DNA_synthesis = Mass_Action_1_222(kdori, ORI);                                                                                                                                                                                  //
    Budding = ksbud * (ebudn2 * CLN2 + ebudn3 * CLN3 + ebudb5 * CLB5);                                                                                                                                                                                     //
    Negative_regulation_of_Cell_budding = Mass_Action_1_222(kdbud, BUD);                                                                                                                                                                                   //
    Spindle_formation = ksspn * CLB2 / (Jspn + CLB2);                                                                                                                                                                                                      //
    Spindle_disassembly = Mass_Action_1_222(kdspn, SPN);                                                                                                                                                                                                   //
    d_BUD_dt = Budding - Negative_regulation_of_Cell_budding;                                                                                                                                                                                              //
    d_C2_dt = (Assoc_of_CLB2_and_SIC1 - Dissoc_of_CLB2SIC1_complex - Phosphorylation_of_C2) + Dephosphorylation_of_C2P - Degradation_of_CLB2_in_C2;                                                                                                        //
    d_C2P_dt = Phosphorylation_of_C2 - Dephosphorylation_of_C2P - Degradation_of_SIC1_in_C2P - Degradation_of_CLB2_in_C2P;                                                                                                                                 //
    d_C5_dt = (Assoc_of_CLB5_and_SIC1 - Dissoc_of_CLB5SIC1 - Phosphorylation_of_C5) + Dephosphorylation_of_C5P - Degradation_of_CLB5_in_C5;                                                                                                                //
    d_C5P_dt = Phosphorylation_of_C5 - Dephosphorylation_of_C5P - Degradation_of_SIC1P_in_C5P_ - Degradation_of_CLB5_in_C5P;                                                                                                                               //
    d_CDC14_dt = ((CDC14_synthesis - CDC14_degradation - Assoc_with_NET1_to_form_RENT) + Dissoc_from_RENT - Assoc_with_NET1P_to_form_RENTP) + Dissoc_from_RENP + Degradation_of_NET1_in_RENT + Degradation_of_NET1P_in_RENTP;                              //
    d_CDC15_dt = CDC15_activation - inactivation_2;                                                                                                                                                                                                        //
    d_CDC20_dt = -Degradation_of_active_CDC20 + Activation_of_CDC20 - Inactivation_2;                                                                                                                                                                      //
    d_CDC20i_dt = (Synthesis_of_inactive_CDC20 - Degradation_of_inactiveCDC20 - Activation_of_CDC20) + Inactivation_2;                                                                                                                                     //
    d_CDC6_dt = (((CDC6_synthesis - Phosphorylation_of_CDC6) + Dephosphorylation_of_CDC6 - CLB2CDC6_complex_formation) + CLB2CDC6_dissociation - CLB5CDC6_complex_formation) + CLB5CDC6_dissociation + CLB2_degradation_in_F2 + CLB5_degradation_in_F5;    //
    d_CDC6P_dt = (Phosphorylation_of_CDC6 - Dephosphorylation_of_CDC6 - Degradation_of_CDC6P) + CLB2_degradation_in_F2P + CLB5_degradation_in_F5P;                                                                                                         //
    d_CDH1_dt = (CDH1_synthesis - CDH1_degradation) + CDH1i_activation - Inactivation_3;                                                                                                                                                                   //
    d_CDH1i_dt = (-CDH1i_degradation - CDH1i_activation) + Inactivation_3;                                                                                                                                                                                 //
    d_CLB2_dt = ((Synthesis_of_CLB2 - Degradation_of_CLB2 - Assoc_of_CLB2_and_SIC1) + Dissoc_of_CLB2SIC1_complex + Degradation_of_SIC1_in_C2P - CLB2CDC6_complex_formation) + CLB2CDC6_dissociation + CDC6_degradation_in_F2P;                             //
    d_CLB5_dt = ((Synthesis_of_CLB5 - Degradation_of_CLB5 - Assoc_of_CLB5_and_SIC1) + Dissoc_of_CLB5SIC1 + Degradation_of_SIC1P_in_C5P_ - CLB5CDC6_complex_formation) + CLB5CDC6_dissociation + CDC6_degradation_in_F5P;                                   //
    d_CLN2_dt = Synthesis_of_CLN2 - Degradation_of_CLN2;                                                                                                                                                                                                   //
    d_ESP1_dt = (Degradation_of_PDS1_in_PE - Assoc_with_ESP1_to_form_PE) + Disso_from_PE;                                                                                                                                                                  //
    d_F2_dt = (CLB2CDC6_complex_formation - CLB2CDC6_dissociation - F2_phosphorylation) + F2P_dephosphorylation - CLB2_degradation_in_F2;                                                                                                                  //
    d_F2P_dt = F2_phosphorylation - F2P_dephosphorylation - CDC6_degradation_in_F2P - CLB2_degradation_in_F2P;                                                                                                                                             //
    d_F5_dt = (CLB5CDC6_complex_formation - CLB5CDC6_dissociation - F5_phosphorylation) + F5P_dephosphorylation - CLB5_degradation_in_F5;                                                                                                                  //
    d_F5P_dt = F5_phosphorylation - F5P_dephosphorylation - CDC6_degradation_in_F5P - CLB5_degradation_in_F5P;                                                                                                                                             //
    d_IEP_dt = Activation_of_IEP - Inactivation_1;                                                                                                                                                                                                         //
    d_MASS_dt = Growth;                                                                                                                                                                                                                                    //
    d_NET1_dt = (-Assoc_with_NET1_to_form_RENT + Dissoc_from_RENT + Net1_synthesis - Net1_degradation - NET1_phosphorylation) + dephosphorylation_1 + Degradation_of_CDC14_in_RENT;                                                                        //
    d_NET1P_dt = ((-Assoc_with_NET1P_to_form_RENTP + Dissoc_from_RENP - Net1P_degradation) + NET1_phosphorylation - dephosphorylation_1) + Degradation_of_CDC14_in_RENTP;                                                                                  //
    d_ORI_dt = DNA_synthesis - Negative_regulation_of_DNA_synthesis;                                                                                                                                                                                       //
    d_PDS1_dt = (PDS1_synthesis - degradation_2 - Assoc_with_ESP1_to_form_PE) + Disso_from_PE;                                                                                                                                                             //
    d_PPX_dt = PPX_synthesis - degradation_1;                                                                                                                                                                                                              //
    d_RENT_dt = (Assoc_with_NET1_to_form_RENT - Dissoc_from_RENT - RENT_phosphorylation) + dephosphorylation_2 - Degradation_of_NET1_in_RENT - Degradation_of_CDC14_in_RENT;                                                                               //
    d_RENTP_dt = (Assoc_with_NET1P_to_form_RENTP - Dissoc_from_RENP) + RENT_phosphorylation - dephosphorylation_2 - Degradation_of_NET1P_in_RENTP - Degradation_of_CDC14_in_RENTP;                                                                         //
    d_SIC1_dt = (((Synthesis_of_SIC1 - Phosphorylation_of_SIC1) + Dephosphorylation_of_SIC1 - Assoc_of_CLB2_and_SIC1) + Dissoc_of_CLB2SIC1_complex - Assoc_of_CLB5_and_SIC1) + Dissoc_of_CLB5SIC1 + Degradation_of_CLB2_in_C2 + Degradation_of_CLB5_in_C5; //
    d_SIC1P_dt = (Phosphorylation_of_SIC1 - Dephosphorylation_of_SIC1 - Fast_Degradation_of_SIC1P) + Degradation_of_CLB2_in_C2P + Degradation_of_CLB5_in_C5P;                                                                                              //
    d_SPN_dt = Spindle_formation - Spindle_disassembly;                                                                                                                                                                                                    //
    d_SWI5_dt = (Synthesis_of_SWI5 - Degradation_of_SWI5) + Activation_of_SWI5 - Inactivation_of_SWI5;                                                                                                                                                     //
    d_SWI5P_dt = (-Degradation_of_SWI5P - Activation_of_SWI5) + Inactivation_of_SWI5;                                                                                                                                                                      //
    d_TEM1GTP_dt = TEM1_activation - inactivation_1;                                                                                                                                                                                                       //

    std::vector<double> derivatives(36);
    derivatives[0] = d_BUD_dt;
    derivatives[1] = d_C2_dt;
    derivatives[2] = d_C2P_dt;
    derivatives[3] = d_C5_dt;
    derivatives[4] = d_C5P_dt;
    derivatives[5] = d_CDC14_dt;
    derivatives[6] = d_CDC15_dt;
    derivatives[7] = d_CDC20_dt;
    derivatives[8] = d_CDC20i_dt;
    derivatives[9] = d_CDC6_dt;
    derivatives[10] = d_CDC6P_dt;
    derivatives[11] = d_CDH1_dt;
    derivatives[12] = d_CDH1i_dt;
    derivatives[13] = d_CLB2_dt;
    derivatives[14] = d_CLB5_dt;
    derivatives[15] = d_CLN2_dt;
    derivatives[16] = d_ESP1_dt;
    derivatives[17] = d_F2_dt;
    derivatives[18] = d_F2P_dt;
    derivatives[19] = d_F5_dt;
    derivatives[20] = d_F5P_dt;
    derivatives[21] = d_IEP_dt;
    derivatives[22] = d_MASS_dt;
    derivatives[23] = d_NET1_dt;
    derivatives[24] = d_NET1P_dt;
    derivatives[25] = d_ORI_dt;
    derivatives[26] = d_PDS1_dt;
    derivatives[27] = d_PPX_dt;
    derivatives[28] = d_RENT_dt;
    derivatives[29] = d_RENTP_dt;
    derivatives[30] = d_SIC1_dt;
    derivatives[31] = d_SIC1P_dt;
    derivatives[32] = d_SPN_dt;
    derivatives[33] = d_SWI5_dt;
    derivatives[34] = d_SWI5P_dt;
    derivatives[35] = d_TEM1GTP_dt;
    return derivatives;
}

// MODEL FUNCTIONS
inline double Chen2004SbmlOdeSystem::BB_218(double A1, double A2, double A3, double A4)
{
    return (A2 - A1) + A3 * A2 + A4 * A1;
}

inline double Chen2004SbmlOdeSystem::GK_219(double A1, double A2, double A3, double A4)
{
    return 2.0 * A4 * A1 / ((A2 - A1) + A3 * A2 + A4 * A1 + sm::root(2.0, std::pow((A2 - A1) + A3 * A2 + A4 * A1, 2.0) - 4.0 * (A2 - A1) * A4 * A1));
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
    this->mDerivedQuantityNames.push_back("cell");
    this->mDerivedQuantityUnits.push_back("non-dim");

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

    this->mDerivedQuantityNames.push_back("D");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("mu");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Vdb5");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Vdb2");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Vasbf");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Visbf");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Vkpc1");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Vkpf6");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Vacdh");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Vicdh");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Vppnet");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Vkpnet");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Vdppx");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Vdpds");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Vaiep");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Vd2c1");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Vd2f6");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Vppc1");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Vppf6");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("F");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Growth");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Synthesis_of_CLN2");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Degradation_of_CLN2");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Synthesis_of_CLB2");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Degradation_of_CLB2");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Synthesis_of_CLB5");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Degradation_of_CLB5");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Synthesis_of_SIC1");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Phosphorylation_of_SIC1");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Dephosphorylation_of_SIC1");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Fast_Degradation_of_SIC1P");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Assoc_of_CLB2_and_SIC1");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Dissoc_of_CLB2SIC1_complex");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Assoc_of_CLB5_and_SIC1");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Dissoc_of_CLB5SIC1");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Phosphorylation_of_C2");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Dephosphorylation_of_C2P");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Phosphorylation_of_C5");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Dephosphorylation_of_C5P");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Degradation_of_CLB2_in_C2");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Degradation_of_CLB5_in_C5");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Degradation_of_SIC1_in_C2P");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Degradation_of_SIC1P_in_C5P_");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Degradation_of_CLB2_in_C2P");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Degradation_of_CLB5_in_C5P");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("CDC6_synthesis");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Phosphorylation_of_CDC6");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Dephosphorylation_of_CDC6");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Degradation_of_CDC6P");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("CLB2CDC6_complex_formation");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("CLB2CDC6_dissociation");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("CLB5CDC6_complex_formation");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("CLB5CDC6_dissociation");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("F2_phosphorylation");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("F2P_dephosphorylation");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("F5_phosphorylation");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("F5P_dephosphorylation");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("CLB2_degradation_in_F2");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("CLB5_degradation_in_F5");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("CDC6_degradation_in_F2P");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("CDC6_degradation_in_F5P");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("CLB2_degradation_in_F2P");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("CLB5_degradation_in_F5P");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Synthesis_of_SWI5");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Degradation_of_SWI5");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Degradation_of_SWI5P");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Activation_of_SWI5");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Inactivation_of_SWI5");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Activation_of_IEP");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Inactivation_1");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Synthesis_of_inactive_CDC20");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Degradation_of_inactiveCDC20");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Degradation_of_active_CDC20");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Activation_of_CDC20");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Inactivation_2");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("CDH1_synthesis");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("CDH1_degradation");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("CDH1i_degradation");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("CDH1i_activation");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Inactivation_3");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("CDC14_synthesis");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("CDC14_degradation");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Assoc_with_NET1_to_form_RENT");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Dissoc_from_RENT");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Assoc_with_NET1P_to_form_RENTP");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Dissoc_from_RENP");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Net1_synthesis");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Net1_degradation");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Net1P_degradation");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("NET1_phosphorylation");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("dephosphorylation_1");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("RENT_phosphorylation");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("dephosphorylation_2");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Degradation_of_NET1_in_RENT");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Degradation_of_NET1P_in_RENTP");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Degradation_of_CDC14_in_RENT");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Degradation_of_CDC14_in_RENTP");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("TEM1_activation");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("inactivation_1");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("CDC15_activation");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("inactivation_2");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("PPX_synthesis");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("degradation_1");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("PDS1_synthesis");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("degradation_2");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Degradation_of_PDS1_in_PE");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Assoc_with_ESP1_to_form_PE");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Disso_from_PE");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("DNA_synthesis");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Negative_regulation_of_DNA_synthesis");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Budding");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Negative_regulation_of_Cell_budding");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Spindle_formation");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Spindle_disassembly");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__BCK2");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__BUB2");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__BUD");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__C2");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__C2P");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__C5");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__C5P");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__CDC14");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__CDC14T");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__CDC15");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__CDC15i");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__CDC20");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__CDC20i");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__CDC6");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__CDC6P");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__CDC6T");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__CDH1");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__CDH1i");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__CKIT");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__CLB2");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__CLB2T");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__CLB5");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__CLB5T");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__CLN2");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__CLN3");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__ESP1");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__F2");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__F2P");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__F5");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__F5P");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__IE");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__IEP");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__LTE1");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__MAD2");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__MASS");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__MCM1");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__NET1");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__NET1P");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__NET1T");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__ORI");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__PDS1");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__PE");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__PPX");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__RENT");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__RENTP");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__SBF");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__SIC1");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__SIC1P");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__SIC1T");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__SPN");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__SWI5");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__SWI5P");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__TEM1GDP");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("conc__TEM1GTP");
    this->mDerivedQuantityUnits.push_back("non-dim");

    // PARAMETERS
    this->mParameterNames.push_back("BUB2");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("LTE1");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("MAD2");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("b0");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("bub2h");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("bub2l");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("C0");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("CDC15T");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Dn3");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ebudb5");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ebudn2");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ebudn3");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ec1b2");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ec1b5");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ec1k2");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ec1n2");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ec1n3");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ef6b2");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ef6b5");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ef6k2");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ef6n2");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ef6n3");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("eicdhb2");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("eicdhb5");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("eicdhn2");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("eicdhn3");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("eorib2");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("eorib5");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("esbfb5");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("esbfn2");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("esbfn3");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ESP1T");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("IET");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("J20ppx");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Jacdh");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Jaiep");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Jamcm");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Jasbf");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Jatem");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Jd2c1");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Jd2f6");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Jicdh");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Jiiep");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Jimcm");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Jisbf");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Jitem");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Jn3");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Jpds");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Jspn");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ka15_p");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ka15_p_p");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ka15p");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ka20_p");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ka20_p_p");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kacdh_p");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kacdh_p_p");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kaiep");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kamcm");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kasb2");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kasb5");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kasbf");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kasesp");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kasf2");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kasf5");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kasrent");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kasrentp");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kaswi");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kd14");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kd1c1");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kd1f6");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kd1pds_p");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kd20");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kd2c1");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kd2f6");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kd2pds_p_p");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kd3c1");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kd3f6");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kd3pds_p_p");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kdb2_p");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kdb2_p_p");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kdb2p");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kdb5_p");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kdb5_p_p");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kdbud");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kdcdh");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kdib2");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kdib5");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kdiesp");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kdif2");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kdif5");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kdirent");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kdirentp");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kdn2");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kdnet");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kdori");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kdppx_p");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kdppx_p_p");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kdspn");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kdswi");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("KEZ");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("KEZ2");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ki15");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kicdh_p");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kicdh_p_p");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kiiep");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kimcm");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kisbf_p");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kisbf_p_p");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kiswi");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kkpnet_p");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kkpnet_p_p");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kppc1");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kppf6");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kppnet_p");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kppnet_p_p");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ks14");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ks1pds_p_p");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ks20_p");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ks20_p_p");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ks2pds_p_p");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ksb2_p");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ksb2_p_p");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ksb5_p");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ksb5_p_p");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ksbud");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ksc1_p");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ksc1_p_p");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kscdh");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ksf6_p");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ksf6_p_p");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ksf6_p_p_p");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ksn2_p");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ksn2_p_p");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ksnet");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ksori");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kspds_p");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ksppx");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ksspn");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ksswi_p");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ksswi_p_p");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("lte1h");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("lte1l");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("mad2h");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("mad2l");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("mdt");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("TEM1T");
    this->mParameterUnits.push_back("non-dim");

    this->mInitialised = true;
}

// Register the ODE system with Boost serialization
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(Chen2004SbmlOdeSystem)