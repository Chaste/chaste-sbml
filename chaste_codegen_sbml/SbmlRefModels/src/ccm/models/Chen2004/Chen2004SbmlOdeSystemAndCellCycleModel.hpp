#ifndef CHEN2004SBMLODESYSTEMANDCELLCYCLEMODEL_HPP_
#define CHEN2004SBMLODESYSTEMANDCELLCYCLEMODEL_HPP_

#include <vector>

#include <boost/serialization/base_object.hpp>

#include "AbstractOdeSystem.hpp"
#include "ChasteSerialization.hpp"

class Chen2004SbmlOdeSystem : public AbstractOdeSystem
{
private:
    // (De-)serialize Chen2004SbmlOdeSystem
    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(AbstractOdeSystem);
    }

    // CONSTANT PARAMETERS
    const double b0 = 0.054; // b0
    const double bub2h = 1.0; // bub2h
    const double C0 = 0.4; // C0
    const double Dn3 = 1.0; // Dn3
    const double ebudb5 = 1.0; // ebudb5
    const double ebudn2 = 0.25; // ebudn2
    const double ebudn3 = 0.05; // ebudn3
    const double ec1b2 = 0.45; // ec1b2
    const double ec1b5 = 0.1; // ec1b5
    const double ec1k2 = 0.03; // ec1k2
    const double ec1n2 = 0.06; // ec1n2
    const double ec1n3 = 0.3; // ec1n3
    const double ef6b2 = 0.55; // ef6b2
    const double ef6b5 = 0.1; // ef6b5
    const double ef6k2 = 0.03; // ef6k2
    const double ef6n2 = 0.06; // ef6n2
    const double ef6n3 = 0.3; // ef6n3
    const double eicdhb2 = 1.2; // eicdhb2
    const double eicdhb5 = 8.0; // eicdhb5
    const double eicdhn2 = 0.4; // eicdhn2
    const double eicdhn3 = 0.25; // eicdhn3
    const double eorib2 = 0.45; // eorib2
    const double eorib5 = 0.9; // eorib5
    const double esbfb5 = 2.0; // esbfb5
    const double esbfn2 = 2.0; // esbfn2
    const double esbfn3 = 10.0; // esbfn3
    const double J20ppx = 0.15; // J20ppx
    const double Jacdh = 0.03; // Jacdh
    const double Jaiep = 0.1; // Jaiep
    const double Jamcm = 0.1; // Jamcm
    const double Jasbf = 0.01; // Jasbf
    const double Jatem = 0.1; // Jatem
    const double Jd2c1 = 0.05; // Jd2c1
    const double Jd2f6 = 0.05; // Jd2f6
    const double Jicdh = 0.03; // Jicdh
    const double Jiiep = 0.1; // Jiiep
    const double Jimcm = 0.1; // Jimcm
    const double Jisbf = 0.01; // Jisbf
    const double Jitem = 0.1; // Jitem
    const double Jn3 = 6.0; // Jn3
    const double Jpds = 0.04; // Jpds
    const double Jspn = 0.14; // Jspn
    const double ka15_p = 0.002; // ka15'
    const double ka15_p_p = 1.0; // ka15''
    const double ka15p = 0.001; // ka15p
    const double ka20_p = 0.05; // ka20'
    const double ka20_p_p = 0.2; // ka20''
    const double kacdh_p = 0.01; // kacdh'
    const double kacdh_p_p = 0.8; // kacdh''
    const double kaiep = 0.1; // kaiep
    const double kamcm = 1.0; // kamcm
    const double kasb2 = 50.0; // kasb2
    const double kasb5 = 50.0; // kasb5
    const double kasbf = 0.38; // kasbf
    const double kasesp = 50.0; // kasesp
    const double kasf2 = 15.0; // kasf2
    const double kasf5 = 0.01; // kasf5
    const double kasrent = 200.0; // kasrent
    const double kasrentp = 1.0; // kasrentp
    const double kaswi = 2.0; // kaswi
    const double kd14 = 0.1; // kd14
    const double kd1c1 = 0.01; // kd1c1
    const double kd1f6 = 0.01; // kd1f6
    const double kd1pds_p = 0.01; // kd1pds'
    const double kd20 = 0.3; // kd20
    const double kd2c1 = 1.0; // kd2c1
    const double kd2f6 = 1.0; // kd2f6
    const double kd2pds_p_p = 0.2; // kd2pds''
    const double kd3c1 = 1.0; // kd3c1
    const double kd3f6 = 1.0; // kd3f6
    const double kd3pds_p_p = 0.04; // kd3pds''
    const double kdb2_p = 0.003; // kdb2'
    const double kdb2_p_p = 0.4; // kdb2''
    const double kdb2p = 0.15; // kdb2p
    const double kdb5_p = 0.01; // kdb5'
    const double kdb5_p_p = 0.16; // kdb5''
    const double kdbud = 0.06; // kdbud
    const double kdcdh = 0.01; // kdcdh
    const double kdib2 = 0.05; // kdib2
    const double kdib5 = 0.06; // kdib5
    const double kdiesp = 0.5; // kdiesp
    const double kdif2 = 0.5; // kdif2
    const double kdif5 = 0.01; // kdif5
    const double kdirent = 1.0; // kdirent
    const double kdirentp = 2.0; // kdirentp
    const double kdn2 = 0.12; // kdn2
    const double kdnet = 0.03; // kdnet
    const double kdori = 0.06; // kdori
    const double kdppx_p = 0.17; // kdppx'
    const double kdppx_p_p = 2.0; // kdppx''
    const double kdspn = 0.06; // kdspn
    const double kdswi = 0.08; // kdswi
    const double ki15 = 0.5; // ki15
    const double kicdh_p = 0.001; // kicdh'
    const double kicdh_p_p = 0.08; // kicdh''
    const double kiiep = 0.15; // kiiep
    const double kimcm = 0.15; // kimcm
    const double kisbf_p = 0.6; // kisbf'
    const double kisbf_p_p = 8.0; // kisbf''
    const double kiswi = 0.05; // kiswi
    const double kkpnet_p = 0.01; // kkpnet'
    const double kkpnet_p_p = 0.6; // kkpnet''
    const double kppc1 = 4.0; // kppc1
    const double kppf6 = 4.0; // kppf6
    const double kppnet_p = 0.05; // kppnet'
    const double kppnet_p_p = 3.0; // kppnet''
    const double ks14 = 0.2; // ks14
    const double ks1pds_p_p = 0.03; // ks1pds''
    const double ks20_p = 0.006; // ks20'
    const double ks20_p_p = 0.6; // ks20''
    const double ks2pds_p_p = 0.055; // ks2pds''
    const double ksb2_p = 0.001; // ksb2'
    const double ksb2_p_p = 0.04; // ksb2''
    const double ksb5_p = 0.0008; // ksb5'
    const double ksb5_p_p = 0.005; // ksb5''
    const double ksbud = 0.2; // ksbud
    const double ksc1_p = 0.012; // ksc1'
    const double ksc1_p_p = 0.12; // ksc1''
    const double kscdh = 0.01; // kscdh
    const double ksf6_p = 0.024; // ksf6'
    const double ksf6_p_p = 0.12; // ksf6''
    const double ksf6_p_p_p = 0.004; // ksf6'''
    const double ksn2_p = 0.0; // ksn2'
    const double ksn2_p_p = 0.15; // ksn2''
    const double ksnet = 0.084; // ksnet
    const double ksori = 2.0; // ksori
    const double kspds_p = 0.0; // kspds'
    const double ksppx = 0.1; // ksppx
    const double ksspn = 0.1; // ksspn
    const double ksswi_p = 0.005; // ksswi'
    const double ksswi_p_p = 0.08; // ksswi''
    const double mad2h = 8.0; // mad2h
    const double mdt = 90.0; // mdt

    // STATE VARIABLES
    double BUD; // BUD
    double C2; // C2
    double C2P; // C2P
    double C5; // C5
    double C5P; // C5P
    double CDC14; // CDC14
    double CDC15; // CDC15
    double CDC20; // CDC20
    double CDC20i; // CDC20i
    double CDC6; // CDC6
    double CDC6P; // CDC6P
    double CDH1; // CDH1
    double CDH1i; // CDH1i
    double CLB2; // CLB2
    double CLB5; // CLB5
    double CLN2; // CLN2
    double ESP1; // ESP1
    double F2; // F2
    double F2P; // F2P
    double F5; // F5
    double F5P; // F5P
    double IEP; // IEP
    double MASS; // MASS
    double NET1; // NET1
    double NET1P; // NET1P
    double ORI; // ORI
    double PDS1; // PDS1
    double PPX; // PPX
    double RENT; // RENT
    double RENTP; // RENTP
    double SIC1; // SIC1
    double SIC1P; // SIC1P
    double SPN; // SPN
    double SWI5; // SWI5
    double SWI5P; // SWI5P
    double TEM1GTP; // TEM1GTP

    // DERIVED QUANTITIES
    double BCK2; // BCK2
    double CDC14T; // CDC14T
    double CDC15i; // CDC15i
    double CDC6T; // CDC6T
    double CKIT; // CKIT
    double CLB2T; // CLB2T
    double CLB5T; // CLB5T
    double CLN3; // CLN3
    double IE; // IE
    double MCM1; // MCM1
    double NET1T; // NET1T
    double PE; // PE
    double SBF; // SBF
    double SIC1T; // SIC1T
    double TEM1GDP; // TEM1GDP

    // VARIABLE PARAMETERS
    double cell; // cell
    double BUB2; // BUB2
    double LTE1; // LTE1
    double MAD2; // MAD2
    double bub2l; // bub2l
    double CDC15T; // CDC15T
    double ESP1T; // ESP1T
    double IET; // IET
    double KEZ; // KEZ
    double KEZ2; // KEZ2
    double lte1h; // lte1h
    double lte1l; // lte1l
    double mad2l; // mad2l
    double TEM1T; // TEM1T

    // RULE-BASED PARAMETERS
    double D; // D
    double mu; // mu
    double Vdb5; // Vdb5
    double Vdb2; // Vdb2
    double Vasbf; // Vasbf
    double Visbf; // Visbf
    double Vkpc1; // Vkpc1
    double Vkpf6; // Vkpf6
    double Vacdh; // Vacdh
    double Vicdh; // Vicdh
    double Vppnet; // Vppnet
    double Vkpnet; // Vkpnet
    double Vdppx; // Vdppx
    double Vdpds; // Vdpds
    double Vaiep; // Vaiep
    double Vd2c1; // Vd2c1
    double Vd2f6; // Vd2f6
    double Vppc1; // Vppc1
    double Vppf6; // Vppf6
    double F; // F

    // REACTIONS
    double Growth; // Growth
    double Synthesis_of_CLN2; // Synthesis of CLN2
    double Degradation_of_CLN2; // Degradation of CLN2
    double Synthesis_of_CLB2; // Synthesis of CLB2
    double Degradation_of_CLB2; // Degradation of CLB2
    double Synthesis_of_CLB5; // Synthesis of CLB5
    double Degradation_of_CLB5; // Degradation of CLB5
    double Synthesis_of_SIC1; // Synthesis of SIC1
    double Phosphorylation_of_SIC1; // Phosphorylation of SIC1
    double Dephosphorylation_of_SIC1; // Dephosphorylation of SIC1
    double Fast_Degradation_of_SIC1P; // Fast Degradation of SIC1P
    double Assoc_of_CLB2_and_SIC1; // Assoc. of CLB2 and SIC1
    double Dissoc_of_CLB2SIC1_complex; // Dissoc. of CLB2/SIC1 complex
    double Assoc_of_CLB5_and_SIC1; // Assoc. of CLB5 and SIC1
    double Dissoc_of_CLB5SIC1; // Dissoc. of CLB5/SIC1
    double Phosphorylation_of_C2; // Phosphorylation of C2
    double Dephosphorylation_of_C2P; // Dephosphorylation of C2P
    double Phosphorylation_of_C5; // Phosphorylation of C5
    double Dephosphorylation_of_C5P; // Dephosphorylation of C5P
    double Degradation_of_CLB2_in_C2; // Degradation of CLB2 in C2
    double Degradation_of_CLB5_in_C5; // Degradation of CLB5 in C5
    double Degradation_of_SIC1_in_C2P; // Degradation of SIC1 in C2P
    double Degradation_of_SIC1P_in_C5P_; // Degradation of SIC1P in C5P
    double Degradation_of_CLB2_in_C2P; // Degradation of CLB2 in C2P
    double Degradation_of_CLB5_in_C5P; // Degradation of CLB5 in C5P
    double CDC6_synthesis; // CDC6 synthesis
    double Phosphorylation_of_CDC6; // Phosphorylation of CDC6
    double Dephosphorylation_of_CDC6; // Dephosphorylation of CDC6
    double Degradation_of_CDC6P; // Degradation of CDC6P
    double CLB2CDC6_complex_formation; // CLB2/CDC6 complex formation
    double CLB2CDC6_dissociation; // CLB2/CDC6 dissociation
    double CLB5CDC6_complex_formation; // CLB5/CDC6 complex formation
    double CLB5CDC6_dissociation; // CLB5/CDC6 dissociation
    double F2_phosphorylation; // F2 phosphorylation
    double F2P_dephosphorylation; // F2P dephosphorylation
    double F5_phosphorylation; // F5 phosphorylation
    double F5P_dephosphorylation; // F5P dephosphorylation
    double CLB2_degradation_in_F2; // CLB2 degradation in F2
    double CLB5_degradation_in_F5; // CLB5 degradation in F5
    double CDC6_degradation_in_F2P; // CDC6 degradation in F2P
    double CDC6_degradation_in_F5P; // CDC6 degradation in F5P
    double CLB2_degradation_in_F2P; // CLB2 degradation in F2P
    double CLB5_degradation_in_F5P; // CLB5 degradation in F5P
    double Synthesis_of_SWI5; // Synthesis of SWI5
    double Degradation_of_SWI5; // Degradation of SWI5
    double Degradation_of_SWI5P; // Degradation of SWI5P
    double Activation_of_SWI5; // Activation of SWI5
    double Inactivation_of_SWI5; // Inactivation of SWI5
    double Activation_of_IEP; // Activation of IEP
    double Inactivation_1; // Inactivation
    double Synthesis_of_inactive_CDC20; // Synthesis of inactive CDC20
    double Degradation_of_inactiveCDC20; // Degradation of inactiveCDC20
    double Degradation_of_active_CDC20; // Degradation of active CDC20
    double Activation_of_CDC20; // Activation of CDC20
    double Inactivation_2; // Inactivation
    double CDH1_synthesis; // CDH1 synthesis
    double CDH1_degradation; // CDH1 degradation
    double CDH1i_degradation; // CDH1i degradation
    double CDH1i_activation; // CDH1i activation
    double Inactivation_3; // Inactivation
    double CDC14_synthesis; // CDC14 synthesis
    double CDC14_degradation; // CDC14 degradation
    double Assoc_with_NET1_to_form_RENT; // Assoc. with NET1 to form RENT
    double Dissoc_from_RENT; // Dissoc. from RENT
    double Assoc_with_NET1P_to_form_RENTP; // Assoc with NET1P to form RENTP
    double Dissoc_from_RENP; // Dissoc. from RENP
    double Net1_synthesis; // Net1 synthesis
    double Net1_degradation; // Net1 degradation
    double Net1P_degradation; // Net1P degradation
    double NET1_phosphorylation; // NET1 phosphorylation
    double dephosphorylation_1; // dephosphorylation
    double RENT_phosphorylation; // RENT phosphorylation
    double dephosphorylation_2; // dephosphorylation
    double Degradation_of_NET1_in_RENT; // Degradation of NET1 in RENT
    double Degradation_of_NET1P_in_RENTP; // Degradation of NET1P in RENTP
    double Degradation_of_CDC14_in_RENT; // Degradation of CDC14 in RENT
    double Degradation_of_CDC14_in_RENTP; // Degradation of CDC14 in RENTP
    double TEM1_activation; // TEM1 activation
    double inactivation_1; // inactivation
    double CDC15_activation; // CDC15 activation
    double inactivation_2; // inactivation
    double PPX_synthesis; // PPX synthesis
    double degradation_1; // degradation
    double PDS1_synthesis; // PDS1 synthesis
    double degradation_2; // degradation
    double Degradation_of_PDS1_in_PE; // Degradation of PDS1 in PE
    double Assoc_with_ESP1_to_form_PE; // Assoc. with ESP1 to form PE
    double Disso_from_PE; // Disso. from PE
    double DNA_synthesis; // DNA synthesis
    double Negative_regulation_of_DNA_synthesis; // Negative regulation of DNA synthesis
    double Budding; // Budding
    double Negative_regulation_of_Cell_budding; // Negative regulation of Cell budding
    double Spindle_formation; // Spindle formation
    double Spindle_disassembly; // Spindle disassembly

    // EVENTS
    std::vector<bool> mEventsSatisfied;
    bool mEventsInitialised;

public:
    Chen2004SbmlOdeSystem(std::vector<double> stateVariables = std::vector<double>());

    ~Chen2004SbmlOdeSystem();

    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY);
    void ProcessRules(double time, const std::vector<double>& rY);

    std::vector<double> ComputeDerivedQuantities(double time, const std::vector<double> &rY);

    double CalculateRootFunction(double time, const std::vector<double>& rY);
    bool CalculateStoppingEvent(double time, const std::vector<double>& rY);
    double ProcessEvents(double time, const std::vector<double>& rY);

    // FUNCTIONS
    inline double BB_218(double A1, double A2, double A3, double A4);
    inline double GK_219(double A1, double A2, double A3, double A4);
    inline double MichaelisMenten_220(double M1, double J1, double k1, double S1);
    inline double Mass_Action_2_221(double k1, double S1, double S2);
    inline double Mass_Action_1_222(double k1, double S1);
};

namespace
{
    namespace serialization
    {
        // Provide constructor for serializing Chen2004SbmlOdeSystem
        template <class Archive>
        inline void save_construct_data(Archive &ar, const Chen2004SbmlOdeSystem * t, const unsigned int version)
        {
            // Save data required to construct instance
            const std::vector<double> state_variables = t->rGetConstStateVariables();
            ar << state_variables;
        }

        // Provide constructor for de-serializing Chen2004SbmlOdeSystem
        template <class Archive>
        inline void load_construct_data(Archive &ar, Chen2004SbmlOdeSystem * t, const unsigned int version)
        {
            // Retrieve data from archive required to construct new instance
            std::vector<double> state_variables;
            ar >> state_variables;

            // Invoke inplace constructor to initialise instance
            ::new (t)Chen2004SbmlOdeSystem(state_variables);
        }
    } // namespace serialization
} // namespace ...

// Define SbmlCellCycleWrapperModel using wrappers
#include "SbmlCellCycleWrapperModel.hpp"
#include "SbmlCellCycleWrapperModel.cpp"

typedef SbmlCellCycleWrapperModel<Chen2004SbmlOdeSystem, 36> Chen2004SbmlCellCycleModel;

// Declare identifiers for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(Chen2004SbmlOdeSystem)
EXPORT_TEMPLATE_CLASS2(SbmlCellCycleWrapperModel, Chen2004SbmlOdeSystem, 36)

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(Chen2004SbmlCellCycleModel)

#endif // CHEN2004SBMLODESYSTEMANDCELLCYCLEMODEL_HPP_