#ifndef CHEN_2004_SBML_ODE_SYSTEM_HPP_
#define CHEN_2004_SBML_ODE_SYSTEM_HPP_

#include <vector>

#include <boost/serialization/base_object.hpp>

#include "AbstractSbmlOdeSystem.hpp"
#include "ChasteSerialization.hpp"
#include "SbmlEventType.hpp"

class Chen2004SbmlOdeSystem : public AbstractSbmlOdeSystem
{
private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Save / load Chen2004SbmlOdeSystem archive
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive& boost::serialization::base_object<AbstractSbmlOdeSystem>(*this);
    }

    // PARAMETERS
    double b0;         // b0
    double bub2h;      // bub2h
    double bub2l;      // bub2l
    double C0;         // C0
    double CDC15T;     // CDC15T
    double Dn3;        // Dn3
    double ebudb5;     // ebudb5
    double ebudn2;     // ebudn2
    double ebudn3;     // ebudn3
    double ec1b2;      // ec1b2
    double ec1b5;      // ec1b5
    double ec1k2;      // ec1k2
    double ec1n2;      // ec1n2
    double ec1n3;      // ec1n3
    double ef6b2;      // ef6b2
    double ef6b5;      // ef6b5
    double ef6k2;      // ef6k2
    double ef6n2;      // ef6n2
    double ef6n3;      // ef6n3
    double eicdhb2;    // eicdhb2
    double eicdhb5;    // eicdhb5
    double eicdhn2;    // eicdhn2
    double eicdhn3;    // eicdhn3
    double eorib2;     // eorib2
    double eorib5;     // eorib5
    double esbfb5;     // esbfb5
    double esbfn2;     // esbfn2
    double esbfn3;     // esbfn3
    double ESP1T;      // ESP1T
    double IET;        // IET
    double J20ppx;     // J20ppx
    double Jacdh;      // Jacdh
    double Jaiep;      // Jaiep
    double Jamcm;      // Jamcm
    double Jasbf;      // Jasbf
    double Jatem;      // Jatem
    double Jd2c1;      // Jd2c1
    double Jd2f6;      // Jd2f6
    double Jicdh;      // Jicdh
    double Jiiep;      // Jiiep
    double Jimcm;      // Jimcm
    double Jisbf;      // Jisbf
    double Jitem;      // Jitem
    double Jn3;        // Jn3
    double Jpds;       // Jpds
    double Jspn;       // Jspn
    double ka15_p;     // ka15'
    double ka15_p_p;   // ka15''
    double ka15p;      // ka15p
    double ka20_p;     // ka20'
    double ka20_p_p;   // ka20''
    double kacdh_p;    // kacdh'
    double kacdh_p_p;  // kacdh''
    double kaiep;      // kaiep
    double kamcm;      // kamcm
    double kasb2;      // kasb2
    double kasb5;      // kasb5
    double kasbf;      // kasbf
    double kasesp;     // kasesp
    double kasf2;      // kasf2
    double kasf5;      // kasf5
    double kasrent;    // kasrent
    double kasrentp;   // kasrentp
    double kaswi;      // kaswi
    double kd14;       // kd14
    double kd1c1;      // kd1c1
    double kd1f6;      // kd1f6
    double kd1pds_p;   // kd1pds'
    double kd20;       // kd20
    double kd2c1;      // kd2c1
    double kd2f6;      // kd2f6
    double kd2pds_p_p; // kd2pds''
    double kd3c1;      // kd3c1
    double kd3f6;      // kd3f6
    double kd3pds_p_p; // kd3pds''
    double kdb2_p;     // kdb2'
    double kdb2_p_p;   // kdb2''
    double kdb2p;      // kdb2p
    double kdb5_p;     // kdb5'
    double kdb5_p_p;   // kdb5''
    double kdbud;      // kdbud
    double kdcdh;      // kdcdh
    double kdib2;      // kdib2
    double kdib5;      // kdib5
    double kdiesp;     // kdiesp
    double kdif2;      // kdif2
    double kdif5;      // kdif5
    double kdirent;    // kdirent
    double kdirentp;   // kdirentp
    double kdn2;       // kdn2
    double kdnet;      // kdnet
    double kdori;      // kdori
    double kdppx_p;    // kdppx'
    double kdppx_p_p;  // kdppx''
    double kdspn;      // kdspn
    double kdswi;      // kdswi
    double KEZ;        // KEZ
    double KEZ2;       // KEZ2
    double ki15;       // ki15
    double kicdh_p;    // kicdh'
    double kicdh_p_p;  // kicdh''
    double kiiep;      // kiiep
    double kimcm;      // kimcm
    double kisbf_p;    // kisbf'
    double kisbf_p_p;  // kisbf''
    double kiswi;      // kiswi
    double kkpnet_p;   // kkpnet'
    double kkpnet_p_p; // kkpnet''
    double kppc1;      // kppc1
    double kppf6;      // kppf6
    double kppnet_p;   // kppnet'
    double kppnet_p_p; // kppnet''
    double ks14;       // ks14
    double ks1pds_p_p; // ks1pds''
    double ks20_p;     // ks20'
    double ks20_p_p;   // ks20''
    double ks2pds_p_p; // ks2pds''
    double ksb2_p;     // ksb2'
    double ksb2_p_p;   // ksb2''
    double ksb5_p;     // ksb5'
    double ksb5_p_p;   // ksb5''
    double ksbud;      // ksbud
    double ksc1_p;     // ksc1'
    double ksc1_p_p;   // ksc1''
    double kscdh;      // kscdh
    double ksf6_p;     // ksf6'
    double ksf6_p_p;   // ksf6''
    double ksf6_p_p_p; // ksf6'''
    double ksn2_p;     // ksn2'
    double ksn2_p_p;   // ksn2''
    double ksnet;      // ksnet
    double ksori;      // ksori
    double kspds_p;    // kspds'
    double ksppx;      // ksppx
    double ksspn;      // ksspn
    double ksswi_p;    // ksswi'
    double ksswi_p_p;  // ksswi''
    double lte1h;      // lte1h
    double lte1l;      // lte1l
    double mad2h;      // mad2h
    double mad2l;      // mad2l
    double mdt;        // mdt
    double TEM1T;      // TEM1T

    // STATE VARIABLES
    double BUD;     // BUD
    double C2;      // C2
    double C2P;     // C2P
    double C5;      // C5
    double C5P;     // C5P
    double CDC14;   // CDC14
    double CDC15;   // CDC15
    double CDC20;   // CDC20
    double CDC20i;  // CDC20i
    double CDC6;    // CDC6
    double CDC6P;   // CDC6P
    double CDH1;    // CDH1
    double CDH1i;   // CDH1i
    double CLB2;    // CLB2
    double CLB5;    // CLB5
    double CLN2;    // CLN2
    double ESP1;    // ESP1
    double F2;      // F2
    double F2P;     // F2P
    double F5;      // F5
    double F5P;     // F5P
    double IEP;     // IEP
    double MASS;    // MASS
    double NET1;    // NET1
    double NET1P;   // NET1P
    double ORI;     // ORI
    double PDS1;    // PDS1
    double PPX;     // PPX
    double RENT;    // RENT
    double RENTP;   // RENTP
    double SIC1;    // SIC1
    double SIC1P;   // SIC1P
    double SPN;     // SPN
    double SWI5;    // SWI5
    double SWI5P;   // SWI5P
    double TEM1GTP; // TEM1GTP

    double d_BUD_dt;
    double d_C2_dt;
    double d_C2P_dt;
    double d_C5_dt;
    double d_C5P_dt;
    double d_CDC14_dt;
    double d_CDC15_dt;
    double d_CDC20_dt;
    double d_CDC20i_dt;
    double d_CDC6_dt;
    double d_CDC6P_dt;
    double d_CDH1_dt;
    double d_CDH1i_dt;
    double d_CLB2_dt;
    double d_CLB5_dt;
    double d_CLN2_dt;
    double d_ESP1_dt;
    double d_F2_dt;
    double d_F2P_dt;
    double d_F5_dt;
    double d_F5P_dt;
    double d_IEP_dt;
    double d_MASS_dt;
    double d_NET1_dt;
    double d_NET1P_dt;
    double d_ORI_dt;
    double d_PDS1_dt;
    double d_PPX_dt;
    double d_RENT_dt;
    double d_RENTP_dt;
    double d_SIC1_dt;
    double d_SIC1P_dt;
    double d_SPN_dt;
    double d_SWI5_dt;
    double d_SWI5P_dt;
    double d_TEM1GTP_dt;

    // DERIVED QUANTITIES
    double cell;    // cell
    double BCK2;    // BCK2
    double BUB2;    // BUB2
    double CDC14T;  // CDC14T
    double CDC15i;  // CDC15i
    double CDC6T;   // CDC6T
    double CKIT;    // CKIT
    double CLB2T;   // CLB2T
    double CLB5T;   // CLB5T
    double CLN3;    // CLN3
    double IE;      // IE
    double LTE1;    // LTE1
    double MAD2;    // MAD2
    double MCM1;    // MCM1
    double NET1T;   // NET1T
    double PE;      // PE
    double SBF;     // SBF
    double SIC1T;   // SIC1T
    double TEM1GDP; // TEM1GDP
    double D;       // D
    double mu;      // mu
    double Vdb5;    // Vdb5
    double Vdb2;    // Vdb2
    double Vasbf;   // Vasbf
    double Visbf;   // Visbf
    double Vkpc1;   // Vkpc1
    double Vkpf6;   // Vkpf6
    double Vacdh;   // Vacdh
    double Vicdh;   // Vicdh
    double Vppnet;  // Vppnet
    double Vkpnet;  // Vkpnet
    double Vdppx;   // Vdppx
    double Vdpds;   // Vdpds
    double Vaiep;   // Vaiep
    double Vd2c1;   // Vd2c1
    double Vd2f6;   // Vd2f6
    double Vppc1;   // Vppc1
    double Vppf6;   // Vppf6
    double F;       // F

    // STOICHIOMETRY VARIABLES

    // REACTIONS
    double Growth;                               // Growth
    double Synthesis_of_CLN2;                    // Synthesis of CLN2
    double Degradation_of_CLN2;                  // Degradation of CLN2
    double Synthesis_of_CLB2;                    // Synthesis of CLB2
    double Degradation_of_CLB2;                  // Degradation of CLB2
    double Synthesis_of_CLB5;                    // Synthesis of CLB5
    double Degradation_of_CLB5;                  // Degradation of CLB5
    double Synthesis_of_SIC1;                    // Synthesis of SIC1
    double Phosphorylation_of_SIC1;              // Phosphorylation of SIC1
    double Dephosphorylation_of_SIC1;            // Dephosphorylation of SIC1
    double Fast_Degradation_of_SIC1P;            // Fast Degradation of SIC1P
    double Assoc_of_CLB2_and_SIC1;               // Assoc. of CLB2 and SIC1
    double Dissoc_of_CLB2SIC1_complex;           // Dissoc. of CLB2/SIC1 complex
    double Assoc_of_CLB5_and_SIC1;               // Assoc. of CLB5 and SIC1
    double Dissoc_of_CLB5SIC1;                   // Dissoc. of CLB5/SIC1
    double Phosphorylation_of_C2;                // Phosphorylation of C2
    double Dephosphorylation_of_C2P;             // Dephosphorylation of C2P
    double Phosphorylation_of_C5;                // Phosphorylation of C5
    double Dephosphorylation_of_C5P;             // Dephosphorylation of C5P
    double Degradation_of_CLB2_in_C2;            // Degradation of CLB2 in C2
    double Degradation_of_CLB5_in_C5;            // Degradation of CLB5 in C5
    double Degradation_of_SIC1_in_C2P;           // Degradation of SIC1 in C2P
    double Degradation_of_SIC1P_in_C5P_;         // Degradation of SIC1P in C5P
    double Degradation_of_CLB2_in_C2P;           // Degradation of CLB2 in C2P
    double Degradation_of_CLB5_in_C5P;           // Degradation of CLB5 in C5P
    double CDC6_synthesis;                       // CDC6 synthesis
    double Phosphorylation_of_CDC6;              // Phosphorylation of CDC6
    double Dephosphorylation_of_CDC6;            // Dephosphorylation of CDC6
    double Degradation_of_CDC6P;                 // Degradation of CDC6P
    double CLB2CDC6_complex_formation;           // CLB2/CDC6 complex formation
    double CLB2CDC6_dissociation;                // CLB2/CDC6 dissociation
    double CLB5CDC6_complex_formation;           // CLB5/CDC6 complex formation
    double CLB5CDC6_dissociation;                // CLB5/CDC6 dissociation
    double F2_phosphorylation;                   // F2 phosphorylation
    double F2P_dephosphorylation;                // F2P dephosphorylation
    double F5_phosphorylation;                   // F5 phosphorylation
    double F5P_dephosphorylation;                // F5P dephosphorylation
    double CLB2_degradation_in_F2;               // CLB2 degradation in F2
    double CLB5_degradation_in_F5;               // CLB5 degradation in F5
    double CDC6_degradation_in_F2P;              // CDC6 degradation in F2P
    double CDC6_degradation_in_F5P;              // CDC6 degradation in F5P
    double CLB2_degradation_in_F2P;              // CLB2 degradation in F2P
    double CLB5_degradation_in_F5P;              // CLB5 degradation in F5P
    double Synthesis_of_SWI5;                    // Synthesis of SWI5
    double Degradation_of_SWI5;                  // Degradation of SWI5
    double Degradation_of_SWI5P;                 // Degradation of SWI5P
    double Activation_of_SWI5;                   // Activation of SWI5
    double Inactivation_of_SWI5;                 // Inactivation of SWI5
    double Activation_of_IEP;                    // Activation of IEP
    double Inactivation_1;                       // Inactivation
    double Synthesis_of_inactive_CDC20;          // Synthesis of inactive CDC20
    double Degradation_of_inactiveCDC20;         // Degradation of inactiveCDC20
    double Degradation_of_active_CDC20;          // Degradation of active CDC20
    double Activation_of_CDC20;                  // Activation of CDC20
    double Inactivation_2;                       // Inactivation
    double CDH1_synthesis;                       // CDH1 synthesis
    double CDH1_degradation;                     // CDH1 degradation
    double CDH1i_degradation;                    // CDH1i degradation
    double CDH1i_activation;                     // CDH1i activation
    double Inactivation_3;                       // Inactivation
    double CDC14_synthesis;                      // CDC14 synthesis
    double CDC14_degradation;                    // CDC14 degradation
    double Assoc_with_NET1_to_form_RENT;         // Assoc. with NET1 to form RENT
    double Dissoc_from_RENT;                     // Dissoc. from RENT
    double Assoc_with_NET1P_to_form_RENTP;       // Assoc with NET1P to form RENTP
    double Dissoc_from_RENP;                     // Dissoc. from RENP
    double Net1_synthesis;                       // Net1 synthesis
    double Net1_degradation;                     // Net1 degradation
    double Net1P_degradation;                    // Net1P degradation
    double NET1_phosphorylation;                 // NET1 phosphorylation
    double dephosphorylation_1;                  // dephosphorylation
    double RENT_phosphorylation;                 // RENT phosphorylation
    double dephosphorylation_2;                  // dephosphorylation
    double Degradation_of_NET1_in_RENT;          // Degradation of NET1 in RENT
    double Degradation_of_NET1P_in_RENTP;        // Degradation of NET1P in RENTP
    double Degradation_of_CDC14_in_RENT;         // Degradation of CDC14 in RENT
    double Degradation_of_CDC14_in_RENTP;        // Degradation of CDC14 in RENTP
    double TEM1_activation;                      // TEM1 activation
    double inactivation_1;                       // inactivation
    double CDC15_activation;                     // CDC15 activation
    double inactivation_2;                       // inactivation
    double PPX_synthesis;                        // PPX synthesis
    double degradation_1;                        // degradation
    double PDS1_synthesis;                       // PDS1 synthesis
    double degradation_2;                        // degradation
    double Degradation_of_PDS1_in_PE;            // Degradation of PDS1 in PE
    double Assoc_with_ESP1_to_form_PE;           // Assoc. with ESP1 to form PE
    double Disso_from_PE;                        // Disso. from PE
    double DNA_synthesis;                        // DNA synthesis
    double Negative_regulation_of_DNA_synthesis; // Negative regulation of DNA synthesis
    double Budding;                              // Budding
    double Negative_regulation_of_Cell_budding;  // Negative regulation of Cell budding
    double Spindle_formation;                    // Spindle formation
    double Spindle_disassembly;                  // Spindle disassembly

    /**
     * Process the events in the model.
     *
     * @param time The current time
     * @param rY The current state variables
     *
     * @return How close we are to the time of the next event
     */
    double ProcessModelEvents(double time, const std::vector<double>& rY) override;

    /**
     * Run the assignment rules to update state.
     *
     * @param time The current time
     * @param rY The current state variables
     */
    void RunAssignmentRules(double time) override;

    /**
     * Run the initial assignments to set initial state.
     *
     * @param time The current time
     */
    void RunInitialAssignments(double time) override;

    /**
     * Run the reactions to update state.
     *
     * @param time The current time
     * @param rY The current state variables
     */
    void RunReactions(double time) override;

    /**
     * Update variable parameters from current ODE system parameter settings.
     *
     * @param time The current time
     */
    void UpdateParameters(double time) override;

    /**
     * Update state variables from the given ODE system state.
     *
     * @param time The current time
     * @param rStateVariables The state variables to use
     */
    void UpdateStateVariables(double time, const std::vector<double>& rStateVariables) override;

    // MODEL FUNCTIONS
    inline double BB_218(double A1, double A2, double A3, double A4);
    inline double GK_219(double A1, double A2, double A3, double A4);
    inline double MichaelisMenten_220(double M1, double J1, double k1, double S1);
    inline double Mass_Action_2_221(double k1, double S1, double S2);
    inline double Mass_Action_1_222(double k1, double S1);

public:
    /**
     * Default constructor
     */
    Chen2004SbmlOdeSystem();

    /**
     * Destructor
     */
    ~Chen2004SbmlOdeSystem();

    /**
     * Compute the derived quantities from the given system state.
     *
     * @param time  the time at which to compute the derived quantities
     * @param rY a vector of values for the state variables
     *
     * @return a vector of derived quantities
     */
    std::vector<double> ComputeDerivedQuantities(double time, const std::vector<double>& rY) override;

    /**
     * Compute the RHS of the ODE system.
     *
     * An ODE solver will call this function repeatedly to solve for y = [y1 ... yn].
     *
     * @param time the time used to evaluate the RHS.
     * @param rY an input solution vector used to evaluate the RHS.
     * @param rDY an output vector to be filled in with the resulting derivatives y' = [y1' ... yn'].
     */
    void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY) override;

    void Initialise(double time = 0.0);

    std::vector<double> RunModelEquations(double time, const std::vector<double>& rStateVariables);
};

// Register the ODE system with Boost serialization
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(Chen2004SbmlOdeSystem)

#endif // CHEN_2004_SBML_ODE_SYSTEM_HPP_