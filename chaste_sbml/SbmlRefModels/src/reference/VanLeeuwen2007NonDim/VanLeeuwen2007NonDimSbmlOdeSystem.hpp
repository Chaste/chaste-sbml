#ifndef VAN_LEEUWEN_2007_NON_DIM_SBML_ODE_SYSTEM_HPP_
#define VAN_LEEUWEN_2007_NON_DIM_SBML_ODE_SYSTEM_HPP_

#include <vector>

#include <boost/serialization/base_object.hpp>

#include "AbstractSbmlOdeSystem.hpp"
#include "ChasteSerialization.hpp"
#include "SbmlEventType.hpp"

class VanLeeuwen2007NonDimSbmlOdeSystem : public AbstractSbmlOdeSystem
{
private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Save / load VanLeeuwen2007NonDimSbmlOdeSystem archive
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
    double K_T;                     // K_T
    double K_C;                     // K_C
    double K_D;                     // K_D
    double p_u;                     // p_u
    double wnt_level;               // wnt_level
    double xi_D;                    // xi_D
    double xi_Dx;                   // xi_Dx
    double xi_X;                    // xi_X
    double xi_C;                    // xi_C
    double d_Dx;                    // d_Dx
    double d_X;                     // d_X
    double d_D;                     // d_D
    double p_c;                     // p_c
    double gamma1;                  // gamma1
    double gamma2;                  // gamma2
    double s_D;                     // s_D
    double s_X;                     // s_X
    double d_u;                     // d_u
    double s_c;                     // s_c
    double d_c;                     // d_c
    double s_CA;                    // s_CA
    double d_CA;                    // d_CA
    double s_CT;                    // s_CT
    double d_CT;                    // d_CT
    double s_A;                     // s_A
    double d_A;                     // d_A
    double s_T;                     // s_T
    double d_T;                     // d_T
    double s_Y;                     // s_Y
    double d_Y;                     // d_Y
    double ComplexTransitThreshold; // ComplexTransitThreshold

    // STATE VARIABLES
    double X;    // X
    double D;    // D
    double C_o;  // C_o
    double C_u;  // C_u
    double C_c;  // C_c
    double A;    // A
    double C_A;  // C_A
    double T;    // T
    double C_oT; // C_oT
    double C_cT; // C_cT
    double Y;    // Y

    double d_X_dt;
    double d_D_dt;
    double d_C_o_dt;
    double d_C_u_dt;
    double d_C_c_dt;
    double d_A_dt;
    double d_C_A_dt;
    double d_T_dt;
    double d_C_oT_dt;
    double d_C_cT_dt;
    double d_Y_dt;

    // DERIVED QUANTITIES
    double cytosolmembraneandnucleus; // cytosolmembraneandnucleus
    double C_F;                       // C_F
    double C_T;                       // C_T
    double drag;                      // drag

    // REACTIONS
    double mwcfbf7716_cc13_473c_979a_033c57a28857; // r1
    double mwab347951_a286_432d_b03b_254dcdba4a25; // r2
    double mwc497befc_0edd_4b0d_8895_77dbfa05a4be; // r22
    double mwfa0f5940_d4e9_44e1_8a0c_379d6fdf9d0f; // r23
    double mwbe108cd1_a5d9_4d55_97c9_ac584df6a005; // r7
    double mwad0ba91a_37fb_48f5_8392_6c043dfdd5dc; // r16
    double mw1d3d836b_77bd_489b_9a37_acc874344652; // r8
    double mw552840aa_6d3c_4188_850d_9a3caafcdaa5; // r3
    double mw5a301d6f_ffa1_4362_98de_0919a68808bc; // r4
    double mw4b073704_f41a_43d1_8ff8_af03f1cb6e1d; // r17
    double mw1c7459d0_b5a9_4a89_8682_79fc8ca4ca6c; // r9
    double mw51e5ed29_82db_47fa_9eed_ab52b7786dd5; // r10
    double mw8704b9b0_f93d_405a_89a1_e000c5b66896; // r11
    double mw322e8f78_68f6_4293_ae2a_ab238ec5cf9b; // r18
    double mwf5440ce3_586c_491f_ba71_bbf1c485027f; // r12
    double mw08c332d6_1aaf_498f_8e1c_13f1b9da85bf; // r19
    double mw3ece1442_a0b3_47eb_a8bf_cb317a46165f; // r15
    double mwfb6f8dfa_2e07_4249_a8fe_af33ca27471e; // r5
    double mwaa4e7692_3290_45ab_8a38_b5ffd49ede87; // r6
    double mw86ea3c28_4745_4a89_9cbd_0f23c603e01b; // r20
    double mwee62535f_931d_41e9_ad53_2cad54b94778; // r21
    double mw1b250e04_290a_4689_9fc5_5f0bf6711e02; // r13
    double mw925599eb_19a0_4434_8be3_67c40721b71d; // r14
    double mw321b3e5e_f6ed_4345_9346_55ffb1ff2137; // r24

    /**
     * Process the events in the model.
     *
     * @param time The current time
     * @param rY The current state variables
     *
     * @return How close we are to the time of the next event
     */
    double ProcessModelEvents(double time, const std::vector<double>& rY) override;

    // MODEL FUNCTIONS

public:
    /**
     * Default constructor
     */
    VanLeeuwen2007NonDimSbmlOdeSystem();

    /**
     * Destructor
     */
    ~VanLeeuwen2007NonDimSbmlOdeSystem();

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
CHASTE_CLASS_EXPORT(VanLeeuwen2007NonDimSbmlOdeSystem)

#endif // VAN_LEEUWEN_2007_NON_DIM_SBML_ODE_SYSTEM_HPP_