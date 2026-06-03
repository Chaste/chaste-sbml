#ifndef VAN_LEEUWEN_2007_SBML_ODE_SYSTEM_HPP_
#define VAN_LEEUWEN_2007_SBML_ODE_SYSTEM_HPP_

#include <vector>

#include <boost/serialization/base_object.hpp>

#include "AbstractSbmlOdeSystem.hpp"
#include "ChasteSerialization.hpp"
#include "SbmlEventType.hpp"

class VanLeeuwen2007SbmlOdeSystem : public AbstractSbmlOdeSystem
{
private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Save / load VanLeeuwen2007SbmlOdeSystem archive
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template <class Archive>
    void serialize(Archive &archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractSbmlOdeSystem>(*this);
    }

    // PARAMETERS
    double K_T; // K_T
    double K_C; // K_C
    double K_D; // K_D
    double p_u; // p_u
    double wnt_level; // wnt_level
    double gamma1; // gamma1
    double gamma2; // gamma2
    double xi_D; // xi_D
    double xi_Dx; // xi_Dx
    double xi_X; // xi_X
    double xi_C; // xi_C
    double s_D; // s_D
    double d_Dx; // d_Dx
    double s_X; // s_X
    double d_X; // d_X
    double d_u; // d_u
    double s_c; // s_c
    double d_c; // d_c
    double s_CA; // s_CA
    double d_CA; // d_CA
    double s_CT; // s_CT
    double d_CT; // d_CT
    double p_c; // p_c
    double s_A; // s_A
    double d_A; // d_A
    double s_T; // s_T
    double d_T; // d_T
    double s_Y; // s_Y
    double d_Y; // d_Y
    double d_D; // d_D
    double ComplexTransitThreshold; // ComplexTransitThreshold

    // STATE VARIABLES
    double X; // X
    double D; // D
    double C_o; // C_o
    double C_u; // C_u
    double C_c; // C_c
    double A; // A
    double C_A; // C_A
    double T; // T
    double C_oT; // C_oT
    double C_cT; // C_cT
    double Y; // Y

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
    double C_F; // C_F
    double C_T; // C_T
    double drag; // drag

    // STOICHIOMETRY VARIABLES

    // REACTIONS
    double mwd6b35759_f098_484c_9c65_e84e7e4b61e4; // r1
    double mweddac6d0_231e_4c92_ba2a_c91edc682ff5; // r2
    double mwdf62dfed_ec88_4d81_bc9d_da0e10f41e4b; // r22
    double mwee9cc998_28e9_4173_a694_f3e278a639b7; // r23
    double mw661e341d_97d1_4e6f_8812_3be7ffc86d42; // r7
    double mw179aa33c_9a7e_43c0_9285_3d8f97719c60; // r16
    double mwff8d34f9_e036_49f1_b3b8_3706ecb98660; // r8
    double mwd4fa317a_a484_4f38_b9b0_6ed404d9adcb; // r3
    double mwcb88a249_a200_4e95_9185_5654bf1ebfc0; // r4
    double mw0aeac2fc_dd2f_4fe7_b60c_97f0b008ed89; // r17
    double mw4b47c66d_37e6_4c33_b043_1f6b3b814449; // r9
    double mw69974db4_8ead_416c_a220_f6dc3be1f3b6; // r10
    double mwccb628c3_76c8_47be_9d6b_6d3f1ae30dcc; // r11
    double mw581d69f1_60b3_4d21_9323_31b05ee89570; // r18
    double mwb17c2c57_279d_4e88_b9cf_896029135cc1; // r12
    double mwe3236fc5_2118_40cb_8db3_ef9da29137cf; // r19
    double mw0be4a28b_e9c6_43da_8f95_d9c564a7caae; // r15
    double mwc360befb_07da_4d19_bbec_523fbef47dc9; // r5
    double mw985b39ad_50ff_4d55_95f7_42e7f1bf6a3d; // r6
    double mwf27b65bb_38d6_44cf_aa66_9b38ae5885eb; // r20
    double mwd7a7bdbb_8f2e_42ed_bf97_48e33e45dafb; // r21
    double mw988a8caf_bd68_462b_86d7_51844c1dcfd3; // r13
    double mw9ab26a4c_bd70_45e0_bacc_f830ab28abca; // r14
    double mw931baf8f_6572_46f6_96eb_cae40ee267b7; // r24

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

public:
    /** 
     * Default constructor
     */
    VanLeeuwen2007SbmlOdeSystem();

    /**
     * Destructor
     */
    ~VanLeeuwen2007SbmlOdeSystem();
    
    /**
     * Compute the derived quantities from the given system state.
     *
     * @param time  the time at which to compute the derived quantities
     * @param rY a vector of values for the state variables
     * 
     * @return a vector of derived quantities
     */
    std::vector<double> ComputeDerivedQuantities(double time, const std::vector<double> &rY) override;

    /**
     * Compute the RHS of the ODE system.
     * 
     * An ODE solver will call this function repeatedly to solve for y = [y1 ... yn].
     *
     * @param time the time used to evaluate the RHS.
     * @param rY an input solution vector used to evaluate the RHS.
     * @param rDY an output vector to be filled in with the resulting derivatives y' = [y1' ... yn'].
     */
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY) override;

    void Initialise(double time = 0.0);

    std::vector<double> RunModelEquations(double time, const std::vector<double>& rStateVariables);
};

// Register the ODE system with Boost serialization
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(VanLeeuwen2007SbmlOdeSystem)

#endif // VAN_LEEUWEN_2007_SBML_ODE_SYSTEM_HPP_