#ifndef VANLEEUWEN2007NONDIMSBMLODESYSTEMANDSRNMODEL_HPP_
#define VANLEEUWEN2007NONDIMSBMLODESYSTEMANDSRNMODEL_HPP_

#include <vector>

#include <boost/serialization/base_object.hpp>

#include "AbstractOdeSystem.hpp"
#include "ChasteSerialization.hpp"

class VanLeeuwen2007NonDimSbmlOdeSystem : public AbstractOdeSystem
{
private:
    // (De-)serialize VanLeeuwen2007NonDimSbmlOdeSystem
    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar& BOOST_SERIALIZATION_BASE_OBJECT_NVP(AbstractOdeSystem);
    }

    // CONSTANT PARAMETERS
    const double K_T = 1500.0;                // K_T
    const double K_C = 6000.0;                // K_C
    const double K_D = 150.0;                 // K_D
    const double p_u = 1.0;                   // p_u
    const double xi_D = 0.00666666666666667;  // xi_D
    const double xi_Dx = 0.00666666666666667; // xi_Dx
    const double xi_X = 0.266666666666667;    // xi_X
    const double xi_C = 0.0;                  // xi_C
    const double d_Dx = 0.00666666666666667;  // d_Dx
    const double d_X = 0.133333333333333;     // d_X
    const double d_D = 0.00666666666666667;   // d_D
    const double p_c = 0.0;                   // p_c
    const double s_D = 0.133333333333333;     // s_D
    const double s_X = 0.0533333333333333;    // s_X
    const double d_u = 0.0666666666666667;    // d_u
    const double s_c = 1.0;                   // s_c
    const double d_c = 0.00133333333333333;   // d_c
    const double s_CA = 0.0111111111111111;   // s_CA
    const double d_CA = 0.466666666666667;    // d_CA
    const double s_CT = 0.00133333333333333;  // s_CT
    const double d_CT = 1.0;                  // d_CT
    const double s_A = 0.8;                   // s_A
    const double d_A = 0.00266666666666667;   // d_A
    const double s_T = 0.4;                   // s_T
    const double d_T = 0.000533333333333333;  // d_T
    const double s_Y = 1.0;                   // s_Y
    const double d_Y = 0.00133333333333333;   // d_Y

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

    // DERIVED QUANTITIES
    double C_F;  // C_F
    double C_T;  // C_T
    double drag; // drag

    // VARIABLE PARAMETERS
    double cytosolmembraneandnucleus; // cytosolmembraneandnucleus
    double wnt_level;                 // wnt_level
    double gamma1;                    // gamma1
    double gamma2;                    // gamma2
    double ComplexTransitThreshold;   // ComplexTransitThreshold

    // RULE-BASED PARAMETERS

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

public:
    VanLeeuwen2007NonDimSbmlOdeSystem(std::vector<double> stateVariables = std::vector<double>());

    ~VanLeeuwen2007NonDimSbmlOdeSystem();

    void AdjustOdeParameters(double time);

    void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY);
    void ProcessRules(double time, const std::vector<double>& rY);

    std::vector<double> ComputeDerivedQuantities(double time, const std::vector<double>& rY);

    // FUNCTIONS
};

namespace
{
namespace serialization
{
    // Provide constructor for serializing VanLeeuwen2007NonDimSbmlOdeSystem
    template <class Archive>
    inline void save_construct_data(Archive& ar, const VanLeeuwen2007NonDimSbmlOdeSystem* t, const unsigned int version)
    {
        // Save data required to construct instance
        const std::vector<double> state_variables = t->rGetConstStateVariables();
        ar << state_variables;
    }

    // Provide constructor for de-serializing VanLeeuwen2007NonDimSbmlOdeSystem
    template <class Archive>
    inline void load_construct_data(Archive& ar, VanLeeuwen2007NonDimSbmlOdeSystem* t, const unsigned int version)
    {
        // Retrieve data from archive required to construct new instance
        std::vector<double> state_variables;
        ar >> state_variables;

        // Invoke inplace constructor to initialise instance
        ::new (t) VanLeeuwen2007NonDimSbmlOdeSystem(state_variables);
    }
} // namespace serialization
} // namespace

// Define SbmlSrnWrapperModel using wrappers
#include "SbmlSrnWrapperModel.cpp"
#include "SbmlSrnWrapperModel.hpp"

typedef SbmlSrnWrapperModel<VanLeeuwen2007NonDimSbmlOdeSystem, 11> VanLeeuwen2007NonDimSbmlSrnModel;

// Declare identifiers for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(VanLeeuwen2007NonDimSbmlOdeSystem)
EXPORT_TEMPLATE_CLASS2(SbmlSrnWrapperModel, VanLeeuwen2007NonDimSbmlOdeSystem, 11)

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(VanLeeuwen2007NonDimSbmlSrnModel)

#endif // VANLEEUWEN2007NONDIMSBMLODESYSTEMANDSRNMODEL_HPP_