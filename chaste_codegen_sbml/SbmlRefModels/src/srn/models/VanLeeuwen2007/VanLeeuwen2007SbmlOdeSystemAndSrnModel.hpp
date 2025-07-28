#ifndef VANLEEUWEN2007SBMLODESYSTEMANDSRNMODEL_HPP_
#define VANLEEUWEN2007SBMLODESYSTEMANDSRNMODEL_HPP_

#include <vector>

#include <boost/serialization/base_object.hpp>

#include "AbstractOdeSystem.hpp"
#include "ChasteSerialization.hpp"

class VanLeeuwen2007SbmlOdeSystem : public AbstractOdeSystem
{
private:
    // (De-)serialize VanLeeuwen2007SbmlOdeSystem
    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar& BOOST_SERIALIZATION_BASE_OBJECT_NVP(AbstractOdeSystem);
    }

    // CONSTANT PARAMETERS
    const double K_T = 50.0;   // K_T
    const double K_C = 200.0;  // K_C
    const double K_D = 5.0;    // K_D
    const double p_u = 100.0;  // p_u
    const double xi_D = 5.0;   // xi_D
    const double xi_Dx = 5.0;  // xi_Dx
    const double xi_X = 200.0; // xi_X
    const double xi_C = 0.0;   // xi_C
    const double s_D = 100.0;  // s_D
    const double d_Dx = 5.0;   // d_Dx
    const double s_X = 10.0;   // s_X
    const double d_X = 100.0;  // d_X
    const double d_u = 50.0;   // d_u
    const double s_c = 25.0;   // s_c
    const double d_c = 1.0;    // d_c
    const double s_CA = 250.0; // s_CA
    const double d_CA = 350.0; // d_CA
    const double s_CT = 30.0;  // s_CT
    const double d_CT = 750.0; // d_CT
    const double p_c = 0.0;    // p_c
    const double s_A = 20.0;   // s_A
    const double d_A = 2.0;    // d_A
    const double s_T = 10.0;   // s_T
    const double d_T = 0.4;    // d_T
    const double s_Y = 10.0;   // s_Y
    const double d_Y = 1.0;    // d_Y
    const double d_D = 5.0;    // d_D

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

public:
    VanLeeuwen2007SbmlOdeSystem(std::vector<double> stateVariables = std::vector<double>());

    ~VanLeeuwen2007SbmlOdeSystem();

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
    // Provide constructor for serializing VanLeeuwen2007SbmlOdeSystem
    template <class Archive>
    inline void save_construct_data(Archive& ar, const VanLeeuwen2007SbmlOdeSystem* t, const unsigned int version)
    {
        // Save data required to construct instance
        const std::vector<double> state_variables = t->rGetConstStateVariables();
        ar << state_variables;
    }

    // Provide constructor for de-serializing VanLeeuwen2007SbmlOdeSystem
    template <class Archive>
    inline void load_construct_data(Archive& ar, VanLeeuwen2007SbmlOdeSystem* t, const unsigned int version)
    {
        // Retrieve data from archive required to construct new instance
        std::vector<double> state_variables;
        ar >> state_variables;

        // Invoke inplace constructor to initialise instance
        ::new (t) VanLeeuwen2007SbmlOdeSystem(state_variables);
    }
} // namespace serialization
} // namespace

// Define SbmlSrnWrapperModel using wrappers
#include "SbmlSrnWrapperModel.cpp"
#include "SbmlSrnWrapperModel.hpp"

typedef SbmlSrnWrapperModel<VanLeeuwen2007SbmlOdeSystem, 11> VanLeeuwen2007SbmlSrnModel;

// Declare identifiers for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(VanLeeuwen2007SbmlOdeSystem)
EXPORT_TEMPLATE_CLASS2(SbmlSrnWrapperModel, VanLeeuwen2007SbmlOdeSystem, 11)

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(VanLeeuwen2007SbmlSrnModel)

#endif // VANLEEUWEN2007SBMLODESYSTEMANDSRNMODEL_HPP_