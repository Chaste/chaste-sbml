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
    void serialize(Archive &ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(AbstractOdeSystem);
    }

    // COMPARTMENTS:
    double cytosolmembraneandnucleus;

    // CONST PARAMETERS:
    const double K_T = 1500.0;
    const double K_C = 6000.0;
    const double K_D = 150.0;
    const double p_u = 1.0;
    const double xi_D = 0.00666666666666667;
    const double xi_Dx = 0.00666666666666667;
    const double xi_X = 0.266666666666667;
    const double xi_C = 0.0;
    const double d_Dx = 0.00666666666666667;
    const double d_X = 0.133333333333333;
    const double d_D = 0.00666666666666667;
    const double p_c = 0.0;
    const double s_D = 0.133333333333333;
    const double s_X = 0.0533333333333333;
    const double d_u = 0.0666666666666667;
    const double s_c = 1.0;
    const double d_c = 0.00133333333333333;
    const double s_CA = 0.0111111111111111;
    const double d_CA = 0.466666666666667;
    const double s_CT = 0.00133333333333333;
    const double d_CT = 1.0;
    const double s_A = 0.8;
    const double d_A = 0.00266666666666667;
    const double s_T = 0.4;
    const double d_T = 0.000533333333333333;
    const double s_Y = 1.0;
    const double d_Y = 0.00133333333333333;

    // STATE VARIABLES:
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
    double C_F; // C_F
    double C_T; // C_T
    double drag; // drag

    // STATE PARAMETERS:

    double wnt_level;
    double gamma1;
    double gamma2;
    double ComplexTransitThreshold;


public:
    VanLeeuwen2007NonDimSbmlOdeSystem(std::vector<double> stateVariables = std::vector<double>());

    ~VanLeeuwen2007NonDimSbmlOdeSystem();

    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY);
    void RefreshState(const std::vector<double> &rY);


    // FUNCTION DEFINITIONS:
};

namespace
{
    namespace serialization
    {
        // Provide constructor for serializing VanLeeuwen2007NonDimSbmlOdeSystem
        template <class Archive>
        inline void save_construct_data(Archive &ar, const VanLeeuwen2007NonDimSbmlOdeSystem * t, const unsigned int version)
        {
            // Save data required to construct instance
            const std::vector<double> state_variables = t->rGetConstStateVariables();
            ar << state_variables;
        }

        // Provide constructor for de-serializing VanLeeuwen2007NonDimSbmlOdeSystem
        template <class Archive>
        inline void load_construct_data(Archive &ar, VanLeeuwen2007NonDimSbmlOdeSystem * t, const unsigned int version)
        {
            // Retrieve data from archive required to construct new instance
            std::vector<double> state_variables;
            ar >> state_variables;

            // Invoke inplace constructor to initialise instance
            ::new (t)VanLeeuwen2007NonDimSbmlOdeSystem(state_variables);
        }
    } // namespace serialization
} // namespace ...

// Define SbmlSrnWrapperModel using wrappers
#include "SbmlSrnWrapperModel.hpp"
#include "SbmlSrnWrapperModel.cpp"

typedef SbmlSrnWrapperModel<VanLeeuwen2007NonDimSbmlOdeSystem, 14> VanLeeuwen2007NonDimSbmlSrnModel;

// Declare identifiers for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(VanLeeuwen2007NonDimSbmlOdeSystem)
EXPORT_TEMPLATE_CLASS2(SbmlSrnWrapperModel, VanLeeuwen2007NonDimSbmlOdeSystem, 14)

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(VanLeeuwen2007NonDimSbmlSrnModel)

#endif // VANLEEUWEN2007NONDIMSBMLODESYSTEMANDSRNMODEL_HPP_