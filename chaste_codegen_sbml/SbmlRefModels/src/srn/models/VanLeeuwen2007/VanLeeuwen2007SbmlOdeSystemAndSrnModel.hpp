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
    void serialize(Archive &ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(AbstractOdeSystem);
    }

    // COMPARTMENTS:
    double cytosolmembraneandnucleus;

    // CONST PARAMETERS:
    const double K_T = 50.0;
    const double K_C = 200.0;
    const double K_D = 5.0;
    const double p_u = 100.0;
    const double xi_D = 5.0;
    const double xi_Dx = 5.0;
    const double xi_X = 200.0;
    const double xi_C = 0.0;
    const double s_D = 100.0;
    const double d_Dx = 5.0;
    const double s_X = 10.0;
    const double d_X = 100.0;
    const double d_u = 50.0;
    const double s_c = 25.0;
    const double d_c = 1.0;
    const double s_CA = 250.0;
    const double d_CA = 350.0;
    const double s_CT = 30.0;
    const double d_CT = 750.0;
    const double p_c = 0.0;
    const double s_A = 20.0;
    const double d_A = 2.0;
    const double s_T = 10.0;
    const double d_T = 0.4;
    const double s_Y = 10.0;
    const double d_Y = 1.0;
    const double d_D = 5.0;

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
    VanLeeuwen2007SbmlOdeSystem(std::vector<double> stateVariables = std::vector<double>());

    ~VanLeeuwen2007SbmlOdeSystem();

    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY);
    void RefreshState(const std::vector<double> &rY);


    // FUNCTION DEFINITIONS:
};

namespace
{
    namespace serialization
    {
        // Provide constructor for serializing VanLeeuwen2007SbmlOdeSystem
        template <class Archive>
        inline void save_construct_data(Archive &ar, const VanLeeuwen2007SbmlOdeSystem * t, const unsigned int version)
        {
            // Save data required to construct instance
            const std::vector<double> state_variables = t->rGetConstStateVariables();
            ar << state_variables;
        }

        // Provide constructor for de-serializing VanLeeuwen2007SbmlOdeSystem
        template <class Archive>
        inline void load_construct_data(Archive &ar, VanLeeuwen2007SbmlOdeSystem * t, const unsigned int version)
        {
            // Retrieve data from archive required to construct new instance
            std::vector<double> state_variables;
            ar >> state_variables;

            // Invoke inplace constructor to initialise instance
            ::new (t)VanLeeuwen2007SbmlOdeSystem(state_variables);
        }
    } // namespace serialization
} // namespace ...

// Define SbmlSrnWrapperModel using wrappers
#include "SbmlSrnWrapperModel.hpp"
#include "SbmlSrnWrapperModel.cpp"

typedef SbmlSrnWrapperModel<VanLeeuwen2007SbmlOdeSystem, 14> VanLeeuwen2007SbmlSrnModel;

// Declare identifiers for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(VanLeeuwen2007SbmlOdeSystem)
EXPORT_TEMPLATE_CLASS2(SbmlSrnWrapperModel, VanLeeuwen2007SbmlOdeSystem, 14)

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(VanLeeuwen2007SbmlSrnModel)

#endif // VANLEEUWEN2007SBMLODESYSTEMANDSRNMODEL_HPP_