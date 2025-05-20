#ifndef TAN2014SBMLODESYSTEMANDSRNMODEL_HPP_
#define TAN2014SBMLODESYSTEMANDSRNMODEL_HPP_

#include <cmath>
#include <iostream>

#include <boost/serialization/base_object.hpp>

#include "AbstractOdeSystem.hpp"
#include "ChasteSerialization.hpp"

class Tan2014SbmlOdeSystem : public AbstractOdeSystem
{
private:
    // (De-)serialize Tan2014SbmlOdeSystem
    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(AbstractOdeSystem);
    }

    // COMPARTMENTS:
    double compartment;
    double CytosolMembrane;
    double nucleus;

    // PARAMETERS:
    double wnt_level;
    double kdegradation;
    double Bsyn;
    double K_n_active_k;
    double kC_k1;
    double kC_k2;
    double kN_k1;
    double kN_k2;
    double kdiffusion_k;
    double K_c_active_k;
    double gamma;
    double ComplexTransitThreshold;


public:
    Tan2014SbmlOdeSystem(std::vector<double> stateVariables = std::vector<double>());

    ~Tan2014SbmlOdeSystem();

    void Init();

    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY);


    // FUNCTION DEFINITIONS:
};

namespace
{
    namespace serialization
    {
        // Provide constructor for serializing Tan2014SbmlOdeSystem
        template <class Archive>
        inline void save_construct_data(Archive &ar, const Tan2014SbmlOdeSystem * t, const unsigned int version)
        {
            // Save data required to construct instance
            const std::vector<double> state_variables = t->rGetConstStateVariables();
            ar << state_variables;
        }

        // Provide constructor for de-serializing Tan2014SbmlOdeSystem
        template <class Archive>
        inline void load_construct_data(Archive &ar, Tan2014SbmlOdeSystem * t, const unsigned int version)
        {
            // Retrieve data from archive required to construct new instance
            std::vector<double> state_variables;
            ar >> state_variables;

            // Invoke inplace constructor to initialise instance
            ::new (t)Tan2014SbmlOdeSystem(state_variables);
        }
    } // namespace serialization
} // namespace ...

// Define SbmlSrnWrapperModel using wrappers
#include "SbmlSrnWrapperModel.hpp"
#include "SbmlSrnWrapperModel.cpp"

typedef SbmlSrnWrapperModel<Tan2014SbmlOdeSystem, 6> Tan2014SbmlSrnModel;

// Declare identifiers for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(Tan2014SbmlOdeSystem)
EXPORT_TEMPLATE_CLASS2(SbmlSrnWrapperModel, Tan2014SbmlOdeSystem, 6)

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(Tan2014SbmlSrnModel)

#endif // TAN2014SBMLODESYSTEMANDSRNMODEL_HPP_