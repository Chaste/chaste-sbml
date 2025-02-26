#ifndef TAN2014ODESYSTEMANDSRNMODEL_HPP_
#define TAN2014ODESYSTEMANDSRNMODEL_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/shared_ptr.hpp>

#include <cmath>
#include <iostream>
#include "AbstractOdeSystem.hpp"

class Tan2014OdeSystem : public AbstractOdeSystem
{
private:
    /* Initialise compartments and values. */
    double compartment;
    double CytosolMembrane;
    double nucleus;

    /* Initialise model parameters. */
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

    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive &archive, const unsigned int version)
    {
        archive &boost::serialization::base_object<AbstractOdeSystem>(*this);
    }

public:
    /* Default constructor. */
    Tan2014OdeSystem(std::vector<double> stateVariables = std::vector<double>());

    /* Destructor. */
    ~Tan2014OdeSystem();

    void Init();

    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY);
};

namespace
{
    namespace serialization
    {
        /* Serialize information required to construct a Tan2014OdeSystem. */
        template <class Archive>
        inline void save_construct_data(
            Archive &ar, const Tan2014OdeSystem *t, const unsigned int file_version)
        {
            const std::vector<double> state_variables = t->rGetConstStateVariables();
            ar & state_variables;
        }
        /* De-serialize constructor parameters and intiialise a Tan2014OdeSystem. */
        template <class Archive>
        inline void load_construct_data(
            Archive &ar, Tan2014OdeSystem *t, const unsigned int file_version)
        {
            std::vector<double> state_variables;
            ar & state_variables;

            // Invoke inplace constructor to initialise instance
            ::new (t) Tan2014OdeSystem(state_variables);
        }
    }
} // namespace ...

/* Define SRN model using Wrappers. */
#include "SbmlSrnWrapperModel.hpp"
#include "SbmlSrnWrapperModel.cpp"

typedef SbmlSrnWrapperModel<Tan2014OdeSystem, 7> Tan2014SrnModel;

// Declare identifiers for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(Tan2014OdeSystem)
EXPORT_TEMPLATE_CLASS2(SbmlSrnWrapperModel, Tan2014OdeSystem, 7)

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(Tan2014SrnModel)

#endif /* TAN2014ODESYSTEMANDSRNMODEL_HPP_ */