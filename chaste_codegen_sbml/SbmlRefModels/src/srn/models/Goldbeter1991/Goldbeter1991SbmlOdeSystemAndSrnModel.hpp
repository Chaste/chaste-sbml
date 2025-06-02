#ifndef GOLDBETER1991SBMLODESYSTEMANDSRNMODEL_HPP_
#define GOLDBETER1991SBMLODESYSTEMANDSRNMODEL_HPP_

#include <cmath>
#include <iostream>

#include <boost/serialization/base_object.hpp>

#include "AbstractOdeSystem.hpp"
#include "ChasteSerialization.hpp"

class Goldbeter1991SbmlOdeSystem : public AbstractOdeSystem
{
private:
    // (De-)serialize Goldbeter1991SbmlOdeSystem
    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(AbstractOdeSystem);
    }

    // COMPARTMENTS:
    double cell;

    // CONST PARAMETERS:
    const double VM1 = 3.0;
    const double VM3 = 1.0;
    const double Kc = 0.5;

    // STATE VARIABLES:
    double C; // Cyclin
    double M; // cdc_2_kinase
    double X; // Cyclin Protease

    // STATE PARAMETERS:

    double V1;
    double V3;


public:
    Goldbeter1991SbmlOdeSystem(std::vector<double> stateVariables = std::vector<double>());

    ~Goldbeter1991SbmlOdeSystem();

    void Init();

    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY);


    // FUNCTION DEFINITIONS:
};

namespace
{
    namespace serialization
    {
        // Provide constructor for serializing Goldbeter1991SbmlOdeSystem
        template <class Archive>
        inline void save_construct_data(Archive &ar, const Goldbeter1991SbmlOdeSystem * t, const unsigned int version)
        {
            // Save data required to construct instance
            const std::vector<double> state_variables = t->rGetConstStateVariables();
            ar << state_variables;
        }

        // Provide constructor for de-serializing Goldbeter1991SbmlOdeSystem
        template <class Archive>
        inline void load_construct_data(Archive &ar, Goldbeter1991SbmlOdeSystem * t, const unsigned int version)
        {
            // Retrieve data from archive required to construct new instance
            std::vector<double> state_variables;
            ar >> state_variables;

            // Invoke inplace constructor to initialise instance
            ::new (t)Goldbeter1991SbmlOdeSystem(state_variables);
        }
    } // namespace serialization
} // namespace ...

// Define SbmlSrnWrapperModel using wrappers
#include "SbmlSrnWrapperModel.hpp"
#include "SbmlSrnWrapperModel.cpp"

typedef SbmlSrnWrapperModel<Goldbeter1991SbmlOdeSystem, 3> Goldbeter1991SbmlSrnModel;

// Declare identifiers for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(Goldbeter1991SbmlOdeSystem)
EXPORT_TEMPLATE_CLASS2(SbmlSrnWrapperModel, Goldbeter1991SbmlOdeSystem, 3)

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(Goldbeter1991SbmlSrnModel)

#endif // GOLDBETER1991SBMLODESYSTEMANDSRNMODEL_HPP_