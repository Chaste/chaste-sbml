#ifndef GOLDBETER1991SBMLODESYSTEMANDSRNMODEL_HPP_
#define GOLDBETER1991SBMLODESYSTEMANDSRNMODEL_HPP_

#include <vector>

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

    // CONSTANT PARAMETERS
    const double VM1 = 3.0; // VM1
    const double VM3 = 1.0; // VM3
    const double Kc = 0.5; // Kc

    // STATE VARIABLES
    double C; // Cyclin
    double M; // cdc_2_kinase
    double X; // Cyclin Protease

    // DERIVED QUANTITIES

    // VARIABLE PARAMETERS
    double cell; // cell

    // RULE-BASED PARAMETERS
    double V1; // V1
    double V3; // V3

    // REACTIONS
    double reaction1; // creation of cyclin
    double reaction2; // default degradation of cyclin
    double reaction3; // cdc2 kinase triggered degration of cyclin
    double reaction4; // activation of cdc2 kinase
    double reaction5; // deactivation of cdc2 kinase
    double reaction6; // activation of cyclin protease
    double reaction7; // deactivation of cyclin protease


public:
    Goldbeter1991SbmlOdeSystem(std::vector<double> stateVariables = std::vector<double>());

    ~Goldbeter1991SbmlOdeSystem();

    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY);
    void ProcessRules(double time, const std::vector<double>& rY);



    // FUNCTIONS
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