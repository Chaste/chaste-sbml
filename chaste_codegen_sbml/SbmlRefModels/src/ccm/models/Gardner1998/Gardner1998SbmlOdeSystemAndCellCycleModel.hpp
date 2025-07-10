#ifndef GARDNER1998SBMLODESYSTEMANDCELLCYCLEMODEL_HPP_
#define GARDNER1998SBMLODESYSTEMANDCELLCYCLEMODEL_HPP_

#include <vector>

#include <boost/serialization/base_object.hpp>

#include "AbstractOdeSystem.hpp"
#include "ChasteSerialization.hpp"

class Gardner1998SbmlOdeSystem : public AbstractOdeSystem
{
private:
    // (De-)serialize Gardner1998SbmlOdeSystem
    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(AbstractOdeSystem);
    }

    // CONSTANT PARAMETERS
    const double K6 = 0.3; // K6
    const double V1p = 0.75; // V1p
    const double V3p = 0.3; // V3p

    // STATE VARIABLES
    double C; // cyclin
    double X; // protease
    double M; // cdc2k
    double Y; // cyclin inhibitor
    double Z; // complex inhibitor-cyclin

    // DERIVED QUANTITIES

    // PARAMETERS
    double Cell; // Cell
    double V1; // V1
    double V3; // V3

    // REACTIONS
    double reaction1; // creation of cyclin
    double reaction2; // cdc2 kinase triggered degration of cyclin
    double reaction3; // default degradation of cyclin
    double reaction4; // activation of cdc2 kinase
    double reaction5; // deactivation of cdc2 kinase
    double reaction6; // activation of cyclin protease
    double reaction7; // deactivation of cyclin protease
    double reaction8; // reaction8
    double reaction9; // reaction9
    double reaction10; // desinhibition of cyclin
    double reaction11; // degradation of inhibited cyclin
    double reaction12; // creation of cyclin inhibitor
    double reaction13; // degradation of cyclin inhibitor


public:
    Gardner1998SbmlOdeSystem(std::vector<double> stateVariables = std::vector<double>());

    ~Gardner1998SbmlOdeSystem();

    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY);
    void ProcessRules(double time, const std::vector<double>& rY);



    // FUNCTIONS
};

namespace
{
    namespace serialization
    {
        // Provide constructor for serializing Gardner1998SbmlOdeSystem
        template <class Archive>
        inline void save_construct_data(Archive &ar, const Gardner1998SbmlOdeSystem * t, const unsigned int version)
        {
            // Save data required to construct instance
            const std::vector<double> state_variables = t->rGetConstStateVariables();
            ar << state_variables;
        }

        // Provide constructor for de-serializing Gardner1998SbmlOdeSystem
        template <class Archive>
        inline void load_construct_data(Archive &ar, Gardner1998SbmlOdeSystem * t, const unsigned int version)
        {
            // Retrieve data from archive required to construct new instance
            std::vector<double> state_variables;
            ar >> state_variables;

            // Invoke inplace constructor to initialise instance
            ::new (t)Gardner1998SbmlOdeSystem(state_variables);
        }
    } // namespace serialization
} // namespace ...

// Define SbmlCellCycleWrapperModel using wrappers
#include "SbmlCellCycleWrapperModel.hpp"
#include "SbmlCellCycleWrapperModel.cpp"

typedef SbmlCellCycleWrapperModel<Gardner1998SbmlOdeSystem, 5> Gardner1998SbmlCellCycleModel;

// Declare identifiers for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(Gardner1998SbmlOdeSystem)
EXPORT_TEMPLATE_CLASS2(SbmlCellCycleWrapperModel, Gardner1998SbmlOdeSystem, 5)

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(Gardner1998SbmlCellCycleModel)

#endif // GARDNER1998SBMLODESYSTEMANDCELLCYCLEMODEL_HPP_