#ifndef GARDNER1998SBMLODESYSTEMANDCELLCYCLEMODEL_HPP_
#define GARDNER1998SBMLODESYSTEMANDCELLCYCLEMODEL_HPP_

#include <cmath>
#include <iostream>

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

    // COMPARTMENTS:
    double Cell;

    // CONST PARAMETERS:
    const double K6 = 0.3;
    const double V1p = 0.75;
    const double V3p = 0.3;

    // STATE VARIABLES:
    double C; // cyclin
    double X; // protease
    double M; // cdc2k
    double Y; // cyclin inhibitor
    double Z; // complex inhibitor-cyclin

    // STATE PARAMETERS:

    double V1;
    double V3;


public:
    Gardner1998SbmlOdeSystem(std::vector<double> stateVariables = std::vector<double>());

    ~Gardner1998SbmlOdeSystem();

    void Init();

    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY);


    // FUNCTION DEFINITIONS:
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