#ifndef GARDNER1998ODESYSTEMANDCELLCYCLEMODEL_HPP_
#define GARDNER1998ODESYSTEMANDCELLCYCLEMODEL_HPP_

#include <cmath>
#include <iostream>

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/shared_ptr.hpp>

#include "AbstractOdeSystem.hpp"
#include "ChasteSerialization.hpp"

class Gardner1998OdeSystem : public AbstractOdeSystem
{
private:

    /* Declare model compartments. */
    double Cell;

    /* Declare model parameters. */
    double V1;
    double K6;
    double V1p;
    double V3;
    double V3p;

    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive &archive, const unsigned int version)
    {
        archive &boost::serialization::base_object<AbstractOdeSystem>(*this);
    }

public:

    /* Default constructor. */
    Gardner1998OdeSystem(std::vector<double> stateVariables = std::vector<double>());

    /* Destructor. */
    ~Gardner1998OdeSystem();

    /* Declare model functions. */

    void Init();

    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY);

};

namespace
{
namespace serialization
{
/* Serialize information required to construct a Gardner1998OdeSystem. */
template <class Archive>
inline void save_construct_data(
    Archive &ar, const Gardner1998OdeSystem *t, const unsigned int file_version)
{
    const std::vector<double> state_variables = t->rGetConstStateVariables();
    ar & state_variables;
}

/* De-serialize constructor parameters and intitialise a Gardner1998OdeSystem. */
template <class Archive>
inline void load_construct_data(
    Archive &ar, Gardner1998OdeSystem *t, const unsigned int file_version)
{
    std::vector<double> state_variables;
    ar & state_variables;

    // Invoke inplace constructor to initialise instance
    ::new (t) Gardner1998OdeSystem(state_variables);
}
}
} // namespace ...

/* Define cell cycle model using wrappers. */
#include "SbmlCellCycleWrapperModel.hpp"
#include "SbmlCellCycleWrapperModel.cpp"

typedef SbmlCellCycleWrapperModel<Gardner1998OdeSystem, 5> Gardner1998CellCycleModel;

/* Declare identifiers for the serializer */
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(Gardner1998OdeSystem)
EXPORT_TEMPLATE_CLASS2(SbmlCellCycleWrapperModel, Gardner1998OdeSystem, 5)

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(Gardner1998CellCycleModel)

#endif // GARDNER1998ODESYSTEMANDCELLCYCLEMODEL_HPP_