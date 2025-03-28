#ifndef GOLDBETER1991ODESYSTEMANDSRNMODEL_HPP_
#define GOLDBETER1991ODESYSTEMANDSRNMODEL_HPP_

#include <cmath>
#include <iostream>

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/shared_ptr.hpp>

#include "AbstractOdeSystem.hpp"
#include "ChasteSerialization.hpp"

class Goldbeter1991OdeSystem : public AbstractOdeSystem
{
private:

    /* Declare model compartments. */
    double cell;

    /* Declare model parameters. */
    double V1;
    double V3;
    double VM1;
    double VM3;
    double Kc;


    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive &archive, const unsigned int version)
    {
        archive &boost::serialization::base_object<AbstractOdeSystem>(*this);
    }

public:

    /* Default constructor. */
    Goldbeter1991OdeSystem(std::vector<double> stateVariables = std::vector<double>());

    /* Destructor. */
    ~Goldbeter1991OdeSystem();

    /* Declare model functions. */

    void Init();

    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY);

};

namespace
{
namespace serialization
{
/* Serialize information required to construct a Goldbeter1991OdeSystem. */
template <class Archive>
inline void save_construct_data(
    Archive &ar, const Goldbeter1991OdeSystem *t, const unsigned int file_version)
{
    const std::vector<double> state_variables = t->rGetConstStateVariables();
    ar & state_variables;
}

/* De-serialize constructor parameters and intitialise a Goldbeter1991OdeSystem. */
template <class Archive>
inline void load_construct_data(
    Archive &ar, Goldbeter1991OdeSystem *t, const unsigned int file_version)
{
    std::vector<double> state_variables;
    ar & state_variables;

    // Invoke inplace constructor to initialise instance
    ::new (t) Goldbeter1991OdeSystem(state_variables);
}
}
} // namespace ...

/* Define SbmlSrnWrapperModel using wrappers. */
#include "SbmlSrnWrapperModel.hpp"
#include "SbmlSrnWrapperModel.cpp"

typedef SbmlSrnWrapperModel<Goldbeter1991OdeSystem, 3> Goldbeter1991SrnModel;

/* Declare identifiers for the serializer */
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(Goldbeter1991OdeSystem)
EXPORT_TEMPLATE_CLASS2(SbmlSrnWrapperModel, Goldbeter1991OdeSystem, 3)

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(Goldbeter1991SrnModel)

#endif // GOLDBETER1991ODESYSTEMANDSRNMODEL_HPP_