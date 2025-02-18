#ifndef GOLDBETER_1991ODESYSTEMANDSRNMODEL_HPP_
#define GOLDBETER_1991ODESYSTEMANDSRNMODEL_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/shared_ptr.hpp>

#include <cmath>
#include <iostream>
#include "AbstractOdeSystem.hpp"

class goldbeter_1991OdeSystem : public AbstractOdeSystem
{
private:

    /* Initialise compartments and values. */
    double cell;

    /* Initialise model parameters. */
    double V1;
    double V3;
    double VM1;
    double VM3;
    double Kc;


    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractOdeSystem>(*this);
    }

public:

    /* Default constructor. */
    goldbeter_1991OdeSystem(std::vector<double> stateVariables=std::vector<double>());

    /* Destructor. */
    ~goldbeter_1991OdeSystem();


    void Init();

    void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY);

};

namespace
{
namespace serialization
{
/* Serialize information required to construct a goldbeter_1991OdeSystem. */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const goldbeter_1991OdeSystem * t, const unsigned int file_version)
{
    const std::vector<double> state_variables = t->rGetConstStateVariables();
    ar & state_variables;
}
/* De-serialize constructor parameters and intiialise a goldbeter_1991OdeSystem. */
template<class Archive>
inline void load_construct_data(
    Archive & ar, goldbeter_1991OdeSystem * t, const unsigned int file_version)
{
    std::vector<double> state_variables;
    ar & state_variables;
    
    // Invoke inplace constructor to initialise instance
    ::new(t)goldbeter_1991OdeSystem(state_variables);
}
}
} // namespace ...

/* Define SRN model using Wrappers. */
#include "SbmlSrnWrapperModel.hpp"
#include "SbmlSrnWrapperModel.cpp"

typedef SbmlSrnWrapperModel<goldbeter_1991OdeSystem,5> goldbeter_1991;

// Declare identifiers for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(goldbeter_1991OdeSystem)
EXPORT_TEMPLATE_CLASS2(SbmlSrnWrapperModel, goldbeter_1991OdeSystem, 5)

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(goldbeter_1991)

#endif /* GOLDBETER_1991ODESYSTEMANDSRNMODEL_HPP_ */