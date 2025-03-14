#ifndef TYSONNOVAK2001ODESYSTEMANDCELLCYCLEMODEL_HPP_
#define TYSONNOVAK2001ODESYSTEMANDCELLCYCLEMODEL_HPP_

#include <cmath>
#include <iostream>

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/shared_ptr.hpp>

#include "AbstractOdeSystem.hpp"
#include "ChasteSerialization.hpp"

class TysonNovak2001OdeSystem : public AbstractOdeSystem
{
private:

    /* Declare model compartments. */
    double cell;

    /* Declare model parameters. */
    double k1;
    double k2p;
    double k2pp;
    double k2ppp;
    double k3p;
    double k3pp;
    double J3;
    double k4;
    double k5p;
    double k5pp;
    double J5;
    double k6;
    double n;
    double k7;
    double J7;
    double k8;
    double J8;
    double k9;
    double k10;
    double mu;
    double k11;
    double k12p;
    double k12pp;
    double mmax;
    double k12ppp;
    double Keq;
    double k13;
    double k14;
    double k15p;
    double k15pp;
    double k16p;
    double k16pp;
    double J15;
    double J16;
    double k4p;
    double J4;
    double TF;

    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive &archive, const unsigned int version)
    {
        archive &boost::serialization::base_object<AbstractOdeSystem>(*this);
    }

public:

    /* Default constructor. */
    TysonNovak2001OdeSystem(std::vector<double> stateVariables = std::vector<double>());

    /* Destructor. */
    ~TysonNovak2001OdeSystem();

    /* Declare model functions. */
    double GK(double A1, double A2, double A3, double A4);

    void Init();

    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY);

    // Stopping event is required for Cell Cycle Models to divide.
    bool CalculateStoppingEvent(double time, const std::vector<double>& rY);

    // Checks if any events have been triggered and updates the system accordingly.
    void CheckAndUpdateEvents(double time, const std::vector<double>& rY);

    // Checks whether all events are satisifed.
    bool AreAllEventsSatisfied(double time, const std::vector<double>& rY);
};

namespace
{
namespace serialization
{
/* Serialize information required to construct a TysonNovak2001OdeSystem. */
template <class Archive>
inline void save_construct_data(
    Archive &ar, const TysonNovak2001OdeSystem *t, const unsigned int file_version)
{
    const std::vector<double> state_variables = t->rGetConstStateVariables();
    ar & state_variables;
}

/* De-serialize constructor parameters and intitialise a TysonNovak2001OdeSystem. */
template <class Archive>
inline void load_construct_data(
    Archive &ar, TysonNovak2001OdeSystem *t, const unsigned int file_version)
{
    std::vector<double> state_variables;
    ar & state_variables;

    // Invoke inplace constructor to initialise instance
    ::new (t) TysonNovak2001OdeSystem(state_variables);
}
}
} // namespace ...

/* Define cell cycle model using wrappers. */
#include "SbmlCellCycleWrapperModel.hpp"
#include "SbmlCellCycleWrapperModel.cpp"

typedef SbmlCellCycleWrapperModel<TysonNovak2001OdeSystem, 11> TysonNovak2001CellCycleModel;

/* Declare identifiers for the serializer */
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(TysonNovak2001OdeSystem)
EXPORT_TEMPLATE_CLASS2(SbmlCellCycleWrapperModel, TysonNovak2001OdeSystem, 11)

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(TysonNovak2001CellCycleModel)

#endif // TYSONNOVAK2001ODESYSTEMANDCELLCYCLEMODEL_HPP_