#ifndef CHEN2004ODESYSTEMANDCELLCYCLEMODEL_HPP_
#define CHEN2004ODESYSTEMANDCELLCYCLEMODEL_HPP_

#include <cmath>
#include <iostream>

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/shared_ptr.hpp>

#include "AbstractOdeSystem.hpp"
#include "ChasteSerialization.hpp"

class Chen2004OdeSystem : public AbstractOdeSystem
{
private:

    /* Declare model compartments. */
    double cell;

    /* Declare model parameters. */
    double b0;
    double bub2h;
    double bub2l;
    double C0;
    double CDC15T;
    double Dn3;
    double ebudb5;
    double ebudn2;
    double ebudn3;
    double ec1b2;
    double ec1b5;
    double ec1k2;
    double ec1n2;
    double ec1n3;
    double ef6b2;
    double ef6b5;
    double ef6k2;
    double ef6n2;
    double ef6n3;
    double eicdhb2;
    double eicdhb5;
    double eicdhn2;
    double eicdhn3;
    double eorib2;
    double eorib5;
    double esbfb5;
    double esbfn2;
    double esbfn3;
    double ESP1T;
    double IET;
    double J20ppx;
    double Jacdh;
    double Jaiep;
    double Jamcm;
    double Jasbf;
    double Jatem;
    double Jd2c1;
    double Jd2f6;
    double Jicdh;
    double Jiiep;
    double Jimcm;
    double Jisbf;
    double Jitem;
    double Jn3;
    double Jpds;
    double Jspn;
    double ka15_p;
    double ka15_p_p;
    double ka15p;
    double ka20_p;
    double ka20_p_p;
    double kacdh_p;
    double kacdh_p_p;
    double kaiep;
    double kamcm;
    double kasb2;
    double kasb5;
    double kasbf;
    double kasesp;
    double kasf2;
    double kasf5;
    double kasrent;
    double kasrentp;
    double kaswi;
    double kd14;
    double kd1c1;
    double kd1f6;
    double kd1pds_p;
    double kd20;
    double kd2c1;
    double kd2f6;
    double kd2pds_p_p;
    double kd3c1;
    double kd3f6;
    double kd3pds_p_p;
    double kdb2_p;
    double kdb2_p_p;
    double kdb2p;
    double kdb5_p;
    double kdb5_p_p;
    double kdbud;
    double kdcdh;
    double kdib2;
    double kdib5;
    double kdiesp;
    double kdif2;
    double kdif5;
    double kdirent;
    double kdirentp;
    double kdn2;
    double kdnet;
    double kdori;
    double kdppx_p;
    double kdppx_p_p;
    double kdspn;
    double kdswi;
    double KEZ;
    double KEZ2;
    double ki15;
    double kicdh_p;
    double kicdh_p_p;
    double kiiep;
    double kimcm;
    double kisbf_p;
    double kisbf_p_p;
    double kiswi;
    double kkpnet_p;
    double kkpnet_p_p;
    double kppc1;
    double kppf6;
    double kppnet_p;
    double kppnet_p_p;
    double ks14;
    double ks1pds_p_p;
    double ks20_p;
    double ks20_p_p;
    double ks2pds_p_p;
    double ksb2_p;
    double ksb2_p_p;
    double ksb5_p;
    double ksb5_p_p;
    double ksbud;
    double ksc1_p;
    double ksc1_p_p;
    double kscdh;
    double ksf6_p;
    double ksf6_p_p;
    double ksf6_p_p_p;
    double ksn2_p;
    double ksn2_p_p;
    double ksnet;
    double ksori;
    double kspds_p;
    double ksppx;
    double ksspn;
    double ksswi_p;
    double ksswi_p_p;
    double lte1h;
    double lte1l;
    double mad2h;
    double mad2l;
    double mdt;
    double TEM1T;
    double D;
    double mu;
    double Vdb5;
    double Vdb2;
    double Vasbf;
    double Visbf;
    double Vkpc1;
    double Vkpf6;
    double Vacdh;
    double Vicdh;
    double Vppnet;
    double Vkpnet;
    double Vdppx;
    double Vdpds;
    double Vaiep;
    double Vd2c1;
    double Vd2f6;
    double Vppc1;
    double Vppf6;
    double F;

    std::vector<bool> eventsSatisfied;

    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive &archive, const unsigned int version)
    {
        archive &boost::serialization::base_object<AbstractOdeSystem>(*this);
    }

public:

    /* Default constructor. */
    Chen2004OdeSystem(std::vector<double> stateVariables = std::vector<double>());

    /* Destructor. */
    ~Chen2004OdeSystem();

    /* Declare model functions. */
    double BB_218(double A1, double A2, double A3, double A4);
    double GK_219(double A1, double A2, double A3, double A4);
    double MichaelisMenten_220(double M1, double J1, double k1, double S1);
    double Mass_Action_2_221(double k1, double S1, double S2);
    double Mass_Action_1_222(double k1, double S1);

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
/* Serialize information required to construct a Chen2004OdeSystem. */
template <class Archive>
inline void save_construct_data(
    Archive &ar, const Chen2004OdeSystem *t, const unsigned int file_version)
{
    const std::vector<double> state_variables = t->rGetConstStateVariables();
    ar & state_variables;
}

/* De-serialize constructor parameters and intitialise a Chen2004OdeSystem. */
template <class Archive>
inline void load_construct_data(
    Archive &ar, Chen2004OdeSystem *t, const unsigned int file_version)
{
    std::vector<double> state_variables;
    ar & state_variables;

    // Invoke inplace constructor to initialise instance
    ::new (t) Chen2004OdeSystem(state_variables);
}
}
} // namespace ...

/* Define SbmlCellCycleWrapperModel using wrappers. */
#include "SbmlCellCycleWrapperModel.hpp"
#include "SbmlCellCycleWrapperModel.cpp"

typedef SbmlCellCycleWrapperModel<Chen2004OdeSystem, 51> Chen2004CellCycleModel;

/* Declare identifiers for the serializer */
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(Chen2004OdeSystem)
EXPORT_TEMPLATE_CLASS2(SbmlCellCycleWrapperModel, Chen2004OdeSystem, 51)

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(Chen2004CellCycleModel)

#endif // CHEN2004ODESYSTEMANDCELLCYCLEMODEL_HPP_