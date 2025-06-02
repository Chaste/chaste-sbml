#ifndef CHEN2004SBMLODESYSTEMANDCELLCYCLEMODEL_HPP_
#define CHEN2004SBMLODESYSTEMANDCELLCYCLEMODEL_HPP_

#include <cmath>
#include <iostream>

#include <boost/serialization/base_object.hpp>

#include "AbstractOdeSystem.hpp"
#include "ChasteSerialization.hpp"

class Chen2004SbmlOdeSystem : public AbstractOdeSystem
{
private:
    // (De-)serialize Chen2004SbmlOdeSystem
    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(AbstractOdeSystem);
    }

    // COMPARTMENTS:
    double cell;

    // CONST PARAMETERS:
    const double b0 = 0.054;
    const double bub2h = 1.0;
    const double C0 = 0.4;
    const double Dn3 = 1.0;
    const double ebudb5 = 1.0;
    const double ebudn2 = 0.25;
    const double ebudn3 = 0.05;
    const double ec1b2 = 0.45;
    const double ec1b5 = 0.1;
    const double ec1k2 = 0.03;
    const double ec1n2 = 0.06;
    const double ec1n3 = 0.3;
    const double ef6b2 = 0.55;
    const double ef6b5 = 0.1;
    const double ef6k2 = 0.03;
    const double ef6n2 = 0.06;
    const double ef6n3 = 0.3;
    const double eicdhb2 = 1.2;
    const double eicdhb5 = 8.0;
    const double eicdhn2 = 0.4;
    const double eicdhn3 = 0.25;
    const double eorib2 = 0.45;
    const double eorib5 = 0.9;
    const double esbfb5 = 2.0;
    const double esbfn2 = 2.0;
    const double esbfn3 = 10.0;
    const double J20ppx = 0.15;
    const double Jacdh = 0.03;
    const double Jaiep = 0.1;
    const double Jamcm = 0.1;
    const double Jasbf = 0.01;
    const double Jatem = 0.1;
    const double Jd2c1 = 0.05;
    const double Jd2f6 = 0.05;
    const double Jicdh = 0.03;
    const double Jiiep = 0.1;
    const double Jimcm = 0.1;
    const double Jisbf = 0.01;
    const double Jitem = 0.1;
    const double Jn3 = 6.0;
    const double Jpds = 0.04;
    const double Jspn = 0.14;
    const double ka15_p = 0.002;
    const double ka15_p_p = 1.0;
    const double ka15p = 0.001;
    const double ka20_p = 0.05;
    const double ka20_p_p = 0.2;
    const double kacdh_p = 0.01;
    const double kacdh_p_p = 0.8;
    const double kaiep = 0.1;
    const double kamcm = 1.0;
    const double kasb2 = 50.0;
    const double kasb5 = 50.0;
    const double kasbf = 0.38;
    const double kasesp = 50.0;
    const double kasf2 = 15.0;
    const double kasf5 = 0.01;
    const double kasrent = 200.0;
    const double kasrentp = 1.0;
    const double kaswi = 2.0;
    const double kd14 = 0.1;
    const double kd1c1 = 0.01;
    const double kd1f6 = 0.01;
    const double kd1pds_p = 0.01;
    const double kd20 = 0.3;
    const double kd2c1 = 1.0;
    const double kd2f6 = 1.0;
    const double kd2pds_p_p = 0.2;
    const double kd3c1 = 1.0;
    const double kd3f6 = 1.0;
    const double kd3pds_p_p = 0.04;
    const double kdb2_p = 0.003;
    const double kdb2_p_p = 0.4;
    const double kdb2p = 0.15;
    const double kdb5_p = 0.01;
    const double kdb5_p_p = 0.16;
    const double kdbud = 0.06;
    const double kdcdh = 0.01;
    const double kdib2 = 0.05;
    const double kdib5 = 0.06;
    const double kdiesp = 0.5;
    const double kdif2 = 0.5;
    const double kdif5 = 0.01;
    const double kdirent = 1.0;
    const double kdirentp = 2.0;
    const double kdn2 = 0.12;
    const double kdnet = 0.03;
    const double kdori = 0.06;
    const double kdppx_p = 0.17;
    const double kdppx_p_p = 2.0;
    const double kdspn = 0.06;
    const double kdswi = 0.08;
    const double ki15 = 0.5;
    const double kicdh_p = 0.001;
    const double kicdh_p_p = 0.08;
    const double kiiep = 0.15;
    const double kimcm = 0.15;
    const double kisbf_p = 0.6;
    const double kisbf_p_p = 8.0;
    const double kiswi = 0.05;
    const double kkpnet_p = 0.01;
    const double kkpnet_p_p = 0.6;
    const double kppc1 = 4.0;
    const double kppf6 = 4.0;
    const double kppnet_p = 0.05;
    const double kppnet_p_p = 3.0;
    const double ks14 = 0.2;
    const double ks1pds_p_p = 0.03;
    const double ks20_p = 0.006;
    const double ks20_p_p = 0.6;
    const double ks2pds_p_p = 0.055;
    const double ksb2_p = 0.001;
    const double ksb2_p_p = 0.04;
    const double ksb5_p = 0.0008;
    const double ksb5_p_p = 0.005;
    const double ksbud = 0.2;
    const double ksc1_p = 0.012;
    const double ksc1_p_p = 0.12;
    const double kscdh = 0.01;
    const double ksf6_p = 0.024;
    const double ksf6_p_p = 0.12;
    const double ksf6_p_p_p = 0.004;
    const double ksn2_p = 0.0;
    const double ksn2_p_p = 0.15;
    const double ksnet = 0.084;
    const double ksori = 2.0;
    const double kspds_p = 0.0;
    const double ksppx = 0.1;
    const double ksspn = 0.1;
    const double ksswi_p = 0.005;
    const double ksswi_p_p = 0.08;
    const double mad2h = 8.0;
    const double mdt = 90.0;

    // STATE VARIABLES:
    double BCK2; // BCK2
    double BUD; // BUD
    double C2; // C2
    double C2P; // C2P
    double C5; // C5
    double C5P; // C5P
    double CDC14; // CDC14
    double CDC14T; // CDC14T
    double CDC15; // CDC15
    double CDC15i; // CDC15i
    double CDC20; // CDC20
    double CDC20i; // CDC20i
    double CDC6; // CDC6
    double CDC6P; // CDC6P
    double CDC6T; // CDC6T
    double CDH1; // CDH1
    double CDH1i; // CDH1i
    double CKIT; // CKIT
    double CLB2; // CLB2
    double CLB2T; // CLB2T
    double CLB5; // CLB5
    double CLB5T; // CLB5T
    double CLN2; // CLN2
    double CLN3; // CLN3
    double ESP1; // ESP1
    double F2; // F2
    double F2P; // F2P
    double F5; // F5
    double F5P; // F5P
    double IE; // IE
    double IEP; // IEP
    double MASS; // MASS
    double MCM1; // MCM1
    double NET1; // NET1
    double NET1P; // NET1P
    double NET1T; // NET1T
    double ORI; // ORI
    double PDS1; // PDS1
    double PE; // PE
    double PPX; // PPX
    double RENT; // RENT
    double RENTP; // RENTP
    double SBF; // SBF
    double SIC1; // SIC1
    double SIC1P; // SIC1P
    double SIC1T; // SIC1T
    double SPN; // SPN
    double SWI5; // SWI5
    double SWI5P; // SWI5P
    double TEM1GDP; // TEM1GDP
    double TEM1GTP; // TEM1GTP

    // STATE PARAMETERS:
    double BUB2;
    double LTE1;
    double MAD2;

    double bub2l;
    double CDC15T;
    double ESP1T;
    double IET;
    double KEZ;
    double KEZ2;
    double lte1h;
    double lte1l;
    double mad2l;
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

    // EVENTS:
    std::vector<bool> eventsSatisfied;

public:
    Chen2004SbmlOdeSystem(std::vector<double> stateVariables = std::vector<double>());

    ~Chen2004SbmlOdeSystem();

    void Init();

    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY);

    bool CalculateStoppingEvent(double time, const std::vector<double>& rY);

    // FUNCTION DEFINITIONS:
    double BB_218(double A1, double A2, double A3, double A4);
    double GK_219(double A1, double A2, double A3, double A4);
    double MichaelisMenten_220(double M1, double J1, double k1, double S1);
    double Mass_Action_2_221(double k1, double S1, double S2);
    double Mass_Action_1_222(double k1, double S1);
};

namespace
{
    namespace serialization
    {
        // Provide constructor for serializing Chen2004SbmlOdeSystem
        template <class Archive>
        inline void save_construct_data(Archive &ar, const Chen2004SbmlOdeSystem * t, const unsigned int version)
        {
            // Save data required to construct instance
            const std::vector<double> state_variables = t->rGetConstStateVariables();
            ar << state_variables;
        }

        // Provide constructor for de-serializing Chen2004SbmlOdeSystem
        template <class Archive>
        inline void load_construct_data(Archive &ar, Chen2004SbmlOdeSystem * t, const unsigned int version)
        {
            // Retrieve data from archive required to construct new instance
            std::vector<double> state_variables;
            ar >> state_variables;

            // Invoke inplace constructor to initialise instance
            ::new (t)Chen2004SbmlOdeSystem(state_variables);
        }
    } // namespace serialization
} // namespace ...

// Define SbmlCellCycleWrapperModel using wrappers
#include "SbmlCellCycleWrapperModel.hpp"
#include "SbmlCellCycleWrapperModel.cpp"

typedef SbmlCellCycleWrapperModel<Chen2004SbmlOdeSystem, 51> Chen2004SbmlCellCycleModel;

// Declare identifiers for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(Chen2004SbmlOdeSystem)
EXPORT_TEMPLATE_CLASS2(SbmlCellCycleWrapperModel, Chen2004SbmlOdeSystem, 51)

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(Chen2004SbmlCellCycleModel)

#endif // CHEN2004SBMLODESYSTEMANDCELLCYCLEMODEL_HPP_