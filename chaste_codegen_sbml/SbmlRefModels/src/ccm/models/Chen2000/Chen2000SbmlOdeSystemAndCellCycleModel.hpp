#ifndef CHEN2000SBMLODESYSTEMANDCELLCYCLEMODEL_HPP_
#define CHEN2000SBMLODESYSTEMANDCELLCYCLEMODEL_HPP_

#include <cmath>
#include <iostream>

#include <boost/serialization/base_object.hpp>

#include "AbstractOdeSystem.hpp"
#include "ChasteSerialization.hpp"

class Chen2000SbmlOdeSystem : public AbstractOdeSystem
{
private:
    // (De-)serialize Chen2000SbmlOdeSystem
    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(AbstractOdeSystem);
    }

    // COMPARTMENTS:
    double COMpartment;

    // CONST PARAMETERS:
    const double ks_n2 = 0.0;
    const double ks_n2_ = 0.05;
    const double kd_n2 = 0.1;
    const double Hct1_T = 1.0;
    const double ks_b2 = 0.002;
    const double ks_b2_ = 0.05;
    const double kd_b2 = 0.01;
    const double kd_b2_ = 2.0;
    const double kd_b2__ = 0.05;
    const double ks_b5 = 0.006;
    const double ks_b5_ = 0.02;
    const double kd_b5 = 0.1;
    const double kd_b5_ = 0.25;
    const double Bck2_0 = 0.0027;
    const double Jn3 = 6.0;
    const double Dn3 = 1.0;
    const double Cln3_max = 0.02;
    const double ks_c1 = 0.02;
    const double ks_c1_ = 0.1;
    const double kas_b2 = 50.0;
    const double kdi_b2 = 0.05;
    const double kas_b5 = 50.0;
    const double kdi_b5 = 0.05;
    const double kd2_c1 = 0.3;
    const double epsilonc1_n3 = 20.0;
    const double epsilonc1_k2 = 2.0;
    const double epsilonc1_b5 = 1.0;
    const double epsilonc1_b2 = 0.067;
    const double ks_20 = 0.005;
    const double ks_20_ = 0.06;
    const double ka_20 = 1.0;
    const double ki_20 = 0.1;
    const double ki_20_ = 10.0;
    const double ka_t1 = 0.04;
    const double ka_t1_ = 2.0;
    const double ki_t1 = 0.0;
    const double ki_t1_ = 0.64;
    const double Ji_t1 = 0.05;
    const double Ja_t1 = 0.05;
    const double epsiloni_t1_n2 = 1.0;
    const double epsiloni_t1_b5 = 0.5;
    const double epsiloni_t1_b2 = 1.0;
    const double mu = 0.005776;
    const double ks_ori = 2.0;
    const double kd_ori = 0.06;
    const double epsilonori_b2 = 0.4;
    const double ks_bud = 0.3;
    const double kd_bud = 0.06;
    const double epsilonbud_b5 = 1.0;
    const double ks_spn = 0.08;
    const double kd_spn = 0.06;
    const double J_spn = 0.2;
    const double ka_sbf = 1.0;
    const double ki_sbf = 0.5;
    const double ki_sbf_ = 6.0;
    const double Ji_sbf = 0.01;
    const double Ja_sbf = 0.01;
    const double epsilonsbf_n3 = 75.0;
    const double epsilonsbf_b5 = 0.5;
    const double ka_mcm = 1.0;
    const double ki_mcm = 0.15;
    const double Ji_mcm = 1.0;
    const double Ja_mcm = 1.0;
    const double ka_swi = 1.0;
    const double ki_swi = 0.3;
    const double ki_swi_ = 0.2;
    const double Ji_swi = 0.1;
    const double Ja_swi = 0.1;
    const double kd1_c1 = 0.01;
    const double kd_20 = 0.08;
    const double Jd2_c1 = 0.05;

    // STATE VARIABLES:

    // STATE PARAMETERS:

    double Cln2;
    double Clb2_T;
    double Vd_b2;
    double Clb2;
    double Clb5;
    double Sic1;
    double Clb5_T;
    double Vd_b5;
    double Bck2;
    double Cln3;
    double Sic1_T;
    double Clb2_Sic1;
    double Clb5_Sic1;
    double Vd2_c1;
    double Cdc20_T;
    double Cdc20;
    double Vi_20;
    double Hct1;
    double Vi_t1;
    double mass;
    double ORI;
    double BUD;
    double SPN;
    double SBF;
    double Va_sbf;
    double MBF;
    double Mcm1;
    double Swi5;


public:
    Chen2000SbmlOdeSystem(std::vector<double> stateVariables = std::vector<double>());

    ~Chen2000SbmlOdeSystem();

    void Init();

    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY);


    // FUNCTION DEFINITIONS:
};

namespace
{
    namespace serialization
    {
        // Provide constructor for serializing Chen2000SbmlOdeSystem
        template <class Archive>
        inline void save_construct_data(Archive &ar, const Chen2000SbmlOdeSystem * t, const unsigned int version)
        {
            // Save data required to construct instance
            const std::vector<double> state_variables = t->rGetConstStateVariables();
            ar << state_variables;
        }

        // Provide constructor for de-serializing Chen2000SbmlOdeSystem
        template <class Archive>
        inline void load_construct_data(Archive &ar, Chen2000SbmlOdeSystem * t, const unsigned int version)
        {
            // Retrieve data from archive required to construct new instance
            std::vector<double> state_variables;
            ar >> state_variables;

            // Invoke inplace constructor to initialise instance
            ::new (t)Chen2000SbmlOdeSystem(state_variables);
        }
    } // namespace serialization
} // namespace ...

// Define SbmlCellCycleWrapperModel using wrappers
#include "SbmlCellCycleWrapperModel.hpp"
#include "SbmlCellCycleWrapperModel.cpp"

typedef SbmlCellCycleWrapperModel<Chen2000SbmlOdeSystem, 0> Chen2000SbmlCellCycleModel;

// Declare identifiers for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(Chen2000SbmlOdeSystem)
EXPORT_TEMPLATE_CLASS2(SbmlCellCycleWrapperModel, Chen2000SbmlOdeSystem, 0)

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(Chen2000SbmlCellCycleModel)

#endif // CHEN2000SBMLODESYSTEMANDCELLCYCLEMODEL_HPP_