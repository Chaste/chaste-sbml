#ifndef CHEN2000ODESYSTEMANDCELLCYCLEMODEL_HPP_
#define CHEN2000ODESYSTEMANDCELLCYCLEMODEL_HPP_

#include <cmath>
#include <iostream>

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/shared_ptr.hpp>

#include "AbstractOdeSystem.hpp"
#include "ChasteSerialization.hpp"

class Chen2000OdeSystem : public AbstractOdeSystem
{
private:

    /* Declare model compartments. */
    double COMpartment;

    /* Declare model parameters. */
    double Cln2;
    double ks_n2;
    double ks_n2_;
    double kd_n2;
    double Clb2_T;
    double Hct1_T;
    double ks_b2;
    double ks_b2_;
    double kd_b2;
    double kd_b2_;
    double kd_b2__;
    double Vd_b2;
    double Clb2;
    double Clb5;
    double Sic1;
    double Clb5_T;
    double ks_b5;
    double ks_b5_;
    double kd_b5;
    double kd_b5_;
    double Vd_b5;
    double Bck2;
    double Bck2_0;
    double Cln3;
    double Jn3;
    double Dn3;
    double Cln3_max;
    double Sic1_T;
    double ks_c1;
    double ks_c1_;
    double Clb2_Sic1;
    double kas_b2;
    double kdi_b2;
    double Clb5_Sic1;
    double kas_b5;
    double kdi_b5;
    double Vd2_c1;
    double kd2_c1;
    double epsilonc1_n3;
    double epsilonc1_k2;
    double epsilonc1_b5;
    double epsilonc1_b2;
    double Cdc20_T;
    double ks_20;
    double ks_20_;
    double Cdc20;
    double ka_20;
    double ki_20;
    double ki_20_;
    double Vi_20;
    double Hct1;
    double ka_t1;
    double ka_t1_;
    double ki_t1;
    double ki_t1_;
    double Vi_t1;
    double Ji_t1;
    double Ja_t1;
    double epsiloni_t1_n2;
    double epsiloni_t1_b5;
    double epsiloni_t1_b2;
    double mass;
    double mu;
    double ORI;
    double ks_ori;
    double kd_ori;
    double epsilonori_b2;
    double BUD;
    double ks_bud;
    double kd_bud;
    double epsilonbud_b5;
    double SPN;
    double ks_spn;
    double kd_spn;
    double J_spn;
    double SBF;
    double ka_sbf;
    double ki_sbf;
    double ki_sbf_;
    double Va_sbf;
    double Ji_sbf;
    double Ja_sbf;
    double epsilonsbf_n3;
    double epsilonsbf_b5;
    double MBF;
    double Mcm1;
    double ka_mcm;
    double ki_mcm;
    double Ji_mcm;
    double Ja_mcm;
    double Swi5;
    double ka_swi;
    double ki_swi;
    double ki_swi_;
    double Ji_swi;
    double Ja_swi;
    double kd1_c1;
    double kd_20;
    double Jd2_c1;

    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive &archive, const unsigned int version)
    {
        archive &boost::serialization::base_object<AbstractOdeSystem>(*this);
    }

public:

    /* Default constructor. */
    Chen2000OdeSystem(std::vector<double> stateVariables = std::vector<double>());

    /* Destructor. */
    ~Chen2000OdeSystem();

    /* Declare model functions. */

    void Init();

    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY);

};

namespace
{
namespace serialization
{
/* Serialize information required to construct a Chen2000OdeSystem. */
template <class Archive>
inline void save_construct_data(
    Archive &ar, const Chen2000OdeSystem *t, const unsigned int file_version)
{
    const std::vector<double> state_variables = t->rGetConstStateVariables();
    ar & state_variables;
}

/* De-serialize constructor parameters and intitialise a Chen2000OdeSystem. */
template <class Archive>
inline void load_construct_data(
    Archive &ar, Chen2000OdeSystem *t, const unsigned int file_version)
{
    std::vector<double> state_variables;
    ar & state_variables;

    // Invoke inplace constructor to initialise instance
    ::new (t) Chen2000OdeSystem(state_variables);
}
}
} // namespace ...

/* Define cell cycle model using wrappers. */
#include "SbmlCellCycleWrapperModel.hpp"
#include "SbmlCellCycleWrapperModel.cpp"

typedef SbmlCellCycleWrapperModel<Chen2000OdeSystem, 0> Chen2000CellCycleModel;

/* Declare identifiers for the serializer */
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(Chen2000OdeSystem)
EXPORT_TEMPLATE_CLASS2(SbmlCellCycleWrapperModel, Chen2000OdeSystem, 0)

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(Chen2000CellCycleModel)

#endif // CHEN2000ODESYSTEMANDCELLCYCLEMODEL_HPP_