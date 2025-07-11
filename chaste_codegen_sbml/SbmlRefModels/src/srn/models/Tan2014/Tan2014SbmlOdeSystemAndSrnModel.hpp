#ifndef TAN2014SBMLODESYSTEMANDSRNMODEL_HPP_
#define TAN2014SBMLODESYSTEMANDSRNMODEL_HPP_

#include <vector>

#include <boost/serialization/base_object.hpp>

#include "AbstractOdeSystem.hpp"
#include "ChasteSerialization.hpp"

class Tan2014SbmlOdeSystem : public AbstractOdeSystem
{
private:
    // (De-)serialize Tan2014SbmlOdeSystem
    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(AbstractOdeSystem);
    }

    // CONSTANT PARAMETERS
    const double kdegradation = 0.0163; // 
    const double Bsyn = 1.306; // 
    const double K_n_active_k = 17.16; // 
    const double kC_k1 = 1e-05; // 
    const double kC_k2 = 0.000647; // 
    const double kN_k1 = 0.0001; // 
    const double kN_k2 = 0.00349; // 
    const double kdiffusion_k = 39.13; // 
    const double K_c_active_k = 4.5; // 
    const double ComplexTransitThreshold = 1.0; // ComplexTransitThreshold

    // STATE VARIABLES
    double bcat_cm; // 
    double ligand_cm; // 
    double complex_cm; // 
    double bcat_nu; // 
    double ligand_nu; // 
    double complex_nu; // 

    // DERIVED QUANTITIES
    double drag; // drag

    // VARIABLE PARAMETERS
    double compartment; // 
    double CytosolMembrane; // 
    double nucleus; // 
    double wnt_level; // wnt_level
    double gamma; // gamma

    // RULE-BASED PARAMETERS

    // REACTIONS
    double Bsynthesis; // 
    double kDegradation; // 
    double kC; // 
    double kN; // 
    double kdiffusion; // 
    double K_c_active; // 
    double K_n_active; // 


public:
    Tan2014SbmlOdeSystem(std::vector<double> stateVariables = std::vector<double>());

    ~Tan2014SbmlOdeSystem();

    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY);
    void ProcessRules(double time, const std::vector<double>& rY);

    std::vector<double> ComputeDerivedQuantities(double time, const std::vector<double> &rY);


    // FUNCTIONS
};

namespace
{
    namespace serialization
    {
        // Provide constructor for serializing Tan2014SbmlOdeSystem
        template <class Archive>
        inline void save_construct_data(Archive &ar, const Tan2014SbmlOdeSystem * t, const unsigned int version)
        {
            // Save data required to construct instance
            const std::vector<double> state_variables = t->rGetConstStateVariables();
            ar << state_variables;
        }

        // Provide constructor for de-serializing Tan2014SbmlOdeSystem
        template <class Archive>
        inline void load_construct_data(Archive &ar, Tan2014SbmlOdeSystem * t, const unsigned int version)
        {
            // Retrieve data from archive required to construct new instance
            std::vector<double> state_variables;
            ar >> state_variables;

            // Invoke inplace constructor to initialise instance
            ::new (t)Tan2014SbmlOdeSystem(state_variables);
        }
    } // namespace serialization
} // namespace ...

// Define SbmlSrnWrapperModel using wrappers
#include "SbmlSrnWrapperModel.hpp"
#include "SbmlSrnWrapperModel.cpp"

typedef SbmlSrnWrapperModel<Tan2014SbmlOdeSystem, 6> Tan2014SbmlSrnModel;

// Declare identifiers for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(Tan2014SbmlOdeSystem)
EXPORT_TEMPLATE_CLASS2(SbmlSrnWrapperModel, Tan2014SbmlOdeSystem, 6)

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(Tan2014SbmlSrnModel)

#endif // TAN2014SBMLODESYSTEMANDSRNMODEL_HPP_