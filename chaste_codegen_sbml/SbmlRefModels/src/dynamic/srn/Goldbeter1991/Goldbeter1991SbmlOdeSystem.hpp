#ifndef GOLDBETER_1991_SBML_ODE_SYSTEM_HPP_
#define GOLDBETER_1991_SBML_ODE_SYSTEM_HPP_

#include <vector>

#include <boost/serialization/base_object.hpp>

#include "AbstractSbmlOdeSystem.hpp"
#include "ChasteSerialization.hpp"
#include "SbmlEventType.hpp"

class Goldbeter1991SbmlOdeSystem : public AbstractSbmlOdeSystem
{
private:
    friend class boost::serialization::access;
    /**
     * Archive / unarchive Goldbeter1991SbmlOdeSystem.
     *
     * @param archive the archive
     * @param version the current version of this class
     */

    template <class Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar& BOOST_SERIALIZATION_BASE_OBJECT_NVP(AbstractSbmlOdeSystem);
    }

    // CONSTANT PARAMETERS
    const double VM1 = 3.0; // VM1
    const double VM3 = 1.0; // VM3
    const double Kc = 0.5;  // Kc

    // STATE VARIABLES
    double C; // Cyclin
    double M; // cdc_2_kinase
    double X; // Cyclin Protease

    // DERIVED QUANTITIES

    // VARIABLE PARAMETERS
    double cell; // cell

    // RULE-BASED PARAMETERS
    double V1; // V1
    double V3; // V3

    // REACTIONS
    double reaction1; // creation of cyclin
    double reaction2; // default degradation of cyclin
    double reaction3; // cdc2 kinase triggered degration of cyclin
    double reaction4; // activation of cdc2 kinase
    double reaction5; // deactivation of cdc2 kinase
    double reaction6; // activation of cyclin protease
    double reaction7; // deactivation of cyclin protease

public:
    /**
     * Default constructor
     */
    Goldbeter1991SbmlOdeSystem();

    /**
     * Destructor
     */
    ~Goldbeter1991SbmlOdeSystem();

    /**
     * Compute the derived quantities from the given system state.
     *
     * @param time  the time at which to compute the derived quantities
     * @param rY a vector of values for the state variables
     *
     * @return a vector of derived quantities
     */
    std::vector<double> ComputeDerivedQuantities(double time, const std::vector<double>& rY) override;

    /**
     * Compute the RHS of the ODE system.
     *
     * An ODE solver will call this function repeatedly to solve for y = [y1 ... yn].
     *
     * @param time the time used to evaluate the RHS.
     * @param rY an input solution vector used to evaluate the RHS.
     * @param rDY an output vector to be filled in with the resulting derivatives y' = [y1' ... yn'].
     */
    void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY) override;

    /**
     * Process the events in the model.
     *
     * @param time The current time
     * @param rY The current state variables
     *
     * @return How close we are to the time of the next event
     */
    double ProcessModelEvents(double time, const std::vector<double>& rY) override;

    /**
     * Run the equations governing the model to update state.
     *
     * @param time The current time
     * @param rY The current state variables
     */
    void RunModelRules(double time, const std::vector<double>& rY) override;

    // MODEL FUNCTIONS
};

// namespace
// {
// namespace serialization
// {
//     // Provide constructor for serializing Goldbeter1991SbmlOdeSystem
//     template <class Archive>
//     inline void save_construct_data(Archive& ar, const Goldbeter1991SbmlOdeSystem* t, const unsigned int version)
//     {
//         // Save data required to construct instance
//         const std::vector<double> state_variables = t->rGetConstStateVariables();
//         ar << state_variables;
//     }

//     // Provide constructor for de-serializing Goldbeter1991SbmlOdeSystem
//     template <class Archive>
//     inline void load_construct_data(Archive& ar, Goldbeter1991SbmlOdeSystem* t, const unsigned int version)
//     {
//         // Retrieve data from archive required to construct new instance
//         std::vector<double> state_variables;
//         ar >> state_variables;

//         // Invoke inplace constructor to initialise instance
//         ::new (t) Goldbeter1991SbmlOdeSystem(state_variables);
//     }
// } // namespace serialization
// } // namespace

// Define SbmlSrnWrapperModel using wrappers
// #include "SbmlSrnWrapperModel.cpp"
// #include "SbmlSrnWrapperModel.hpp"

// typedef SbmlSrnWrapperModel<Goldbeter1991SbmlOdeSystem, 3> Goldbeter1991SbmlSrnModel;

// Declare identifiers for the serializer
// #include "SerializationExportWrapper.hpp"
// CHASTE_CLASS_EXPORT(Goldbeter1991SbmlOdeSystem)
// EXPORT_TEMPLATE_CLASS2(SbmlSrnWrapperModel, Goldbeter1991SbmlOdeSystem, 3)

// #include "CellCycleModelOdeSolverExportWrapper.hpp"
// EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(Goldbeter1991SbmlSrnModel)

#endif // GOLDBETER_1991_SBML_ODE_SYSTEM_HPP_