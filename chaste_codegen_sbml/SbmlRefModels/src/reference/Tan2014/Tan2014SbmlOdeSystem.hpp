#ifndef TAN_2014_SBML_ODE_SYSTEM_HPP_
#define TAN_2014_SBML_ODE_SYSTEM_HPP_

#include <vector>

#include <boost/serialization/base_object.hpp>

#include "AbstractSbmlOdeSystem.hpp"
#include "ChasteSerialization.hpp"
#include "SbmlEventType.hpp"

class Tan2014SbmlOdeSystem : public AbstractSbmlOdeSystem
{
private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Save / load Tan2014SbmlOdeSystem archive
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive& boost::serialization::base_object<AbstractSbmlOdeSystem>(*this);
    }

    // FIXED PARAMETERS
    const double kdegradation = 0.0163;         //
    const double Bsyn = 1.306;                  //
    const double K_n_active_k = 17.16;          //
    const double kC_k1 = 1e-05;                 //
    const double kC_k2 = 0.000647;              //
    const double kN_k1 = 0.0001;                //
    const double kN_k2 = 0.00349;               //
    const double kdiffusion_k = 39.13;          //
    const double K_c_active_k = 4.5;            //
    const double ComplexTransitThreshold = 1.0; // ComplexTransitThreshold

    // STATE VARIABLES
    double bcat_cm;    //
    double ligand_cm;  //
    double complex_cm; //
    double bcat_nu;    //
    double ligand_nu;  //
    double complex_nu; //

    // DERIVED QUANTITIES
    double drag; // Convert drag amount to concentration

    // VARIABLE PARAMETERS
    double compartment;     //
    double CytosolMembrane; //
    double nucleus;         //
    double wnt_level;       // wnt_level
    double gamma;           // gamma

    // REACTIONS
    double Bsynthesis;   //
    double kDegradation; //
    double kC;           //
    double kN;           //
    double kdiffusion;   //
    double K_c_active;   //
    double K_n_active;   //

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
     * Run the assignment rules to update state.
     *
     * @param time The current time
     * @param rY The current state variables
     */
    void RunAssignmentRules(double time) override;

    /**
     * Run the reactions to update state.
     *
     * @param time The current time
     * @param rY The current state variables
     */
    void RunReactions(double time) override;

    /**
     * Update variable parameters from current ODE system parameter settings.
     *
     * @param time The current time
     */
    void UpdateParameters(double time) override;

    /**
     * Update state variables from the given ODE system state.
     *
     * @param time The current time
     * @param rStateVariables The state variables to use
     */
    void UpdateStateVariables(double time, const std::vector<double>& rStateVariables) override;

    // MODEL FUNCTIONS

public:
    /**
     * Default constructor
     */
    Tan2014SbmlOdeSystem();

    /**
     * Destructor
     */
    ~Tan2014SbmlOdeSystem();

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
};

// Register the ODE system with Boost serialization
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(Tan2014SbmlOdeSystem)

#endif // TAN_2014_SBML_ODE_SYSTEM_HPP_