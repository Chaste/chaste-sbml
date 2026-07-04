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

    // PARAMETERS
    double wnt_level;               // wnt_level
    double kdegradation;            //
    double Bsyn;                    //
    double K_n_active_k;            //
    double kC_k1;                   //
    double kC_k2;                   //
    double kN_k1;                   //
    double kN_k2;                   //
    double kdiffusion_k;            //
    double K_c_active_k;            //
    double gamma;                   // gamma
    double ComplexTransitThreshold; // ComplexTransitThreshold

    // STATE VARIABLES
    double bcat_cm;    //
    double ligand_cm;  //
    double complex_cm; //
    double bcat_nu;    //
    double ligand_nu;  //
    double complex_nu; //

    double d_bcat_cm_dt;
    double d_ligand_cm_dt;
    double d_complex_cm_dt;
    double d_bcat_nu_dt;
    double d_ligand_nu_dt;
    double d_complex_nu_dt;

    // DERIVED QUANTITIES
    double compartment;     //
    double CytosolMembrane; //
    double nucleus;         //
    double drag;            // drag

    // STOICHIOMETRY VARIABLES

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

    void Initialise(double time = 0.0);

    std::vector<double> RunModelEquations(double time, const std::vector<double>& rStateVariables);
};

// Register the ODE system with Boost serialization
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(Tan2014SbmlOdeSystem)

#endif // TAN_2014_SBML_ODE_SYSTEM_HPP_