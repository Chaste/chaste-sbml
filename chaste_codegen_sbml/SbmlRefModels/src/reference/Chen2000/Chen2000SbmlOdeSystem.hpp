#ifndef CHEN_2000_SBML_ODE_SYSTEM_HPP_
#define CHEN_2000_SBML_ODE_SYSTEM_HPP_

#include <vector>

#include <boost/serialization/base_object.hpp>

#include "AbstractSbmlOdeSystem.hpp"
#include "ChasteSerialization.hpp"
#include "SbmlEventType.hpp"

class Chen2000SbmlOdeSystem : public AbstractSbmlOdeSystem
{
private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Save / load Chen2000SbmlOdeSystem archive
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
    const double ks_n2 = 0.0;          // ks_n2
    const double ks_n2_ = 0.05;        // ks_n2_
    const double kd_n2 = 0.1;          // kd_n2
    const double Hct1_T = 1.0;         // Hct1_T
    const double ks_b2 = 0.002;        // ks_b2
    const double ks_b2_ = 0.05;        // ks_b2_
    const double kd_b2 = 0.01;         // kd_b2
    const double kd_b2_ = 2.0;         // kd_b2_
    const double kd_b2__ = 0.05;       // kd_b2__
    const double ks_b5 = 0.006;        // ks_b5
    const double ks_b5_ = 0.02;        // ks_b5_
    const double kd_b5 = 0.1;          // kd_b5
    const double kd_b5_ = 0.25;        // kd_b5_
    const double Bck2_0 = 0.0027;      // Bck2_0
    const double Jn3 = 6.0;            // Jn3
    const double Dn3 = 1.0;            // Dn3
    const double Cln3_max = 0.02;      // Cln3_max
    const double ks_c1 = 0.02;         // ks_c1
    const double ks_c1_ = 0.1;         // ks_c1_
    const double kas_b2 = 50.0;        // kas_b2
    const double kdi_b2 = 0.05;        // kdi_b2
    const double kas_b5 = 50.0;        // kas_b5
    const double kdi_b5 = 0.05;        // kdi_b5
    const double kd2_c1 = 0.3;         // kd2_c1
    const double epsilonc1_n3 = 20.0;  // epsilonc1_n3
    const double epsilonc1_k2 = 2.0;   // epsilonc1_k2
    const double epsilonc1_b5 = 1.0;   // epsilonc1_b5
    const double epsilonc1_b2 = 0.067; // epsilonc1_b2
    const double ks_20 = 0.005;        // ks_20
    const double ks_20_ = 0.06;        // ks_20_
    const double ka_20 = 1.0;          // ka_20
    const double ki_20 = 0.1;          // ki_20
    const double ki_20_ = 10.0;        // ki_20_
    const double ka_t1 = 0.04;         // ka_t1
    const double ka_t1_ = 2.0;         // ka_t1_
    const double ki_t1 = 0.0;          // ki_t1
    const double ki_t1_ = 0.64;        // ki_t1_
    const double Ji_t1 = 0.05;         // Ji_t1
    const double Ja_t1 = 0.05;         // Ja_t1
    const double epsiloni_t1_n2 = 1.0; // epsiloni_t1_n2
    const double epsiloni_t1_b5 = 0.5; // epsiloni_t1_b5
    const double epsiloni_t1_b2 = 1.0; // epsiloni_t1_b2
    const double mu = 0.005776;        // mu
    const double ks_ori = 2.0;         // ks_ori
    const double kd_ori = 0.06;        // kd_ori
    const double epsilonori_b2 = 0.4;  // epsilonori_b2
    const double ks_bud = 0.3;         // ks_bud
    const double kd_bud = 0.06;        // kd_bud
    const double epsilonbud_b5 = 1.0;  // epsilonbud_b5
    const double ks_spn = 0.08;        // ks_spn
    const double kd_spn = 0.06;        // kd_spn
    const double J_spn = 0.2;          // J_spn
    const double ka_sbf = 1.0;         // ka_sbf
    const double ki_sbf = 0.5;         // ki_sbf
    const double ki_sbf_ = 6.0;        // ki_sbf_
    const double Ji_sbf = 0.01;        // Ji_sbf
    const double Ja_sbf = 0.01;        // Ja_sbf
    const double epsilonsbf_n3 = 75.0; // epsilonsbf_n3
    const double epsilonsbf_b5 = 0.5;  // epsilonsbf_b5
    const double ka_mcm = 1.0;         // ka_mcm
    const double ki_mcm = 0.15;        // ki_mcm
    const double Ji_mcm = 1.0;         // Ji_mcm
    const double Ja_mcm = 1.0;         // Ja_mcm
    const double ka_swi = 1.0;         // ka_swi
    const double ki_swi = 0.3;         // ki_swi
    const double ki_swi_ = 0.2;        // ki_swi_
    const double Ji_swi = 0.1;         // Ji_swi
    const double Ja_swi = 0.1;         // Ja_swi
    const double kd1_c1 = 0.01;        // kd1_c1
    const double kd_20 = 0.08;         // kd_20
    const double Jd2_c1 = 0.05;        // Jd2_c1

    // STATE VARIABLES
    double Cln2;      // Cln2
    double Clb2_T;    // Clb2_T
    double Clb5_T;    // Clb5_T
    double Sic1_T;    // Sic1_T
    double Clb2_Sic1; // Clb2_Sic1
    double Clb5_Sic1; // Clb5_Sic1
    double Cdc20_T;   // Cdc20_T
    double Cdc20;     // Cdc20
    double Hct1;      // Hct1
    double mass;      // mass
    double ORI;       // ORI
    double BUD;       // BUD
    double SPN;       // SPN

    // DERIVED QUANTITIES
    double COMpartment; //
    double Vd_b2;       // Vd_b2
    double Clb2;        // Clb2
    double Clb5;        // Clb5
    double Sic1;        // Sic1
    double Vd_b5;       // Vd_b5
    double Bck2;        // Bck2
    double Cln3;        // Cln3
    double Vd2_c1;      // Vd2_c1
    double Vi_20;       // Vi_20
    double Vi_t1;       // Vi_t1
    double SBF;         // SBF
    double Va_sbf;      // Va_sbf
    double MBF;         // MBF
    double Mcm1;        // Mcm1
    double Swi5;        // Swi5

    // VARIABLE PARAMETERS

    // REACTIONS

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
    Chen2000SbmlOdeSystem();

    /**
     * Destructor
     */
    ~Chen2000SbmlOdeSystem();

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
CHASTE_CLASS_EXPORT(Chen2000SbmlOdeSystem)

#endif // CHEN_2000_SBML_ODE_SYSTEM_HPP_