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

    // PARAMETERS
    double ks_n2;          // ks_n2
    double ks_n2_;         // ks_n2_
    double kd_n2;          // kd_n2
    double Hct1_T;         // Hct1_T
    double ks_b2;          // ks_b2
    double ks_b2_;         // ks_b2_
    double kd_b2;          // kd_b2
    double kd_b2_;         // kd_b2_
    double kd_b2__;        // kd_b2__
    double ks_b5;          // ks_b5
    double ks_b5_;         // ks_b5_
    double kd_b5;          // kd_b5
    double kd_b5_;         // kd_b5_
    double Bck2_0;         // Bck2_0
    double Jn3;            // Jn3
    double Dn3;            // Dn3
    double Cln3_max;       // Cln3_max
    double ks_c1;          // ks_c1
    double ks_c1_;         // ks_c1_
    double kas_b2;         // kas_b2
    double kdi_b2;         // kdi_b2
    double kas_b5;         // kas_b5
    double kdi_b5;         // kdi_b5
    double kd2_c1;         // kd2_c1
    double epsilonc1_n3;   // epsilonc1_n3
    double epsilonc1_k2;   // epsilonc1_k2
    double epsilonc1_b5;   // epsilonc1_b5
    double epsilonc1_b2;   // epsilonc1_b2
    double ks_20;          // ks_20
    double ks_20_;         // ks_20_
    double ka_20;          // ka_20
    double ki_20;          // ki_20
    double ki_20_;         // ki_20_
    double ka_t1;          // ka_t1
    double ka_t1_;         // ka_t1_
    double ki_t1;          // ki_t1
    double ki_t1_;         // ki_t1_
    double Ji_t1;          // Ji_t1
    double Ja_t1;          // Ja_t1
    double epsiloni_t1_n2; // epsiloni_t1_n2
    double epsiloni_t1_b5; // epsiloni_t1_b5
    double epsiloni_t1_b2; // epsiloni_t1_b2
    double mu;             // mu
    double ks_ori;         // ks_ori
    double kd_ori;         // kd_ori
    double epsilonori_b2;  // epsilonori_b2
    double ks_bud;         // ks_bud
    double kd_bud;         // kd_bud
    double epsilonbud_b5;  // epsilonbud_b5
    double ks_spn;         // ks_spn
    double kd_spn;         // kd_spn
    double J_spn;          // J_spn
    double ka_sbf;         // ka_sbf
    double ki_sbf;         // ki_sbf
    double ki_sbf_;        // ki_sbf_
    double Ji_sbf;         // Ji_sbf
    double Ja_sbf;         // Ja_sbf
    double epsilonsbf_n3;  // epsilonsbf_n3
    double epsilonsbf_b5;  // epsilonsbf_b5
    double ka_mcm;         // ka_mcm
    double ki_mcm;         // ki_mcm
    double Ji_mcm;         // Ji_mcm
    double Ja_mcm;         // Ja_mcm
    double ka_swi;         // ka_swi
    double ki_swi;         // ki_swi
    double ki_swi_;        // ki_swi_
    double Ji_swi;         // Ji_swi
    double Ja_swi;         // Ja_swi
    double kd1_c1;         // kd1_c1
    double kd_20;          // kd_20
    double Jd2_c1;         // Jd2_c1

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

    double d_Cln2_dt;
    double d_Clb2_T_dt;
    double d_Clb5_T_dt;
    double d_Sic1_T_dt;
    double d_Clb2_Sic1_dt;
    double d_Clb5_Sic1_dt;
    double d_Cdc20_T_dt;
    double d_Cdc20_dt;
    double d_Hct1_dt;
    double d_mass_dt;
    double d_ORI_dt;
    double d_BUD_dt;
    double d_SPN_dt;

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

    // STOICHIOMETRY VARIABLES

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

    void Initialise(double time = 0.0);

    std::vector<double> RunModelEquations(double time, const std::vector<double>& rStateVariables);
};

// Register the ODE system with Boost serialization
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(Chen2000SbmlOdeSystem)

#endif // CHEN_2000_SBML_ODE_SYSTEM_HPP_