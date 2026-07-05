#ifndef TYSON_NOVAK_2001_SBML_ODE_SYSTEM_HPP_
#define TYSON_NOVAK_2001_SBML_ODE_SYSTEM_HPP_

#include <vector>

#include <boost/serialization/base_object.hpp>

#include "AbstractSbmlOdeSystem.hpp"
#include "ChasteSerialization.hpp"
#include "SbmlEventType.hpp"

class TysonNovak2001SbmlOdeSystem : public AbstractSbmlOdeSystem
{
private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Save / load TysonNovak2001SbmlOdeSystem archive
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
    double k1;     // k1
    double k2p;    // k2p
    double k2pp;   // k2pp
    double k2ppp;  // k2ppp
    double k3p;    // k3p
    double k3pp;   // k3pp
    double J3;     // J3
    double k4;     // k4
    double k5p;    // k5p
    double k5pp;   // k5pp
    double J5;     // J5
    double k6;     // k6
    double n;      // n
    double k7;     // k7
    double J7;     // J7
    double k8;     // k8
    double J8;     // J8
    double k9;     // k9
    double k10;    // k10
    double mu;     // mu
    double k11;    // k11
    double k12p;   // k12p
    double k12pp;  // k12pp
    double mmax;   // mmax
    double k12ppp; // k12ppp
    double Keq;    // Keq
    double k13;    // k13
    double k14;    // k14
    double k15p;   // k15p
    double k15pp;  // k15pp
    double k16p;   // k16p
    double k16pp;  // k16pp
    double J15;    // J15
    double J16;    // J16
    double k4p;    // k4p
    double J4;     // J4

    // STATE VARIABLES
    double CycBt;  // CycBt
    double Cdc20a; // Cdc20a
    double Cdh1;   // Cdh1
    double m;      // m
    double Cdc20t; // Cdc20t
    double IEP;    // IEP
    double CKIt;   // CKIt
    double SK;     // SK

    double d_CycBt_dt;
    double d_Cdc20a_dt;
    double d_Cdh1_dt;
    double d_m_dt;
    double d_Cdc20t_dt;
    double d_IEP_dt;
    double d_CKIt_dt;
    double d_SK_dt;

    // DERIVED QUANTITIES
    double cell;   // cell
    double CycB;   // CycB
    double Trimer; // Trimer
    double Mad;    // Mad
    double TF;     // TF

    // STOICHIOMETRY VARIABLES

    // REACTIONS
    double CycBt_synthesis;           // CycBt synthesis
    double CycBdegradation;           // CycBt degradation
    double CycBdegradationviaCdh1;    // CycBt degradation via Cdh1
    double CycBtdegradationviaCdc20a; // CycBt degradation via Cdc20a
    double Cdh1synthesis;             // Cdh1 synthesis
    double Cdh1degradation;           // Cdh1 degradation
    double Cdc20tsynthesis;           // Cdc20t synthesis
    double Cdc20t_deg;                // Cdc20t degradation
    double Cdc20activation;           // Cdc20 activation
    double Cdc20ainhibition;          // Cdc20a inhibition
    double Cdc20adegradation;         // Cdc20a degradation
    double IEPsynthesis;              // IEP synthesis
    double IEPdegradation;            // IEP degradation
    double growth;                    // growth
    double CKItsynthesis;             // CKIt synthesis
    double CKIdegradation;            // CKIt degradation
    double CKItphosphorilationviaSK;  // CKIt phosphorilation via SK
    double eq_7;                      // CKIt Trimer sequestred
    double SKsynthesis;               // SK synthesis
    double SKdegradation;             // SK degradation

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
    inline double GK(double A1, double A2, double A3, double A4);

public:
    /**
     * Default constructor
     */
    TysonNovak2001SbmlOdeSystem();

    /**
     * Destructor
     */
    ~TysonNovak2001SbmlOdeSystem();

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
CHASTE_CLASS_EXPORT(TysonNovak2001SbmlOdeSystem)

#endif // TYSON_NOVAK_2001_SBML_ODE_SYSTEM_HPP_