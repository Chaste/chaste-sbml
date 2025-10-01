#ifndef TYSONNOVAK2001SBMLODESYSTEM_HPP_
#define TYSONNOVAK2001SBMLODESYSTEM_HPP_

#include <vector>

#include <boost/serialization/base_object.hpp>

#include "AbstractSbmlOdeSystem.hpp"
#include "ChasteSerialization.hpp"
#include "SbmlEventType.hpp"

class TysonNovak2001SbmlOdeSystem : public AbstractSbmlOdeSystem
{
private:
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

    // CONSTANT PARAMETERS
    const double k1 = 0.04;      // k1
    const double k2p = 0.04;     // k2p
    const double k2pp = 1.0;     // k2pp
    const double k2ppp = 1.0;    // k2ppp
    const double k3p = 1.0;      // k3p
    const double k3pp = 10.0;    // k3pp
    const double J3 = 0.04;      // J3
    const double k4 = 35.0;      // k4
    const double k5p = 0.005;    // k5p
    const double k5pp = 0.2;     // k5pp
    const double J5 = 0.3;       // J5
    const double k6 = 0.1;       // k6
    const double n = 4.0;        // n
    const double k7 = 1.0;       // k7
    const double J7 = 0.001;     // J7
    const double k8 = 0.5;       // k8
    const double J8 = 0.001;     // J8
    const double k9 = 0.1;       // k9
    const double k10 = 0.02;     // k10
    const double mu = 0.005;     // mu
    const double k11 = 1.0;      // k11
    const double k12p = 0.2;     // k12p
    const double k12pp = 50.0;   // k12pp
    const double mmax = 10.0;    // mmax
    const double k12ppp = 100.0; // k12ppp
    const double Keq = 1000.0;   // Keq
    const double k13 = 1.0;      // k13
    const double k14 = 1.0;      // k14
    const double k15p = 1.5;     // k15p
    const double k15pp = 0.05;   // k15pp
    const double k16p = 1.0;     // k16p
    const double k16pp = 3.0;    // k16pp
    const double J15 = 0.01;     // J15
    const double J16 = 0.01;     // J16
    const double k4p = 2.0;      // k4p
    const double J4 = 0.04;      // J4

    // STATE VARIABLES
    double CycBt;  // CycBt
    double Cdc20a; // Cdc20a
    double Cdh1;   // Cdh1
    double m;      // m
    double Cdc20t; // Cdc20t
    double IEP;    // IEP
    double CKIt;   // CKIt
    double SK;     // SK

    // DERIVED QUANTITIES
    double CycB;   // CycB
    double Trimer; // Trimer
    double Mad;    // Mad

    // VARIABLE PARAMETERS
    double cell; // cell

    // RULE-BASED PARAMETERS
    double TF; // TF

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
    inline double GK(double A1, double A2, double A3, double A4);
};

// Register the ODE system with Boost serialization
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(TysonNovak2001SbmlOdeSystem)

#endif // TYSONNOVAK2001SBMLODESYSTEM_HPP_