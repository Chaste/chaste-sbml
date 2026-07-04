#include <cmath>
#include <limits>
#include <vector>

#include "CellwiseOdeSystemInformation.hpp"
#include "SbmlEventType.hpp"
#include "SbmlMath.hpp"

#include "Tan2014SbmlOdeSystem.hpp"

namespace sm = sbmlmath;

Tan2014SbmlOdeSystem::Tan2014SbmlOdeSystem()
        : AbstractSbmlOdeSystem(6, 12, 0)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<Tan2014SbmlOdeSystem>);

    Initialise();

    // EVENTS
}

Tan2014SbmlOdeSystem::~Tan2014SbmlOdeSystem()
{
}

std::vector<double> Tan2014SbmlOdeSystem::ComputeDerivedQuantities(double time, const std::vector<double>& rY)
{
    std::vector<double> dqs;
    RunModelEquations(time, rY);

    // AMOUNT / CONCENTRATION CONVERSIONS
    double amt__bcat_cm = bcat_cm * CytosolMembrane;       //
    double amt__ligand_cm = ligand_cm * CytosolMembrane;   //
    double amt__complex_cm = complex_cm * CytosolMembrane; //
    double amt__bcat_nu = bcat_nu * nucleus;               //
    double amt__ligand_nu = ligand_nu * nucleus;           //
    double amt__complex_nu = complex_nu * nucleus;         //
    double amt__drag = drag * CytosolMembrane;             //

    dqs.push_back(compartment);
    dqs.push_back(CytosolMembrane);
    dqs.push_back(nucleus);
    dqs.push_back(drag);
    dqs.push_back(Bsynthesis);
    dqs.push_back(kDegradation);
    dqs.push_back(kC);
    dqs.push_back(kN);
    dqs.push_back(kdiffusion);
    dqs.push_back(K_c_active);
    dqs.push_back(K_n_active);
    dqs.push_back(amt__bcat_cm);
    dqs.push_back(amt__ligand_cm);
    dqs.push_back(amt__complex_cm);
    dqs.push_back(amt__bcat_nu);
    dqs.push_back(amt__ligand_nu);
    dqs.push_back(amt__complex_nu);
    dqs.push_back(amt__drag);

    return dqs;
}

void Tan2014SbmlOdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    std::vector<double> derivatives = RunModelEquations(time, rY);
    for (unsigned i = 0; i < rDY.size(); ++i)
    {
        rDY[i] = derivatives[i];
    }

    // TODO: Scale time appropriately
}

void Tan2014SbmlOdeSystem::Initialise(double time)
{
    compartment = 1.0;                                                                                          //
    CytosolMembrane = 1.16;                                                                                     //
    nucleus = 0.65;                                                                                             //
    bcat_cm = 46.6;                                                                                             //
    ligand_cm = 581.1;                                                                                          //
    complex_cm = 418.9;                                                                                         //
    bcat_nu = 32.6;                                                                                             //
    ligand_nu = 516.8;                                                                                          //
    complex_nu = 483.2;                                                                                         //
    drag = 1.0;                                                                                                 //
    wnt_level = 0.0;                                                                                            //
    kdegradation = 0.0163;                                                                                      //
    Bsyn = 1.306;                                                                                               //
    K_n_active_k = 17.16;                                                                                       //
    kC_k1 = 1e-5;                                                                                               //
    kC_k2 = 0.000647;                                                                                           //
    kN_k1 = 0.0001;                                                                                             //
    kN_k2 = 0.00349;                                                                                            //
    kdiffusion_k = 39.13;                                                                                       //
    K_c_active_k = 4.5;                                                                                         //
    gamma = 1.0;                                                                                                //
    ComplexTransitThreshold = 1.0;                                                                              //
    drag = sm::piecewise((complex_cm - 700.0) / 10.0, ((complex_cm - 700.0) / 10.0) > 1.0, 1.0);                //
    Bsynthesis = Bsyn * CytosolMembrane;                                                                        //
    kDegradation = CytosolMembrane * kdegradation * gamma * bcat_cm * (1.0 - 0.5 * wnt_level);                  //
    kC = CytosolMembrane * (kC_k1 * bcat_cm * ligand_cm - kC_k2 * complex_cm);                                  //
    kN = nucleus * (kN_k1 * bcat_nu * ligand_nu - kN_k2 * complex_nu);                                          //
    kdiffusion = kdiffusion_k * (bcat_cm - bcat_nu);                                                            //
    K_c_active = K_c_active_k * bcat_cm;                                                                        //
    K_n_active = K_n_active_k * bcat_nu;                                                                        //
    d_bcat_cm_dt = ((Bsynthesis - kDegradation - kC - kdiffusion - K_c_active) + K_n_active) / CytosolMembrane; //
    d_ligand_cm_dt = -kC / CytosolMembrane;                                                                     //
    d_complex_cm_dt = kC / CytosolMembrane;                                                                     //
    d_bcat_nu_dt = (-kN + kdiffusion + K_c_active - K_n_active) / nucleus;                                      //
    d_ligand_nu_dt = -kN / nucleus;                                                                             //
    d_complex_nu_dt = kN / nucleus;                                                                             //

    mStateVariables.push_back(bcat_cm);
    mStateVariables.push_back(ligand_cm);
    mStateVariables.push_back(complex_cm);
    mStateVariables.push_back(bcat_nu);
    mStateVariables.push_back(ligand_nu);
    mStateVariables.push_back(complex_nu);

    SetDefaultInitialCondition(0, bcat_cm);
    SetDefaultInitialCondition(1, ligand_cm);
    SetDefaultInitialCondition(2, complex_cm);
    SetDefaultInitialCondition(3, bcat_nu);
    SetDefaultInitialCondition(4, ligand_nu);
    SetDefaultInitialCondition(5, complex_nu);

    mParameters.push_back(wnt_level);
    mParameters.push_back(kdegradation);
    mParameters.push_back(Bsyn);
    mParameters.push_back(K_n_active_k);
    mParameters.push_back(kC_k1);
    mParameters.push_back(kC_k2);
    mParameters.push_back(kN_k1);
    mParameters.push_back(kN_k2);
    mParameters.push_back(kdiffusion_k);
    mParameters.push_back(K_c_active_k);
    mParameters.push_back(gamma);
    mParameters.push_back(ComplexTransitThreshold);
}

double Tan2014SbmlOdeSystem::ProcessModelEvents(double time, const std::vector<double>& rY)
{
    // Ensure all member variables (state vars, parameters, derived quantities) reflect
    // the rY passed in. Without this, event triggers and assignments would use stale
    // values from the last EvaluateYDerivatives call, which may differ from rY when
    // called from CalculateRootFunction or CalculateStoppingEvent with a different state.
    RunModelEquations(time, rY);

    // Do NOT clear mEventAdjustedStateVars/Parameters here. Once set by an event fire,
    // they must persist across all CVODE bisection calls until AdjustParameters applies
    // them. Clearing here would erase the stored assignment when a later bisection call
    // lands in the clamped state (mEventSatisfied=true), causing the halving to be lost.
    // CalculateStoppingEvent (BackwardEuler path) clears these itself before calling.

    // Root function for CVODE: the maximum signed event distance, where each distance is
    // positive exactly when its event's trigger condition holds. Taking the MAXIMUM (not the
    // minimum absolute value) means the combined function crosses zero the moment ANY event
    // becomes triggered, and cannot be masked by another event that happens to sit just below
    // its own boundary (a small negative distance). A min-abs combination misses an event
    // whose rising edge coincides with another event re-arming near its threshold.
    double max_dist = -std::numeric_limits<double>::max();

    return max_dist; // Signed distance of the event closest to triggering
}

std::vector<double> Tan2014SbmlOdeSystem::RunModelEquations(double time, const std::vector<double>& rStateVariables)
{
    bcat_cm = rStateVariables[0];
    ligand_cm = rStateVariables[1];
    complex_cm = rStateVariables[2];
    bcat_nu = rStateVariables[3];
    ligand_nu = rStateVariables[4];
    complex_nu = rStateVariables[5];

    wnt_level = GetParameter(0);
    kdegradation = GetParameter(1);
    Bsyn = GetParameter(2);
    K_n_active_k = GetParameter(3);
    kC_k1 = GetParameter(4);
    kC_k2 = GetParameter(5);
    kN_k1 = GetParameter(6);
    kN_k2 = GetParameter(7);
    kdiffusion_k = GetParameter(8);
    K_c_active_k = GetParameter(9);
    gamma = GetParameter(10);
    ComplexTransitThreshold = GetParameter(11);

    drag = sm::piecewise((complex_cm - 700.0) / 10.0, ((complex_cm - 700.0) / 10.0) > 1.0, 1.0);                //
    Bsynthesis = Bsyn * CytosolMembrane;                                                                        //
    kDegradation = CytosolMembrane * kdegradation * gamma * bcat_cm * (1.0 - 0.5 * wnt_level);                  //
    kC = CytosolMembrane * (kC_k1 * bcat_cm * ligand_cm - kC_k2 * complex_cm);                                  //
    kN = nucleus * (kN_k1 * bcat_nu * ligand_nu - kN_k2 * complex_nu);                                          //
    kdiffusion = kdiffusion_k * (bcat_cm - bcat_nu);                                                            //
    K_c_active = K_c_active_k * bcat_cm;                                                                        //
    K_n_active = K_n_active_k * bcat_nu;                                                                        //
    d_bcat_cm_dt = ((Bsynthesis - kDegradation - kC - kdiffusion - K_c_active) + K_n_active) / CytosolMembrane; //
    d_ligand_cm_dt = -kC / CytosolMembrane;                                                                     //
    d_complex_cm_dt = kC / CytosolMembrane;                                                                     //
    d_bcat_nu_dt = (-kN + kdiffusion + K_c_active - K_n_active) / nucleus;                                      //
    d_ligand_nu_dt = -kN / nucleus;                                                                             //
    d_complex_nu_dt = kN / nucleus;                                                                             //

    std::vector<double> derivatives(6);
    derivatives[0] = d_bcat_cm_dt;
    derivatives[1] = d_ligand_cm_dt;
    derivatives[2] = d_complex_cm_dt;
    derivatives[3] = d_bcat_nu_dt;
    derivatives[4] = d_ligand_nu_dt;
    derivatives[5] = d_complex_nu_dt;
    return derivatives;
}

// MODEL FUNCTIONS

template <>
void CellwiseOdeSystemInformation<Tan2014SbmlOdeSystem>::Initialise()
{
    // STATE VARIABLES
    this->mVariableNames.push_back("bcat_cm");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(46.6);

    this->mVariableNames.push_back("ligand_cm");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(581.1);

    this->mVariableNames.push_back("complex_cm");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(418.9);

    this->mVariableNames.push_back("bcat_nu");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(32.6);

    this->mVariableNames.push_back("ligand_nu");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(516.8);

    this->mVariableNames.push_back("complex_nu");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(483.2);

    // DERIVED QUANTITIES
    this->mDerivedQuantityNames.push_back("compartment");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("CytosolMembrane");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("nucleus");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("drag");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("Bsynthesis");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("kDegradation");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("kC");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("kN");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("kdiffusion");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("K_c_active");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("K_n_active");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("amt__bcat_cm");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("amt__ligand_cm");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("amt__complex_cm");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("amt__bcat_nu");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("amt__ligand_nu");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("amt__complex_nu");
    this->mDerivedQuantityUnits.push_back("non-dim");

    this->mDerivedQuantityNames.push_back("amt__drag");
    this->mDerivedQuantityUnits.push_back("non-dim");

    // PARAMETERS
    this->mParameterNames.push_back("wnt_level");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kdegradation");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Bsyn");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("K_n_active_k");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kC_k1");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kC_k2");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kN_k1");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kN_k2");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("kdiffusion_k");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("K_c_active_k");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("gamma");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ComplexTransitThreshold");
    this->mParameterUnits.push_back("non-dim");

    this->mInitialised = true;
}

// Register the ODE system with Boost serialization
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(Tan2014SbmlOdeSystem)