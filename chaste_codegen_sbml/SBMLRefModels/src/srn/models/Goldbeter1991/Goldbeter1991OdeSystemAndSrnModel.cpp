#include "CellwiseOdeSystemInformation.hpp"
#include "SbmlMath.hpp"

#include "Goldbeter1991OdeSystemAndSrnModel.hpp"

/* SBML ODE System */
Goldbeter1991OdeSystem::Goldbeter1991OdeSystem(std::vector<double> stateVariables)
    : AbstractOdeSystem(3)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<Goldbeter1991OdeSystem>);

    Init();

    SetDefaultInitialCondition(0, 0.01); // Cyclin
    SetDefaultInitialCondition(1, 0.01); // cdc_2_kinase
    SetDefaultInitialCondition(2, 0.01); // Cyclin Protease


    if (stateVariables != std::vector<double>())
    {
        SetStateVariables(stateVariables);
    }
}

Goldbeter1991OdeSystem::~Goldbeter1991OdeSystem()
{
}


void Goldbeter1991OdeSystem::Init()
{
    /* Initialise model compartments. */
    cell = 1.0;

    /* Initialise model parameters. */
    V1 = 0.0;
    V3 = 0.0;
    VM1 = 3.0;
    VM3 = 1.0;
    Kc = 0.5;

}

void Goldbeter1991OdeSystem::EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY)
{
    /* Define state variables */
    double C = rY[0]; // Cyclin
    double M = rY[1]; // cdc_2_kinase
    double X = rY[2]; // Cyclin Protease

    /* Define state parameters. */


    /* Define algebraic rules. */
    V1 = this->GetStateVariable("C") * VM1 * pow(this->GetStateVariable("C") + Kc, -1);
    V3 = this->GetStateVariable("M") * VM3;

    /* Define the reactions in this model. */
    // creation of cyclin
    double vi = 0.025; // vi
    double reaction1 = cell * vi;

    // default degradation of cyclin
<<<<<<< Updated upstream
    double kd = 0.01; // kd
    double reaction2 = C * cell * kd;

    // cdc2 kinase triggered degration of cyclin
    double vd = 0.25; // vd
    double Kd = 0.02; // Kd
    double reaction3 = C * cell * vd * X * pow(C + Kd, -1);

    // activation of cdc2 kinase
    double K1 = 0.005; // K1
    double reaction4 = cell * (1 + -1 * M) * V1 * pow(K1 + -1 * M + 1, -1);

    // deactivation of cdc2 kinase
    double V2 = 1.5; // V2
    double K2 = 0.005; // K2
    double reaction5 = cell * M * V2 * pow(K2 + M, -1);

    // activation of cyclin protease
    double K3 = 0.005; // K3
    double reaction6 = cell * V3 * (1 + -1 * X) * pow(K3 + -1 * X + 1, -1);

    // deactivation of cyclin protease
    double K4 = 0.005; // K4
    double V4 = 0.5; // V4
    double reaction7 = cell * V4 * X * pow(K4 + X, -1);
=======
    double reaction2 = 0.0;
    {
        double kd = 0.01;
        reaction2 = this->GetStateVariable("C") * cell * kd;
    }

    // cdc2 kinase triggered degration of cyclin
    double reaction3 = 0.0;
    {
        double vd = 0.25;
        double Kd = 0.02;
        reaction3 = this->GetStateVariable("C") * cell * vd * this->GetStateVariable("X") * pow(this->GetStateVariable("C") + Kd, -1);
    }

    // activation of cdc2 kinase
    double reaction4 = 0.0;
    {
        double K1 = 0.005;
        reaction4 = cell * (1 + -1 * this->GetStateVariable("M")) * V1 * pow(K1 + -1 * this->GetStateVariable("M") + 1, -1);
    }

    // deactivation of cdc2 kinase
    double reaction5 = 0.0;
    {
        double V2 = 1.5;
        double K2 = 0.005;
        reaction5 = cell * this->GetStateVariable("M") * V2 * pow(K2 + this->GetStateVariable("M"), -1);
    }

    // activation of cyclin protease
    double reaction6 = 0.0;
    {
        double K3 = 0.005;
        reaction6 = cell * V3 * (1 + -1 * this->GetStateVariable("X")) * pow(K3 + -1 * this->GetStateVariable("X") + 1, -1);
    }

    // deactivation of cyclin protease
    double reaction7 = 0.0;
    {
        double K4 = 0.005;
        double V4 = 0.5;
        reaction7 = cell * V4 * this->GetStateVariable("X") * pow(K4 + this->GetStateVariable("X"), -1);
    }
>>>>>>> Stashed changes


    rDY[0] = (reaction1 - reaction2 - reaction3) / cell; // dCyclin/dt
    rDY[1] = (reaction4 - reaction5) / cell; // dcdc_2_kinase/dt
    rDY[2] = (reaction6 - reaction7) / cell; // dCyclin Protease/dt

    /* Account for the differences in timescales. */
}


template <>
void CellwiseOdeSystemInformation<Goldbeter1991OdeSystem>::Initialise()
{
    this->mVariableNames.push_back("Cyclin");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.01);

    this->mVariableNames.push_back("cdc_2_kinase");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.01);

    this->mVariableNames.push_back("Cyclin Protease");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.01);


    /* Define state parameters. */
    // Parameters without set values must be externally defined
    this->mInitialised = true;
}

/* Define SRN model using Wrappers. */
#include "SbmlSrnWrapperModel.hpp"
#include "SbmlSrnWrapperModel.cpp"

typedef SbmlSrnWrapperModel<Goldbeter1991OdeSystem, 3> Goldbeter1991SrnModel;

// Declare identifiers for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(Goldbeter1991OdeSystem)
EXPORT_TEMPLATE_CLASS2(SbmlSrnWrapperModel, Goldbeter1991OdeSystem, 3)

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(Goldbeter1991SrnModel)