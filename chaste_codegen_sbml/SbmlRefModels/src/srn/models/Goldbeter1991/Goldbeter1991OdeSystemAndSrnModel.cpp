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
    /* Define algebraic rules. */
    double V1 = rY[0] * VM1 * std::pow(rY[0] + Kc, -1);
    double V3 = rY[1] * VM3;

    /* Define the reactions in this model. */
    // creation of cyclin
    double reaction1 = 0.0;
    {
        double vi = 0.025;
        reaction1 = cell * vi;
    }

    // default degradation of cyclin
    double reaction2 = 0.0;
    {
        double kd = 0.01;
        reaction2 = rY[0] * cell * kd;
    }

    // cdc2 kinase triggered degration of cyclin
    double reaction3 = 0.0;
    {
        double vd = 0.25;
        double Kd = 0.02;
        reaction3 = rY[0] * cell * vd * rY[2] * std::pow(rY[0] + Kd, -1);
    }

    // activation of cdc2 kinase
    double reaction4 = 0.0;
    {
        double K1 = 0.005;
        reaction4 = cell * (1 + -1 * rY[1]) * V1 * std::pow(K1 + -1 * rY[1] + 1, -1);
    }

    // deactivation of cdc2 kinase
    double reaction5 = 0.0;
    {
        double V2 = 1.5;
        double K2 = 0.005;
        reaction5 = cell * rY[1] * V2 * std::pow(K2 + rY[1], -1);
    }

    // activation of cyclin protease
    double reaction6 = 0.0;
    {
        double K3 = 0.005;
        reaction6 = cell * V3 * (1 + -1 * rY[2]) * std::pow(K3 + -1 * rY[2] + 1, -1);
    }

    // deactivation of cyclin protease
    double reaction7 = 0.0;
    {
        double K4 = 0.005;
        double V4 = 0.5;
        reaction7 = cell * V4 * rY[2] * std::pow(K4 + rY[2], -1);
    }


    rDY[0] = (reaction1 - reaction2 - reaction3) / cell; // dCyclin/dt
    rDY[1] = (reaction4 - reaction5) / cell; // dcdc_2_kinase/dt
    rDY[2] = (reaction6 - reaction7) / cell; // dCyclin Protease/dt

    /* Account for the differences in timescales. */
}


template <>
void CellwiseOdeSystemInformation<Goldbeter1991OdeSystem>::Initialise()
{
    this->mVariableNames.push_back("C");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.01);

    this->mVariableNames.push_back("M");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.01);

    this->mVariableNames.push_back("X");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.01);


    /* Define state parameters. */
    // Parameters without set values must be externally defined
    this->mInitialised = true;
}

/* Define SbmlSrnWrapperModel using wrappers. */
#include "SbmlSrnWrapperModel.hpp"
#include "SbmlSrnWrapperModel.cpp"

typedef SbmlSrnWrapperModel<Goldbeter1991OdeSystem, 3> Goldbeter1991SrnModel;

// Declare identifiers for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(Goldbeter1991OdeSystem)
EXPORT_TEMPLATE_CLASS2(SbmlSrnWrapperModel, Goldbeter1991OdeSystem, 3)

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(Goldbeter1991SrnModel)