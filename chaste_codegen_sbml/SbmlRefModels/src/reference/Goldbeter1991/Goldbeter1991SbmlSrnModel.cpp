#include "Goldbeter1991SbmlSrnModel.hpp"

Goldbeter1991SbmlSrnModel::Goldbeter1991SbmlSrnModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
    : AbstractSbmlSrnModel(3, pOdeSolver)
{
}

Goldbeter1991SbmlSrnModel::Goldbeter1991SbmlSrnModel(const Goldbeter1991SbmlSrnModel& rModel)
        : AbstractSbmlSrnModel(rModel)
{
    /*
     * Set each member variable of the new SRN model that inherits
     * its value from the parent.
     *
     * Note 1: some of the new SRN model's member variables
     * will already have been correctly initialized in its constructor.
     *
     * Note 2: one or more of the new SRN model's member variables
     * may be set/overwritten as soon as InitialiseDaughterCell() is called on
     * the new SRN model.
     *
     * Note 3: Only set the variables defined in this class. Variables defined
     * in parent classes will be defined there.
     */
    Goldbeter1991SbmlOdeSystem* p_ode_system = static_cast<Goldbeter1991SbmlOdeSystem*>(rModel.GetOdeSystem());
    assert(p_ode_system != nullptr);
    this->SetOdeSystem(new Goldbeter1991SbmlOdeSystem(*p_ode_system));
}

AbstractSrnModel* Goldbeter1991SbmlSrnModel::CreateSrnModel()
{
    return new Goldbeter1991SbmlSrnModel(*this);
}

void Goldbeter1991SbmlSrnModel::Initialise()
{
    assert(mpOdeSystem == nullptr);
    AbstractSbmlSrnModel::Initialise(new Goldbeter1991SbmlOdeSystem);
}

void Goldbeter1991SbmlSrnModel::OutputSrnModelParameters(out_stream& rParamsFile)
{
    AbstractSbmlSrnModel::OutputSrnModelParameters(rParamsFile);
}

// Register the SRN model with Boost serialization
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(Goldbeter1991SbmlSrnModel)

// Register the CellCycleModel<OdeSolver> classes with Boost serialization
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(Goldbeter1991SbmlSrnModel)