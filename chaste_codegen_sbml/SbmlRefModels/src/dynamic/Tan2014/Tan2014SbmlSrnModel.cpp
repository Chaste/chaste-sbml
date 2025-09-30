#include "Tan2014SbmlSrnModel.hpp"

Tan2014SbmlSrnModel::Tan2014SbmlSrnModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
        : AbstractSbmlSrnModel(6, pOdeSolver)
{
}

Tan2014SbmlSrnModel::Tan2014SbmlSrnModel(const Tan2014SbmlSrnModel& rModel)
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
    assert(rModel.GetOdeSystem());
    SetOdeSystem(new Tan2014SbmlOdeSystem(*rModel.GetOdeSystem()));
}

Tan2014SbmlSrnModel* Tan2014SbmlSrnModel::CreateSrnModel()
{
    return new Tan2014SbmlSrnModel(*this);
}

void Tan2014SbmlSrnModel::Initialise()
{
    AbstractSbmlSrnModel::Initialise(new Tan2014SbmlOdeSystem);
}

Tan2014SbmlOdeSystem* Tan2014SbmlSrnModel::GetOdeSystem() const
{
    return static_cast<Tan2014SbmlOdeSystem*>(mpOdeSystem);
}

void Tan2014SbmlSrnModel::OutputSrnModelParameters(out_stream& rParamsFile)
{
    AbstractSbmlSrnModel::OutputSrnModelParameters(rParamsFile);
}

// Register the SRN model with Boost serialization
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(Tan2014SbmlSrnModel)

// Register the CellCycleModel<OdeSolver> classes with Boost serialization
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(Tan2014SbmlSrnModel)