#include "Chen2000SbmlCellCycleModel.hpp"

Chen2000SbmlCellCycleModel::Chen2000SbmlCellCycleModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
    : AbstractSbmlCellCycleModel(pOdeSolver)
{
}

Chen2000SbmlCellCycleModel::Chen2000SbmlCellCycleModel(const Chen2000SbmlCellCycleModel& rModel)
        : AbstractSbmlCellCycleModel(rModel)
{
    /*
     * Set each member variable of the new cell-cycle model that inherits
     * its value from the parent.
     *
     * Note 1: some of the new cell-cycle model's member variables
     * will already have been correctly initialized in its constructor.
     *
     * Note 2: one or more of the new cell-cycle model's member variables
     * may be set/overwritten as soon as InitialiseDaughterCell() is called on
     * the new cell-cycle model.
     *
     * Note 3: Only set the variables defined in this class. Variables defined
     * in parent classes will be defined there.
     */
    Chen2000SbmlOdeSystem* p_ode_system = static_cast<Chen2000SbmlOdeSystem*>(rModel.GetOdeSystem());
    assert(p_ode_system != nullptr);
    this->SetOdeSystem(new Chen2000SbmlOdeSystem(*p_ode_system));
}

AbstractCellCycleModel* Chen2000SbmlCellCycleModel::CreateCellCycleModel()
{
    return new Chen2000SbmlCellCycleModel(*this);
}

void Chen2000SbmlCellCycleModel::Initialise()
{
    assert(mpOdeSystem == nullptr);
    mpOdeSystem = new Chen2000SbmlOdeSystem;
    AbstractSbmlCellCycleModel::Initialise();
}

void Chen2000SbmlCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    AbstractSbmlCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Register the cell-cycle model with Boost serialization
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(Chen2000SbmlCellCycleModel)

// Register the CellCycleModel<OdeSolver> classes with Boost serialization
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(Chen2000SbmlCellCycleModel)