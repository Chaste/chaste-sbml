#include "Gardner1998SbmlCellCycleModel.hpp"

Gardner1998SbmlCellCycleModel::Gardner1998SbmlCellCycleModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
        : AbstractSbmlCellCycleModel(pOdeSolver)
{
}

Gardner1998SbmlCellCycleModel::Gardner1998SbmlCellCycleModel(const Gardner1998SbmlCellCycleModel& rModel)
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
    Gardner1998SbmlOdeSystem* p_ode_system = static_cast<Gardner1998SbmlOdeSystem*>(rModel.GetOdeSystem());
    assert(p_ode_system != nullptr);
    this->SetOdeSystem(new Gardner1998SbmlOdeSystem(*p_ode_system));
}

Gardner1998SbmlCellCycleModel* Gardner1998SbmlCellCycleModel::CreateCellCycleModel()
{
    return new Gardner1998SbmlCellCycleModel(*this);
}

void Gardner1998SbmlCellCycleModel::Initialise()
{
    assert(mpOdeSystem == nullptr);
    mpOdeSystem = new Gardner1998SbmlOdeSystem;
    AbstractSbmlCellCycleModel::Initialise();
}

void Gardner1998SbmlCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    AbstractSbmlCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Register the cell-cycle model with Boost serialization
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(Gardner1998SbmlCellCycleModel)

// Register the CellCycleModel<OdeSolver> classes with Boost serialization
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(Gardner1998SbmlCellCycleModel)