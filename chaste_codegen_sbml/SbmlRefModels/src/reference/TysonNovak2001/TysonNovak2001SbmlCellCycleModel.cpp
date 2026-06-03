#include "TysonNovak2001SbmlCellCycleModel.hpp"

TysonNovak2001SbmlCellCycleModel::TysonNovak2001SbmlCellCycleModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
        : AbstractSbmlCellCycleModel(pOdeSolver)
{
}

TysonNovak2001SbmlCellCycleModel::TysonNovak2001SbmlCellCycleModel(const TysonNovak2001SbmlCellCycleModel& rModel)
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
    TysonNovak2001SbmlOdeSystem* p_ode_system = static_cast<TysonNovak2001SbmlOdeSystem*>(rModel.GetOdeSystem());
    assert(p_ode_system != nullptr);
    this->SetOdeSystem(new TysonNovak2001SbmlOdeSystem(*p_ode_system));
}

AbstractCellCycleModel* TysonNovak2001SbmlCellCycleModel::CreateCellCycleModel()
{
    return new TysonNovak2001SbmlCellCycleModel(*this);
}

void TysonNovak2001SbmlCellCycleModel::Initialise()
{
    assert(mpOdeSystem == nullptr);
    mpOdeSystem = new TysonNovak2001SbmlOdeSystem;
    AbstractSbmlCellCycleModel::Initialise();
}

void TysonNovak2001SbmlCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    AbstractSbmlCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Register the cell-cycle model with Boost serialization
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(TysonNovak2001SbmlCellCycleModel)

// Register the CellCycleModel<OdeSolver> classes with Boost serialization
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(TysonNovak2001SbmlCellCycleModel)