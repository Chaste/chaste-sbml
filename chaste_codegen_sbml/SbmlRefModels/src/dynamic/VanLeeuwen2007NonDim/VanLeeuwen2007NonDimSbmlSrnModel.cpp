#include "VanLeeuwen2007NonDimSbmlSrnModel.hpp"

VanLeeuwen2007NonDimSbmlSrnModel::VanLeeuwen2007NonDimSbmlSrnModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
        : AbstractSbmlSrnModel(11, pOdeSolver)
{
}

VanLeeuwen2007NonDimSbmlSrnModel::VanLeeuwen2007NonDimSbmlSrnModel(const VanLeeuwen2007NonDimSbmlSrnModel& rModel)
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
    SetOdeSystem(new VanLeeuwen2007NonDimSbmlOdeSystem(*rModel.GetOdeSystem()));
}

VanLeeuwen2007NonDimSbmlSrnModel* VanLeeuwen2007NonDimSbmlSrnModel::CreateSrnModel()
{
    return new VanLeeuwen2007NonDimSbmlSrnModel(*this);
}

void VanLeeuwen2007NonDimSbmlSrnModel::Initialise()
{
    AbstractSbmlSrnModel::Initialise(new VanLeeuwen2007NonDimSbmlOdeSystem);
}

VanLeeuwen2007NonDimSbmlOdeSystem* VanLeeuwen2007NonDimSbmlSrnModel::GetOdeSystem() const
{
    return static_cast<VanLeeuwen2007NonDimSbmlOdeSystem*>(mpOdeSystem);
}

void VanLeeuwen2007NonDimSbmlSrnModel::OutputSrnModelParameters(out_stream& rParamsFile)
{
    AbstractSbmlSrnModel::OutputSrnModelParameters(rParamsFile);
}

// Register the SRN model with Boost serialization
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(VanLeeuwen2007NonDimSbmlSrnModel)

// Register the CellCycleModel<OdeSolver> classes with Boost serialization
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(VanLeeuwen2007NonDimSbmlSrnModel)