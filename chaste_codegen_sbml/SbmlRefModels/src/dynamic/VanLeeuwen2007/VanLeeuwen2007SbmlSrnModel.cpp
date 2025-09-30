#include "VanLeeuwen2007SbmlSrnModel.hpp"

VanLeeuwen2007SbmlSrnModel::VanLeeuwen2007SbmlSrnModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
        : AbstractSbmlSrnModel(11, pOdeSolver)
{
}

VanLeeuwen2007SbmlSrnModel::VanLeeuwen2007SbmlSrnModel(const VanLeeuwen2007SbmlSrnModel& rModel)
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
    SetOdeSystem(new VanLeeuwen2007SbmlOdeSystem(*rModel.GetOdeSystem()));
}

VanLeeuwen2007SbmlSrnModel* VanLeeuwen2007SbmlSrnModel::CreateSrnModel()
{
    return new VanLeeuwen2007SbmlSrnModel(*this);
}

void VanLeeuwen2007SbmlSrnModel::Initialise()
{
    AbstractSbmlSrnModel::Initialise(new VanLeeuwen2007SbmlOdeSystem);
}

VanLeeuwen2007SbmlOdeSystem* VanLeeuwen2007SbmlSrnModel::GetOdeSystem() const
{
    return static_cast<VanLeeuwen2007SbmlOdeSystem*>(mpOdeSystem);
}

void VanLeeuwen2007SbmlSrnModel::OutputSrnModelParameters(out_stream& rParamsFile)
{
    AbstractSbmlSrnModel::OutputSrnModelParameters(rParamsFile);
}

// Register the SRN model with Boost serialization
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(VanLeeuwen2007SbmlSrnModel)

// Register the CellCycleModel<OdeSolver> classes with Boost serialization
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(VanLeeuwen2007SbmlSrnModel)