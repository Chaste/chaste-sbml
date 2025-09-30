#include "{{ srn_hpp_file }}"

{{ srn_class_name }}::{{ srn_class_name }}(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
    : AbstractSbmlSrnModel({{ state_variables|length }}, pOdeSolver)
{
}

{{ srn_class_name }}::{{ srn_class_name }}(const {{ srn_class_name }}& rModel)
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
    SetOdeSystem(new {{ ode_class_name }}(*rModel.GetOdeSystem()));
}

{{ srn_class_name }}* {{ srn_class_name }}::CreateSrnModel()
{
    return new {{ srn_class_name }}(*this);
}

void {{ srn_class_name }}::Initialise()
{
    assert(mpOdeSystem == nullptr);
    AbstractSbmlSrnModel::Initialise(new {{ ode_class_name }});
}

{{ ode_class_name }}* {{ srn_class_name }}::GetOdeSystem() const
{
    return static_cast<{{ ode_class_name }}*>(mpOdeSystem);
}

void {{ srn_class_name }}::OutputSrnModelParameters(out_stream& rParamsFile)
{
    AbstractSbmlSrnModel::OutputSrnModelParameters(rParamsFile);
}

// Register the SRN model with Boost serialization
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT({{ srn_class_name }})

// Register the CellCycleModel<OdeSolver> classes with Boost serialization
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER({{ srn_class_name }})
