#include "{{ cell_cycle_hpp_file }}"

{{ cell_cycle_class_name }}::{{ cell_cycle_class_name }}(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
    : AbstractSbmlCellCycleModel(pOdeSolver)
{
}

{{ cell_cycle_class_name }}::{{ cell_cycle_class_name }}(const {{ cell_cycle_class_name }}& rModel)
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
    {{ ode_class_name }}* p_ode_system = static_cast<{{ ode_class_name }}*>(rModel.GetOdeSystem());
    assert(p_ode_system != nullptr);
    this->SetOdeSystem(new {{ ode_class_name }}(*p_ode_system));
}

{{ cell_cycle_class_name }}* {{ cell_cycle_class_name }}::CreateCellCycleModel()
{
    return new {{ cell_cycle_class_name }}(*this);
}

void {{ cell_cycle_class_name }}::Initialise()
{
    assert(mpOdeSystem == nullptr);
    mpOdeSystem = new {{ ode_class_name }};
    AbstractSbmlCellCycleModel::Initialise();
}

void {{ cell_cycle_class_name }}::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    AbstractSbmlCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Register the cell-cycle model with Boost serialization
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT({{ cell_cycle_class_name }})

// Register the CellCycleModel<OdeSolver> classes with Boost serialization
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER({{ cell_cycle_class_name }})
