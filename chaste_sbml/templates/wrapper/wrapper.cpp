#include "{{ wrapper_hpp_file }}"

{{ wrapper_class_name }}::{{ wrapper_class_name }}(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
    : {{ wrapper_base_class }}({% if model_type == ModelType.SRN %}{{ state_variables|length }}, {% endif %}pOdeSolver)
{
}

{{ wrapper_class_name }}::{{ wrapper_class_name }}(const {{ wrapper_class_name }}& rModel)
        : {{ wrapper_base_class }}(rModel)
{
    /*
     * Set each member variable of the new {{ wrapper_model_noun }} that inherits
     * its value from the parent.
     *
     * Note 1: some of the new {{ wrapper_model_noun }}'s member variables
     * will already have been correctly initialized in its constructor.
     *
     * Note 2: one or more of the new {{ wrapper_model_noun }}'s member variables
     * may be set/overwritten as soon as InitialiseDaughterCell() is called on
     * the new {{ wrapper_model_noun }}.
     *
     * Note 3: Only set the variables defined in this class. Variables defined
     * in parent classes will be defined there.
     */
    {{ ode_class_name }}* p_ode_system = static_cast<{{ ode_class_name }}*>(rModel.GetOdeSystem());
    assert(p_ode_system != nullptr);
    this->SetOdeSystem(new {{ ode_class_name }}(*p_ode_system));
}

{{ wrapper_abstract_class }}* {{ wrapper_class_name }}::{{ wrapper_create_method }}()
{
    return new {{ wrapper_class_name }}(*this);
}

void {{ wrapper_class_name }}::Initialise()
{
    assert(mpOdeSystem == nullptr);
{% if model_type == ModelType.SRN %}
    {{ wrapper_base_class }}::Initialise(new {{ ode_class_name }});
{% else %}
    mpOdeSystem = new {{ ode_class_name }};
    {{ wrapper_base_class }}::Initialise();
{% endif %}
}

void {{ wrapper_class_name }}::{{ wrapper_output_method }}(out_stream& rParamsFile)
{
    {{ wrapper_base_class }}::{{ wrapper_output_method }}(rParamsFile);
}

// Register the {{ wrapper_model_noun }} with Boost serialization
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT({{ wrapper_class_name }})

// Register the CellCycleModel<OdeSolver> classes with Boost serialization
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER({{ wrapper_class_name }})
