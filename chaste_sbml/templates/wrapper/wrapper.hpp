#ifndef {{ wrapper_header_guard }}
#define {{ wrapper_header_guard }}

#include <boost/serialization/base_object.hpp>

#include "{{ wrapper_base_class }}.hpp"
#include "ChasteSerialization.hpp"

#include "{{ ode_hpp_file }}"

class {{ wrapper_class_name }} : public {{ wrapper_base_class }}
{
private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Save / load {{ wrapper_class_name }} archive
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template <class Archive>
    void serialize(Archive &archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<{{ wrapper_base_class }}>(*this);
    }

protected:
    /**
     * Protected copy-constructor for use by {{ wrapper_create_method }}().
     * The only way to copy an instance of a subclass of {{ wrapper_abstract_class }} is
     * by calling {{ wrapper_create_method }}(), which ensures that the instance is copied
     * correctly.
     *
     * This copy-constructor helps subclasses of {{ wrapper_abstract_class }} to
     * ensure that all their members are copied over correctly. It is primarily
     * used during cell division to set member variables for a daughter cell.
     * Note that the {{ wrapper_model_noun }} of the parent cell will have run ResetForDivision()
     * just before calling {{ wrapper_create_method }}(), so performing an exact copy of the
     * parent cell's {{ wrapper_model_noun }} is suitable behaviour. Any further initialisation
     * specific to the daughter cell can be completed via InitialiseDaughterCell().
     *
     * @param rModel the {{ wrapper_model_noun }} to copy.
     */
    {{ wrapper_class_name }}(const {{ wrapper_class_name }}& rModel);

public:
    /**
     * Default constructor.
     *
     * @param pOdeSolver Optional pointer to a CellCycleModelOdeSolver, which allows using different solvers)
     *
     */
    {{ wrapper_class_name }}(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver = boost::shared_ptr<AbstractCellCycleModelOdeSolver>());

    /**
     * Overridden builder method to create new copies of this {{ wrapper_model_noun }}.
     *
     * @return Returns a copy of the current {{ wrapper_model_noun }}.
     */
    {{ wrapper_abstract_class }}* {{ wrapper_create_method }}() override;

    /**
     * Initialise the {{ wrapper_model_noun }} at the start of a simulation.
     *
     * This overridden method sets up a new ODE system.
     */
    void Initialise() override;

    /**
     * Output {{ wrapper_model_noun }} parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void {{ wrapper_output_method }}(out_stream& rParamsFile) override;
};

// Register the {{ wrapper_model_noun }} with Boost serialization
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT({{ wrapper_class_name }})

// Register the CellCycleModel<OdeSolver> classes with Boost serialization
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER({{ wrapper_class_name }})

#endif // {{ wrapper_header_guard }}
