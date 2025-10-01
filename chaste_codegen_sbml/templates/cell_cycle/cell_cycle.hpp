#ifndef {{ cell_cycle_header_guard }}
#define {{ cell_cycle_header_guard }}

#include <boost/serialization/base_object.hpp>

#include "AbstractSbmlCellCycleModel.hpp"
#include "ChasteSerialization.hpp"

#include "{{ ode_hpp_file }}"

class {{ cell_cycle_class_name }} : public AbstractSbmlCellCycleModel
{
private:
    friend class boost::serialization::access;
    /**
     * Save / load {{ cell_cycle_class_name }} archive
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template <class Archive>
    void serialize(Archive &archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractSbmlCellCycleModel>(*this);
    }

protected:
    /**
     * Protected copy-constructor for use by CreateCellCycleModel().
     * The only way to copy an instance of a subclass of AbstractCellCycleModel is
     * by calling CreateCellCycleModel(), which ensures that the instance is copied 
     * correctly.
     *
     * This copy-constructor helps subclasses of AbstractCellCycleModel to
     * ensure that all their members are copied over correctly. It is primarily
     * used during cell division to set member variables for a daughter cell.
     * Note that the cell-cycle model of the parent cell will have run ResetForDivision() 
     * just before calling CreateCellCycleModel(), so performing an exact copy of the
     * parent cell's cell-cycle model is suitable behaviour. Any further initialisation
     * specific to the daughter cell can be completed via InitialiseDaughterCell().
     *
     * @param rModel the cell-cycle model to copy.
     */
    {{ cell_cycle_class_name }}(const {{ cell_cycle_class_name }}& rModel);

public:
    /**
     * Default constructor.
     *
     * @param pOdeSolver Optional pointer to a CellCycleModelOdeSolver, which allows using different solvers)
     *
     */
    {{ cell_cycle_class_name }}(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver = boost::shared_ptr<AbstractCellCycleModelOdeSolver>());

    /**
     * Overridden builder method to create new copies of this SRN model.
     *
     * @return Returns a copy of the current SRN model.
     */
    {{ cell_cycle_class_name }}* CreateCellCycleModel();

    /**
     * @return #mpOdeSystem.
     */
    {{ ode_class_name }}* GetOdeSystem() const;

    /**
     * Initialise the SRN model at the start of a simulation.
     *
     * This overridden method sets up a new ODE system.
     */
    void Initialise() override; // override

    /**
     * Output SRN model parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellCycleModelParameters(out_stream& rParamsFile) override;
};

// Register the cell-cycle model with Boost serialization
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT({{ cell_cycle_class_name }})

// Register the CellCycleModel<OdeSolver> classes with Boost serialization
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER({{ cell_cycle_class_name }})

#endif // {{ cell_cycle_header_guard }}
