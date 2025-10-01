#ifndef {{ srn_header_guard }}
#define {{ srn_header_guard }}

#include <boost/serialization/base_object.hpp>

#include "AbstractSbmlSrnModel.hpp"
#include "ChasteSerialization.hpp"

#include "{{ ode_hpp_file }}"

class {{ srn_class_name }} : public AbstractSbmlSrnModel
{
private:
    friend class boost::serialization::access;
    /**
     * Save / load {{ srn_class_name }} archive
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template <class Archive>
    void serialize(Archive &archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractSbmlSrnModel>(*this);
    }

protected:
    /**
     * Protected copy-constructor for use by CreateSrnModel().
     * The only way to copy an instance of a subclass of AbstractSrnModel is
     * by calling CreateSrnModel(), which ensures that the instance is copied 
     * correctly.
     *
     * This copy-constructor helps subclasses of AbstractSrnModel to
     * ensure that all their members are copied over correctly. It is primarily
     * used during cell division to set member variables for a daughter cell.
     * Note that the SRN model of the parent cell will have run ResetForDivision() 
     * just before calling CreateSrnModel(), so performing an exact copy of the
     * parent cell's SRN model is suitable behaviour. Any further initialisation
     * specific to the daughter cell can be completed via InitialiseDaughterCell().
     *
     * @param rModel the SRN model to copy.
     */
    {{ srn_class_name }}(const {{ srn_class_name }}& rModel);

public:
    /**
     * Default constructor.
     *
     * @param pOdeSolver Optional pointer to a CellCycleModelOdeSolver, which allows using different solvers)
     *
     */
    {{ srn_class_name }}(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver = boost::shared_ptr<AbstractCellCycleModelOdeSolver>());

    /**
     * Overridden builder method to create new copies of this SRN model.
     *
     * @return Returns a copy of the current SRN model.
     */
    {{ srn_class_name }}* CreateSrnModel() override;

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
    void OutputSrnModelParameters(out_stream& rParamsFile) override;
};

// Register the SRN model with Boost serialization
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT({{ srn_class_name }})

// Register the CellCycleModel<OdeSolver> classes with Boost serialization
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER({{ srn_class_name }})

#endif // {{ srn_header_guard }}
