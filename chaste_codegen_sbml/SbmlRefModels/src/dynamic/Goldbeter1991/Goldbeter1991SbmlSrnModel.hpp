#ifndef GOLDBETER1991SBMLSRNMODEL_HPP_
#define GOLDBETER1991SBMLSRNMODEL_HPP_

#include <boost/serialization/base_object.hpp>

#include "AbstractSbmlSrnModel.hpp"
#include "ChasteSerialization.hpp"

#include "Goldbeter1991SbmlOdeSystem.hpp"

class Goldbeter1991SbmlSrnModel : public AbstractSbmlSrnModel
{
private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the cell-cycle model and member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive& BOOST_SERIALIZATION_BASE_OBJECT_NVP(AbstractSbmlSrnModel);
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
     * Note that the SrnModel of the parent cell will have run ResetForDivision() 
     * just before calling CreateSrnModel(), so performing an exact copy of the
     * parent cell's SrnModel is suitable behaviour. Any further initialisation
     * specific to the daughter cell can be completed via InitialiseDaughterCell().
     *
     * @param rModel the SRN model to copy.
     */
    Goldbeter1991SbmlSrnModel(const Goldbeter1991SbmlSrnModel& rModel);

public:
    /**
     * Default constructor.
     *
     * @param pOdeSolver Optional pointer to a CellCycleModelOdeSolver, which allows using different solvers)
     *
     */
    Goldbeter1991SbmlSrnModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver = boost::shared_ptr<AbstractCellCycleModelOdeSolver>());

    /**
     * Overridden builder method to create new copies of this srn model.
     *
     * @return Returns a copy of the current srn model.
     */
    Goldbeter1991SbmlSrnModel* CreateSrnModel();

    /**
     * @return #mpOdeSystem.
     */
    Goldbeter1991SbmlOdeSystem* GetOdeSystem() const;

    /**
     * Initialise the cell-cycle model at the start of a simulation.
     *
     * This overridden method sets up a new Ode system.
     */
    void Initialise() override; // override

    /**
     * Outputs SRN model parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSrnModelParameters(out_stream& rParamsFile) override;
};

// Register the SRN model with Boost serialization
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(Goldbeter1991SbmlSrnModel)

// Register the CellCycleModel<OdeSolver> classes with Boost serialization
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(Goldbeter1991SbmlSrnModel)

#endif // GOLDBETER1991SBMLSRNMODEL_HPP_
