#ifndef CHEN2004SBMLCELLCYCLEMODEL_HPP_
#define CHEN2004SBMLCELLCYCLEMODEL_HPP_

#include <boost/serialization/base_object.hpp>

#include "AbstractSbmlCellCycleModel.hpp"
#include "ChasteSerialization.hpp"

#include "Chen2004SbmlOdeSystem.hpp"

class Chen2004SbmlCellCycleModel : public AbstractSbmlCellCycleModel
{
private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Save / load Chen2004SbmlCellCycleModel archive
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive& boost::serialization::base_object<AbstractSbmlCellCycleModel>(*this);
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
    Chen2004SbmlCellCycleModel(const Chen2004SbmlCellCycleModel& rModel);

public:
    /**
     * Default constructor.
     *
     * @param pOdeSolver Optional pointer to a CellCycleModelOdeSolver, which allows using different solvers)
     *
     */
    Chen2004SbmlCellCycleModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver = boost::shared_ptr<AbstractCellCycleModelOdeSolver>());

    /**
     * Overridden builder method to create new copies of this SRN model.
     *
     * @return Returns a copy of the current SRN model.
     */
    AbstractCellCycleModel* CreateCellCycleModel() override;

    /**
     * Initialise the SRN model at the start of a simulation.
     *
     * This overridden method sets up a new ODE system.
     */
    void Initialise() override;

    /**
     * Output SRN model parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellCycleModelParameters(out_stream& rParamsFile) override;
};

// Register the cell-cycle model with Boost serialization
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(Chen2004SbmlCellCycleModel)

// Register the CellCycleModel<OdeSolver> classes with Boost serialization
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(Chen2004SbmlCellCycleModel)

#endif // CHEN2004SBMLCELLCYCLEMODEL_HPP_