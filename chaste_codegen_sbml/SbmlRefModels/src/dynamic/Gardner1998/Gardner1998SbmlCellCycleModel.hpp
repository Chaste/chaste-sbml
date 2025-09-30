#ifndef GARDNER1998SBMLCELLCYCLEMODEL_HPP_
#define GARDNER1998SBMLCELLCYCLEMODEL_HPP_

#include <boost/serialization/base_object.hpp>

#include "AbstractSbmlCellCycleModel.hpp"
#include "ChasteSerialization.hpp"

#include "Gardner1998SbmlOdeSystem.hpp"

class Gardner1998SbmlCellCycleModel : public AbstractSbmlCellCycleModel
{
private:
    friend class boost::serialization::access;
    /**
     * Save / load Gardner1998SbmlCellCycleModel archive
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive& BOOST_SERIALIZATION_BASE_OBJECT_NVP(AbstractSbmlCellCycleModel);
    }

protected:
    /**
     * Protected copy-constructor for use by CreateCellCycleModel().
     * 
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
    Gardner1998SbmlCellCycleModel(const Gardner1998SbmlCellCycleModel& rModel);

public:
    /**
     * Default constructor.
     *
     * @param pOdeSolver Optional pointer to a CellCycleModelOdeSolver, which allows using different solvers)
     *
     */
    Gardner1998SbmlCellCycleModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver = boost::shared_ptr<AbstractCellCycleModelOdeSolver>());

    /**
     * Overridden builder method to create new copies of this SRN model.
     *
     * @return Returns a copy of the current SRN model.
     */
    Gardner1998SbmlCellCycleModel* CreateCellCycleModel();

    /**
     * @return #mpOdeSystem.
     */
    Gardner1998SbmlOdeSystem* GetOdeSystem() const;

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

// Register the SRN model with Boost serialization
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(Gardner1998SbmlCellCycleModel)

// Register the CellCycleModel<OdeSolver> classes with Boost serialization
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(Gardner1998SbmlCellCycleModel)

#endif // GARDNER1998SBMLCELLCYCLEMODEL_HPP_