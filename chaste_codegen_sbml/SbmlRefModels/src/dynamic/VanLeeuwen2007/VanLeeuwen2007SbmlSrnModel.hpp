#ifndef VANLEEUWEN2007SBMLSRNMODEL_HPP_
#define VANLEEUWEN2007SBMLSRNMODEL_HPP_

#include <boost/serialization/base_object.hpp>

#include "AbstractSbmlSrnModel.hpp"
#include "ChasteSerialization.hpp"

#include "VanLeeuwen2007SbmlOdeSystem.hpp"

class VanLeeuwen2007SbmlSrnModel : public AbstractSbmlSrnModel
{
private:
    friend class boost::serialization::access;
    /**
     * Save / load VanLeeuwen2007SbmlSrnModel archive
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive& boost::serialization::base_object<AbstractSbmlSrnModel>(*this);
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
    VanLeeuwen2007SbmlSrnModel(const VanLeeuwen2007SbmlSrnModel& rModel);

public:
    /**
     * Default constructor.
     *
     * @param pOdeSolver Optional pointer to a CellCycleModelOdeSolver, which allows using different solvers)
     *
     */
    VanLeeuwen2007SbmlSrnModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver = boost::shared_ptr<AbstractCellCycleModelOdeSolver>());

    /**
     * Overridden builder method to create new copies of this SRN model.
     *
     * @return Returns a copy of the current SRN model.
     */
    VanLeeuwen2007SbmlSrnModel* CreateSrnModel();

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
CHASTE_CLASS_EXPORT(VanLeeuwen2007SbmlSrnModel)

// Register the CellCycleModel<OdeSolver> classes with Boost serialization
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(VanLeeuwen2007SbmlSrnModel)

#endif // VANLEEUWEN2007SBMLSRNMODEL_HPP_