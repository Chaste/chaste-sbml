#ifndef GOLDBETER_1991_SBML_SRN_MODEL_HPP_
#define GOLDBETER_1991_SBML_SRN_MODEL_HPP_

#include <boost/serialization/base_object.hpp>

#include "AbstractSbmlSrnModel.hpp"
#include "ChasteSerialization.hpp"

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
        archive& boost::serialization::base_object<AbstractSbmlSrnModel>(*this);
    }

protected:
    /**
     * Protected copy-constructor for use by CreateSrnModel.  The only way for external code to create a copy of a SRN model
     * is by calling that method, to ensure that a model of the correct subclass is created.
     * This copy-constructor helps subclasses to ensure that all member variables are correctly copied when this happens.
     *
     * This method is called by child classes to set member variables for a daughter cell upon cell division.
     * Note that the parent SRN model will have had ResetForDivision() called just before CreateSrnModel() is called,
     * so performing an exact copy of the parent is suitable behaviour. Any daughter-cell-specific initialisation
     * can be done in InitialiseDaughterCell().
     *
     * @param rModel the SRN model to copy.
     */
    Goldbeter1991SbmlSrnModel(const Goldbeter1991SbmlSrnModel& rModel) = default;

public:
    /**
     * Default constructor.
     *
     * @param pOdeSolver Optional pointer to a CellCycleModelOdeSolver, which allows using different solvers)
     * .
     */
    Goldbeter1991SbmlSrnModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver = boost::shared_ptr<AbstractCellCycleModelOdeSolver>());

    /**
     * Overridden builder method to create new copies of
     * this srn model.
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
    void Initialise(); // override
};

#endif // GOLDBETER_1991_SBML_SRN_MODEL_HPP_
