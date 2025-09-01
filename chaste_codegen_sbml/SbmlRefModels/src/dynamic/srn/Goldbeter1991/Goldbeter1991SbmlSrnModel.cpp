
#include "AbstractSbmlSrnModel.hpp"
#include "Goldbeter1991SbmlOdeSystem.hpp"

#include "Goldbeter1991SbmlSrnModel.hpp"

Goldbeter1991SbmlSrnModel::Goldbeter1991SbmlSrnModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
        : AbstractSbmlSrnModel(3, pOdeSolver)
{
}

// New method for copy constructor
// Goldbeter1991SbmlSrnModel::Goldbeter1991SbmlSrnModel(const Goldbeter1991SbmlSrnModel& rModel)
//         : AbstractSbmlSrnModel(rModel)
// {
//     /*
//      * Set each member variable of the new SRN model that inherits
//      * its value from the parent.
//      *
//      * Note 1: some of the new SRN model's member variables
//      * will already have been correctly initialized in its constructor.
//      *
//      * Note 2: one or more of the new SRN model's member variables
//      * may be set/overwritten as soon as InitialiseDaughterCell() is called on
//      * the new SRN model.
//      *
//      * Note 3: Only set the variables defined in this class. Variables defined
//      * in parent classes will be defined there.
//      */

//     assert(rModel.GetOdeSystem());
//     // SetOdeSystem(new Goldbeter1991SbmlOdeSystem(*rModel.GetOdeSystem()));
//     // SetOdeSystem(new Goldbeter1991SbmlOdeSystem(static_cast<Goldbeter1991SbmlOdeSystem&>(*rModel.GetOdeSystem())));
//     SetOdeSystem(new Goldbeter1991SbmlOdeSystem(*rModel.GetOdeSystem()));
// }

Goldbeter1991SbmlSrnModel* Goldbeter1991SbmlSrnModel::CreateSrnModel()
{
    return new Goldbeter1991SbmlSrnModel(*this);
}

Goldbeter1991SbmlOdeSystem* Goldbeter1991SbmlSrnModel::GetOdeSystem() const
{
    return static_cast<Goldbeter1991SbmlOdeSystem*>(mpOdeSystem);
}

void Goldbeter1991SbmlSrnModel::Initialise()
{
    AbstractSbmlSrnModel::Initialise(new Goldbeter1991SbmlOdeSystem);
}
