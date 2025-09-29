/*

Copyright (c) 2005-2025, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef ABSTRACT_SBML_SRN_MODEL_HPP_
#define ABSTRACT_SBML_SRN_MODEL_HPP_

#include <boost/serialization/base_object.hpp>

#include "AbstractOdeSrnModel.hpp"
#include "AbstractSbmlOdeSystem.hpp"
#include "ChasteSerialization.hpp"

/**
 * A base class for SRN models generated from SBML
 */

class AbstractSbmlSrnModel : public AbstractOdeSrnModel
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
        archive& BOOST_SERIALIZATION_BASE_OBJECT_NVP(AbstractOdeSrnModel);
    }

public:
    /**
     * Default constructor calls base class.
     *
     * @param pOdeSolver An optional pointer to a cell-cycle model ODE solver object (allows the use of different ODE solvers)
     */
    AbstractSbmlSrnModel(unsigned stateSize, boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver = boost::shared_ptr<AbstractCellCycleModelOdeSolver>());

    virtual ~AbstractSbmlSrnModel() = default;

    // void Initialise() override;

    using AbstractOdeSrnModel::Initialise;
    // /**
    //  * Overridden Initialise() method to set up the ODE system.
    //  *
    //  * @param pOdeSystem pointer to an ODE system
    //  */
    void Initialise(AbstractSbmlOdeSystem* pOdeSystem);

    /**
     * Outputs cell-cycle model parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSrnModelParameters(out_stream& rParamsFile);

    /**
     * Overridden SimulateToCurrentTime() method for custom behaviour
     */
    void SimulateToCurrentTime();

    /**
     * @return the value of a given state variable.
     *
     * @param rName the name of the state variable
     */
    double GetStateVariable(const std::string& rName);
};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(AbstractSbmlSrnModel)

// Provide a unique identifier for serialization
CLASS_IS_ABSTRACT(AbstractSbmlSrnModel)

#endif // ABSTRACT_SBML_SRN_MODEL_HPP_
