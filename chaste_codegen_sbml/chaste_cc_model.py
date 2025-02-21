import os
import shutil
import pathlib

from libsbml import SBMLReader

from . import filewriters

ODESYSTEMANDCELLCYCLEMODEL_HPP_

class ChasteSRNModel(ChasteModel):

    def __init__(self, sbml_file: str, model_name: str = None) -> None:
        pass

    def generate_hpp(self):
        """Construct the Chaste header file from the SBML file."""

        hpp = []

        # Add header guards
        guard = self.model_name.upper() + "ODESYSTEMANDSRNMODEL_HPP_"
        hpp.append(f"#ifndef {guard}")
        hpp.append(f"#define {guard}")
        hpp.append("")

        # Add includes
        hpp.append('''
#include ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>"
#include <boost/serialization/shared_ptr.hpp>"
                   
#include <cmath>"
#include <iostream>"
                   
#include AbstractOdeSystem.hpp"
''')

        # Define the ODE System class
        compartments = self.get_compartments()
        compartment_decl = "\n    ".join(
            [f"double {c};"] for c in compartments
        )

        parameters = self.get_parameters()
        parameter_decl = "\n    ".join(
            [f"double {p};"] for p in parameters
        )

        ode_class_def = f"""
class {self.model_name}OdeSystem : public AbstractOdeSystem
{{
private:
    // Compartments
    {compartment_decl}

    // Parameters
    {parameter_decl}

    // Vector to check if all events have been triggered
    std::vector<bool> eventsSatisfied;
 
public:


    bool CalculateStoppingEvent(double time, const std::vector<double>& rY);
    void CheckAndUpdateEvents(double time, const std::vector<double>& rY);
    bool AreAllEventsSatisfied(double time, const std::vector<double>& rY);



    // Define Cell Cycle model using Wrappers.
    #include "SbmlCcmWrapperModel.hpp"
    #include "SbmlCcmWrapperModel.cpp"

}};
"""
    def generate_code(self, output_directory=None):
        pass
