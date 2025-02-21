import os
import shutil
import pathlib

from libsbml import SBMLReader

from . import filewriters
from .chaste_model import ChasteModel

HEADER_GUARD_SUFFIX = "ODESYSTEMANDSRNMODEL_HPP_"

ODE_SUFFIX = "OdeSystem"
SRN_SUFFIX = "SrnModel"

TAB = "    "


class ChasteSRNModel(ChasteModel):

    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)

        self.ode_name = self.model_name + ODE_SUFFIX
        self.srn_name = self.model_name + SRN_SUFFIX

        # Size is defined as the number of ODEs
        odes = self.get_odes()
        rules = self.get_rules()
        self.num_species = len(odes) + len(rules)

    def generate_hpp(self) -> str:
        """Generate the Chaste header for an SRN model from SBML data.

        return: The generated header file as a string.
        """

        # Get inputs for the header file template
        header_guard = self.model_name.upper() + HEADER_GUARD_SUFFIX

        compartments = self.get_compartments()
        compartment_decls = f"\n{TAB}".join([f"double {var};" for var in compartments])

        parameters = self.get_parameters()
        parameter_decls = f"\n{TAB}".join([f"double {var};" for var in parameters])

        functions = self.get_function_definitions()
        function_decls = f"\n{TAB}".join(
            [f"double {fn} ({', '.join(args)});" for fn, args in functions.items()]
        )

        hpp = f"""
#ifndef {header_guard}
#define {header_guard}

#include ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>"
#include <boost/serialization/shared_ptr.hpp>"
                   
#include <cmath>"
#include <iostream>"
                   
#include AbstractOdeSystem.hpp"

class {self.ode_name} : public AbstractOdeSystem
{{
private:
    // Model compartments
    {compartment_decls}

    // Model parameters
    {parameter_decls}

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {{
        archive & boost::serialization::base_object<AbstractOdeSystem>(*this);"
    }}

public:
    // Default constructor.
    {self.model_name}(std::vector<double> stateVariables=std::vector<double>());

    // Destructor.
    ~{self.model_name}();

    // Model function definitions
    {function_decls}

    void Init();

    void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY);

}};

namespace
{{
namespace serialization
{{
// Serialize information required to construct a {self.model_name}.
template<class Archive>
inline void save_construct_data(
    Archive & ar, const {self.model_name} * t, const unsigned int file_version)
{{
    const std::vector<double> state_variables = t->rGetConstStateVariables();
    ar & state_variables;
}}

// De-serialize constructor parameters and intiialise a {self.model_name}.
template<class Archive>
inline void load_construct_data(
    Archive & ar, {self.model_name} * t, const unsigned int file_version)
{{
    std::vector<double> state_variables;
    ar & state_variables;
    
    // Invoke inplace constructor to initialise instance
    ::new(t){self.model_name}(state_variables);
}}
}}
}} // namespace ...

// Define SRN model using Wrappers.
#include "SbmlSrnWrapperModel.hpp"
#include "SbmlSrnWrapperModel.cpp"

typedef SbmlSrnWrapperModel<{self.ode_name}, {self.num_species}> {self.srn_name};

// Declare identifiers for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT({self.ode_name})
EXPORT_TEMPLATE_CLASS2(SbmlSrnWrapperModel, {self.ode_name}, {self.num_species})

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER({self.srn_name})

#endif // {header_guard}
"""
        return hpp

    #########################################
    def generate_cpp(self) -> str:
        """Generate the Chaste source for an SRN model from SBML data.

        return: The generated source file as a string.
        """

        species = self.get_species()
        species_init = f"\n{TAB}".join(
            [f"SetDefaultInitialCondition({i}, {val});" for i, val in species.items()]
        )

        params = self.get_parameters()
        params_init = f"\n{TAB}".join(
            [f"this->mParameters.push_back({val});" for val in params.values()]
        )

        cpp = f"""
#include "{self.srn_name}.hpp"
#include "CellwiseOdeSystemInformation.hpp"

// SBML ODE System
{self.ode_name}::{self.ode_name} (std::vector<double> stateVariables)
    : AbstractOdeSystem({self.num_species})
{{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<{self.ode_name}>);

    Init();

    {species_init}
    {params_init}

    if (stateVariables != std::vector<double>())
    {{
        SetStateVariables(stateVariables);
    }}
}}

{self.ode_name}::~{self.ode_name}()
{{
}}

{function_definitions}

void {self.ode_name}::Init()
 {{
    /* Initialise the parameters. */
    // cell = 1.0;
    // V1 = 0.0;
    // V3 = 0.0;
    // VM1 = 3.0;
    // VM3 = 1.0;
    // Kc = 0.5;
}}

void {self.ode_name}::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{{
    
    /* Define state variables */
    double C = rY[0]; // Cyclin
    double M = rY[1]; // cdc_2_kinase
    double X = rY[2]; // Cyclin_Protease

    /* Define state parameters */

     /* Define algebraic rules. */
    V1 = C * VM1 * pow(C + Kc, -1);
    V3 = M * VM3;

    /* Define the reactions in this model. */
    // creation of cyclin
    double vi = 0.025;
    double reaction1 = cell * vi;

    // default degradation of cyclin
    double kd = 0.01;
    double reaction2 = C * cell * kd;

    // cdc2 kinase triggered degration of cyclin
    double vd = 0.25;
    double Kd = 0.02;
    double reaction3 = C * cell * vd * X * pow(C + Kd, -1);

    // activation of cdc2 kinase
    double K1 = 0.005;
    double reaction4 = cell * (1 + -1 * M) * V1 * pow(K1 + -1 * M + 1, -1);

    // deactivation of cdc2 kinase
    double V2 = 1.5;
    double K2 = 0.005;
    double reaction5 = cell * M * V2 * pow(K2 + M, -1);

    // activation of cyclin protease
    double K3 = 0.005;
    double reaction6 = cell * V3 * (1 + -1 * X) * pow(K3 + -1 * X + 1, -1);

    // deactivation of cyclin protease
    double K4 = 0.005;
    double V4 = 0.5;
    double reaction7 = cell * V4 * X * pow(K4 + X, -1);


    rDY[0] = (reaction1 - reaction2 - reaction3) / cell; // dCyclin/dt
    rDY[1] = (reaction4 - reaction5) / cell; // dcdc_2_kinase/dt
    rDY[2] = (reaction6 - reaction7) / cell; // dCyclin Protease/dt

    /* Account for the differences in timescales. */
    // rDY[0] *= 3600.0;
    // rDY[1] *= 3600.0;
    // rDY[2] *= 3600.0;

    // std::cout << rDY[0] << ", " << rDY[1] << ", " << rDY[2] << ", " << std::endl;
}}

template<>
void CellwiseOdeSystemInformation<{self.ode_name}>::Initialise()
{{
    this->mVariableNames.push_back("Cyclin");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.01);

    this->mVariableNames.push_back("cdc_2_kinase");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.01);

    this->mVariableNames.push_back("Cyclin Protease");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.01);


    this->mInitialised = true;
}}

/* Define SRN model using Wrappers. */
#include "SbmlSrnWrapperModel.hpp"
#include "SbmlSrnWrapperModel.cpp"

typedef SbmlSrnWrapperModel<{self.ode_name}, {self.num_species}> {self.srn_name};

// Declare identifiers for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT({self.ode_name})
EXPORT_TEMPLATE_CLASS2(SbmlSrnWrapperModel, {self.ode_name}, {self.num_species})

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER({self.srn_name})
"""

        # Functiond efinitions
        funct_defn_str = GetFunctionDefinitionsForSource(filename, model)
        source_file.write(funct_defn_str)

        # Initialise parameters
        init_fn = GetInitForSource(filename, model)
        source_file.write(init_fn)

        # Get the derivative function
        derivs_fn = GetEvaluateYDerivativesVoidString(filename, model)
        source_file.write(derivs_fn)

        # Initialise function
        initialise_fn = GetInitialiseString(filename, model)
        source_file.write(initialise_fn)

        # Define SRN Model
        srn_model_defn = GetModelDefinitionString(filename, model, False)
        source_file.write(srn_model_defn)

    def generate_code(self, output_directory=None):
        filewriters.WriteSrnModelToFile(self.model_name, self.model)
