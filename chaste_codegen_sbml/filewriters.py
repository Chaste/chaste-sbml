from . import translator

############################################################################################################
#####################                 SRN-specific functiions                ###############################
############################################################################################################


def WriteSourceFileForSrnModel(filename, model):
    """Function write source file for SrnModel."""

    # Initialise function
    initialise_fn = GetInitialiseString(filename, model)
    source_file.write(initialise_fn)

    # Define SRN Model
    srn_model_defn = GetModelDefinitionString(filename, model, False)
    source_file.write(srn_model_defn)

    source_file.close()

    print(srn_model_name + ".cpp written!\n")


############################################################################################################
#####################                 CCM-specific functiions                ###############################
############################################################################################################
def GetStoppingEventBooleanString(filename, model):
    """String to define Stopping Event boolean required in Chaste for CCMs to divide."""
    ode_name = GetOdeSystemName(filename)

    stopping_event_str = (
        "bool "
        + ode_name
        + "::CalculateStoppingEvent(double time, const std::vector<double> & rY)\n"  # Define function
        + "{\n"
        + translator.GetBlockCommentDefinition(
            1, "Return true if all events have been triggered.", True
        )
        + translator.AddTabs(1)
        + "return AreAllEventsSatisfied(time, rY);\n"  # Define stopping event
        + "}\n\n"
    )

    return stopping_event_str


def GetCheckAndUpdateEventsVoidString(filename, model):
    """Construct the void that checks if any events have been triggered and
    updates the system accordingly."""
    ode_name = GetOdeSystemName(filename)

    update_events_str = (
        "void "
        + ode_name
        + "::CheckAndUpdateEvents(double time, const std::vector<double> & rY)\n"  # Define function
        + "{\n"
        + translator.AddTabs(1)
        + "std::vector<double> dy(rY.size());\n"  # Initialise derivatives vector
        + translator.AddTabs(1)
        + "EvaluateYDerivatives(time, rY, dy);\n\n"  # Evaluate derivatives
        + translator.GetAllEventsString(model)
        + "}\n\n"
    )

    return update_events_str


def GetAreAllEventsSatisfiedBooleanString(filename):
    """Boolean to check whether or not all events are satisifed."""
    ode_name = GetOdeSystemName(filename)

    are_events_satisfied_str = (
        "bool "
        + ode_name
        + "::AreAllEventsSatisfied(double time, const std::vector<double>& rY)\n"  # Define function
        + "{\n"
        + translator.AddTabs(1)
        + "CheckAndUpdateEvents(time, rY);\n\n"
        + translator.AddTabs(1)
        + "bool events_satisfied = true;\n\n"
        + translator.AddTabs(1)
        + "if ( std::find(eventsSatisfied.begin(), eventsSatisfied.end(), false) != eventsSatisfied.end() )\n"
        + translator.AddTabs(1)
        + "{\n"
        + translator.AddTabs(2)
        + "events_satisfied = false;\n"
        + translator.AddTabs(1)
        + "}\n\n"
        + translator.AddTabs(1)
        + "if (events_satisfied) "
        + translator.GetCommentDefinition(0, "Reset events vector if division is triggered", True)
        + translator.AddTabs(1)
        + "{\n"
        + translator.AddTabs(2)
        + "std::fill(eventsSatisfied.begin(), eventsSatisfied.end(), false);\n"
        + translator.AddTabs(1)
        + "}\n\n"
        + translator.AddTabs(1)
        + "return events_satisfied;\n"
        + "}\n\n"
    )

    return are_events_satisfied_str


def WriteHeaderFileForCcmModel(filename, model):
    """Construct the Chaste header file from the SBML file."""

    ccm_model_name = GetModelName(filename, model)  # Get the name of the file we will write

    # Open to file to write
    header_file = open(ccm_model_name + ".hpp", "w")

    # Define the header files
    header_file_defn = GetHeaderFileDefinitionString(filename, model)
    header_file.write(header_file_defn)

    # Include the appropriate files
    include_files = GetIncludedFilesForHeaderString()
    header_file.write(include_files)

    # Define the ODE System class
    ode_class = GetOdeClassDefinitionString(filename, model)
    header_file.write(ode_class)

    # Define the serialization
    serialization = GetSerializationInformationString(filename)
    header_file.write(serialization)

    # Define the SRN model
    srn_model_defn = GetModelDefinitionString(filename, model, True)
    header_file.write(srn_model_defn)

    # Close the file
    header_close = GetHeaderFileClosingString(filename, model)
    header_file.write(header_close)

    header_file.close()

    print(ccm_model_name + ".hpp written!\n")


def WriteSourceFileForCcmModel(filename, model):
    """Function write source file for SrnModel."""
    ccm_model_name = GetModelName(filename, model)  # Get the name of the file we will write

    # Open to file to write
    source_file = open(ccm_model_name + ".cpp", "w")

    # Include header files
    included_files = GetIncludedFilesForSourceString(filename, model)
    source_file.write(included_files)

    # Initialise class
    class_def = GetClassDefinition(filename, model)
    source_file.write(class_def)

    # Constructor for system
    constructor = GetClassConstructor(filename)
    source_file.write(constructor)

    # Function definitions
    funct_defn_str = GetFunctionDefinitionsForSource(filename, model)
    source_file.write(funct_defn_str)

    # Initialise parameters
    init_fn = GetInitForSource(filename, model)
    source_file.write(init_fn)

    # Get the derivative function
    derivs_fn = GetEvaluateYDerivativesVoidString(filename, model)
    source_file.write(derivs_fn)

    # Get the stopping event function
    stopping_event_fn = GetStoppingEventBooleanString(filename, model)
    source_file.write(stopping_event_fn)

    # Get the void to check and update SBML events
    events_fn = GetCheckAndUpdateEventsVoidString(filename, model)
    source_file.write(events_fn)

    # Get the void to check and update SBML events
    events_satisfied_fn = GetAreAllEventsSatisfiedBooleanString(filename)
    source_file.write(events_satisfied_fn)

    # Initialise function
    initialise_fn = GetInitialiseString(filename, model)
    source_file.write(initialise_fn)

    # Define SRN Model
    srn_model_defn = GetModelDefinitionString(filename, model, False)
    source_file.write(srn_model_defn)

    source_file.close()

    print(ccm_model_name + ".cpp written!\n")


def WriteCcmModelToFile(filename, model):
    """Script that calls upon functions to write the header
    and source files needed for CCM models in Chaste."""
    # Write the .hpp file
    WriteHeaderFileForCcmModel(filename, model)

    # Write the .cpp fil
    WriteSourceFileForCcmModel(filename, model)
