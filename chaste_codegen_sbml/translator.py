def GetOdeDefinitionFromAssignment(num_tabs, n, rhs, comp, newline):
    """Get a string defining the ODE in C++."""
    ode_def = AddTabs(num_tabs) + "rDY[" + str(n) + "] = (" + rhs + ") / " + comp + ";"

    if newline:
        ode_def += "\n"

    return ode_def


def GetCppEquivalent(ast_name):
    """Function to get the equivalent of an AST node type. As the list of these operators used is not too large,
    we can essentially hardcode the dictionary of operators and their equivalents."""
    operators = {
        "and": " && ",
        "or": " || ",
        "not": " ! ",
        "eq": " == ",
        "neq": " != ",
        "gt": " > ",
        "geq": " >= ",
        "lt": " < ",
        "leq": " <= ",
    }

    if ast_name in operators:
        return operators[ast_name]
    else:
        return ast_name


def ConvertVariableNamesIntoChasteEquivalent(string_with_variables, model):
    """Function to convert a string with variable names into the Chaste equivalent, replacing
    the variable names with rY[INDEX]."""

    # Get the species ids and their Chaste equivalents
    species_ids_dict = GetSpeciesIndices(model)

    # Create the argument to be passed into the holder function
    argument_list = [string_with_variables]

    for key, value in species_ids_dict.iteritems():
        argument_list.append((key, value))

    string_with_chaste_names = reduce(lambda s, p: s.replace(p[0], p[1]), argument_list)

    return string_with_chaste_names


def GetVariableString(node):
    """Function to get the string equivalent of the node type, whether that is an operator,
    variable name or value."""
    variable = ""

    if node.isNumber():  # If the node is a number, we convert the value into a string
        node_value = node.getValue()
        variable += str(node_value) + " )"
    else:  # We can just grab the name
        node_name = node.getName()
        if node.isLogical():  # If the node represents a logical operator
            node_operator = GetCppEquivalent(node_name)
            variable += node_operator
        else:
            if node.isName():  # If we have a variable name
                variable += "( " + node_name
            else:  # Else we have a relational operator, e.g. <, <= etc
                node_relator = GetCppEquivalent(node_name)
                variable += node_relator

    return variable


## Constant parameter functions
def GetConstantParameterString(model):
    """Get all the constant parameters used in the reactions. These include
    local reaction parameters and compartment volumes (assuming single compartment)
    for now."""
    num_parameters = model.getNumParameters()

    const_param_string = GetBlockCommentDefinition(
        1, "Define constant parameters and their values", True
    )

    for i in range(num_parameters):
        parameter = model.getParameter(i)
        if parameter.getConstant():
            param_id = parameter.getId()
            param_value = parameter.getValue()
            param_def = GetConstDoubleDefinition(1, param_id, param_value, True)
            const_param_string += param_def

    return const_param_string


## Non-constant parameter functions
def GetNonConstantParameterString(model):
    """Get all the non-constant parameters used in the reactions. These include
    local parameters used to flag division events (e.g Tyson and Novak (2001))."""
    num_parameters = model.getNumParameters()

    non_const_param_string = GetBlockCommentDefinition(
        1, "Define non-constant parameters and their values", True
    )

    for i in range(num_parameters):
        parameter = model.getParameter(i)
        if not parameter.getConstant():
            param_id = parameter.getId()
            param_value = parameter.getValue()
            param_def = GetConstDoubleDefinition(1, param_id, param_value, True)
            non_const_param_string += param_def

    return non_const_param_string


###################################################################################################
######################        SRN AND Cell Cycle Model Functions         ##########################
###################################################################################################
# Parameter functions
def GetParameter(n, model):
    """Get the ID of the parameter"""
    parameter = model.getParameter(n)
    parameter_id = parameter.getId()
    return parameter_id


def GetParameterValue(n, model):
    """Get the value of the parameter"""
    parameter = model.getParameter(n)
    parameter_id = parameter.getValue()
    return parameter_id


def GetParameterName(n, model):
    """Get the ID of the parameter"""
    parameter = model.getParameter(n)
    parameter_name = parameter.getName()
    return parameter_name


def IsParameterConstant(n, model):
    """Boolean to check if parameter is constant"""
    parameter = model.getParameter(n)

    if parameter.getConstant():
        return True
    else:
        return False


#############################################################################################

## State variables functions


def GetSpeciesIndices(model):
    """Function to return a dictionary of species IDs and their
    Chaste variable name equivalent."""
    odes_dict = GetOdesDictionary(model)

    species_dict = {}  # Initialise
    num_species = model.getNumSpecies()

    species_count = 0
    for i in range(num_species):
        species = model.getSpecies(i)
        species_id = species.getId()
        if species_id in odes_dict:
            species_dict[" " + species_id + " "] = " rY[" + str(species_count) + "] "
            species_count += 1

    return species_dict


## Compartment functions


def GetCompartmentVariableString(model):
    """Extract compartments and their sizes to produce string that is written to
    C++ code."""
    num_compartments = model.getNumCompartments()

    compartments_string = GetBlockCommentDefinition(
        1, "Define compartments and their sizes", True
    )

    for i in range(num_compartments):
        compartment_id = GetCompartment(i, model)
        ################# EDITED ##################### use name rather than id
        compartment_name = GetCompartmentName(i, model)
        compartment_val = GetCompartmentValue(i, model)
        compartment_def = GetConstDoubleDefinition(
            1, compartment_name, compartment_val, True
        )
        compartments_string += compartment_def

    return compartments_string


## Functions to return information needed to initialise ODEs
def GetInitialInformationForASpecies(species, model):
    """Get the string describing the initial information for a species; that is,
    its name, units and initial concentrations."""

    # Define the string for the name
    species_name = species.getName()
    species_id = species.getId()

    initial_information_string = ""

    if IsSpeciesDefinedAsOde(species_id, model):
        if species_name == "":  # If no name, set the name to be the species ID
            name_string = (
                AddTabs(1) + 'this->mVariableNames.push_back("' + species_id + '");\n'
            )
        else:
            name_string = (
                AddTabs(1) + 'this->mVariableNames.push_back("' + species_name + '");\n'
            )

        # Define the string for the units
        species_compartment = species.getCompartment()
        species_units = species.getSubstanceUnits()

        if (
            species_compartment
        ):  # if there's a compartment, then we would have normalised the ODE, so declare it as non-dimensional
            units_string = AddTabs(1) + 'this->mVariableUnits.push_back("non-dim");\n'
        else:
            units_string = (
                AddTabs(1)
                + 'this->mVariableUnits.push_back("'
                + species_units
                + '");\n'
            )

        # Define the string for the initial condition
        if species.isSetInitialAmount():
            species_init_conc = species.getInitialAmount()
        else:
            species_init_conc = species.getInitialConcentration()
        init_conc_string = (
            AddTabs(1)
            + "this->mInitialConditions.push_back("
            + str(species_init_conc)
            + ");\n"
        )

        initial_information_string += (
            name_string + units_string + init_conc_string + "\n"
        )
    elif IsSpeciesDefinedAsParameter(species_id, model):
        if species_name == "":  # If no name, set the name to be the species ID
            name_string = (
                AddTabs(1) + 'this->mParameterNames.push_back("' + species_id + '");\n'
            )
        else:
            name_string = (
                AddTabs(1)
                + 'this->mParameterNames.push_back("'
                + species_name
                + '");\n'
            )

        # Define the string for the units
        species_compartment = species.getCompartment()
        species_units = species.getSubstanceUnits()

        if (
            species_compartment
        ):  # if there's a compartment, then we would have normalised the ODE, so declare it as non-dimensional
            units_string = AddTabs(1) + 'this->mParameterUnits.push_back("non-dim");\n'
        else:
            units_string = (
                AddTabs(1)
                + 'this->mParameterUnits.push_back("'
                + species_units
                + '");\n'
            )

    elif IsVariableDefinedAsRule(
        species_id, model
    ):  ############################ EDITED ##################################### if species in rule from simbio
        if species_name == "":
            name_string = (
                AddTabs(1) + 'this->mVariableNames.push_back("' + species_id + '");\n'
            )

        else:
            name_string = (
                AddTabs(1) + 'this->mVariableNames.push_back("' + species_name + '");\n'
            )

        # Define the string for the units
        species_compartment = species.getCompartment()
        species_units = species.getSubstanceUnits()

        if (
            species_compartment
        ):  # if there's a compartment, then we would have normalised the ODE, so declare it as non-dimensional
            units_string = AddTabs(1) + 'this->mVariableUnits.push_back("non-dim");\n'
        else:
            units_string = (
                AddTabs(1)
                + 'this->mVariableUnits.push_back("'
                + species_units
                + '");\n'
            )

        # Define the string for the initial condition
        if species.isSetInitialAmount():
            species_init_conc = species.getInitialAmount()
        else:
            species_init_conc = species.getInitialConcentration()
        init_conc_string = (
            AddTabs(1)
            + "this->mInitialConditions.push_back("
            + str(species_init_conc)
            + ");\n"
        )

        # Define string (should be three lines)
        initial_information_string += (
            name_string + units_string + init_conc_string + "\n"
        )

    return initial_information_string


def GetInitialInformationForAParameter(parameter, model):
    """String describing initial information for a parameter; that is, its name and units. This function should only be called if
    the parameter does not have a set value (and thus must be externally defined)."""
    parameter_name = parameter.getName()
    parameter_id = parameter.getId()

    initial_information_string = ""

    if parameter_name == "":
        name_string = (
            AddTabs(1) + 'this->mParameterNames.push_back("' + parameter_id + '");\n'
        )
    else:
        name_string = (
            AddTabs(1) + 'this->mParameterNames.push_back("' + parameter_name + '");\n'
        )

    # Define the string for the units
    if parameter.isSetUnits():
        parameter_units = species.getUnits()
        units_string = (
            AddTabs(1) + 'this->mParameterUnits.push_back("' + parameter_units + '");\n'
        )
    else:
        units_string = AddTabs(1) + 'this->mParameterUnits.push_back("non-dim");\n'

    initial_information_string += name_string + units_string + "\n"

    return initial_information_string


def GetInitialInformationString(model):
    """Get all the information required to intiialise the system as a string."""
    initial_information_string = ""  # Initialise string

    num_species = model.getNumSpecies()

    for i in range(num_species):
        species = model.getSpecies(i)
        species_string = GetInitialInformationForASpecies(species, model)
        initial_information_string += species_string

    num_parameters = model.getNumParameters()

    for i in range(num_parameters):
        parameter = model.getParameter(i)
        parameter_id = parameter.getId()
        parameter_name = parameter.getName()

        if (not IsParameterValueSet(i, model)) & (
            not IsVariableDefinedAsRule(parameter_id, model)
        ):
            parameter_string = GetInitialInformationForAParameter(parameter, model)
            initial_information_string += parameter_string
        ############### EDITED ###################
        elif parameter_name:
            if (
                ("wnt" in parameter_name)
                or ("gamma" in parameter_name)
                or ("ComplexTransit") in parameter_name
            ):
                parameter_string = GetInitialInformationForAParameter(parameter, model)
                initial_information_string += parameter_string

    return initial_information_string


###################################################################################################
#####################            Cell Cycle Model Functions         ###############################
###################################################################################################


def GetEventTriggerString(event):
    """Function to get the condition required to trigger an Event. Output is a boolean string
    that can go straight into C++."""
    trigger_string = ""  # Initialise string

    # Get the AST describing the trigger
    event_trigger = event.getTrigger()
    trigger_ast = event_trigger.getMath()

    preorder_list = (
        trigger_ast.getListOfNodes()
    )  # Get the list of nodes in preorder (DFS)

    # Get the list of nodes in-order, so that we can readily generate the condition
    inorder_list = []
    inorder_list = GetInOrderTraversal(inorder_list, trigger_ast, preorder_list)

    # Generate the string using in-order list
    for i in range(len(inorder_list)):
        ast_node = preorder_list.get(inorder_list[i])
        variable_string = GetVariableString(ast_node)
        trigger_string += variable_string

    # We need to replace the species variable names with their equivalents in Chaste.
    original_model = event.getModel()

    trigger_string = ConvertVariableNamesIntoChasteEquivalent(
        trigger_string, original_model
    )
    return trigger_string


def GetEventAssignmentString(event_assignment):
    """Get string equivalent for an event assignment."""
    event_assignment_string = ""
    original_model = event_assignment.getModel()

    event_assignment_variable = (
        " " + event_assignment.getVariable() + " "
    )  # Get the variable name
    event_assignment_variable = ConvertVariableNamesIntoChasteEquivalent(
        event_assignment_variable, original_model
    )  # Convert into Chaste code
    event_assignment_variable = event_assignment_variable.replace(
        "rY", "this->rGetStateVariables()"
    )  # Replace rY with modifiable state
    event_assignment_variable = event_assignment_variable.replace(
        " ", ""
    )  # Remove whitespace

    event_assignment_math = event_assignment.getMath()  # Get the AST for the maths
    assignment_formula = " " + formulaToString(event_assignment_math)
    assignment_formula = ConvertVariableNamesIntoChasteEquivalent(
        assignment_formula, original_model
    )  # Convert into Chaste code

    event_assignment_string = (
        AddTabs(2)
        + event_assignment_variable
        + " = double ( "
        + assignment_formula.strip()
        + " );\n"
    )

    return event_assignment_string


def GetEventAssignments(event):
    """Get all the event assignments that occur during an event."""
    assignments_string = ""

    for i in range(event.getNumEventAssignments()):
        event_assignment = event.getEventAssignment(i)
        event_assignment_string = GetEventAssignmentString(event_assignment)
        assignments_string += event_assignment_string

    return assignments_string


def GetTriggerAndAssignmentString(n, event):
    """Get the trigger and assignment for an event."""

    # Get the condition for the trigger
    trigger_string = GetEventTriggerString(event)

    # Get the assignment upon the trigger
    event_assignment = GetEventAssignments(event)

    condition_string = (
        AddTabs(1)
        + "if ( "
        + trigger_string
        + " )\n"
        + AddTabs(1)
        + "{\n"
        + event_assignment
        + AddTabs(2)
        + "eventsSatisfied["
        + str(n)
        + "] = true;\n"
        + AddTabs(1)
        + "}\n"
    )

    return condition_string


def GetAllEventsString(model):
    """Generate string to describe all the events, their triggers and assignments."""
    num_events = model.getNumEvents()
    all_events = ""

    for i in range(num_events):
        event = model.getEvent(i)
        trigger_and_assignment = GetTriggerAndAssignmentString(i, event)
        all_events += trigger_and_assignment + "\n"

    return all_events  #!
