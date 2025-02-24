srn_hpp_template = """
#ifndef {header_guard_str}
#define {header_guard_str}

#include ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>"
#include <boost/serialization/shared_ptr.hpp>"
                   
#include <cmath>"
#include <iostream>"
                   
#include AbstractOdeSystem.hpp"

class {ode_name} : public AbstractOdeSystem
{{
private:
    // Model compartments
    {compartment_decls_str}

    // Model parameters
    {parameter_decls_str}

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {{
        archive & boost::serialization::base_object<AbstractOdeSystem>(*this);"
    }}

public:
    // Default constructor.
    {model_name}(std::vector<double> stateVariables=std::vector<double>());

    // Destructor.
    ~{model_name}();

    // Model function definitions
    {function_decls_str}

    void Init();

    void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY);

}};

namespace
{{
namespace serialization
{{
// Serialize information required to construct a {model_name}.
template<class Archive>
inline void save_construct_data(
    Archive & ar, const {model_name} * t, const unsigned int file_version)
{{
    const std::vector<double> state_variables = t->rGetConstStateVariables();
    ar & state_variables;
}}

// De-serialize constructor parameters and intiialise a {model_name}.
template<class Archive>
inline void load_construct_data(
    Archive & ar, {model_name} * t, const unsigned int file_version)
{{
    std::vector<double> state_variables;
    ar & state_variables;
    
    // Invoke inplace constructor to initialise instance
    ::new(t){model_name}(state_variables);
}}
}}
}} // namespace ...

// Define SRN model using wrappers.
#include "SbmlSrnWrapperModel.hpp"
#include "SbmlSrnWrapperModel.cpp"

typedef SbmlSrnWrapperModel<{ode_name}, {size}> {srn_name};

// Declare identifiers for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT({ode_name})
EXPORT_TEMPLATE_CLASS2(SbmlSrnWrapperModel, {ode_name}, {size})

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER({srn_name})

#endif // {header_guard_str}
"""
