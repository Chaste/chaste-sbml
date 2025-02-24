function_impl_template = """
    double {ode_name}::{fn_id}({fn_args})
    {{
        return {fn_body_cpp};
    }}
"""

state_param_template = "double {par_id} = this->mParameters[{par_num}]; // {par_name}"
