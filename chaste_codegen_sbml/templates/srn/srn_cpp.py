function_impl_template = """
    double {ode_name}::{fn_id}({fn_args})
    {{
        return {fn_body_cpp};
    }}
"""
