#ifndef {{ test_header_guard }}
#define {{ test_header_guard }}

#include <cxxtest/TestSuite.h>
{% if model_type == ModelType.SRN or model_type == ModelType.CELL_CYCLE %}

#include "AbstractCellBasedTestSuite.hpp"
{% endif %}

#include "{{ ode_hpp_file }}"
{% if model_type == ModelType.SRN %}
#include "{{ srn_hpp_file }}"
{% elif model_type == ModelType.CELL_CYCLE %}
#include "{{ cell_cycle_hpp_file }}"
{% endif %}

// This test is never run in parallel
#include "FakePetscSetup.hpp"

class Test{{ ode_class_name }} : public CxxTest::TestSuite
{
    // TODO: Add tests
public:
    void TestOdeSystem()
    {
        TS_ASSERT_THROWS_NOTHING({{ ode_class_name }} ode_system);
    }
};
{% if model_type == ModelType.SRN %}

class Test{{ srn_class_name }} : public AbstractCellBasedTestSuite
{
    // TODO: Add tests
public:
    void TestSrnModel()
    {
        TS_ASSERT_THROWS_NOTHING({{ srn_class_name }} srn_model);
    }
};
{% elif model_type == ModelType.CELL_CYCLE %}

class Test{{ cell_cycle_class_name }} : public AbstractCellBasedTestSuite
{
    // TODO: Add tests
public:
    void TestCellCycleModel()
    {
        TS_ASSERT_THROWS_NOTHING({{ cell_cycle_class_name }} cell_cycle_model);
    }
};
{% endif %}

#endif // {{ test_header_guard }}
