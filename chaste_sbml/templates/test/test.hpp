#ifndef {{ test_header_guard }}
#define {{ test_header_guard }}

#include <cxxtest/TestSuite.h>
{% if model_type != ModelType.GENERIC %}

#include "AbstractCellBasedTestSuite.hpp"
{% endif %}

#include "{{ ode_hpp_file }}"
{% if model_type != ModelType.GENERIC %}
#include "{{ wrapper_hpp_file }}"
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
{% if model_type != ModelType.GENERIC %}

class Test{{ wrapper_class_name }} : public AbstractCellBasedTestSuite
{
    // TODO: Add tests
public:
    void TestModel()
    {
        TS_ASSERT_THROWS_NOTHING({{ wrapper_class_name }} model);
    }
};
{% endif %}

#endif // {{ test_header_guard }}
