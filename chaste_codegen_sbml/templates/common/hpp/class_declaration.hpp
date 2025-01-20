class {{class_name}} : public AbstractOdeSystem
{
private:
  /* Initialise compartments and values. */
  double cell;

  /* Initialise model parameters. */
  {%for state_var in state_vars -%}
  double {{state_var}};
  {%endfor -%};

  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive &archive, const unsigned int version)
  {
    archive &boost::serialization::base_object<AbstractOdeSystem>(*this);
  }

public:
  /* Default constructor. */
  {{class_name}}(std::vector<double> stateVariables = std::vector<double>());

  /* Destructor. */
  ~{{class_name}}();

  void Init();

  void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY);
};
