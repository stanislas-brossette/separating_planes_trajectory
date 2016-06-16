#include <cube-stacks/utils/ProblemConfig.hh>

namespace cubestacks
{
ProblemConfig::CustomString::CustomString() {}
ProblemConfig::CustomString::CustomString(const CustomString& cS)
    : std::string(cS), v(cS.v)
{
}
ProblemConfig::CustomString::CustomString(std::string s)
    : std::string(s), v(Eigen::VectorXd::Zero(0))
{
}
ProblemConfig::CustomString::CustomString(std::string* s)
    : std::string(*s), v(Eigen::VectorXd::Zero(0))
{
}

ProblemConfig::CustomString::CustomString(Eigen::VectorXd& vect, std::string s)
    : std::string(s), v(vect)
{
}

double ProblemConfig::CustomString::asDouble() const
{
  try
  {
    return std::stod(this->c_str());
  }
  catch (std::exception& e)
  {
    return 0;
  }
}
int ProblemConfig::CustomString::asInt() const
{
  try
  {
    return std::atoi(this->c_str());
  }
  catch (std::exception& e)
  {
    return 0;
  }
}
bool ProblemConfig::CustomString::asBool() const { return *this == "true"; }
Eigen::VectorXd ProblemConfig::CustomString::asVectorXd() const { return v; }
Eigen::Vector3d ProblemConfig::CustomString::asVector3d() const
{
  assert(v.size() == 3);
  return Eigen::Vector3d(v);
}

ProblemConfig::CustomString::operator double() const { return asDouble(); }
ProblemConfig::CustomString::operator int() const { return asInt(); }
ProblemConfig::CustomString::operator bool() const { return asBool(); }
ProblemConfig::CustomString::operator Eigen::VectorXd() const
{
  return asVectorXd();
}
ProblemConfig::CustomString::operator Eigen::Vector3d() const
{
  return asVector3d();
}

ProblemConfig::ProblemConfig(std::string configFile)
    : ProblemConfig(configFile,
                    std::make_shared<std::map<std::string, CustomString>>())
{
}

ProblemConfig::ProblemConfig(
    std::string configFile,
    std::shared_ptr<std::map<std::string, CustomString>> p)
    : propPtr(p), prop(*propPtr), globalCategory("")
{
  loadFile(configFile);
}

void ProblemConfig::loadFile(char* string) { loadFile(std::string(string)); }

void ProblemConfig::loadFile(std::string configFile)
{
  if (configFile.length() <= 0)
  {
    return;
  }

  YAML::Node node = YAML::LoadFile(configFile);

  std::function<void(const YAML::Node&, std::string)> getAsCustomString =
      [this, &getAsCustomString](const YAML::Node& tree, std::string prefix)
  {
    if (tree.Type() == YAML::NodeType::Map)
    {
      for (YAML::const_iterator it = tree.begin(); it != tree.end(); ++it)
      {
        if (it->second.Type() == YAML::NodeType::Scalar)
        {
          prop[globalCategory + prefix + it->first.as<std::string>()] =
              it->second.as<std::string>();
          // std::cout << prefix + it->first.as<std::string>()
          //<< ": " << std::flush;
          // std::cout << it->second.as<std::string>()
          //<< std::endl;
        }
        else
        {
          getAsCustomString(it->second,
                            prefix + it->first.as<std::string>() + ".");
        }
      }
    }
    if (tree.Type() == YAML::NodeType::Sequence)
    {
      Eigen::VectorXd vect(static_cast<long>(tree.size()));
      std::string acc("");
      int i = 0;
      for (YAML::const_iterator it = tree.begin(); it != tree.end(); ++it)
      {
        acc += it->as<std::string>();
        vect(i++) = it->as<double>();
      }
      prop[globalCategory + prefix.substr(0, prefix.length() - 1)] =
          CustomString(vect, acc);
    }
  };

  getAsCustomString(node, "");
}

const ProblemConfig::CustomString ProblemConfig::get(std::string key)
{
  if (prop.find(globalCategory + key) == prop.end())
  {
    throw std::out_of_range("No value was loaded for key " + globalCategory +
                            key);
  }

  return prop[globalCategory + key];
}

const ProblemConfig::CustomString ProblemConfig::operator[](std::string key)
{
  return get(key);
}

ProblemConfig ProblemConfig::subSection(std::string section)
{
  ProblemConfig pC("", propPtr);

  pC.globalCategory = globalCategory + section + ".";

  return std::move(pC);
}

ProblemConfig ProblemConfig::operator()(std::string section)
{
  return std::move(subSection(section));
}

ProblemConfig ProblemConfig::loadUserConfig(std::string envVar)
{
  ProblemConfig config;
  char* userConfigPath;
  if ((userConfigPath = std::getenv("PROBLEM_GENERATOR_USER_CONFIG")) !=
      nullptr)
  {
    std::cout << "loading config from:" << userConfigPath << std::endl;
    config.loadFile(userConfigPath);
  }
  else
  {
    std::cerr
        << "\033[31mWARNING\033[0m: no user config file found (environment "
           "variable \033[34m" << envVar
        << "\033[0m not set). Robot models and other user-dependent parameters"
           "may not be set correctly." << std::endl;
  }
  return config;
}
}  // end of namespace inertialidentification
