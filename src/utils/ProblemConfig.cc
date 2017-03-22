#include <sstream>
#include <stdexcept>
#include <boost/filesystem.hpp>

#include <feet-trajectory/utils/ProblemConfig.hh>

namespace feettrajectory
{
ProblemConfig::CustomString::CustomString() {}
ProblemConfig::CustomString::CustomString(const CustomString& cS)
    : std::string(cS),
      vEigen(cS.vEigen),
      vBox(cS.vBox),
      vFixedPlan(cS.vFixedPlan),
      vString(cS.vString)
{
}
ProblemConfig::CustomString::CustomString(std::string s) : std::string(s) {}
ProblemConfig::CustomString::CustomString(std::string* s) : std::string(*s) {}

ProblemConfig::CustomString::CustomString(Eigen::VectorXd& vect, std::string s)
    : std::string(s), vEigen(vect)
{
}

ProblemConfig::CustomString::CustomString(std::vector<FixedPlan>& vect,
                                          std::string s)
    : std::string(s), vFixedPlan(vect)
{
}
ProblemConfig::CustomString::CustomString(std::vector<Box>& vect, std::string s)
    : std::string(s), vBox(vect)
{
}

ProblemConfig::CustomString::CustomString(std::vector<std::string>& v,
                                          std::string s)
    : std::string(s), vString(v)
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
size_t ProblemConfig::CustomString::asSize_t() const
{
  int res(std::stoi(this->c_str()));
  assert(
      res >= 0 &&
      "Trying to cast a negative number into a size_t, BAD things will happen");
  try
  {
    return static_cast<size_t>(res);
  }
  catch (std::exception& e)
  {
    return 0;
  }
}

bool ProblemConfig::CustomString::asBool() const { return *this == "true"; }

Eigen::VectorXd ProblemConfig::CustomString::asVectorXd() const
{
  return vEigen;
}

Eigen::Vector3d ProblemConfig::CustomString::asVector3d() const
{
  assert(vEigen.size() == 3);
  return Eigen::Vector3d(vEigen);
}

std::vector<Box> ProblemConfig::CustomString::asVecBox() const { return vBox; }
std::vector<FixedPlan> ProblemConfig::CustomString::asVecFixedPlan() const
{
  return vFixedPlan;
}

std::vector<std::string> ProblemConfig::CustomString::asVecString() const
{
  return vString;
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
ProblemConfig::CustomString::operator std::vector<Box>() const
{
  return asVecBox();
}
ProblemConfig::CustomString::operator std::vector<FixedPlan>() const
{
  return asVecFixedPlan();
}
ProblemConfig::CustomString::operator std::vector<std::string>() const
{
  return asVecString();
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
  if (!boost::filesystem::exists(configFile))
  {
    std::stringstream errorMsg;
    errorMsg << "\nError loading YAML File:\n" << configFile
             << "\nDoes not exist!\n" << std::endl;
    throw std::runtime_error(errorMsg.str());
  }
  YAML::Node node = YAML::LoadFile(configFile);

  std::function<void(const YAML::Node&, std::string)> getAsCustomString =
      [this, &getAsCustomString](const YAML::Node& node, std::string prefix)
  {
    if (node.IsMap())
    {
      for (auto subNode : node)
      {
        if (subNode.second.IsScalar())
        {
          prop[globalCategory + prefix + subNode.first.as<std::string>()] =
              subNode.second.as<std::string>();
        }
        else if (subNode.second.IsSequence())
        {
          bool isBoxMap = false;
          bool isFixedPlanMap = false;
          for (auto subSubNode : subNode.second)
          {
            bool hasNormal = false;
            bool hasD= false;
            bool hasCenter = false;
            bool hasSize = false;
            if (subSubNode.IsMap() && subSubNode.size() == 2)
            {
              for (auto i : subSubNode)
              {
                if (std::string("normal").compare(i.first.as<std::string>()) ==
                    0)
                  hasNormal = true;
                else if (std::string("d")
                             .compare(i.first.as<std::string>()) == 0)
                  hasD= true;
                else if (std::string("center")
                             .compare(i.first.as<std::string>()) == 0)
                  hasCenter = true;
                else if (std::string("size")
                             .compare(i.first.as<std::string>()) == 0)
                  hasSize = true;
              }
            }
            if (hasNormal && hasD)
            {
              isFixedPlanMap = true;
            }
            else if (hasCenter && hasSize)
            {
              isBoxMap = true;
            }
            else
            {
              isBoxMap = false;
              isFixedPlanMap = false;
              getAsCustomString(subNode.second,
                                prefix + subNode.first.as<std::string>() + ".");
              break;
            }
          }
          std::string entry = prefix + subNode.first.as<std::string>();
          if (isBoxMap)
          {
            std::vector<Box> vecC;
            int index(0);
            std::string acc("");
            for (auto subSubNode : subNode.second)
            {
              Eigen::Vector3d size, center;
              for (int i = 0; i < 3; i++)
              {
                size[i] = subSubNode["size"][i].as<double>();
                center[i] = subSubNode["center"][i].as<double>();
              }
              vecC.push_back(Box(index, size, center));
              index++;
            }
            prop[entry] = CustomString(vecC, acc);
          }
          else if (isFixedPlanMap)
          {
            std::vector<FixedPlan> vecC;
            std::string acc("");
            for (auto subSubNode : subNode.second)
            {
              double d = subSubNode["d"].as<double>();
              Eigen::Vector3d normal;
              for (int i = 0; i < 3; i++)
              {
                normal[i] = subSubNode["normal"][i].as<double>();
              }
              vecC.push_back(FixedPlan(d, normal));
            }
            prop[entry] = CustomString(vecC, acc);
          }
        }
        else
        {
          getAsCustomString(subNode.second,
                            prefix + subNode.first.as<std::string>() + ".");
        }
      }
    }
    if (node.IsSequence())
    {
      if (node.size() && node[0].IsScalar())
      {
        bool is_double_sequence = true;
        try
        {
          double d = node[0].as<double>();
          (void)(d);
        }
        catch (std::runtime_error&)
        {
          is_double_sequence = false;
        }
        if (is_double_sequence)
        {
          Eigen::VectorXd vect(static_cast<long>(node.size()));
          std::string acc("");
          int i = 0;
          for (auto subNode : node)
          {
            acc += subNode.as<std::string>();
            vect(i++) = subNode.as<double>();
          }
          prop[globalCategory + prefix.substr(0, prefix.length() - 1)] =
              CustomString(vect, acc);
        }
        else
        {
          std::vector<std::string> v(static_cast<size_t>(node.size()));
          std::string acc("");
          size_t i = 0;
          for (const auto& subNode : node)
          {
            acc += ";";
            acc += subNode.as<std::string>();
            v[i++] = subNode.as<std::string>();
          }
          prop[globalCategory + prefix.substr(0, prefix.length() - 1)] =
              CustomString(v, acc);
        }
      }
    }
  };

  getAsCustomString(node, "");
}

const ProblemConfig::CustomString ProblemConfig::get(
    const std::string& key) const
{
  if (prop.find(globalCategory + key) == prop.end())
  {
    throw std::out_of_range("No value was loaded for key " + globalCategory +
                            key);
  }

  return prop.at(globalCategory + key);
}

bool ProblemConfig::has(const std::string& key) const
{
  if (prop.find(globalCategory + key) == prop.end())
    return false;
  else
    return true;
}

const ProblemConfig::CustomString ProblemConfig::operator[](
    const std::string& key) const
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

// ProblemConfig ProblemConfig::loadUserConfig()
//{
// ProblemConfig config;
// char* userConfigPath;
// if ((userConfigPath = std::getenv("PROBLEM_GENERATOR_USER_CONFIG")) !=
// nullptr)
//{
// config.loadFile(userConfigPath);
//}
// else
//{
// std::cerr
//<< "\033[31mWARNING\033[0m: no user config file found (environment "
//"variable \033[34m" << "PROBLEM_GENERATOR_USER_CONFIG"
//<< "\033[0m not set)."
//<< " Will default to the default file: " << PG_DEFAULT_CONFIG
//<< " Robot models and other user-dependent parameters"
//"may not be set correctly." << std::endl;
// config.loadFile(PG_DEFAULT_CONFIG);
//}
// return config;
//}
}  // end of namespace feettrajectory
