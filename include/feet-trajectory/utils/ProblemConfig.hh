#pragma once

#include <yaml-cpp/yaml.h>

#include <iostream>
#include <map>
#include <string>
#include <functional>
#include <memory>
#include <exception>
#include <tuple>

#include <Eigen/Core>
#include <feet-trajectory/utils/Box.hh>
#include <feet-trajectory/utils/FixedPlan.hh>

namespace feettrajectory
{

struct ProblemConfig
{
  class CustomString : public std::string
  {
   public:
    Eigen::VectorXd vEigen = Eigen::VectorXd::Zero(0);
    std::vector<Box> vBox = {};
    std::vector<FixedPlan> vFixedPlan = {};
    std::vector<std::string> vString = {};

    CustomString();
    CustomString(const CustomString& cS);
    CustomString(std::string s);
    CustomString(std::string* s);
    CustomString(Eigen::VectorXd& vect, std::string s = "");
    CustomString(std::vector<Box>& vect, std::string s = "");
    CustomString(std::vector<FixedPlan>& vect, std::string s = "");
    CustomString(std::vector<std::string>& v, std::string s = "");

    double asDouble() const;
    int asInt() const;
    size_t asSize_t() const;
    bool asBool() const;
    Eigen::VectorXd asVectorXd() const;
    Eigen::Vector3d asVector3d() const;
    std::vector<Box> asVecBox() const;
    std::vector<FixedPlan> asVecFixedPlan() const;
    std::vector<std::string> asVecString() const;

    operator double() const;
    operator int() const;
    operator size_t() const;
    operator bool() const;
    operator Eigen::VectorXd() const;
    operator Eigen::Vector3d() const;
    operator std::vector<Box>() const;
    operator std::vector<FixedPlan>() const;
    operator std::vector<std::string>() const;
  };

  ProblemConfig(std::string configFile = "");
  ProblemConfig(std::string configFile,
                std::shared_ptr<std::map<std::string, CustomString>> p);

  std::shared_ptr<std::map<std::string, CustomString>> propPtr;
  std::map<std::string, CustomString>& prop;
  void loadFile(char* string);
  void loadFile(std::string configFile);

  /// @brief Checks if key exists
  bool has(std::string key);

  const CustomString get(std::string key);
  const CustomString operator[](std::string key);
  ProblemConfig subSection(std::string section);
  ProblemConfig operator()(std::string section);
  //static ProblemConfig loadUserConfig();

 private:
  std::string globalCategory;
};
} // end of namespace feettrajectory
