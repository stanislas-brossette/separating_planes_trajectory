#pragma once

#include <yaml-cpp/yaml.h>

#include <iostream>
#include <map>
#include <string>
#include <functional>
#include <memory>
#include <exception>

#include <Eigen/Core>

namespace feettrajectory
{
struct ProblemConfig
{
  class CustomString : public std::string
  {
   public:
    Eigen::VectorXd v;

    CustomString();
    CustomString(const CustomString& cS);
    CustomString(std::string s);
    CustomString(std::string* s);
    CustomString(Eigen::VectorXd& vect, std::string s = "");

    double asDouble() const;
    int asInt() const;
    size_t asSize_t() const;
    bool asBool() const;
    Eigen::VectorXd asVectorXd() const;
    Eigen::Vector3d asVector3d() const;
    Eigen::Vector4d asVector4d() const;

    operator double() const;
    operator int() const;
    operator bool() const;
    operator Eigen::VectorXd() const;
    operator Eigen::Vector3d() const;
  };

  ProblemConfig(std::string configFile = "");
  ProblemConfig(std::string configFile,
                std::shared_ptr<std::map<std::string, CustomString>> p);

  std::shared_ptr<std::map<std::string, CustomString>> propPtr;
  std::map<std::string, CustomString>& prop;
  void loadFile(char* string);
  void loadFile(std::string configFile);
  const CustomString get(std::string key);
  const CustomString operator[](std::string key);
  ProblemConfig subSection(std::string section);
  ProblemConfig operator()(std::string section);
  bool has(std::string key);
  static ProblemConfig loadUserConfig(
      std::string envVar = "PROBLEM_GENERATOR_USER_CONFIG");

 private:
  std::string globalCategory;
};

}  // end of namespace feettrajectory
