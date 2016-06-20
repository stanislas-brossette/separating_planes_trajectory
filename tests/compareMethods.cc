#include <cstdlib>
#include <fstream>
#include <ctime>

#include <pgsolver/utils/finiteDiff.h>
#include <pgsolver/utils/loadFromYaml.h>
#include <pgsolver/solver/Problem.h>
#include <pgsolver/solver/Results.h>
#include <pgsolver/solver/SolverTrustRegionFilter.h>

#include <manifolds/CartesianProduct.h>
#include <manifolds/Point.h>

#include <cube-stacks/CubeStackProblemOnManifold.hh>
#include <cube-stacks/CubeStackProblemOnR.hh>
#include <cube-stacks/utils/ProblemConfig.hh>

using namespace cubestacks;

struct logResult
{
  logResult(size_t nCubes, size_t nPlans, long dim, long repDim, std::string manifold,
             long nCstr, int status, int iteration, double time)
      : nCubes_(nCubes),
        nPlans_(nPlans),
        dim_(dim),
        repDim_(repDim),
        manifold_(manifold),
        nCstr_(nCstr),
        status_(status),
        iteration_(iteration),
        time_(time)
  {
  };
  size_t nCubes_;
  size_t nPlans_;
  long dim_;
  long repDim_;
  std::string manifold_;
  long nCstr_;
  int status_;
  int iteration_;
  double time_;

  std::string print() const
  {
    std::stringstream s;
    s << nCubes_ << ", ";
    s << nPlans_ << ", ";
    s << dim_ << ", ";
    s << repDim_ << ", ";
    //s << manifold_ << ", ";
    s << nCstr_ << ", ";
    s << status_ << ", ";
    s << iteration_ << ", ";
    s << time_;
    return s.str();
  }
};

inline std::ostream& operator<<(std::ostream& os, const logResult& m)
{
  os << m.print();
  return os;
}

std::string printRes(const std::vector<logResult> res)
{
  std::stringstream ss;
  ss << "nCubes_, nPlans_, dim_, repDim_, nCstr_, status_, iteration_, time" << std::endl;
  double averageIter = 0;
  double timePerIter;
  double averageTimePerIter = 0;
  double successRate = 0;
  for (auto s : res)
  {
    ss << s << std::endl;
    if(s.status_==1)
    {
      successRate ++;
      averageIter += s.iteration_;
      timePerIter = s.time_/s.iteration_;
      averageTimePerIter += timePerIter;
    }
  }
  successRate = successRate / res.size();
  averageIter = averageIter/res.size();
  averageTimePerIter = averageTimePerIter/res.size();
  averageTimePerIter = averageTimePerIter/CLOCKS_PER_SEC;
  ss << "successRate: " << 100*successRate << "%" << std::endl;
  ss << "averageIter: " << averageIter << std::endl;
  ss << "averageTimePerIter: " << averageTimePerIter << std::endl;
  return ss.str();
}

inline std::ostream& operator<<(std::ostream& os, const std::vector<logResult>& m)
{
  os << printRes(m);
  return os;
}

logResult solveOnManifold(const mnf::CartesianProduct& M, const Eigen::VectorXd& v0, const std::string& ymlPath)
{
  CubeStackProblemOnManifold myProb(M, ymlPath);
  //pgs::utils::finiteDiffCheck(myProb);
  mnf::Point x0 = myProb.M().createPoint(v0);
  pgs::SolverTrustRegionFilter mySolver;

  pgs::utils::loadOptionsFromYaml(mySolver.opt_, ymlPath);

  mySolver.init(myProb, x0);
  
  auto t0 = std::clock();
  pgs::Results res = mySolver.solve();
  auto tf = std::clock();
  auto solveTime = tf - t0;

  myProb.fileForMatlab(
      "/media/stanislas/Data/virtual_box/shared_data_windows_7/Matlab/tmp/"
      "stackCubesR.m",
      res.x_star);
  logResult out(myProb.nCubes_, myProb.nPlans_, M.dim(), M.representationDim(), M.name(), myProb.totalCstrDim(), res.status, res.iterations, solveTime);
  return out;
}

logResult solveOnR(const mnf::CartesianProduct& R, const Eigen::VectorXd& v0, const std::string& ymlPath)
{
  CubeStackProblemOnR myProb(R, ymlPath);
  //pgs::utils::finiteDiffCheck(myProb);
  mnf::Point x0 = myProb.M().createPoint(v0);
  pgs::SolverTrustRegionFilter mySolver;

  pgs::utils::loadOptionsFromYaml(mySolver.opt_, ymlPath);

  mySolver.init(myProb, x0);
  auto t0 = std::clock();
  pgs::Results res = mySolver.solve();
  auto tf = std::clock();
  auto solveTime = tf - t0;
  myProb.fileForMatlab(
      "/media/stanislas/Data/virtual_box/shared_data_windows_7/Matlab/tmp/"
      "stackCubesM.m",
      res.x_star);
  logResult out(myProb.nCubes_, myProb.nPlans_, R.dim(), R.representationDim(), R.name(), myProb.totalCstrDim(), res.status, res.iterations, solveTime);
  return out;
}

int main(void)
{
  std::ofstream myFile;
  myFile.open("log/CompareMethods.txt");
  std::string testDir = TESTS_DATA_DIR;
  std::string ymlPath = testDir + "/stackCubes.yml";
  ProblemConfig config(ymlPath);

  for (int nCubes = config["minCubes"].asInt(); nCubes < config["maxCubes"].asInt(); nCubes++) 
  {
    myFile << "=============== " << nCubes << " cubes ==============" << std::endl;
    mnf::CartesianProduct* Mptr = CubeStackProblemOnManifold::buildManifold(nCubes);
    mnf::CartesianProduct* Rptr = CubeStackProblemOnR::buildManifold(nCubes);
    myFile << Mptr->name() << std::endl;
    myFile << Rptr->name() << std::endl;

    Mptr->display();
    Rptr->display();

    CubeStackProblemOnManifold tmpProb(*Mptr, ymlPath);

    std::vector<logResult> resManifold, resRealSpace;
    for (int i = 0; i < config["numberOfTests"].asInt(); i++) 
    {
      Eigen::VectorXd v0(Mptr->representationDim());
      v0 = tmpProb.findInitPoint();

      resManifold.push_back(solveOnManifold(*Mptr, v0, ymlPath));
      resRealSpace.push_back(solveOnR(*Rptr, v0, ymlPath));
    }
    
    myFile << "resManifold:\n" << resManifold << std::endl;
    myFile << "resRealSpace:\n" << resRealSpace << std::endl;
  }

  return 0;
}
