#include <pgsolver/utils/loadFromYaml.h>
#include <pgsolver/utils/defs.h>
#include <pgsolver/solver/Problem.h>
#include <pgsolver/solver/Results.h>
#include <pgsolver/solver/SolverTrustRegionFilter.h>

#include <manifolds/CartesianProduct.h>
#include <manifolds/Point.h>

#include <cube-stacks/CubeStackProblemOnManifold.hh>
#include <cube-stacks/CubeStackProblemOnSO3noS2.hh>
#include <cube-stacks/CubeStackProblemOnR.hh>
#include <cube-stacks/utils/ProblemConfig.hh>

using namespace cubestacks;

struct logResult
{
  logResult(const size_t& nCubes, const size_t& nPlans, const long& dim,
            const long& repDim, const std::string& manifold, const long& nCstr,
            const int& status, const int& iteration, const double& time,
            const double& obj_star, const Eigen::VectorXd& v0, std::shared_ptr<pgs::TimeLogger> tLog)
      : nCubes_(nCubes),
        nPlans_(nPlans),
        dim_(dim),
        repDim_(repDim),
        manifold_(manifold),
        nCstr_(nCstr),
        status_(status),
        iteration_(iteration),
        time_(time),
        obj_star_(obj_star),
        v0_(v0),
        tLog_(tLog)
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
  double obj_star_;
  Eigen::VectorXd v0_;
  std::shared_ptr<pgs::TimeLogger> tLog_;
  Eigen::IOFormat fmt_;

  std::ostream& print(std::ostream& o) const
  {
    o << nCubes_ << ", ";
    o << nPlans_ << ", ";
    o << dim_ << ", ";
    o << repDim_ << ", ";
    o << nCstr_ << ", ";
    o << status_ << ", ";
    o << iteration_ << ", ";
    o << time_ << ", ";
    o << obj_star_ << "\n";
    o << v0_.transpose().format(fmt::matlabVector);
    return o;
  }

};

inline std::ostream& operator<<(std::ostream& os, const logResult& m)
{
  return m.print(os);
}

std::string printSummary(const std::vector<logResult> res)
{
  std::stringstream ss;
  ss << "nCubes, nPlans, dim, repDim, nCstr, status, iteration, time, obj_star, v0" << std::endl;
  double averageIter = 0;
  double timePerIter;
  double averageTimePerIter = 0;
  double successRate = 0;
  for (auto s : res)
  {
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

std::string printRes(const std::vector<logResult>& res)
{
  std::stringstream ss;
  ss << "nCubes, nPlans, dim, repDim, nCstr, status, iteration, time, obj_star, v0" << std::endl;
  double averageIter = 0;
  double timePerIter;
  double averageTimePerIter = 0;
  double successRate = 0;
  for (size_t i = 0; i < res.size(); i++)
  {
    ss << res[i] << "\n";
    if(res[i].status_==1)
    {
      successRate ++;
      averageIter += res[i].iteration_;
      timePerIter = res[i].time_/res[i].iteration_;
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

std::string avTimePerIter(const std::vector<pgs::VecTimer>& timer)
{
  std::stringstream ss;
  long maxIter = timer[0].vec().size();
  Eigen::VectorXd perIter(maxIter);
  Eigen::VectorXd count(maxIter);
  perIter.setZero();
  count.setZero();
  for (long i = 0; i < maxIter; i++)
  {
    for (auto s : timer)
    {
      if(s.vec()(i) != 0)
      {
        perIter(i) += s.vec()(i);
        count(i) ++;
      }
    }
  }
  for (long i = 0; i < maxIter; i++)
  {
    perIter(i) = perIter(i)/count(i);
  }
  ss << perIter.transpose().format(fmt::matlabVector);
  return ss.str();
}

std::string printAverageTimesPerIter(const std::vector<logResult>& res, const std::string& prefix)
{
  std::vector<pgs::VecTimer> iter, hessian, qp, regularisation,
      regularisationRestoration, problemUpdate;
  for (auto i : res)
  {
    iter.push_back(i.tLog_->iter);
    hessian.push_back(i.tLog_->hessian);
    qp.push_back(i.tLog_->qp);
    regularisation.push_back(i.tLog_->regularisation);
    regularisationRestoration.push_back(i.tLog_->regularisationRestoration);
    problemUpdate.push_back(i.tLog_->problemUpdate);
  }

  std::stringstream ss;
  ss << prefix << "iterTimePerIter = " << avTimePerIter(iter) << "\n";
  ss << prefix << "hessianTimePerIter = " << avTimePerIter(hessian) << "\n";
  ss << prefix << "qpTimePerIter = " << avTimePerIter(qp) << "\n";
  ss << prefix << "regularisationTimePerIter = " << avTimePerIter(regularisation) << "\n";
  ss << prefix << "regularisationRestorationTimePerIter = " << avTimePerIter(regularisationRestoration) << "\n";
  ss << prefix << "problemUpdateTimePerIter = " << avTimePerIter(problemUpdate) << "\n";
  return ss.str();
}

std::string printTimes(const std::vector<logResult> res, const std::string& prefix)
{
  std::stringstream ss;
  ss << prefix << "initTime = [";
  for (auto s : res)
    ss << s.tLog_->initializationTime() << ",";
  ss << "];\n";

  ss << prefix << "totalTime = [";
  for (auto s : res)
    ss << s.tLog_->totalTime() << ",";
  ss << "];\n";

  ss << prefix << "iter = [";
  for (auto s : res)
    ss << s.tLog_->iter << "\n";
  ss << "];\n";

  ss << prefix << "hessian = [";
  for (auto s : res)
    ss << s.tLog_->hessian << "\n";
  ss << "];\n";

  ss << prefix << "qp = [";
  for (auto s : res)
    ss << s.tLog_->qp << "\n";
  ss << "];\n";

  ss << prefix << "regularisation = [";
  for (auto s : res)
    ss << s.tLog_->regularisation << "\n";
  ss << "];\n";

  ss << prefix << "regularisationRestoration = [";
  for (auto s : res)
    ss << s.tLog_->regularisationRestoration
       << "\n";
  ss << "];\n";

  ss << prefix << "problemUpdate = [";
  for (auto s : res)
    ss << s.tLog_->problemUpdate << "\n";
  ss << "];\n";
  return ss.str();
}


logResult solveOnManifold(const mnf::CartesianProduct& M, const Eigen::VectorXd& v0, const std::string& ymlPath, const std::string& pbName)
{
  CubeStackProblemOnManifold myProb(M, ymlPath);
  myProb.name(pbName);
  //pgs::utils::finiteDiffCheck(myProb);
  mnf::Point x0 = myProb.M().createPoint(v0);
  pgs::SolverTrustRegionFilter mySolver;

  pgs::utils::loadOptionsFromYaml(mySolver.opt_, ymlPath);

  mySolver.init(myProb, x0);

  auto t0 = std::clock();
  pgs::Results res = mySolver.solve();
  auto tf = std::clock();
  auto solveTime = tf - t0;

  //myProb.fileForMatlab(
      //"/media/stanislas/Data/virtual_box/shared_data_windows_7/Matlab/tmp/"
      //"stackCubesR.m",
      //res.x_star);
  logResult out(myProb.nCubes_, myProb.nPlans_, M.dim(), M.representationDim(), M.name(), myProb.totalCstrDim(), res.status, res.iterations, solveTime, res.obj_star, v0, mySolver.getTimeLogger());
  return out;
}

logResult solveOnSO3noS2(const mnf::CartesianProduct& M, const Eigen::VectorXd& v0, const std::string& ymlPath, const std::string& pbName)
{
  CubeStackProblemOnSO3noS2 myProb(M, ymlPath);
  myProb.name(pbName);
  //pgs::utils::finiteDiffCheck(myProb);
  mnf::Point x0 = myProb.M().createPoint(v0);
  pgs::SolverTrustRegionFilter mySolver;

  pgs::utils::loadOptionsFromYaml(mySolver.opt_, ymlPath);

  mySolver.init(myProb, x0);
  
  auto t0 = std::clock();
  pgs::Results res = mySolver.solve();
  auto tf = std::clock();
  auto solveTime = tf - t0;

  //myProb.fileForMatlab(
      //"/media/stanislas/Data/virtual_box/shared_data_windows_7/Matlab/tmp/"
      //"stackCubesSO3noS2.m",
      //res.x_star);
  logResult out(myProb.nCubes_, myProb.nPlans_, M.dim(), M.representationDim(), M.name(), myProb.totalCstrDim(), res.status, res.iterations, solveTime, res.obj_star, v0, mySolver.getTimeLogger());
  return out;
}

logResult solveOnR(const mnf::CartesianProduct& R, const Eigen::VectorXd& v0, const std::string& ymlPath, const std::string& pbName)
{
  CubeStackProblemOnR myProb(R, ymlPath);
  myProb.name(pbName);
  //pgs::utils::finiteDiffCheck(myProb);
  mnf::Point x0 = myProb.M().createPoint(v0);
  pgs::SolverTrustRegionFilter mySolver;

  pgs::utils::loadOptionsFromYaml(mySolver.opt_, ymlPath);

  mySolver.init(myProb, x0);
  auto t0 = std::clock();
  pgs::Results res = mySolver.solve();
  auto tf = std::clock();
  auto solveTime = tf - t0;
  //myProb.fileForMatlab(
      //"/media/stanislas/Data/virtual_box/shared_data_windows_7/Matlab/tmp/"
      //"stackCubesM.m",
      //res.x_star);
  logResult out(myProb.nCubes_, myProb.nPlans_, R.dim(), R.representationDim(), R.name(), myProb.totalCstrDim(), res.status, res.iterations, solveTime, res.obj_star, v0, mySolver.getTimeLogger());
  return out;
}

std::string currentTime()
{
  time_t rawtime;
  struct tm * timeinfo;
  char buffer[80];

  time (&rawtime);
  timeinfo = localtime(&rawtime);

  strftime(buffer,80,"%Y-%m-%d.%H:%M:%S",timeinfo);
  std::string str(buffer);
  return str;
}

int main(void)
{
  std::cout << "CLOCKS_PER_SEC: " << (double) CLOCKS_PER_SEC << std::endl;
  std::ofstream myFile, mySumFile, myTimeFile, myTimePerIterFile;
  std::string testDir = TESTS_DATA_DIR;
  std::string ymlPath = testDir + "/stackCubes.yml";
  ProblemConfig config(ymlPath);

  if(config["timeLog"])
  {
    std::string baseName("log/" + config["logName"]/* + "_" + currentTime()*/);
    std::string fileName(baseName + ".m");
    std::string summaryFileName(baseName + "_summary.m");
    std::string timeFileName(baseName + "_time.m");
    std::string timePerIterFileName(baseName + "_timePerIter.m");
    myFile.open(fileName);
    mySumFile.open(summaryFileName);
    myTimeFile.open(timeFileName);
    myTimePerIterFile.open(timePerIterFileName);
  }

  for (int nCubes = config["minCubes"].asInt(); nCubes <= config["maxCubes"].asInt(); nCubes++)
  {
    std::cerr << "\n\n\n((((------ "<< nCubes << " cubes ------))))" << std::endl;

    mnf::CartesianProduct* Mptr = CubeStackProblemOnManifold::buildManifold(nCubes);
    mnf::CartesianProduct* MSO3noS2ptr = CubeStackProblemOnSO3noS2::buildManifold(nCubes);
    mnf::CartesianProduct* Rptr = CubeStackProblemOnR::buildManifold(nCubes);

    Mptr->display();
    MSO3noS2ptr->display();
    Rptr->display();

    CubeStackProblemOnManifold tmpProb(*Mptr, ymlPath);

    std::vector<logResult> resManifold, resSO3noS2, resRealSpace;
    for (int i = 0; i < config["numberOfTests"].asInt(); i++) 
    {
      std::cerr << "i: " << i << std::endl;
      Eigen::VectorXd v0(Mptr->representationDim());
      v0 = tmpProb.findInitPoint();

      std::string manifoldPbName("BFGS_manifoldPb" + std::to_string(nCubes) + "Cubes" + std::to_string(i));
      std::string SO3noS2PbName("BFGS_SO3noS2Pb" + std::to_string(nCubes) + "Cubes" + std::to_string(i));
      std::string RealSpacePbName("BFGS_RealSpacePb" + std::to_string(nCubes) + "Cubes" + std::to_string(i));
      resManifold.push_back(solveOnManifold(*Mptr, v0, ymlPath, manifoldPbName));
      resSO3noS2.push_back(solveOnSO3noS2(*MSO3noS2ptr, v0, ymlPath, SO3noS2PbName));
      resRealSpace.push_back(solveOnR(*Rptr, v0, ymlPath, RealSpacePbName));
    }
    
    if(config["timeLog"])
    {
      myFile << "%=============== " << nCubes << " cubes ==============" << std::endl;
      myFile << "%" << Mptr->name() << std::endl;
      myFile << "%" << MSO3noS2ptr->name() << std::endl;
      myFile << "%" << Rptr->name() << std::endl;
      myFile << "%resManifold:\n" << printRes(resManifold) << std::endl;
      myFile << "%resSO3noS2:\n" << printRes(resSO3noS2) << std::endl;
      myFile << "%resRealSpace:\n" << printRes(resRealSpace) << std::endl;

      mySumFile << "%=============== " << nCubes << " cubes ==============" << std::endl;
      mySumFile << "%resManifold:\n" << printSummary(resManifold) << std::endl;
      mySumFile << "%resSO3noS2:\n" << printSummary(resSO3noS2) << std::endl;
      mySumFile << "%resRealSpace:\n" << printSummary(resRealSpace) << std::endl;

      myTimeFile << "%=============== " << nCubes << " cubes ==============" << std::endl;
      myTimeFile << "%" << Mptr->name() << std::endl;
      myTimeFile << "%" << MSO3noS2ptr->name() << std::endl;
      myTimeFile << "%" << Rptr->name() << std::endl;
      myTimeFile << "%resManifold:\n" << printTimes(resManifold, "manifold_" + std::to_string(nCubes) + "cubes.") << std::endl;
      myTimeFile << "%resSO3noS2:\n" << printTimes(resSO3noS2, "SO3noS2_" + std::to_string(nCubes) + "cubes.") << std::endl;
      myTimeFile << "%resRealSpace:\n" << printTimes(resRealSpace, "RealSpace_" + std::to_string(nCubes) + "cubes.") << std::endl;

      myTimePerIterFile << "%=============== " << nCubes << " cubes ==============" << std::endl;
      myTimePerIterFile << "%resManifold:\n" << printAverageTimesPerIter(resManifold, "manifold_" + std::to_string(nCubes) + "cubes.") << std::endl;
      myTimePerIterFile << "%resSO3noS2:\n" << printAverageTimesPerIter(resSO3noS2, "SO3noS2_" + std::to_string(nCubes) + "cubes.") << std::endl;
      myTimePerIterFile << "%resRealSpace:\n" << printAverageTimesPerIter(resRealSpace, "RealSpace_" + std::to_string(nCubes) + "cubes.") << std::endl;
    }
    else
    {
      std::cout << "=============== " << nCubes << " cubes ==============" << std::endl;
      std::cout << Mptr->name() << std::endl;
      std::cout << MSO3noS2ptr->name() << std::endl;
      std::cout << Rptr->name() << std::endl;
      std::cout << "resManifold:\n" << printRes(resManifold) << std::endl;
      std::cout << "resSO3noS2:\n" << printRes(resSO3noS2) << std::endl;
      std::cout << "resRealSpace:\n" << printRes(resRealSpace) << std::endl;
    }
  }

  return 0;
}
