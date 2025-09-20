// RunPoissonSolver.h

#ifndef RUN_POISSON_SOLVER_H
#define RUN_POISSON_SOLVER_H

#include "Core_Export.h"
#include "ReadInputData.h"
#include "Timer.h"
#include <filesystem>
#include <memory>

namespace Poisson{

class RunPoissonSolver{
public:
    POISSONCORE_API explicit RunPoissonSolver(const std::filesystem::path& jsonPath);
    
    RunPoissonSolver(const RunPoissonSolver&) = delete;
    RunPoissonSolver& operator=(const RunPoissonSolver&) = delete;
    ~RunPoissonSolver() = default;

    POISSONCORE_API bool readFromJson();
    POISSONCORE_API bool simulate();

private:
    std::filesystem::path m_jsonPath;
    std::unique_ptr<ReadInputData> m_readInputData;
    Timer m_timer;
};

}  // namespace Poisson

#endif  //RUN_POISSON_SOLVER_H