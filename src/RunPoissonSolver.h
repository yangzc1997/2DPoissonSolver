// RunPoissonSolver.h

#ifndef RUN_POISSON_SOLVER_H
#define RUN_POISSON_SOLVER_H

#include "Core_Export.h"
#include "Help.h"
#include "Timer.h"
#include "Mesh.h"
#include "BoundaryCondition.h"
#include "PoissonSolver.h"
#include "ReadInputData.h"
#include <filesystem>
#include <memory>

namespace Poisson{

class POISSONCORE_API RunPoissonSolver{
public:
    explicit RunPoissonSolver(const std::filesystem::path& jsonPath);
    
    RunPoissonSolver(const RunPoissonSolver&) = delete;
    RunPoissonSolver& operator=(const RunPoissonSolver&) = delete;
    ~RunPoissonSolver() = default;

    bool readFromJson();
    bool simulate();

private:
    std::filesystem::path m_jsonPath;
    std::unique_ptr<ReadInputData> m_readInputData;
    std::unique_ptr<Mesh> m_mesh;
    std::unique_ptr<BoundaryCondition> m_boundaryCondition;
    std::unique_ptr<PoissonSolver> m_solver;
    Timer m_timer;
};

}  // namespace Poisson

#endif  //RUN_POISSON_SOLVER_H