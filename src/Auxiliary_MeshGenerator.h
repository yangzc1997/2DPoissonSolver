// Auxiliary_MeshGenerator.h
#ifndef AUXILIARY_MESH_GENERATOR_H
#define AUXILIARY_MESH_GENERATOR_H

#include <string>
#include "Mesh.h"

namespace Poisson {

namespace Auxiliary {

// 网格生成模块
namespace MeshGenerator {
    Mesh generate_mesh(double lx, double ly, int Nx, int Ny, const std::string& mesh_type);
}

} // namespace Auxiliary

} // namespace Poisson

#endif // AUXILIARY_MESH_GENERATOR_H