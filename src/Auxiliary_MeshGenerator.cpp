// Auxiliary_MeshGenerator.cpp
#include "Auxiliary_MeshGenerator.h"
#include <vector>

namespace Poisson {

namespace Auxiliary {

namespace MeshGenerator{

    Mesh generate_mesh(double lx, double ly, int Nx, int Ny, const std::string& mesh_type)
    {
        Mesh mesh;

        // 生成节点
        const int total_nodes = (Nx + 1) * (Ny + 1);
        mesh.nodes.reserve(total_nodes);
        
        double dx = lx / Nx;
        double dy = ly / Ny;
        for (int j = 0; j <= Ny; j++) {
            for (int i = 0; i <= Nx; i++) {
                int id = j * (Nx + 1) + i;
                mesh.nodes.emplace_back(i * dx, j * dy, id);
            }
        }
        
        // 生成单元
        if (mesh_type == "4" || mesh_type == "rectangle") {
            // 生成四边形单元
            mesh.elements.reserve(Nx * Ny);  // 预分配单元内存
            for (int j = 0; j < Ny; j++) {
                for (int i = 0; i < Nx; i++) {
                    int n1 = j * (Nx + 1) + i;
                    int n2 = n1 + 1;
                    int n3 = n1 + (Nx + 1);
                    int n4 = n3 + 1;

                    mesh.elements.emplace_back(std::vector<Node*>{
                        &mesh.nodes[n1], &mesh.nodes[n2], &mesh.nodes[n4], &mesh.nodes[n3]
                    });
                }
            }
        } else {
            // 生成三角形单元
            mesh.elements.reserve(2 * Nx * Ny);
            for (int j = 0; j < Ny; j++) {
                for (int i = 0; i < Nx; i++) {
                    int n1 = j * (Nx + 1) + i;
                    int n2 = n1 + 1;
                    int n3 = n1 + (Nx + 1);
                    int n4 = n3 + 1;

                    // 一个矩形可以划分为两个三角形
                    mesh.elements.emplace_back(std::vector<Node*>{
                        &mesh.nodes[n1], &mesh.nodes[n2], &mesh.nodes[n3]
                    });
                    mesh.elements.emplace_back(std::vector<Node*>{
                        &mesh.nodes[n2], &mesh.nodes[n4], &mesh.nodes[n3]
                    });
                }
            }
        }

        return mesh;
    }

}  // namespace MeshGenerator

} // namespace Auxiliary

} // namespace Poisson
