// Mesh.cpp

#include "Mesh.h"
#include <iostream>

namespace Poisson{

void Mesh::print_mesh_information(bool nodesInf) const {
    std::cout << "网格信息:\n";
    std::cout << "   节点数: " << get_num_nodes() << "\n";
    std::cout << "   单元数: " << get_num_elements() << "\n";
    
    std::cout << "  【单元信息】: " << "\n";
    int n_triangle = 0;
    for (const auto& element: elements) {
        if (element.get_num_nodes() == 3) {
             n_triangle++;
        };
    }
    std::cout << "      三角形单元数：" << n_triangle << "\n";
    std::cout << "      四边形单元数：" << elements.size()-n_triangle << "\n";
    
    if (nodesInf){
        std::cout << "  【节点信息】: " << "\n";
        for (const auto& node: nodes) {
            std::cout << "          Node Id: " << node.id  << " NodePosition: (" << node.x << ", " << node.y << ")"  << std::endl;
        }
    }
}

} // namespace Poisson
