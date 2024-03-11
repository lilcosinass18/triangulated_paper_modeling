#include "utils.hpp"

void createOBJfromTXT(const std::string& txtFilename, const std::string& objFilename) {
    std::ifstream inputFile(txtFilename);
    if (!inputFile.is_open()) {
        std::cerr << "Failed to open input file." << std::endl;
        return;
    }

    std::vector<Point> vertices;
    std::vector<std::pair<int, int>> edges;
    std::vector<std::pair<int, int>> creases;

    std::string line;
    while (std::getline(inputFile, line)) {
        std::istringstream iss(line);
        std::string type;
        iss >> type;
        if (type == "Vertex:") {
            Point vertex;
            iss >> vertex.x >> vertex.y >> vertex.z;
            vertices.push_back(vertex);
        } else if (type == "Edge:") {
            std::pair<int, int> edge;
            iss >> edge.first >> edge.second;
            edges.push_back(edge);
        } else if (type == "Crease:") {
            std::pair<int, int> crease;
            iss >> crease.first >> crease.second;
            creases.push_back(crease);
            // Добавляем сгибы также как рёбра
            edges.push_back(crease);
        }
    }

    inputFile.close();

    // Создание файла OBJ
    std::ofstream objFile(objFilename);
    if (!objFile.is_open()) {
        std::cerr << "Failed to create output file." << std::endl;
        return;
    }

    // Запись вершин в файл OBJ
    for (const auto& vertex : vertices) {
        objFile << "v " << vertex.x << " " << vertex.y << " " << vertex.z << "\n";
    }

    // Запись рёбер в файл OBJ
    for (const auto& edge : edges) {
        objFile << "l " << edge.first + 1 << " " << edge.second + 1 << "\n"; // +1 для приведения к индексации с 1
    }

    objFile.close();

    std::cout << "OBJ file created: " << objFilename << std::endl;
}
