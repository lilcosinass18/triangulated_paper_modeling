#include "paper.hpp"
#include <unordered_set>
#include <queue>
#include <fstream>

void Paper::addEdge(int startIdx, int endIdx) {
    edges.emplace_back(startIdx, endIdx);
}

bool Paper::creaseIntersection(std::pair<int, int> firstCrease, std::pair<int, int> secondCrease) {

    if (!doIntersect(vertices[firstCrease.first], vertices[firstCrease.second],
                     vertices[secondCrease.first], vertices[secondCrease.second])) {
        return false;
    }

    return true;
}

bool Paper::isPointOnEdge(const Point &point, const Point &edgeStart, const Point &edgeEnd) {
    Point vec1 = point - edgeStart;
    Point vec2 = point - edgeEnd;

    Point cross = crossProduct(vec1, vec2);

    if (cross.x != 0 || cross.y != 0 || cross.z != 0) {
        return false;
    }

    // Проверка нахождения точки в пределах отрезка
    Point edgeVec = edgeEnd - edgeStart;

    double dotProduct1 = dotProduct(vec1, edgeVec);
    double dotProduct2 = dotProduct(vec2, edgeVec);

    if (dotProduct1 * dotProduct2 > 0) {
        return false;
    }

    return true;
}


std::vector<Point> Paper::pointsNeedToRotate(const int indexStart, const int indexEnd, const int index) {
    std::unordered_set<int> pointsSet; // Сюда будем сохранять индексы вершин на пути
    std::queue<int> q;
    q.push(indexStart);

    // Создаем вектор предков для восстановления пути
    std::vector<int> parent(vertices.size(), -1);

    while (!q.empty()) {
        int current = q.front();
        q.pop();

        if (current == indexEnd) {
            while (current != indexStart) {
                pointsSet.insert(current);
                current = parent[current];
            }
            break;
        }

        for (const auto& edge : edges) {
            if (edge.first == current || edge.second == current) {
                int next = (edge.first == current) ? edge.second : edge.first;
                if (parent[next] == -1) {
                    parent[next] = current;
                    q.push(next);
                }
            }
        }
    }

    pointsSet.erase(indexStart);
    pointsSet.erase(indexEnd);

    // Проверяем, содержит ли множество вершину с индексом index
    if (pointsSet.find(index) == pointsSet.end()) {
        std::vector<Point> points;
        for (int i = 0; i < vertices.size(); ++i) {
            if (i != indexStart && i != indexEnd && pointsSet.find(i) == pointsSet.end()) {
                points.push_back(vertices[i]);
            }
        }
        return points;
    }



    std::vector<Point> points;
    for (int v : pointsSet) {
        points.push_back(vertices[v]);
    }
    return points;
}


// Для проверки принадлежности lineStart и lineEnd к одному из рёбер
bool Paper::pointOnEdges(const Point &point) {
    for (const auto &edge : edges) {
        const Point &edgeStart = vertices[edge.first];
        const Point &edgeEnd = vertices[edge.second];

        if (isPointOnEdge(point, edgeStart, edgeEnd)) {
            return true;
        }
    }
    return false;
}

// Для проверки принадлежности lineStart и lineEnd к одному из рёбер
bool Paper::checkPointsOnEdge(const Point &lineStart, const Point &lineEnd) {
    for (const auto &edge : edges) {
        const Point &edgeStart = vertices[edge.first];
        const Point &edgeEnd = vertices[edge.second];

        if (isPointOnEdge(lineStart, edgeStart, edgeEnd) && isPointOnEdge(lineEnd, edgeStart, edgeEnd)) {
            return true; // Если точки лежат на одном из ребер
        }
    }
    return false; // Точки не лежат на одом ребре
}

// Для проверки принадлежности lineStart и lineEnd к одному из рёбер
bool Paper::checkPointsOnEdgeAndCrease(const Point &lineStart, const Point &lineEnd) {
    std::vector<std::pair<int, int>> allEdges;
    for (const auto& edge : edges) {
        allEdges.push_back(edge);
    }
    for (const auto& crease : creases) {
        allEdges.push_back(crease);
    }
    for (const auto &edge : allEdges) {
        const Point &edgeStart = vertices[edge.first];
        const Point &edgeEnd = vertices[edge.second];

        if (isPointOnEdge(lineStart, edgeStart, edgeEnd) && isPointOnEdge(lineEnd, edgeStart, edgeEnd)) {
            return true; // Если точки лежат на одном из ребер
        }
    }
    return false; // Точки не лежат на одом ребре
}


void Paper::addCrease(const Point &start, const Point &end) {
    // проверка на то, что вершины принадлежат ребрам
    if (!pointOnEdges(start) || !pointOnEdges(end)) {
        throw std::logic_error("Can't add a crease, because one of points not in edges");
    }

    //проверка, что вершины не лежат на одном ребре
    if (checkPointsOnEdge(start, end)) {
        throw std::invalid_argument("Points does lie on the same edge");
    }

    std::vector<std::pair<int, int>> edgesToRemove = {};
    std::vector<std::pair<int, int>> edgesToAdd = {};

    std::pair<int, int> creaseCandidate = {};

    // тут идет логика обновления вершин и ребер
    if ((std::find(vertices.begin(), vertices.end(), start) == vertices.end())
        && (std::find(vertices.begin(), vertices.end(), end) == vertices.end())) {

        creaseCandidate = std::make_pair(vertices.size(), vertices.size() + 1);
        vertices.push_back(start);
        vertices.push_back(end);

        for (int i = 0; i < edges.size(); ++i) {
            const auto &edge = edges[i];
            const Point &edgeStart = vertices[edge.first];
            const Point &edgeEnd = vertices[edge.second];

            if (isPointOnEdge(start, edgeStart, edgeEnd)) {
                edgesToRemove.push_back(edge);
                edgesToAdd.push_back(std::make_pair(findIndex(vertices, edgeStart), findIndex(vertices, start)));
                edgesToAdd.push_back(std::make_pair(findIndex(vertices, start), findIndex(vertices, edgeEnd)));
            }
            if (isPointOnEdge(end, edgeStart, edgeEnd)) {
                edgesToRemove.push_back(edge);
                edgesToAdd.push_back(std::make_pair(findIndex(vertices, edgeStart), findIndex(vertices, end)));
                edgesToAdd.push_back(std::make_pair(findIndex(vertices, end), findIndex(vertices, edgeEnd)));
            }
        }
    } else if ((std::find(vertices.begin(), vertices.end(), start) != vertices.end())
               && (std::find(vertices.begin(), vertices.end(), end) == vertices.end())) {
        creaseCandidate = std::make_pair(findIndex(vertices,start), vertices.size());
        vertices.push_back(end);
        for (int i = 0; i < edges.size(); ++i) {
            const auto &edge = edges[i];
            const Point &edgeStart = vertices[edge.first];
            const Point &edgeEnd = vertices[edge.second];
            if (isPointOnEdge(end, edgeStart, edgeEnd)) {
                edgesToRemove.push_back(edge);
                edgesToAdd.push_back(std::make_pair(findIndex(vertices, edgeStart), findIndex(vertices, end)));
                edgesToAdd.push_back(std::make_pair(findIndex(vertices, end), findIndex(vertices, edgeEnd)));
            }
        }
    } else if ((std::find(vertices.begin(), vertices.end(), start) == vertices.end())
               && (std::find(vertices.begin(), vertices.end(), end) != vertices.end())) {
        creaseCandidate = std::make_pair(findIndex(vertices,end), vertices.size());
        vertices.push_back(start);
        for (int i = 0; i < edges.size(); ++i) {
            const auto &edge = edges[i];
            const Point &edgeStart = vertices[edge.first];
            const Point &edgeEnd = vertices[edge.second];
            if (isPointOnEdge(start, edgeStart, edgeEnd)) {
                edgesToRemove.push_back(edge);
                edgesToAdd.push_back(std::make_pair(findIndex(vertices, edgeStart), findIndex(vertices, start)));
                edgesToAdd.push_back(std::make_pair(findIndex(vertices, start), findIndex(vertices, edgeEnd)));
            }
        }
    } else if ((std::find(vertices.begin(), vertices.end(), start) != vertices.end())
               && (std::find(vertices.begin(), vertices.end(), end) != vertices.end())) {
        creaseCandidate = std::make_pair(findIndex(vertices,start), findIndex(vertices, end));
    }

    std::pair<int, int> crease_tmp = {};

    if (creases.empty()) {
        creases.push_back(creaseCandidate);
    } else {
        for (const auto tmp: creases) {
            if (!creaseIntersection(tmp, creaseCandidate)) {
                crease_tmp = creaseCandidate;
            } else {
                throw std::logic_error("Crease does intersected existing set of creases");
            }
        }
        creases.push_back(crease_tmp); // добавляем ребро, если оно не пересекает ничего
    }

    // Удаление пересекающих рёбер
    for (const auto &edgeToRemove: edgesToRemove) {
        auto it = std::find(edges.begin(), edges.end(), edgeToRemove);
        if (it != edges.end()) edges.erase(it);
    }

    // Добавление новых рёбер
    edges.insert(edges.end(), edgesToAdd.begin(), edgesToAdd.end());
}



void Paper::rotateVertices(int indexCrease, int indexVertex, double angleDegrees) {
    std::vector<Point> pointsNeedToRotate = this->pointsNeedToRotate(creases[indexCrease].first,
                                                                creases[indexCrease].second, indexVertex);
    std::vector<int> indexNeedToRotate;
    for (auto &vertex : pointsNeedToRotate) {
        indexNeedToRotate.push_back(findIndex(vertices, vertex));
    }
    for (auto &index : indexNeedToRotate) {
        for (auto &vertex : vertices) {
            if (findIndex(vertices,vertex) == index) {
                vertex = rotatePointAroundAxis(vertex, vertices[creases[indexCrease].first],
                                               vertices[creases[indexCrease].second], angleDegrees);
            }
        }
    }
}

void Paper::writePaperToFile(const std::string& filename) {
    std::ofstream outputFile(filename);
    if (outputFile.is_open()) {
        // Записываем информацию о вершинах
        for (size_t i = 0; i < vertices.size(); ++i) {
            const auto& vertex = vertices[i];
            outputFile << "Vertex: " << vertex.x << " " << vertex.y << " " << vertex.z << "\n";
        }
        // Записываем информацию о рёбрах
        for (const auto& edge : edges) {
            outputFile << "Edge: " << edge.first << " " << edge.second << "\n";
        }
        // Записываем информацию о сгибах
        for (const auto& crease : creases) {
            outputFile << "Crease: " << crease.first << " " << crease.second << "\n";
        }
        outputFile.close();
    } else {
        std::cerr << "Unable to open file " << filename << std::endl;
    }
}

std::vector<std::vector<int>> Paper::findQuadrangles() {
    std::vector<std::pair<int, int>> allEdges;
    for (const auto& edge : edges) {
        allEdges.push_back(edge);
    }
    for (const auto& crease : creases) {
        allEdges.push_back(crease);
    }
    std::vector<std::vector<int>> candidateQuadrangles;
    for (int i = 0; i < vertices.size() - 3; ++i) {
        for (int j = i + 1; j < vertices.size() - 2; ++j) {
            for (int k = j + 1; k < vertices.size() - 1; ++k) {
                for (int l = k + 1; l < vertices.size(); ++l) {
                    candidateQuadrangles.push_back({i, j, k, l});
                }
            }
        }
    }

    std::vector<std::vector<int>> quadrangles;
    for (const auto& candidate : candidateQuadrangles) {
        if (isPointsQuadrangle(candidate[0], candidate[1], candidate[2], candidate[3], allEdges)) {
            quadrangles.push_back(candidate);
        }
    }
    return quadrangles;
}

std::vector<std::vector<int>> Paper::findPentagons() {
    std::vector<std::pair<int, int>> allEdges;
    for (const auto& edge : edges) {
        allEdges.push_back(edge);
    }
    for (const auto& crease : creases) {
        allEdges.push_back(crease);
    }
    std::vector<std::vector<int>> candidatePentagons;
    for (int i = 0; i < vertices.size() - 4; ++i) {
        for (int j = i + 1; j < vertices.size() - 3; ++j) {
            for (int k = j + 1; k < vertices.size() - 2; ++k) {
                for (int l = k + 1; l < vertices.size() - 1; ++l) {
                    for (int m = l + 1; m < vertices.size(); ++m) {
                        candidatePentagons.push_back({i, j, k, l, m});
                    }
                }
            }
        }
    }
    std::vector<std::vector<int>> pentagons;
    for (const auto& candidate : candidatePentagons) {
        if (isPointsPentagon(candidate[0], candidate[1], candidate[2], candidate[3], candidate[4], allEdges)) {
            pentagons.push_back(candidate);
        }
    }
    return pentagons;
}

std::vector<std::vector<int>> Paper::findHexagons() {
    std::vector<std::pair<int, int>> allEdges;
    for (const auto& edge : edges) {
        allEdges.push_back(edge);
    }
    for (const auto& crease : creases) {
        allEdges.push_back(crease);
    }
    std::vector<std::vector<int>> candidateHexagons;
    for (int i = 0; i < vertices.size() - 5; ++i) {
        for (int j = i + 1; j < vertices.size() - 4; ++j) {
            for (int k = j + 1; k < vertices.size() - 3; ++k) {
                for (int l = k + 1; l < vertices.size() - 2; ++l) {
                    for (int m = l + 1; m < vertices.size() - 1; ++m) {
                        for (int n = m + 1; n < vertices.size(); ++n) {
                            candidateHexagons.push_back({i, j, k, l, m, n});
                        }
                    }
                }
            }
        }
    }
    std::vector<std::vector<int>> hexagons;
    for (const auto& candidate : candidateHexagons) {
        if (isPointsHexagon(candidate[0], candidate[1], candidate[2], candidate[3], candidate[4], candidate[5], allEdges)) {
            hexagons.push_back(candidate);
        }
    }
    return hexagons;
}

std::vector<int> Paper::findOppositeIndices(std::vector<int> indices, int startIndex) {
    std::vector<int> oppositeIndices;
    for (const auto& index : indices) {
        if (index != startIndex) {
            if (!checkPointsOnEdgeAndCrease(vertices[startIndex], vertices[index])) {
                oppositeIndices.push_back(index);
            }
        }
    }
    return oppositeIndices;
}

void Paper::triangulate(std::vector<std::vector<int>> triangulatedArea) {
    for (const auto& area : triangulatedArea) {
        int startIndex = area[0];
        std::vector<int> oppositeIndices = findOppositeIndices(area, startIndex);
        for (const auto& oppositeIndex : oppositeIndices) {
            addEdge(startIndex, oppositeIndex);
        }
    }
}


void processPaperFromJson(const std::string& filename) {
    // Открываем JSON файл
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open JSON file." << std::endl;
        return;
    }

    // Парсим JSON
    json jsonData;
    try {
        file >> jsonData;
    } catch (json::parse_error& e) {
        std::cerr << "Failed to parse JSON: " << e.what() << std::endl;
        file.close();
        return;
    }

    // Проверяем наличие необходимых полей в JSON
    if (!jsonData.contains("width") || !jsonData.contains("length") || !jsonData.contains("creases") || !jsonData.contains("angles")) {
        std::cerr << "Invalid JSON format. Missing required fields." << std::endl;
        file.close();
        return;
    }

    // Получаем данные из JSON
    double width = jsonData["width"];
    double length = jsonData["length"];
    std::vector<std::pair<Point, Point>> creases;
    std::vector<double> angles;

    //Делаем все в относительных координатах
    for (const auto& creaseData : jsonData["creases"]) {
        Point start(creaseData["start"]["x"], creaseData["start"]["y"], 0);
        Point end(creaseData["end"]["x"], creaseData["end"]["y"], 0);
        start.x = start.x * length;
        start.y = start.y * width;
        end.x = end.x * length;
        end.y = end.y * width;
        creases.push_back(std::make_pair(start, end));
    }

    for (const auto& angle : jsonData["angles"]) {
        angles.push_back(angle);
    }

    // Создаём и обрабатываем бумажку
    Paper paper(width, length);

    for (const auto& crease : creases) {
        paper.addCrease(crease.first, crease.second);
    }

    paper.writePaperToFile("../txt_info/paper_info2D.txt");

    for (int i = 0; i < angles.size(); ++i) {
        paper.rotateVertices(i, 0 ,angles[i]);
    }

    std::vector<std::vector<int>> quadrangles = paper.findQuadrangles();
    std::vector<std::vector<int>> pentagons = paper.findPentagons();
    std::vector<std::vector<int>> hexagons = paper.findHexagons();

    paper.triangulate(quadrangles);
    paper.triangulate(pentagons);
    paper.triangulate(hexagons);

    // Записываем результат в файл
    paper.writePaperToFile("../txt_info/paper_info3D.txt");

    // Закрываем файл
    file.close();

    createOBJfromTXT("../txt_info/paper_info3D.txt", "../output/paper.obj");
}