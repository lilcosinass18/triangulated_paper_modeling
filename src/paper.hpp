#ifndef PAPER_MODELING_PAPER_HPP
#define PAPER_MODELING_PAPER_HPP

#include "geometry.hpp"
#include "utils.hpp"

class Paper {

public:
    Paper() = default;

    Paper(double width, double length) {
        vertices.push_back(Point(0, 0, 0));
        vertices.push_back(Point(length, 0, 0));
        vertices.push_back(Point(length, width, 0));
        vertices.push_back(Point(0, width, 0));
        edges.push_back(std::make_pair(0,1));
        edges.push_back(std::make_pair(1,2));
        edges.push_back(std::make_pair(2,3));
        edges.push_back(std::make_pair(3,0));
    }

    ~Paper() = default;

    void addCrease(const Point &start, const Point &end);
    void rotateVertices(int indexCrease, int indexVertex, double angleDegrees);
    void writePaperToFile(const std::string& filename);

private:
    void addEdge(int startIdx, int endIdx);
    bool creaseIntersection(std::pair<int, int> firstCrease, std::pair<int, int> secondCrease);
    bool isPointOnEdge(const Point &point, const Point &edgeStart, const Point &edgeEnd);
    std::vector<Point> pointsNeedToRotate(const int indexStart, const int indexEnd, const int index);
    bool pointOnEdges(const Point &point);
    bool checkPointsOnEdge(const Point &lineStart, const Point &lineEnd);


    std::vector<Point> vertices;
    std::vector<std::pair<int, int>> edges;
    std::vector<std::pair<int, int>> creases;
};

void processPaperFromJson(const std::string& filename);

#endif //PAPER_MODELING_PAPER_HPP
