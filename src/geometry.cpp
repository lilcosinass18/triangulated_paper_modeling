#include "geometry.hpp"

Point::Point() {
    Point(0, 0, 0);
};

Point Point::get_mid(const Point& another, double t) const {
    return Point(x + (another.x - x)*t, y + (another.y - y)*t, z + (another.z - z)*t);
}

double Point::get_dist(const Point& another) const {
    return sqrt(pow(x - another.x, 2) + pow(y - another.y, 2) + pow(z - another.z, 2));
}

double Point::get_angle(const Point& first, const Point& second) const {
    double a = first.get_dist(second);
    double b = get_dist(second);
    double c = get_dist(first);
    double cos = (pow(b , 2) + pow(c, 2) - pow(a, 2))/(2 * b * c);
    if (cos > 1) {
        cos = 1;
    }
    else if (cos < -1) {
        cos = -1;
    }
    return acos(cos) * 180 / M_PI;
}

Point crossProduct(const Point &first, const Point &second) {
    return Point(first.y * second.z - first.z * second.y,
                 first.z * second.x - first.x * second.z,
                 first.x * second.y - first.y * second.x);
}

double dotProduct(const Point &first, const Point &second) {
    return first.x * second.x + first.y * second.y + first.z * second.z;
}

// Проверка, лежит ли точка q на отрезке pr
bool onSegment(const Point& p, const Point& q, const Point& r)
{
    if (q.x <= std::max(p.x, r.x) && q.x >= std::min(p.x, r.x) &&
        q.y <= std::max(p.y, r.y) && q.y >= std::min(p.y, r.y) &&
        q.z <= std::max(p.z, r.z) && q.z >= std::min(p.z, r.z))
        return true;

    return false;
}

// Определение ориентации трех точек (p, q, r)
// 0 - коллинеарны
// 1 - по часовой стрелке
// 2 - против часовой стрелки
int orientation(const Point& p, const Point& q, const Point& r)
{
    int val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);

    if (val == 0) return 0;

    return (val > 0) ? 1 : 2;
}

// Проверка пересечения двух отрезков 'p1q1' и 'p2q2'
bool doIntersect(const Point& p1, const Point& q1, const Point& p2, const Point& q2)
{
    int o1 = orientation(p1, q1, p2);
    int o2 = orientation(p1, q1, q2);
    int o3 = orientation(p2, q2, p1);
    int o4 = orientation(p2, q2, q1);

    if (o1 != o2 && o3 != o4) {
        // Отрезки пересекаются, но нужно проверить, что точки пересечения не являются концами отрезков
        if ((o1 == 0 && onSegment(p1, p2, q1)) ||
            (o2 == 0 && onSegment(p1, q2, q1)) ||
            (o3 == 0 && onSegment(p2, p1, q2)) ||
            (o4 == 0 && onSegment(p2, q1, q2))) {
            return false; // Отрезки пересекаются только своими концами
        }
        return true;
    }

    if (o1 == 0 && onSegment(p1, p2, q1)) return true;

    if (o2 == 0 && onSegment(p1, q2, q1)) return true;

    if (o3 == 0 && onSegment(p2, p1, q2)) return true;

    if (o4 == 0 && onSegment(p2, q1, q2)) return true;

    return false;
}

double distance(const Point& first, const Point& second) {
    return sqrt(pow(first.x - second.x, 2) + pow(first.y - second.y, 2) + pow(first.z - second.z, 2));
}

// Функция вращения точки вокруг оси
Point rotatePointAroundAxis(const Point &point,
                            const Point &lineStart,
                            const Point &lineEnd,
                            double angleDegrees) {

    double angleRadians = angleDegrees * M_PI / 180.0;

    // Вычисление вектора оси вращения
    double axisX = lineEnd.x - lineStart.x;
    double axisY = lineEnd.y - lineStart.y;
    double axisZ = lineEnd.z - lineStart.z;
    double axisLength = sqrt(axisX * axisX + axisY * axisY + axisZ * axisZ);
    axisX /= axisLength;
    axisY /= axisLength;
    axisZ /= axisLength;

    // Координаты точки относительно начала координат (переносим начало координат)
    double x = point.x - lineStart.x;
    double y = point.y - lineStart.y;
    double z = point.z - lineStart.z;

    // Вращение вокруг оси с использованием матрицы Родрига
    double cosTheta = cos(angleRadians);
    double sinTheta = sin(angleRadians);

    double newX = (cosTheta + (1 - cosTheta) * axisX * axisX) * x;
    newX += ((1 - cosTheta) * axisX * axisY - axisZ * sinTheta) * y;
    newX += ((1 - cosTheta) * axisX * axisZ + axisY * sinTheta) * z;

    double newY = ((1 - cosTheta) * axisY * axisX + axisZ * sinTheta) * x;
    newY += (cosTheta + (1 - cosTheta) * axisY * axisY) * y;
    newY += ((1 - cosTheta) * axisY * axisZ - axisX * sinTheta) * z;

    double newZ = ((1 - cosTheta) * axisZ * axisX - axisY * sinTheta) * x;
    newZ += ((1 - cosTheta) * axisZ * axisY + axisX * sinTheta) * y;
    newZ += (cosTheta + (1 - cosTheta) * axisZ * axisZ) * z;

    // Возвращаем новую точку, возвращая начало координат на место
    return Point(newX + lineStart.x, newY + lineStart.y, newZ + lineStart.z);
}