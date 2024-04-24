#ifndef PAPER_MODELING_GEOMETRY_HPP
#define PAPER_MODELING_GEOMETRY_HPP

#define _USE_MATH_DEFINES

#include <vector>
#include <cstdlib>
#include <cmath>
#include "permution.hpp"

/*! @brief Структура точка
 *
 * Класс реализует точку.
 */
struct Point {
    /*!
     * @param x - координата по оси х
     * @param y - координата по оси у
     * @param z - координата по оси z
     * @note Создает точку по координатам
     */
    Point(double x, double y, double z): x(x), y(y), z(z) {};

    /*!
     * @param another - другая точка
     * @note Создает точку по другой точке
     */
    Point(const Point& another): x(another.x), y(another.y), z(another.z) {};

    /*!
     * @note Дефолтное значение точки
     */
    Point();

    /*!
     * @note Деструктор
     */
    ~Point() = default;

    /*!
     * @param another - другая точка
     * @note Оператор +
     */
    Point operator+(const Point& another) const {
        return Point(x+another.x, y + another.y, z + another.z);
    }

    /*!
     * @param another - другая точка
     * @note Оператор -
     */
    Point operator-(const Point& another) const {
        return Point(x-another.x, y - another.y, z - another.z);
    }

    /*!
     * @param coeff - коеффициент умножения
     * @note Оператор *
     */
    Point operator*(double coeff) const{
        return Point(x * coeff, y * coeff, z * coeff);
    }

    /*!
     * @param another - другая точка
     * @note Оператор ==
     */
    bool operator==(const Point& another) const {
        return fabs(x - another.x) < 1e-6 && fabs(y - another.y) < 1e-6 && fabs(z - another.z) < 1e-6;
    }

    /*!
     * @param another - другая точка
     * @note Оператор !=
     */
    bool operator!=(const Point& another) const {
        return fabs(x - another.x) > 1e-6 || fabs(y - another.y) > 1e-6 || fabs(z - another.z) > 1e-6;
    }

    /*!
     * @param another - точка другого конца отрезка
     * @param t - отношение в котором делится отрезок
     * @note Функция, возвращающая точку, делящую отрезок в отношении t
     */
    Point get_mid(const Point& another, double t = 0.5) const;

    /*!
     * @param another - другая точка
     * @note Функция, находящая расстояние между точками
     */
    double get_dist(const Point& another) const;

    /*!
     * @param first - другая точка
     * @param second - другая точка
     * @note Функция, находящая угол, между лучами, проходящими через переданные точки
     */
    double get_angle(const Point& first, const Point& second) const;

    /*!
     * координата по оси х
     */
    double x;
    /*!
     * координата по оси у
     */
    double y;
    /*!
     * координата по оси z
     */
    double z;
};

/*!
 * @param first - координаты первого ветора
 * @param second - координаты второго вектора
 * @note Функция, находящая векторное произведение двух векторов
 */
Point crossProduct(const Point& first, const Point& second);

/*!
 * @param first - координаты первого ветора
 * @param second - координаты второго вектора
 * @note Функция, находящая скалярное произведение двух векторов
 */
double dotProduct(const Point& first, const Point& second);

/*!
 * @param p - координаты точки p
 * @param q - координаты точки q
 * @param r - координаты точки r
 * @note Функция, которая проверяет лежит ли точка q на отрезке pr
 */
bool onSegment(const Point& p, const Point& q, const Point& r);

/*!
 * @param p - координаты точки p
 * @param q - координаты точки q
 * @param r - координаты точки r
 * @note Функция, которая определяет ориентацию точек p, q, r
 */
int orientation(const Point& p, const Point& q, const Point& r);

/*!
 * @param p1 - координаты точки p1
 * @param q1 - координаты точки q1
 * @param p2 - координаты точки p2
 * @param q2 - координаты точки q2
 * @note Функция, которая определяет пересекаются ли отрезки p1q1 и p2q2
 */
bool doIntersect(const Point& p1, const Point& q1, const Point& p2, const Point& q2);

/*!
 * @param first - координаты первой точки
 * @param second - координаты второй точки
 * @note Функция, которая находит L2-метрику от точек first и second
 */
double distance(const Point& first, const Point& second);

/*!
 * @param point - координаты точки, которую будем вращать относительно оси
 * @param lineStart - координаты первой точки, задающей ось
 * @param lineEnd - координаты второй точки, задающей ось
 * @param angleDegrees - угол в градусах, на который нужно повернуть точку point относительно оси вращения
 * @note Функция, которая вращает точку point на угол angleDegrees отеносительно оси заданной точками lineStart и lineEnd
 */
Point rotatePointAroundAxis(const Point &point,
                            const Point &lineStart,
                            const Point &lineEnd,
                            double angleDegrees);

bool isPointsQuadrangle(int a, int b, int c, int d, std::vector<std::pair<int, int>> edges);

bool isPointsPentagon(int a, int b, int c, int d, int e, std::vector<std::pair<int, int>> edges);

bool isPointsHexagon(int a, int b, int c, int d, int e, int f,std::vector<std::pair<int, int>> edges);

#endif //PAPER_MODELING_GEOMETRY_HPP
