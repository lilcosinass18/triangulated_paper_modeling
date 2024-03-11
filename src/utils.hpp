#ifndef PAPER_MODELING_UTILS_HPP
#define PAPER_MODELING_UTILS_HPP

#include <iostream>
#include <vector>
#include "geometry.hpp"
#include <sstream>
#include <fstream>
#include <json.hpp>

using json = nlohmann::json;

/*!
 * @param vec - вектор шаблонного типа
 * @param value - значение шаблоного типа
 * @note Функция, которая ищет элемент с значением value в векторе vec
 */
template <typename T>
int findIndex(const std::vector<T>& vec, const T& value) {
    auto it = std::find(vec.begin(), vec.end(), value);
    if (it != vec.end()) {
        return std::distance(vec.begin(), it);
    } else {
        return -1;
    }
}

void createOBJfromTXT(const std::string& txtFilename, const std::string& objFilename);

#endif //PAPER_MODELING_UTILS_HPP
