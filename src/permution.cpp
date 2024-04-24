#include "permution.hpp"

bool isUniqueCycle(const std::vector<int>& permutation, const std::vector<std::vector<int>>& cycles) {
    for (const std::vector<int>& cycle : cycles) {
        // Проверяем циклические сдвиги
        for (size_t i = 0; i < permutation.size(); ++i) {
            std::vector<int> shiftedPermutation(permutation.begin() + i, permutation.end());
            shiftedPermutation.insert(shiftedPermutation.end(), permutation.begin(), permutation.begin() + i);
            if (shiftedPermutation == cycle)
                return false;
        }
    }
    return true;
}

std::vector<int> mirrorPermutation(const std::vector<int>& permutation) {
    std::vector<int> mirrored = permutation;
    std::reverse(mirrored.begin() + 1, mirrored.end());
    return mirrored;
}

std::vector<std::vector<int>> generatePermutations(std::vector<int>& digits) {
    std::vector<std::vector<int>> permutations;
    std::vector<std::vector<int>> cycles;

    do {
        if (isUniqueCycle(digits, cycles)) {
            permutations.push_back(digits);
            cycles.push_back(digits);
            cycles.push_back(mirrorPermutation(digits)); // Добавляем зеркальное отображение
        }
    } while (std::next_permutation(digits.begin(), digits.end()));

    return permutations;
}