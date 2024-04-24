#ifndef PAPER_MODELING_PERMUTION_HPP
#define PAPER_MODELING_PERMUTION_HPP

#include <iostream>
#include <vector>
#include <algorithm>
#include <string>

bool isUniqueCycle(const std::vector<int>& permutation, const std::vector<std::vector<int>>& cycles);

std::vector<int> mirrorPermutation(const std::vector<int>& permutation);

std::vector<std::vector<int>> generatePermutations(std::vector<int>& digits);

#endif //PAPER_MODELING_PERMUTION_HPP
