#ifndef INCLUDED_util_hpp_
#define INCLUDED_util_hpp_

// useful macro
#define LOOP for(int i = 0; i < NL; i++) for(int j = 0; j < NL; j++) for(int k = 0; k < NL; k++)
#define LOOPLONG for(int i = 0; i < NLnoise; i++) for(int j = 0; j < NLnoise; j++) for(int k = 0; k < NLnoise; k++)


std::vector<std::vector<double>> transpose(const std::vector<std::vector<double>>& matrix) {
  if (matrix.empty() || matrix[0].empty()) return {}; // 空の場合の処理

  size_t rows = matrix.size();
  size_t cols = matrix[0].size();
  std::vector<std::vector<double>> result(cols, std::vector<double>(rows));

  for (size_t i = 0; i < rows; i++) {
    for (size_t j = 0; j < cols; j++) {
      result[j][i] = matrix[i][j];
    }
  }
  return result;
}

#endif