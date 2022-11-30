//
// Created by Hao-Xin on 2022/9/21.
//

#ifndef HONEYCOMBHEISENBERG_MONTECARLO_SRC_COMMON_H_
#define HONEYCOMBHEISENBERG_MONTECARLO_SRC_COMMON_H_


#include <numeric>
#include <fstream>
#include <algorithm>

template<typename DataType>
void DumpData(
    const std::string &filename,
    const std::vector<DataType> &data
) {
  std::ofstream ofs(filename, std::ofstream::binary);
  for (auto datum: data) {
    ofs << datum << '\n';
  }
  ofs << std::endl;
  ofs.close();
}

template<typename T>
T Mean(const std::vector<T> data) {
  if (data.empty()) {
    return T(0);
  }
  auto const count = static_cast<T>(data.size());
  return std::reduce(data.begin(), data.end()) / count;
}

/// Note the definition
template<typename T>
T Variance(const std::vector<T> data,
           const T &mean) {
  size_t data_size = data.size();
  std::vector<T> diff(data_size);
  std::transform(data.begin(), data.end(), diff.begin(), [mean](double x) { return x - mean; });
  T sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
  auto const count = static_cast<T>(data_size);
  return sq_sum / count;
}

template<typename T>
T Variance(const std::vector<T> data) {
  return Variance(data, Mean(data));
}




#endif //HONEYCOMBHEISENBERG_MONTECARLO_SRC_COMMON_H_
