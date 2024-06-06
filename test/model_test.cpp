//
// Created by 16058 on 2024/5/23.
//

#include "../Model/Expression.hpp"
#include "../Matrix/Matrix.hpp"

int main() {
    Matrix::Matrix<Model::Expression<double>> var = {{1.1, 2., 3., 4.},
                                                     {2.1, 3., 1.1, 4.4}};

    std::cout << var << std::endl;

    return 0;
}
