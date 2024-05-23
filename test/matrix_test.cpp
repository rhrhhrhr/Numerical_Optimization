//
// Created by 16058 on 2024/5/23.
//

#include "../Matrix/Matrix.hpp"

int main() {
    Matrix::Matrix<int> vector1 = {1, 2, 3};
    Matrix::Matrix<double> vector2 = {1.1, 1.2, 1.3};

    std::cout << "vector1:" << std::endl;
    std::cout << vector1 << std::endl << std::endl;
    std::cout << "vector2:" << std::endl;
    std::cout << vector2 << std::endl << std::endl;

    auto res1 = vector1 + vector2;
    auto res2 = vector1.transpose() * vector2;
    auto res3 = vector1 * 3.1;
    auto res4 = Matrix::norm1(vector1);
    auto res5 = Matrix::norm2(vector1);

    std::cout << "vector1 * vector2:" << std::endl;
    std::cout << res1 << std::endl << std::endl;
    std::cout << "vector1.T * vector2:" << std::endl;
    std::cout << res2 << std::endl << std::endl;
    std::cout << "vector1 * 3.1:" << std::endl;
    std::cout << res3 << std::endl << std::endl;
    std::cout << "norm1(vector1) = " << res4 << std::endl << std::endl;
    std::cout << "norm2(vector1) = " << res5 << std::endl << std::endl;

    Matrix::Matrix<int> matrix1 = {{1, 2, 3},
                                   {4, 5, 6},
                                   {7, 8, 7}};
    Matrix::Matrix<double> matrix2 = {{2.1, 1.2, 6.3},
                                      {2.2, 3.2, 5.7}};
    Matrix::Matrix<int> matrix3;
    Matrix::Matrix<int> matrix4(3, 6);
    Matrix::Matrix<int> matrix5(5, 8, 1);

    std::cout << "matrix1:" << std::endl;
    std::cout << matrix1 << std::endl << std::endl;
    std::cout << "matrix2:" << std::endl;
    std::cout << matrix2 << std::endl << std::endl;
    std::cout << "matrix3:" << std::endl;
    std::cout << matrix3 << std::endl << std::endl;
    std::cout << "matrix4:" << std::endl;
    std::cout << matrix4 << std::endl << std::endl;
    std::cout << "matrix5:" << std::endl;
    std::cout << matrix5 << std::endl << std::endl;

    auto res6 = Matrix::identity<int>(3);
    auto res7 = static_cast<Matrix::Matrix<int>>(matrix2);
    auto res8 = matrix1.inverse();
    auto res9 = matrix1 * vector2;
    auto res10 = matrix2 * matrix1;
    auto res11 = Matrix::norm1(matrix2);
    auto res12 = Matrix::norm2(matrix2);

    std::cout << "Matrix::identity(3):" << std::endl;
    std::cout << res6 << std::endl << std::endl;
    std::cout << "static_cast<Matrix::Matrix<int>>(matrix2):" << std::endl;
    std::cout << res7 << std::endl << std::endl;
    std::cout << "matrix1^-1:" << std::endl;
    std::cout << res8 << std::endl << std::endl;
    std::cout << "matrix1 * vector2:" << std::endl;
    std::cout << res9 << std::endl << std::endl;
    std::cout << "matrix2 * matrix1:" << std::endl;
    std::cout << res10 << std::endl << std::endl;
    std::cout << "norm1(matrix1) = " << res11 << std::endl << std::endl;
    std::cout << "norm2(matrix1) = " << res12 << std::endl << std::endl;

    return 0;
}
