//
//  matrixFunctions.hpp
//  linearAlgerbaEigen
//
//  Created by Mihai Dunareanu on 10.12.2017.
//  Copyright © 2017 Mihai Dunareanu. All rights reserved.
//

#ifndef matrixFunctions_hpp
#define matrixFunctions_hpp

#include "include/eigen3/Eigen/Dense"
#include <map>
#include <set>
#include <iostream>

typedef Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> eigenMatrix;
typedef Eigen::Matrix<float, Eigen::Dynamic, 1> eigenVector;
typedef std::pair<eigenMatrix, eigenVector> MyPair;

struct matrixCompare {
    bool operator() (const eigenMatrix& a, const eigenMatrix& b) const;
};

class linearRegression {
public:
    linearRegression();
    ~linearRegression();
    std::map<std::string, MyPair>::iterator begin();
    std::map<std::string, MyPair>::iterator end();
    long unsigned int getSize();
    
    std::map<std::string, MyPair> map_;
    std::map<std::string, MyPair>::iterator it;
};

class polynomialHelper {
public:
    polynomialHelper();
    ~polynomialHelper();

    template <typename Derived> Eigen::Matrix<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic> buildCompanionMatrix (const Eigen::MatrixBase<Derived>& in_monicPolynomialCoeff);
    template <typename Derived> void printPolynome(const Eigen::MatrixBase<Derived>& in_monicPolynomialCoeff);
    
    template <typename Derived> void computeEigenValues(const Eigen::MatrixBase<Derived>& in_Matrix);
    Eigen::Matrix<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic> getEigenVectors();
    Eigen::Matrix<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic> getEigenValues();
    long unsigned int getSize();
    std::set<eigenMatrix, matrixCompare>::iterator begin();
    std::set<eigenMatrix, matrixCompare>::iterator end();
    
    Eigen::ComplexEigenSolver<Eigen::Matrix<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic>> eigensolver;
    std::set<eigenMatrix, matrixCompare>::iterator it;
    std::set<eigenMatrix, matrixCompare> set_;
};

template <typename Derived> Eigen::Matrix<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic> polynomialHelper::buildCompanionMatrix (const Eigen::MatrixBase<Derived>& in_monicPolynomialCoeff) {
    
    /*
     P(X) = c0+c1*x+c2*x^2+...+cn-1*x^n-1+x^n
     
     then the frobenius companion matrix is:
     
     [0     ...   -c0  ]
     [0 1   ...   -c1  ]
C =  [0 0 1 ...   -c2  ]
     [0 0 0 ... 1 -cn-1]
     */
    
    //create companion matrix, square matrix with the same rows and columns as the polynomial has coefficients
    Eigen::Matrix<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic> companionMatrix(in_monicPolynomialCoeff.rows(), in_monicPolynomialCoeff.rows());
    //populate In in the bottom left corner, excluding first row and last column
    companionMatrix.bottomLeftCorner(in_monicPolynomialCoeff.rows()-1, in_monicPolynomialCoeff.rows()-1)  = Eigen::Matrix<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic>::Identity(in_monicPolynomialCoeff.rows()-1, in_monicPolynomialCoeff.rows()-1);
    //fill first row with zeroes
    companionMatrix.topLeftCorner(1, in_monicPolynomialCoeff.rows()) = Eigen::Matrix<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic>::Zero(1, in_monicPolynomialCoeff.rows());
    //fill last column with coefficient column multipled by -1
    companionMatrix.rightCols(1) = in_monicPolynomialCoeff*-1;
    
    return companionMatrix;
}

template <typename Derived> void polynomialHelper::printPolynome(const Eigen::MatrixBase<Derived>& in_monicPolynomialCoeff) {
    
    std::cout << "The Polynomial in monic form: " << std::endl;
    std::cout << "x^" << in_monicPolynomialCoeff.rows();
    
    for (long int i=in_monicPolynomialCoeff.rows()-1; i>0; i--) {
        if (in_monicPolynomialCoeff(i, 0) > 0) std::cout << "+";
        std::cout << in_monicPolynomialCoeff(i, 0) << "x^" << i;
    }
    if (in_monicPolynomialCoeff(0, 0) > 0) std::cout << "+";
    std::cout << in_monicPolynomialCoeff(0, 0) << std::endl;
}

template <typename Derived> void polynomialHelper::computeEigenValues(const Eigen::MatrixBase<Derived>& in_Matrix) {
    this->eigensolver.compute(in_Matrix);
    if (this->eigensolver.info() != Eigen::Success)
        throw eigensolver.info();
}

#endif /* matrixFunctions_hpp */
