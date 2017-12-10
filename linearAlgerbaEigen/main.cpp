//
//  main.cpp
//  vectorPracticeEigen
//
//  Created by Mihai Dunareanu on 30.11.2017.
//  Copyright © 2017 Mihai Dunareanu. All rights reserved.
//

#include <iostream>
#include "include/eigen3/Eigen/Dense"
#include "include/eigen3/Eigen/Eigenvalues"
#include <map>
#include <set>
#include <chrono>

using namespace Eigen;

typedef Eigen::Matrix<float, Dynamic, Dynamic> eigenMatrix;
typedef Eigen::Matrix<float, Dynamic, 1> eigenVector;
typedef std::pair<eigenMatrix, eigenVector> MyPair;

class linearRegression {
public:
    std::map<std::string, MyPair>::iterator begin();
    std::map<std::string, MyPair>::iterator end();
    
    std::map<std::string, MyPair> map_;
    std::map<std::string, MyPair>::iterator it;
    
};

std::map<std::string, MyPair>::iterator linearRegression::begin() {
    return this->map_.begin();
}

std::map<std::string, MyPair>::iterator linearRegression::end() {
    return this->map_.end();
}

class polynomialHelper {
public:
    polynomialHelper() {};
    struct matrixCompare {
        bool operator() (const eigenMatrix& a, const eigenMatrix& b) const {
            if (a.rows() == b.rows() && a.cols() == b.cols())
                return not(a.isApprox(b));
            else return true;
        }
    };
    template <typename Derived> Eigen::Matrix<std::complex<float>, Dynamic, Dynamic> buildCompanionMatrix (const Eigen::MatrixBase<Derived>& in_monicPolynomialCoeff);
    template <typename Derived> void printPolynome(const Eigen::MatrixBase<Derived>& in_monicPolynomialCoeff);
    
    template <typename Derived> void computeEigenValues(const Eigen::MatrixBase<Derived>& in_Matrix);
    Eigen::Matrix<std::complex<float>, Dynamic, Dynamic> getEigenVectors();
    Eigen::Matrix<std::complex<float>, Dynamic, Dynamic> getEigenValues();
    std::set<eigenMatrix, matrixCompare>::iterator begin();
    std::set<eigenMatrix, matrixCompare>::iterator end();
    
    ComplexEigenSolver<Eigen::Matrix<std::complex<float>, Dynamic, Dynamic>> eigensolver;
    std::set<eigenMatrix, matrixCompare>::iterator it;
    std::set<eigenMatrix, matrixCompare> set_;
};

template <typename Derived> Eigen::Matrix<std::complex<float>, Dynamic, Dynamic> polynomialHelper::buildCompanionMatrix (const Eigen::MatrixBase<Derived>& in_monicPolynomialCoeff) {
    
    /*
     P(X) = c0+c1*x+c2*x^2+...+cn-1*x^n-1+x^n
     
     then the frobenius companion matrix is:
     
          [0     ...   -c0  ]
          [0 1   ...   -c1  ]
     C =  [0 0 1 ...   -c2  ]
          [0 0 0 ... 1 -cn-1]
     */
    
    //create companion matrix, square matrix with the same rows and columns as the polynomial has coefficients
    Eigen::Matrix<std::complex<float>, Dynamic, Dynamic> companionMatrix(in_monicPolynomialCoeff.rows(), in_monicPolynomialCoeff.rows());
    //populate In in the bottom left corner, excluding first row and last column
    companionMatrix.bottomLeftCorner(in_monicPolynomialCoeff.rows()-1, in_monicPolynomialCoeff.rows()-1)  = Eigen::Matrix<std::complex<float>, Dynamic, Dynamic>::Identity(in_monicPolynomialCoeff.rows()-1, in_monicPolynomialCoeff.rows()-1);
    //fill first row with zeroes
    companionMatrix.topLeftCorner(1, in_monicPolynomialCoeff.rows()) = Eigen::Matrix<std::complex<float>, Dynamic, Dynamic>::Zero(1, in_monicPolynomialCoeff.rows());
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
    if (this->eigensolver.info() != Success)
        throw eigensolver.info();
}

Eigen::Matrix<std::complex<float>, Dynamic, Dynamic> polynomialHelper::getEigenVectors() {
    if (this->eigensolver.info() != Success)
        throw eigensolver.info();
    else
        return this->eigensolver.eigenvectors();
}

Eigen::Matrix<std::complex<float>, Dynamic, Dynamic> polynomialHelper::getEigenValues() {
    if (this->eigensolver.info() != Success)
        throw eigensolver.info();
    else
    return this->eigensolver.eigenvalues();
}

std::set<eigenMatrix, polynomialHelper::matrixCompare>::iterator polynomialHelper::begin() {
    return this->set_.begin();
}

std::set<eigenMatrix, polynomialHelper::matrixCompare>::iterator polynomialHelper::end() {
    return this->set_.end();
}

polynomialHelper mySet;
linearRegression myMap;
polynomialHelper PolynomialCoeffVectors;

int main(int argc, const char * argv[])
{
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    
    eigenMatrix A(2,2);
    A(0,0) = 3;
    A(1,0) = 2.5;
    A(0,1) = -1;
    A(1,1) = A(1,0) + A(0,1);
    
    eigenMatrix B(2,2);
    //B << 1, 2, 3, 4;
    
    eigenVector v1(2);
    v1(0) = 3;
    v1(1) = -4;
    
    eigenVector v2(2);
    v2 << 1, 2;
    
    eigenVector v3(7);
    v3 << 1, 4, -5, 2, 2.7, 15, 28.1;
    
    eigenVector v4(100);
    v4 = eigenVector::Random(100);
    
    std::chrono::high_resolution_clock::time_point timePoint1 = std::chrono::high_resolution_clock::now();
    
    std::cout << "Extracting the roots of a polynomial" << std::endl;
    PolynomialCoeffVectors.set_.insert(v1);
    PolynomialCoeffVectors.set_.insert(v2);
    PolynomialCoeffVectors.set_.insert(v3);
    PolynomialCoeffVectors.set_.insert(v4);
    std::cout << "Number of equations to solve: " << PolynomialCoeffVectors.set_.size() << std::endl;
    for (PolynomialCoeffVectors.it = PolynomialCoeffVectors.begin(); PolynomialCoeffVectors.it != PolynomialCoeffVectors.end(); PolynomialCoeffVectors.it++) {
        PolynomialCoeffVectors.printPolynome(*PolynomialCoeffVectors.it);
        std::cout << "The associated companion matrix:\n" << PolynomialCoeffVectors.buildCompanionMatrix(*PolynomialCoeffVectors.it) << std::endl;
    PolynomialCoeffVectors.computeEigenValues(PolynomialCoeffVectors.buildCompanionMatrix(*PolynomialCoeffVectors.it));
        std::cout << "The roots of the Polynomial are:\n" <<  PolynomialCoeffVectors.getEigenValues() << std::endl;
        std::cout << "Here's a matrix whose columns are eigenvectors of the matrix corresponding to these eigenvalues:\n" << PolynomialCoeffVectors.getEigenVectors() << std::endl;
    }
    std::cout << std::endl;
    
    std::chrono::high_resolution_clock::time_point timePoint2 = std::chrono::high_resolution_clock::now();
    
    std::cout << "Playing around with a set of matrices" << std::endl;
    mySet.set_.insert(A);
    mySet.set_.insert(B);
    std::cout << "Set size is: " << mySet.set_.size() << std::endl;
    for (mySet.it = mySet.begin(); mySet.it != mySet.end(); mySet.it++) {
        std::cout << "Set matrix: \n" << *mySet.it << std::endl;
        mySet.computeEigenValues(*mySet.it);
        std::cout << "The eigenvalues of the matrix are:\n" << mySet.getEigenValues() << std::endl;
        std::cout << "Here's a matrix whose columns are eigenvectors of the matrix corresponding to these eigenvalues:\n" << mySet.getEigenVectors() << std::endl;
        std::cout << "The determinant of the matrix is " << mySet.it->determinant() << std::endl;
        std::cout << "The inverse of the matrix is:\n" << mySet.it->inverse() << std::endl;
    }
    std::cout << std::endl;
    
    std::chrono::high_resolution_clock::time_point timePoint3 = std::chrono::high_resolution_clock::now();
    
    std::cout << "Playing around with a map of <string, std::pair<Eigen::MatrixXf, Eigen::VectorXf>>" << std::endl;
    myMap.map_.insert(std::make_pair("vector1", std::make_pair(A, v1)));
    myMap.map_.insert(std::make_pair("vector2", std::make_pair(B, v2)));
    std::cout << "Map size is: "  << myMap.map_.size() << std::endl;
    for ( myMap.it = myMap.begin(); myMap.it != myMap.end(); ++myMap.it) {
        std::cout << "Map key: "<< myMap.it->first << std::endl;
        std::cout << "First element in pair, the matrix: " << std::endl << myMap.it->second.first << std::endl;
        std::cout << "Second element in pair, the vector: "<< std::endl << myMap.it->second.second << std::endl;
        std::cout << "In a Ax=b equation, x=A.inverse*b and has solutions if det(A) <> 0" << std::endl;
        std::cout << "The determinant of the matrix is " << myMap.it->second.first.determinant() << std::endl;
        if ( myMap.it->second.first.determinant() !=0 /*|| myMap.it->second.first.determinant().imag() != 0*/) {
            std::cout << "The solution is the vector x: " << myMap.it->second.first.inverse()*myMap.it->second.second << std::endl;
        }
        else std::cout << "There are no solutions" << std::endl;
    }
    std::cout << std::endl;
    
    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
    std::cout << "Initialisation took "
    << std::chrono::duration_cast<std::chrono::microseconds>(timePoint1 - start).count()
    << "us." << std::endl;
    std::cout << "Polynomial root extraction took "
    << std::chrono::duration_cast<std::chrono::microseconds>(timePoint2 - timePoint1).count()
    << "us." << std::endl;
    std::cout << "Matrices in set took "
    << std::chrono::duration_cast<std::chrono::microseconds>(timePoint3 - timePoint2).count()
    << "us." << std::endl;
    std::cout << "Matrices in map took "
    << std::chrono::duration_cast<std::chrono::microseconds>(end - timePoint3).count()
    << "us." << std::endl;
}
