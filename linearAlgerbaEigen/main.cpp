//
//  main.cpp
//  vectorPracticeEigen
//
//  Created by Mihai Dunareanu on 30.11.2017.
//  Copyright © 2017 Mihai Dunareanu. All rights reserved.
//

#include "matrixFunctions.hpp"
#include <chrono>

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
    B << 1, 2, 3, 4;
    
    eigenVector v1(2);
    v1(0) = 3;
    v1(1) = -4;
    
    eigenVector v2(2);
    v2 << 1, 2;
    
    eigenVector v3(7);
    v3 << 1, 4, -5, 2, 2.7, 15, 28.1;
    
    eigenVector v4(10);
    v4 = eigenVector::Random(10);
    
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
