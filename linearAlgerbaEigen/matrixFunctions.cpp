//
//  matrixFunctions.cpp
//  linearAlgerbaEigen
//
//  Created by Mihai Dunareanu on 10.12.2017.
//  Copyright © 2017 Mihai Dunareanu. All rights reserved.
//

#include "matrixFunctions.hpp"

bool matrixCompare::operator() (const eigenMatrix& a, const eigenMatrix& b) const {
    if (a.rows() == b.rows() && a.cols() == b.cols())
        return not(a.isApprox(b));
    else return true;
}

linearRegression::linearRegression() {
    std::cout << "Calling Constructor" << std::endl;
}

linearRegression::~linearRegression() {
    std::cout << "Calling Destructor" << std::endl;
}

std::map<std::string, MyPair>::iterator linearRegression::begin() {
    return this->map_.begin();
}

std::map<std::string, MyPair>::iterator linearRegression::end() {
    return this->map_.end();
}

long unsigned int linearRegression::getSize() {
    return this->map_.size();
}

polynomialHelper::polynomialHelper () {
    std::cout << "Calling Constructor" << std::endl;
}

polynomialHelper::~polynomialHelper () {
    std::cout << "Calling Destructor" << std::endl;
}

Eigen::Matrix<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic> polynomialHelper::getEigenVectors() {
    if (this->eigensolver.info() != Eigen::Success)
        throw eigensolver.info();
    else
        return this->eigensolver.eigenvectors();
}

Eigen::Matrix<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic> polynomialHelper::getEigenValues() {
    if (this->eigensolver.info() != Eigen::Success)
        throw eigensolver.info();
    else
        return this->eigensolver.eigenvalues();
}

std::set<eigenMatrix, matrixCompare>::iterator polynomialHelper::begin() {
    return this->set_.begin();
}

std::set<eigenMatrix, matrixCompare>::iterator polynomialHelper::end() {
    return this->set_.end();
}

long unsigned int polynomialHelper::getSize() {
    return this->set_.size();
}
