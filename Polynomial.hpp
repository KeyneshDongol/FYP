//
//  Polynomial.hpp
//  Polynoimial Operator
//
//  Created by Keynesh on 15/2/23.
//

#ifndef Polynomial_h
#define Polynomial_h

//#define EIGEN_USE_BLAS
//#define EIGEN_USE_LAPACKE //TODO: Try implemting Intel's MLK library extension and compare it with Eigen's BLAS implementation
//#define EIGEN_USE_MKL_ALL

#include <iostream>
#include <Eigen/Dense>
#include <cassert>

template<int ndim> class Polynomial;

template<> class Polynomial<1> {
// Polynomial of degree 4 of order 1
public:
    template <typename Derived>  Polynomial(const Eigen::DenseBase<Derived>& coeffpass): coeff(coeffpass) { // CoeffPass should be a 5x1 Eigen Vector
        assert(coeffpass.rows() == 5);  assert( coeffpass.cols() == 1);
    }
    explicit Polynomial(int degree){
        assert(degree>=0 & degree <5);
        coeff = Eigen::Matrix <double, 5, 1>::Zero();
        coeff(degree) = 1.;
    }
    
    //********************************
    //* Arithmetic
    //********************************
    Polynomial<1>& operator+=(const Polynomial<1>& other){ coeff += other.coeff; return *this; }
    Polynomial<1>& operator-=(const Polynomial<1>& other){ coeff -= other.coeff; return *this; }
    Polynomial<1>& operator+=(const double scalar){ coeff(0) += scalar; return *this; }
    Polynomial<1>& operator-=(const double scalar){ coeff(0) -= scalar; return *this; }
    Polynomial<1>& operator*=(const double scalar){ coeff *= scalar; return *this; }
    Polynomial<1>& operator/=(const double scalar){ coeff /= scalar; return *this; }

    friend std::ostream &operator<<(std::ostream &t_os, Polynomial<1> const& t_poly);
    
// =============================
// IMPLEMENTATION
// =============================
private:
    // The polynomial is expressed as
    // 0th order a0, 1th order a1, etc
    // a0, a1, etc are scalars
//Data
    Eigen::Matrix<double, 5, 1> coeff; // Polyn coefficients of incresing degree
//    double b0, b1, b2, b3, b4;

public:

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    friend Polynomial<1> operator+(Polynomial<1> lhs, const Polynomial<1>& rhs){ lhs += rhs; return lhs;}
    friend Polynomial<1> operator-(Polynomial<1> lhs, const Polynomial<1>& rhs){ lhs -= rhs; return lhs;}
    friend Polynomial<1> operator+(Polynomial<1> lhs, const double rhs){ lhs += rhs; return lhs;}
    friend Polynomial<1> operator+(const double lhs, Polynomial<1> rhs){ rhs += lhs; return rhs;}
    friend Polynomial<1> operator-(Polynomial<1> lhs, const double rhs){ lhs -= rhs; return lhs;}
    friend Polynomial<1> operator-(const double lhs, Polynomial<1> rhs){ rhs -= lhs; rhs *=  -1.; return rhs;}
    friend Polynomial<1> operator*(Polynomial<1> lhs, const double rhs){ lhs *= rhs; return lhs;}
    friend Polynomial<1> operator*(const double lhs, Polynomial<1> rhs){ rhs *= lhs; return rhs;}
    friend Polynomial<1> operator/(Polynomial<1> lhs, const double rhs){ lhs /= rhs; return lhs;}
    
    double evalPoint (double point) const {
        return coeff(0) +
                coeff(1) * point +
                coeff(2) * point * point +
                coeff(3) * point * point * point +
                coeff(4) * point * point * point * point;
    }

    template <typename Derived> auto evalAll(const Eigen::DenseBase<Derived>& points) const {
        assert(points.rows() == 1);
        
        return (coeff(0) +
                coeff(1) * points.array() +
                coeff(2) * points.array()* points.array() +
                coeff(3) * points.array()* points.array()* points.array() +
                coeff(4) * points.array()* points.array()* points.array()* points.array()).matrix();
    }

private:
    //Algorithm to find real roots from companion matrix
    template <typename Derived> Eigen::Matrix<double,1,Eigen::Dynamic> EigenRealRoots(const Eigen::MatrixBase<Derived>& companionMatrix){
        auto roots = companionMatrix.eigenvalues().eval();
        Eigen::Matrix<double, 1, Eigen::Dynamic> realRoots;
        realRoots.resize(Eigen::NoChange, (roots.imag().array() == 0).template cast<int>().sum());
        int currentRootIndex = 0;
        for(int i=0; i< roots.size(); ++i){
            if( roots(i).imag() == 0 ){
                realRoots(currentRootIndex++) = roots(i).real();
            }
        }
        return realRoots;
    }

public:
    //Contruction and solving companion matrix
    Eigen::Matrix<double, 1, Eigen::Dynamic> findRealRoots() {
        if (coeff(4) != 0) {
            Eigen::Matrix<double, 4, 4> companionMatrix (Eigen::Matrix<double,4,4>::Zero());
            companionMatrix.block<3,3>(1,0).setIdentity();
            companionMatrix.col(3) = - coeff.head(4) / coeff(4);
            
//            companionMatrix(0,3) = -(b0)/(b4); //Try with other test cases
//            companionMatrix(1,3) = -(b1)/(b4);
//            companionMatrix(2,3) = -(b2)/(b4);
//            companionMatrix(3,3) = -(b3)/(b4);
            
            return EigenRealRoots(companionMatrix);

        } else if ( coeff(3) != 0 ) {
            Eigen::Matrix<double, 3, 3> companionMatrix (Eigen::Matrix<double,3,3>::Zero());
            companionMatrix.block<2,2>(1,0).setIdentity();
            companionMatrix.col(2) = - coeff.head(3) / coeff(3);
            
//            companionMatrix(0,2) = -(b0)/(b3);
//            companionMatrix(1,2) = -(b1)/(b3);
//            companionMatrix(2,2) = -(b2)/(b3);
            
            return EigenRealRoots(companionMatrix);

        }else if ( coeff(2) != 0) {
            
            Eigen::Matrix<double, 2, 2> companionMatrix (Eigen::Matrix<double,2,2>::Zero());
            companionMatrix.block<2,2>(1,0).setIdentity();
            companionMatrix.col(1) = - coeff.head(2) / coeff(2);
            
//            companionMatrix(0,2) = -(poly.b0)/(poly.b2);
//            companionMatrix(1,2) = -(poly.b1)/(poly.b2);
            
            return EigenRealRoots(companionMatrix);

        } else if ( coeff(1) != 0) {
            Eigen::Matrix<double,1,Eigen::Dynamic> realRoot; //Discuss with prof to improve this
            realRoot.resize(Eigen::NoChange, 1);
            realRoot(0) = -(coeff(0))/(coeff(1));
//            Eigen::Matrix<double,1,Eigen::Dynamic> roots = {{-(coeff(0))/(coeff(1))}};
//            return Eigen::Matrix<double,1,Eigen::Dynamic>({{-(coeff(0))/(coeff(1))}});
            return realRoot;
        }
        return Eigen::Matrix<double,1,Eigen::Dynamic>();
    }
};

 std::ostream &operator<<(std::ostream &t_os, Polynomial<1> const& t_poly){
     t_os << t_poly.coeff(0) ;
     if (t_poly.coeff(1)!= 0) {
         if (t_poly.coeff(1) > 0) {
             t_os << " + " << t_poly.coeff(1) << " x ";
         } else {
             t_os << " - " << -t_poly.coeff(1) << " x ";
         }
     }
        
     return t_os;
 }//TODO: Print up the rest of the polynomial terms





template<> class Polynomial<2> {//TODO: Operator overload for 2D case plus write out a better way to cout the monomials
// Polynomial of degree 4 of order 2

public:
   template <typename Derived> explicit Polynomial(const Eigen::DenseBase<Derived>& coeffpass); // CoeffPass should be a 15x1 Eigen Vector
   explicit Polynomial(int degreeX, int degreeY);


//    double evalPoint(Eigen::Matrix <double, 2, 1> point) const;
//    Eigen::Matrix <double, 1, Eigen::Dynamic> evalAll(Eigen::Matrix <double, 2, Eigen::Dynamic> points)

   explicit Polynomial(Eigen::Matrix <double, 15, 1> coeffpass) {
       a0 = coeffpass(0);
       a1 = coeffpass.block<2,1>(1,0);
       a2 = Eigen::Matrix <double, 2, 2>::Zero();
       a2(0,0) = coeffpass(3,0);
       a2.block<2,1>(0,1) = coeffpass.block<2,1>(4,0);
       a3.block<2,1>(0,0) = coeffpass.block<2,1>(6,0);
       a3.block<2,1>(0,1) = coeffpass.block<2,1>(8,0);
       a4.block<2,1>(0,0) = coeffpass.block<2,1>(10,0);
       a4.block<2,1>(0,1) = coeffpass.block<2,1>(12,0);
       a4m = coeffpass(14);


   }

   //product of row values
   template <typename Derived> Eigen::Matrix <double, 1, Eigen::Dynamic> prod(const Eigen::DenseBase<Derived>& vector) {// vector should be a 2xD Eigen Matrix

       assert(vector.rows() == 2);//TODO: delete assert upon confirming it works

       Eigen::Matrix <double, 1, Eigen::Dynamic> temp;
       temp.resize(Eigen::NoChange, vector.cols());

       temp.topRows(1) = (vector.topRows(1).array() * vector.bottomRows(1).array()).matrix();

       return temp;
   }
   double evalPoint(Eigen::Matrix <double, 2, 1> point) {//TODO: Try expressing this function in terms of an <typename Derived> template
    return a0 +
            a1 * point +
            point.transpose() * a2 * point +
            ((point.array() *point.array()).matrix()).transpose() * a3 * point +
            ((point.array()*point.array()*point.array()).matrix()).transpose() * a4 * point + a4m * (point.prod())*(point.prod());
}

private:

   // The polynomial is expressed as
   // 0th order     a0 +
   // 1st order     a1 * P +
   // 2nd order     P.transpose * a2 * P
   // 3rd order     P^2.transpose * a3 * P
   // 4th order     P^3.transpose * a4 * P + a4m * (x * y)^2

   // a0 is scaler
   // a1 is [1,Ndim]
   // a2 is [Ndim, Ndim] Upper Triangular
   // a3 is [Ndim, Ndim]
   // a4 is [Ndim,Ndim]
   // a4m is scalar
//Data

   double                          a0;
   Eigen::Matrix <double, 1, 2>    a1;
   Eigen::Matrix <double, 2, 2>    a2;
   Eigen::Matrix <double, 2, 2>    a3;
   Eigen::Matrix <double, 2, 2>    a4;
   double                          a4m;

public:

   EIGEN_MAKE_ALIGNED_OPERATOR_NEW

   //values are respective points
       double zeroth(){
           return a0;
       } //TODO: Figure out what exactly tmemplate<typename Derived> does
    
    template <typename Derived> Eigen::Matrix <double, 1, Eigen::Dynamic> first(const Eigen::MatrixBase<Derived>& points){//points should be a 2xD Eigen Matrix

        assert(points.rows() == 2);//TODO: delete assert upon confirming it works

        return a1 * points;
       }
    template <typename Derived> Eigen::Matrix <double, 1, Eigen::Dynamic> second(const Eigen::MatrixBase<Derived>& points) {//points should be a 2xD Eigen Matrix

        assert(points.rows() == 2);//TODO: delete assert upon confirming it works

           return (points.transpose() * a2 * points).diagonal().array();
       }
    template <typename Derived> Eigen::Matrix <double, 1, Eigen::Dynamic> third(const Eigen::MatrixBase<Derived>& points) {//points should be a 2xD Eigen Matrix

        assert(points.rows() == 2);//TODO: delete assert upon confirming it works

           return (((points.array() * points.array()).matrix()).transpose() * a3 * points).diagonal().array();
       }
    template <typename Derived> Eigen::Matrix <double, 1, Eigen::Dynamic> forth_first(const Eigen::MatrixBase<Derived>& points) {//points should be a 2xD Eigen Matrix

        assert(points.rows() == 2);//TODO: delete assert upon confirming it works

           return (((points.array() *points.array()*points.array()).matrix()).transpose() * a4 * points).diagonal().array();
       }
    template <typename Derived> Eigen::Matrix <double, 1, Eigen::Dynamic> forth_second(const Eigen::MatrixBase<Derived>& points) {//points should be a 2xD Eigen Matrix

        assert(points.rows() == 2);//TODO: delete assert upon confirming it works

           return (a4m * prod((points.array() * points.array()).matrix())).array();
       }

        Eigen::Matrix <double, 1, Eigen::Dynamic> evalAll(Eigen::Matrix <double, 2, Eigen::Dynamic> points) {
            //TODO: Check if you can do the same template <typename Derived> trick here after ensuring the above works
        return (zeroth() +
                first(points).array() +
                second(points).array() +
                third(points).array() +
                forth_first(points).array() +
                forth_second(points).array()).matrix();
    }

};


template<> class Polynomial<3> {
// Polynomial of degree 4 of order 3
// The polynomial is expressed as
// 0th degree   c0 +
// 1st degree   c1 * P +
// 2nd degree   P.transpose * c2 * P +
// 3rd degree   P^2.transpose * c3 * P + c3m * x * y * z
// 4th degree   P^3.transpose * c4 * P + c4m1 * P * x * y * z + c4m2 * (P * RotateLeft(P))^2

// c0 is scalar
// c1 is [1,Ndim]
// c2 is [Ndim,Ndim]  Upper Triangular
// c3 is [Ndim,Ndim]
// c3m is scalar
// c4 is [Ndim,Ndim]
// c4m1 is [1,Ndim]
// c4m2 is [1,Ndim]

public:
//Data
   double                          c0;
   Eigen::Matrix <double, 1, 3>    c1;
   Eigen::Matrix <double, 3, 3>    c2;
   Eigen::Matrix <double, 3, 3>    c3;
   double                          c3m;
   Eigen::Matrix <double, 3, 3>    c4;
   Eigen::Matrix <double, 1, 3>    c4m1;
   Eigen::Matrix <double, 1, 3>    c4m2;


   EIGEN_MAKE_ALIGNED_OPERATOR_NEW


   explicit Polynomial( Eigen::Matrix<double,35,1> coeffpass) {
       c0 = coeffpass(0);
       c1 = coeffpass.block<3,1>(1,0);
       c2 = Eigen::Matrix <double, 3, 3>::Zero();
       c2(0,0) = coeffpass(4,0);
       c2.block<2,1>(0,1) = coeffpass.block<2,1>(5,0);
       c2.block<3,1>(0,2) = coeffpass.block<3,1>(7,0);
       c3.block<3,1>(0,0) = coeffpass.block<3,1>(10,0);
       c3.block<3,1>(0,1) = coeffpass.block<3,1>(13,0);
       c3.block<3,1>(0,2) = coeffpass.block<3,1>(16,0);
       c3m = coeffpass(19);
       c4.block<3,1>(0,0) = coeffpass.block<3,1>(20,0);
       c4.block<3,1>(0,1) = coeffpass.block<3,1>(23,0);
       c4.block<3,1>(0,2) = coeffpass.block<3,1>(26,0);
       c4m1 = coeffpass.block<3,1>(29,0);
       c4m2 = coeffpass.block<3,1>(32,0);
   }


   //rotate left
   template <typename Derived> Eigen::Matrix <double, 3, Eigen::Dynamic> leftRotated(const Eigen::DenseBase<Derived>& vector) {//vector should be a 3xD Eigen Matrix

       assert(vector.rows() == 3);


       Eigen::Matrix <double, 3, Eigen::Dynamic> new_;
       new_.resize(Eigen::NoChange, vector.cols());

       new_.topRows(2) = (vector.topRows(2).array() * vector.bottomRows(2).array() * vector.topRows(2).array() * vector.bottomRows(2).array()).matrix();
       new_.bottomRows(1) = (vector.bottomRows(1).array() * vector.topRows(1).array() * vector.bottomRows(1).array() * vector.topRows(1).array()).matrix();

       return new_;
   }

   //product of row values
   static Eigen::Matrix <double, 1, Eigen::Dynamic> prod(Eigen::Matrix <double, 3, Eigen::Dynamic> vector) {


       Eigen::Matrix <double, 1, Eigen::Dynamic> temp;
       temp.resize(Eigen::NoChange, vector.cols());

       temp.topRows(1) = (vector.topRows(1).array() * vector.bottomRows(1).array() * vector.row(1).array()).matrix();

       return temp;
   }

   //product of row values to return scalar
   static double prodScalar(Eigen::Matrix <double, 3, 1> vector) {

       double newDle = vector(0)*vector(1)*vector(2);


       return newDle;
   }


   //value at specific point
   double evalPoint(Eigen::Matrix <double, 3, 1> point) {
       return c0 +
               c1 * point +
               point.transpose() * c2 * point +
               (point.array() * point.array()).matrix().transpose() * c3 * point + c3m * point.prod() +
               (point.array() * point.array() * point.array()).matrix().transpose() * c4 * point + c4m1 * (point.array() * (point.array()).prod()).matrix()
                   + (c4m2 * leftRotated(point))(0);
   }


   //values at respective points
       double zeroth() {
           return c0;
       }
       Eigen::Array <double, 1, Eigen::Dynamic> first(Eigen::Matrix <double, 3, Eigen::Dynamic> points) {
           return c1 * points;
       }
       Eigen::Array <double, 1, Eigen::Dynamic> second(Eigen::Matrix <double, 3, Eigen::Dynamic> points) {
           return (points.transpose() * c2 * points).diagonal().array();
       }
       Eigen::Array <double, 1, Eigen::Dynamic> third_first(Eigen::Matrix<double, 3, Eigen::Dynamic> points) {
           return ((points.array() * points.array()).matrix().transpose() * c3 * points).diagonal().transpose().array();
       }
       Eigen::Array <double, 1, Eigen::Dynamic> third_second(Eigen::Matrix <double, 3, Eigen::Dynamic> points) {
           return (c3m * prod(points)).array();
       }
       Eigen::Array <double, 1, Eigen::Dynamic> forth_first(Eigen::Matrix <double, 3, Eigen::Dynamic> points) {
           return ((points.array() *points.array() *points.array()).matrix().transpose() * c4 * points).diagonal().transpose().array();
       }
       Eigen::Array <double, 1, Eigen::Dynamic> forth_second(Eigen::Matrix <double, 3, Eigen::Dynamic> points) {
           return (c4m1 * points).array() * prod((points)).array();
       }
       Eigen::Array <double, 1, Eigen::Dynamic> forth_third(Eigen::Matrix <double, 3, Eigen::Dynamic> points) {
           return (c4m2 * leftRotated(points)).array();
       }

   Eigen::Matrix <double, 1, Eigen::Dynamic> evalAll(Eigen::Matrix<double, 3, Eigen::Dynamic> points) {
       return (zeroth() +
               first(points) +
               second(points) +
               third_first(points) +
               third_second(points) +
               forth_first(points) +
               forth_second(points) +
               forth_third(points)).matrix();
   }

//    Polynomial<1> reductionOperator(Eigen::Matrix<double,3,1> A_,Eigen::Matrix<double,3,1> a_){
//        auto                            A = A_(0);
//        auto                            B = A_(1);
//        auto                            C = A_(2);
//        auto                            a = a_(0);
//        auto                            b = a_(1);
//        auto                            c = a_(2);
//        auto                            ALR = leftRotated(A_);
//        auto                            aLR = leftRotated(a_);
//        auto                            Atwo = (A_.array()*A_.array()).matrix();
//        auto                            atwo = (a_.array()*a_.array()).matrix();
//        auto                            Athree = (A_.array()*A_.array()*A_.array()).matrix();
//        auto                            athree = (a_.array()*a_.array()*a_.array()).matrix();
//        auto                            aLRtwo = (aLR.array()*aLR.array()).matrix();
//        auto                            ALRtwo = (ALR.array()*ALR.array()).matrix();
//
//
//        // Eigen expression -> double
//        auto w_0 = (c0
//                        + c1*a_
//                        + a_.transpose()*c2*a_
//                        + atwo.transpose()*c3*a_
//                        + c3m*(a*b*c)
//                        + athree.transpose()*c4*a_
//                        + c4m1*(a_*(a*b*c))
//                        + c4m2*(atwo.array()*aLRtwo.array()).matrix());
//
//        auto w_1 = (c1*A_
//                    + A_.transpose()*c*a_
//                    + a_.transpose()*c2*A_
//                    + aLRtwo.transpose()*c3*A_
//                    + (2*A_.array()*a_.array()).matrix().transpose()*c3*a_
//                    + (3*A_.array()*aLRtwo.array()).matrix().transpose()*c4*a_
//                    + (3*athree).transpose()*c4*A_
//                    + c4m1*(A_*a*b*c) + c4m1*(a_*a*b*C) + c4m1*(a_*a*c*B) + c4m1*(A_*b*c)
//                    + c4m2*(a_.array()*aLRtwo.array()*A_.array()).matrix().value()).value()
//                        + (c3m * (a*b*C + a*c*B + b*c*A));
//
//        auto w_2 = (A_.transpose()*c2*A_
//                                    + (2*A_.array()*a_.array()).matrix().transpose()*c3*A_
//                                    + Atwo.transpose()*c3*a_
//                                    + 3*(Atwo.array()*a_.array()).matrix().transpose()*c4*a_
//                                    + 3*(A_.array()*atwo.array()).matrix().transpose()*c4*A_
//                                    + c4m1*(A_(a*b*c)) + c4m1*(A_(a*c*B)) + c4m1*(A_(b*c*A)) + c4m1*(a_(a*B*C)) + c4m1*(a_(b*A*C)) + c4m1*(a_(c*A*B)) + c4m1*(a_*(A*B*C))
//                                    + c4m2*(aLRtwo.array()*ALRtwo.array()).matrix()
//                                    + c4m2*(atwo.array()*ALRtwo.array()).matrix()
//                                    + c4m2*(4*a_.array()*aLR.array()*A_.array()*ALR.array()).matrix()).value()
//                                        + c3m*(a*B*C + b*A*C + c*A*B);
//
//        auto w_3 = Atwo.transpose()*c3*A_
//                            + c3m*A*B*C
//                            + Athree.transpose()*c4*a_
//                            + 3*(Atwo.array()*a_.array()).matrix().transpose()*c4*A_
//                            + c4m1*(A_(a*B*C)) + c4m1*(A_(b*A*C)) + c4m1*(A_*(c*A*B))
//                            + c4m2*(2*aLR.array()*Atwo.array()*ALR.array()).matrix()
//                            + c4m2*(2*A_.array()*Atwo.array()).matrix();
//
//        auto w_4 = Athree.transpose()*c4*A_
//                                        + (c4m1*(A_*(A*B*C))).value()
//                                        + c4m2*(Atwo.array()*ALRtwo.array()).matrix();
//
//        return Polynomial<1>(Eigen::Matrix <double, 5, 1>{w_0 , w_1 , w_2 , w_3 , w_4});
//    }


};








#endif /* Polynomial_h */
