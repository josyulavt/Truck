#include <iostream> 
#include <vector>
#include <cmath>
#include <utility>
// #include "matplotlibcpp.h"

using namespace std; 
// namespace plt = matplotlibcpp;

class MotionModel{
    public: 
        float theta0; 
        float theta1; 
        float delta;
        float v; 

        MotionModel(float theta0, float theta1, float delta){
            this->theta0 = theta0; 
            this->theta1 = theta1; 
            this->delta = delta; 
        };

        void kinematics(float v, float delta, float theta0, float theta1){
            float d1 = 1.0;
            float x1dot = v*cos(theta0-theta1)*cos(theta1);
            float y1dot = v*cos(theta0-theta1)*sin(theta1);
            float theta0dot = (v/d1)*tan(delta); 
            float theta1dot = (v/d1)*sin(theta0-theta1); 

        };

        vector<double> calculateAlpha(double x_A, double x_B, double eta_1, double eta_2, double eta_3, double eta_4, double eta_5,   double eta_6, double eta_7, double eta_8,  double theta_A, double theta_B, double kappa_A, double kappa_B, double dot_kappa_A, double dd_kappa_A, double dot_kappa_B, double dd_kappa_B) {
            double alpha_0 = x_A;
            double alpha_1 = eta_1 * cos(theta_A);
            double alpha_2 = (1.0 / 2.0) * eta_3 * cos(theta_A) - (1.0 / 2.0) * pow(eta_1, 2) * kappa_A * sin(theta_A);
            double alpha_3 = -((1.0 / 2.0) * eta_1 * eta_3 * kappa_A + (1.0 / 6.0) * pow(eta_1, 3) * dot_kappa_A) * sin(theta_A) + (1.0 / 6.0) * eta_5 * cos(theta_A);
            double alpha_4 = -((1.0 / 6.0) * eta_1 * eta_5 * kappa_A + (1.0 / 4.0) * pow(eta_1, 2) * eta_3 * dot_kappa_A + (1.0 / 8.0) * pow(eta_1, 4) * pow(kappa_A, 3) + (1.0 / 24.0) * pow(eta_1, 4) * dd_kappa_A) * sin(theta_A) - (1.0 / 8.0) * pow(eta_3, 2) * kappa_A * sin(theta_A) + (1.0 / 24.0) * eta_7 * cos(theta_A);
            double alpha_5 = 126.0 * (x_B - x_A) - (70.0 * eta_1 + (35.0 / 2.0) * eta_3 + (5.0 / 2.0) * eta_5 + (5.0 / 24.0) * eta_7) * cos(theta_A) + ((35.0 / 2.0) * pow(eta_1, 2) * kappa_A + 
                            (5.0 / 8.0) * pow(eta_1, 4) * pow(kappa_A, 3) + (5.0 / 2.0) * pow(eta_1, 3) * dot_kappa_A + (5.0 / 24.0) * pow(eta_1, 4) * dd_kappa_A) * sin(theta_A) +
                             ((15.0 / 2.0) * eta_1 * eta_3 * kappa_A + (5.0 / 6.0) * eta_1 * eta_5 * kappa_A + (5.0 / 4.0) * pow(eta_1, 2) * eta_3 * dot_kappa_A + (5.0 / 8.0) * pow(eta_3, 2) * kappa_A) * sin(theta_A) 
                             - (56.0 * eta_2 - (21.0 / 2.0) * eta_4 + eta_6 - (1.0 / 24.0) * eta_8) * cos(theta_B) - ((21.0 / 2.0) * pow(eta_2, 2) * kappa_B + (1.0 / 8.0) * pow(eta_2, 4) * pow(kappa_B, 3) 
                             - eta_2 * pow(kappa_B, 3) + (1.0 / 24.0) * pow(eta_2, 4) * dd_kappa_B) * sin(theta_B) + (3.0 * eta_2 * eta_4 * kappa_B - (1.0 / 6.0) * eta_2 * eta_6 * kappa_B - (1.0 / 4.0) * pow(eta_2, 2) * eta_4 * dot_kappa_B 
                             - (1.0 / 8.0) * pow(eta_4, 2) * kappa_B) * sin(theta_B);
            double alpha_6 = -420.0 * (x_B - x_A) +
                            (224.0 * eta_1 + (105.0 / 2.0) * eta_3 + (20.0 / 3.0) * eta_5 + (5.0 / 12.0) * eta_7) * cos(theta_A) -
                            ((105.0 / 2.0) * pow(eta_1, 2) * kappa_A + (5.0 / 4.0) * pow(eta_1, 4) * pow(kappa_A, 3) + (20.0 / 3.0) * pow(eta_1, 3) * dot_kappa_A + (5.0 / 12.0) * pow(eta_1, 4) * pow(kappa_A, 4)) * sin(theta_A) -
                            (20.0 * eta_1 * eta_3 * kappa_A + (5.0 / 3.0) * eta_1 * eta_5 * kappa_A + (5.0 / 2.0) * pow(eta_1, 2) * eta_3 * dot_kappa_A + (5.0 / 4.0) * pow(eta_3, 2) * kappa_A) * sin(theta_A) +
                            (196.0 * eta_2 - (77.0 / 2.0) * eta_4 + (23.0 / 6.0) * eta_6 - (1.0 / 6.0) * eta_8) * cos(theta_B) +
                            ((77.0 / 2.0) * pow(eta_2, 2) * kappa_B + (1.0 / 2.0) * pow(eta_2, 4) * pow(kappa_B, 3) - (23.0 / 6.0) * pow(eta_2, 3) * dot_kappa_B + (1.0 / 6.0) * pow(eta_2, 4) * dd_kappa_B) * sin(theta_B) -
                            ((23.0 / 2.0) * eta_2 * eta_4 * kappa_B - (2.0 / 3.0) * eta_2 * eta_6 * kappa_B - pow(eta_2, 2) * eta_4 * dot_kappa_B - (1.0 / 2.0) * pow(eta_4, 2) * kappa_B);

            double alpha_7 = 540.0 * (x_B - x_A) -
                            (280.0 * eta_1 + 63.0 * eta_3 + (15.0 / 2.0) * eta_5 + (5.0 / 12.0) * eta_7) * cos(theta_A) +
                            (63.0 * pow(eta_1, 2) * kappa_A + (5.0 / 4.0) * pow(eta_1, 4) * pow(kappa_A, 3) + (15.0 / 2.0) * pow(eta_1, 3) * dot_kappa_A + (5.0 / 12.0) * pow(eta_1, 4) * dot_kappa_A) * sin(theta_A) +
                            ((45.0 / 2.0) * eta_1 * eta_3 * kappa_A + (5.0 / 3.0) * eta_1 * eta_5 * kappa_A + (5.0 / 2.0) * pow(eta_1, 2) * eta_3 * dot_kappa_A + (5.0 / 4.0) * pow(eta_3, 2) * kappa_A) * sin(theta_A) -
                            (260.0 * eta_2 - 53.0 * eta_4 + (11.0 / 2.0) * eta_6 - (1.0 / 4.0) * eta_8) * cos(theta_B) -
                            (53.0 * pow(eta_2, 2) * kappa_B + (3.0 / 4.0) * pow(eta_2, 4) * pow(kappa_B, 3) - (11.0 / 2.0) * pow(eta_2, 3) * dot_kappa_B + (1.0 / 4.0) * pow(eta_2, 4) * dd_kappa_B) * sin(theta_B) +
                            ((33.0 / 2.0) * eta_2 * eta_4 * kappa_B - eta_2 * eta_6 * kappa_B - (3.0 / 2.0) * pow(eta_2, 2) * eta_4 * dot_kappa_B - (3.0 / 4.0) * pow(eta_4, 2) * kappa_B);

            double alpha_8 = -315.0 * (x_B - x_A) +
                            (160.0 * eta_1 + 35.0 * eta_3 + 4.0 * eta_5 + (5.0 / 24.0) * eta_7) * cos(theta_A) -
                            (35.0 * pow(eta_1, 2) * kappa_A + (5.0 / 8.0) * pow(eta_1, 4) * pow(kappa_A, 3));

            vector<double> alpha = {alpha_0, alpha_1, alpha_2, alpha_3, alpha_4, alpha_5, alpha_6, alpha_7, alpha_8};

            return alpha;
        };

        vector<double> calculateBeta(double y_A, double y_B,  double eta_1, double eta_2, double eta_3, double eta_4, double eta_5, double eta_6, double eta_7, double eta_8,  double theta_A, double theta_B, double kappa_A, double kappa_B, double dot_kappa_A, double dd_kappa_A, double dot_kappa_B, double dd_kappa_B){

            double beta_0 = y_A;
            double beta_1 = eta_1 * sin(theta_A);
            double beta_2 = (1.0 / 2.0) * eta_3 * sin(theta_A) + (1.0 / 2.0) * pow(eta_1, 2) * kappa_A * cos(theta_A);
            double beta_3 = ((1.0 / 2.0) * eta_1 * eta_3 * kappa_A + (1.0 / 6.0) * pow(eta_1, 3) * dot_kappa_A) * cos(theta_A) + (1.0 / 6.0) * eta_5 * sin(theta_A);
            
            double beta_4 = ((1.0/6.0)*eta_1*eta_5*kappa_A+(1.0/4.0)*pow(eta_1,2)*eta_3*dot_kappa_A+(1.0/8.0)*pow(eta_1,4)*pow(kappa_A,3)+(1.0/24.0)*pow(eta_1,4)*dd_kappa_A)*cos(theta_A)+(1.0/8.0)*pow(eta_7,2)*kappa_A*cos(theta_A)+(1.0/24.0)*eta_7*sin(theta_A);
            
            double beta_5 = 126.0*(y_B-y_A)-(70.0*eta_1+(35.0/2.0)*eta_3+(5.0/2.0)*eta_5+(5.0/24.0)*eta_7)*sin(theta_A)-((35.0/2.0)*pow(eta_1,2)*kappa_A+(5.0/8.0)*pow(eta_1,4)*pow(kappa_A,3)+(5.0/2.0)*pow(eta_1,3)*dot_kappa_A+(5.0/24.0)*pow(eta_1,4)*dd_kappa_A)*cos(theta_A)-((15.0/2.0)*eta_1*eta_3*kappa_A+(5.0/6.0)*eta_1*eta_5*kappa_A+(5.0/4.0)*pow(eta_1,2)*eta_3*dot_kappa_A+(5.0/8.0)*pow(eta_3,2)*kappa_A)*cos(theta_A)-(56.0*eta_2-(21.0/2.0)*eta_4+eta_6-(1.0/24.0)*eta_8)*sin(theta_B)+((21.0/2.0)*pow(eta_2,2)*kappa_A+(1.0/8.0)*pow(eta_2,4)*pow(kappa_B,3)-pow(eta_2,3)*dot_kappa_B+(1.0/24.0)*pow(eta_2,4)*dd_kappa_B)*cos(theta_B)-(3.0*eta_2*eta_4*kappa_B-(1.0/6.0)*eta_2*eta_6*kappa_B-(1.0/4.0)*pow(eta_2,2)*eta_4*dot_kappa_B-(1.0/8.0)*pow(eta_4,2)*kappa_B)*cos(theta_B);

            double beta_6 = -420.0*(y_B-y_A)+(224.0*eta_1+(105.0/2.0)*eta_3+(20.0/3.0)*eta_5+(5.0/12.0)*eta_7)*sin(theta_A)+((105.0/2.0)*pow(eta_1,2)*kappa_A+(5.0/4.0)*pow(eta_1,4)*pow(kappa_A,3)+(20.0/3.0)*pow(eta_1,3)*dot_kappa_A+(5.0/12.0)*pow(eta_1,4)*dd_kappa_A)*cos(theta_A)+(20.0*eta_1*eta_3*kappa_A+(5.0/3.0)*eta_1*eta_5*kappa_A+(5.0/2.0)*pow(eta_1,2)*eta_3*dot_kappa_A+(5.0/4.0)*pow(eta_3,2)*kappa_A)*cos(theta_A)+(196.0*eta_2-(77.0/2.0)*eta_4+(23.0/6.0)*eta_6-(1.0/6.0)*eta_8)*sin(theta_B)-((77.0/2.0)*pow(eta_2,2)*kappa_B+(1.0/2.0)*pow(eta_2,4)*pow(kappa_B,3)-(23.0/6.0)*pow(eta_2,3)*dot_kappa_B+(1.0/6.0)*pow(eta_2,4)*dd_kappa_B)*cos(theta_B)+((23.0/2.0)*eta_2*eta_4*kappa_B-(2.0/3.0)*eta_2*eta_6*kappa_B-pow(eta_2,2)*eta_4*dot_kappa_B-(1.0/2.0)*pow(eta_4,2)*kappa_B)*cos(theta_B);

            double beta_7 = 540.0*(y_B-y_A)-(280.0*eta_1+63.0*eta_3+(15.0/2.0)*eta_5+(5.0/12.0)*eta_7)*sin(theta_A)-(63.0*pow(eta_1,2)*kappa_A+(5.0/4.0)*pow(eta_1,4)*pow(kappa_A,3)+(15.0/2.0)*pow(eta_1,3)*dot_kappa_A+(5.0/12.0)*pow(eta_1,4)*dd_kappa_A)*cos(theta_A)-((45.0/2.0)*eta_1*eta_3*kappa_A+(5.0/3.0)*eta_1*eta_5*kappa_A+(5.0/2.0)*pow(eta_1,2)*eta_3*dot_kappa_A+(5.0/4.0)*pow(eta_3,2)*kappa_A)*cos(theta_A)-(260.0*eta_2-53.0*eta_4+(11.0/2.0)*eta_6-(1.0/4.0)*eta_8)*sin(theta_B)+(53.0*pow(eta_2,2)*kappa_B+(3.0/4.0)*pow(eta_2,4)*pow(kappa_B,3)-(11.0/2.0)*pow(eta_2,3)*dot_kappa_B+(1.0/4.0)*pow(eta_2,4)*dd_kappa_B)*cos(theta_B)-((33.0/2.0)*eta_2*eta_4*kappa_B-eta_2*eta_6*kappa_B-(3.0/2.0)*pow(eta_2,2)*eta_4*dot_kappa_B-(3.0/4.0)*pow(eta_4,2)*kappa_B);

            double beta_8 = -315.0*(y_B-y_A)+ (160.0*eta_1+35.0*eta_3+4.0*eta_5+(5.0/24.0)*eta_7)*sin(theta_A)+ (35.0*pow(eta_1,2)*kappa_A+(5.0/8.0)*pow(eta_1,4)*pow(kappa_A,3)+ 4.0*pow(eta_1,3)*dot_kappa_A+(5.0/24.0)*pow(eta_1,4)*dd_kappa_A)*cos(theta_A)+ (12.0*eta_1*eta_3*kappa_A+(5.0/6.0)*eta_1*eta_5*kappa_A+ (5.0/4.0)*pow(eta_1,2)*eta_3*dot_kappa_A+(5.0/8.0)*pow(eta_3,2)*kappa_A)*cos(theta_A)+ (155.0*eta_2-(65.0/2.0)*eta_4+(7.0/2.0)*eta_6-(1.0/6.0)*eta_8)*sin(theta_B)- ((65.0/2.0)*pow(eta_2,2)*kappa_B+(1.0/2.0)*pow(eta_2,4)*pow(kappa_B,3)- (7.0/2.0)*pow(eta_2,3)*dot_kappa_B+(1.0/6.0)*pow(eta_2,4)*dd_kappa_B)*cos(theta_B)+ ((21.0/2.0)*eta_2*eta_4*kappa_B-(2.0/3.0)*eta_2*eta_6*kappa_B- pow(eta_2,2)*eta_4*dot_kappa_B-(1.0/2.0)*pow(eta_4,2)*kappa_B);

            double beta_9 = 70.0*(y_B-y_A)-(35.0*eta_1+(15.0/2.0)*eta_3+(5.0/6.0)*eta_5+(1.0/24.0)*eta_7)*sin(theta_A)-((15.0/2.0)*pow(eta_1,2)*kappa_A+(1.0/8.0)*pow(eta_1,4)*pow(kappa_A,3)+(5.0/6.0)*pow(eta_1,3)*dot_kappa_A+(1.0/24.0)*pow(eta_1,4)*dd_kappa_A)*cos(theta_A)-((5.0/2.0)*eta_1*eta_3*kappa_A+(1.0/6.0)*eta_1*eta_5*kappa_A+(1.0/4.0)*pow(eta_1,2)*eta_3*dot_kappa_A+(1.0/8.0)*pow(eta_3,2)*kappa_A)*cos(theta_A)-(35.0*eta_2-(15.0/2.0)*eta_4+(5.0/6.0)*eta_6-(1.0/24.0)*eta_8)*sin(theta_B)+((15.0/2.0)*pow(eta_2,2)*kappa_B+(1.0/8.0)*pow(eta_2,4)*pow(kappa_B,3)-(5.0/6.0)*pow(eta_2,3)*dot_kappa_B+(1.0/24.0)*pow(eta_2,4)*dd_kappa_B)*cos(theta_B)-((5.0/2.0)*eta_2*eta_4*kappa_B-(1.0/6.0)*eta_2*eta_6*kappa_B-(1.0/4.0)*pow(eta_2,2)*eta_4*dot_kappa_B-(1.0/8.0)*pow(eta_4,2)*kappa_B)*cos(theta_B);
            
            vector<double> beta = {beta_0, beta_1, beta_2, beta_3, beta_4, beta_5, beta_6, beta_7, beta_8, beta_9};
            return beta;

        };

        double calX(double u, vector<double> alpha){
            double x = 0;
            for(int i = 0; i < 10; i++){
                x += alpha[i] * pow(u, i);
            }
            return x;
        };

        double calY(double u, vector<double> beta){
            double y = 0;
            for(int i = 0; i < 10; i++){
                y += beta[i] * pow(u, i);
            }
            return y;
        };
};


template<typename T>
std::vector<T> arange(T start, T stop, T step = 1) {
    std::vector<T> values;
    for (T value = start; value < stop; value += step)
        values.push_back(value);
    return values;
}


int main(int argc, char* argv[]){
    MotionModel kin(0,0,0);
    
    double x_A = 0.; 
    double x_B = 0.;

    double y_A = 5; 
    double y_B = 4;

    double eta_1 = 1; 
    double eta_2 = 1; 
    double eta_3 = 1;
    double eta_4 = 1;
    double eta_5 = 1; 
    double eta_6 = 1;
    double eta_7 = 1; 
    double eta_8 = 1; 

    double theta_A = 1;
    double theta_B = 1;

    double kappa_A = 1;
    double kappa_B = 1;
    double dot_kappa_A = 1;
    double dd_kappa_A = 1;
    double dot_kappa_B = 1;
    double dd_kappa_B = 1;

    //test trajectory generation 
    vector<double> alpha = kin.calculateAlpha(x_A, x_B, eta_1, eta_2, eta_3, eta_4, eta_5, eta_6, eta_7, eta_8, theta_A, theta_B, kappa_A, kappa_B, dot_kappa_A, dd_kappa_A, dot_kappa_B, dd_kappa_B);
    vector<double> beta = kin.calculateBeta(y_A, y_B, eta_1, eta_2, eta_3, eta_4, eta_5, eta_6, eta_7, eta_8, theta_A, theta_B, kappa_A, kappa_B, dot_kappa_A, dd_kappa_A, dot_kappa_B, dd_kappa_B);

    auto u = arange(0., 1. ,0.1);

    vector<double> x,y; 
    std::vector<std::pair<double,double>> data;

    for(double i: u){
        x.push_back(kin.calX(i, alpha));
        y.push_back(kin.calY(i, beta));
        // double x1 = kin.calX(i, alpha);
        // double y1 = kin.calY(i, beta);
        // data.emplace_back(x1, y1);
    }
    cout<<"x"<<endl;
    for(double x1: x){
        cout<<x1<<endl;
    }
    cout<<"y"<<endl;
    for(double x1: y){
        cout<<x1<<endl;
    }
    // Gnuplot g1("lines");
    // g1 << "plot data'\n";
    // g1.send1d(data);
    // plt::plot(x,y);
    // plt::show();
    return 0;

}
