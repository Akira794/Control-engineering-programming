/* RungeKutta test */
#include <Eigen/Dense>
#include <Eigen/Core>
#include <iostream>
#include <fstream>

using namespace Eigen;
using namespace std;

void RungeKutta(MatrixXd dX, MatrixXd &x, MatrixXd u, double tt, double dt, MatrixXd A, MatrixXd b, MatrixXd c, MatrixXd d)
{
	MatrixXd k1 = A*x + b*u;
	MatrixXd k2 = A*(x + 0.5*k1*dt) + b*u;
	MatrixXd k3 = A*(x + 0.5*k2*dt) + b*u;
	MatrixXd k4 = A*(x + k3*dt) + b*u;
	MatrixXd k = (k1 + 2.0*k2 + 2.0*k3 + k4)*dt / 6.0;
	x = x + k;
}

int main(int argc, char* argv[])
{
	double k = 40.0;
	double m = 250.0;
	double c = 60.0;

	double dt = 0.01;
	double tt = 0.0;

	MatrixXd A(3,3);
	A << 1.0, -1.0,  1.0,
		 1.0, -k/m, -c/m,
		 1.0,  1.0,  1.0;
	
	MatrixXd B(3,1);
	B << 0.0, -1.0/m, 1.0;

	MatrixXd C(1,3);
	C << 1.0, 0.0, 1.0;

	MatrixXd D(1,3);
	D << 1.0, 0.0, 0.0;

	MatrixXd K(1,3); /* FeedBack Gain */
	K << 15.80864382, -3.63298165, 7.85453193;

	MatrixXd X(3,1);
	X << 10.0, 0.0, 0.0;

	MatrixXd dX(3,1);
	dX << 0.0, 0.0, 0.0;

	MatrixXd u(1,1);
	u << 0.0;

	MatrixXd Y(1,1);
	Y << 0.0;

	ofstream ofs("outdata.csv");
	for(int i = 0; i < 1000; i++)
	{
		RungeKutta(dX, X, u, tt, dt, A, B, C, D);
		Y = C * X;
		u = -K*X;

		ofs << tt << "," << Y(0,0) << endl;
		tt += dt;
	}
	
	return 0;
}
