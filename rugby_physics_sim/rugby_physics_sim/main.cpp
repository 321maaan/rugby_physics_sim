
#include <Eigen/Dense>
#include <cmath>
#include <fstream>
#include <iostream>

using Vec3 = Eigen::Vector3d;
using Mat3 = Eigen::Matrix3d;

double Get_rotation_kinetic_energy(const Vec3 &omega,const Mat3 &I)
{
    // 0.5 * w*Iw
    double e = 0.5*(omega.dot(I*omega));
    return e;
}
double Get_center_of_mass_energy(const double &m, const Vec3 &v)
{
    return 0.5 * m * (v.dot(v));
}

Vec3 Get_omega_derivative(const Mat3 &Iinv, const  Vec3 &omega, const Mat3 &I, const Vec3 &tau)
{
    // calculate base on the formula:
    // domega/dt = I(-1) * (-omega(x)I*omega + N(tau))

    Vec3 Iomega = I * omega;
    Vec3 insideSquareBrackts = -omega.cross(Iomega) + tau;
    return Iinv * insideSquareBrackts;

}

Vec3 Get_position_by_time(const double t, const Vec3& x0, const Vec3& v0, const Vec3 &acceleration)
{
    //calculate by the formula
    return x0 + t * v0 + 0.5 * acceleration * (t * t);
}

inline void WriteDataToFile(std::ofstream& out, const double time, const Vec3& pos, const Vec3& omega, const Vec3& omega_derivative, const double &re, const double &cme)
{
    out << time << ','
        << pos.x() << ',' << pos.y() << ',' << pos.z() << ','
        << omega.x() << ',' << omega.y() << ',' << omega.z() << ','
        << omega_derivative.x() << ',' << omega_derivative.y() << ',' << omega_derivative.z() << ','
        << re << ',' << cme
        << '\n';
}

int main()
{

    const double m = 0.44;       // kg
    const double a = 0.14;       // חצי אורך, m
    const double b = 0.095;      // חצי קוטר מקסימלי, m
    const double I3 = 0.4 * m * b * b;            // סביב ציר האורך
    const double I1 = 0.2 * m * (a * a + b * b);  // סביב ציר רוחבי
    const double I2 = I1;
    Mat3 I = Vec3(I1, I2, I3).asDiagonal(); // make I matrix
    std::cout << "I1: " << I1 << ". I2: " << I2 << ". I3: " << I3 << std::endl;
    Mat3 Iinv = I.inverse(); // the inverse of I
    const double g = 9.81;



    // start condition:
    Vec3 x(0, 0, 1.0), v(20, 0, 15), acceleration(0,0,-g);
    Vec3 w(2.0, 0.5, 30.0);          // ω_body, ספין בעיקר סביב ציר האורך + קצת הטיה
    Mat3 R = Mat3::Identity();       // אוריינטציה התחלתית (תחליף לפי הצורך)
    Mat3 Rinv = R.inverse();
    Vec3 tau = Vec3::Zero(); // the force/moment on the ball (N)


    const double dt = 1e-3;
    const double T = 3; // 3 seconds

    std::ofstream out("sim.csv");
    out << "t,x,y,z,wx,wy,wz,dwx,dwy,dwz,re,cme,\n";

    // the simulator loop
    Vec3 lastOmega = w;
    Vec3 lastPosition = x;
    Vec3 lastVelocity = v;

    int N = static_cast<int>(T / dt);
    std::cout << "Starting sim!" << std::endl;
    for (int i = 0; i < N; i++)
    {
        double t = i * dt;

        // derivative at the current state
        Vec3 OmegaDerivative = Get_omega_derivative(Iinv, lastOmega, I, tau);

        // Calculate the energy
        double RotationEnergy = Get_rotation_kinetic_energy(lastOmega, I);
        double CenterOfMassEnergy = Get_center_of_mass_energy(m, lastVelocity);

        // write the current state (row 0 = initial conditions)
        WriteDataToFile(out, t, lastPosition, lastOmega, OmegaDerivative, RotationEnergy, CenterOfMassEnergy);

        // advance the state by one step
        Vec3 currentOmega = lastOmega + OmegaDerivative * dt;
        Vec3 currentAcceleration = acceleration; // for now no change
        Vec3 currentVelocity = lastVelocity + currentAcceleration * dt;
        Vec3 currentPosition = lastPosition + currentVelocity * dt;

        lastOmega = currentOmega;
        lastVelocity = currentVelocity;
        lastPosition = currentPosition;
    }
    std::cout << "Finished sim!" << std::endl << "File named sim.csv has been created";

	return 0;
}