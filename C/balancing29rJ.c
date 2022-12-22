////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////           ---- NOTES ----                  //////////////////////////////
//1. -- Make returntoinitialconfiguration function to return gimbals to initial-zero configuration
//2. change h0 in h matrix to J*w
//3. check the axiserror sign. propably is correct but check it anyway
//4. -- FIX printing the matrices/values after finishing the manoeuvre / maybe insert delay in prints
////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////                                 START IMU                     ///////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
#include <Wire.h>
#include <SPI.h>
#include <SparkFunLSM9DS1.h>
#include <BasicLinearAlgebra.h>
using namespace BLA;
LSM9DS1 imu;
#define PRINT_RAW
#define PRINT_SPEED 50 //50 // 250 ms between prints
static unsigned long lastPrint = 0; // Keep track of print time
#define DECLINATION 4.56
void printGyro();
void printAccel();
void printMag();
int ledpin = 13;
double rollc;
double pitchc;
double headingc;
double headingbias = 0;
double rollcprev = 0;
double pitchcprev = 0;
double headingcprev = 0;
double aroll, apitch, ahead, factor;
double headingcfprev = 0;
double gxprev = 0;
double gyprev = 0;
double gzprev = 0;
double gxglobal = 0;
double gyglobal = 0;
double gzglobal = 0;
double axprev = 0;
double ayprev = 0;
double azprev = 0;
int print_attitude = 0;
int print_gyro = 0;
int print_mag = 0;
Matrix<3, 1> RPY;
Matrix<3, 1> RPYprev;


Matrix<3, 1> printAttitude(float ax, float ay, float az, float gx, float gy, float gz, float mx, float my, float mz)
{
    mx = mx -  (5368.7);
    my = my - (-7339.5);
    mz = mz - (-3015.9);
    float M11 =   1.000;
    float M12 =  0.0875;
    float M13 = -0.0135;
    float M21 =  0.0875;
    float M22 =  0.8884;
    float M23 =  0.0264;
    float M31 = -0.0135;
    float M32 =  0.0264;
    float M33 =  0.9382;
    mx = mx * M11 + my * M21 + mz * M31;
    my = mx * M12 + my * M22 + mz * M32;
    mz = mx * M13 + my * M23 + mz * M33;


    float  fact_gyro = 0.3;
    gx = (1 - fact_gyro) * gxprev + fact_gyro * gx;
    gy = (1 - fact_gyro) * gyprev + fact_gyro * gy;
    gz = (1 - fact_gyro) * gzprev + fact_gyro * gz;
    gxprev = gx;
    gyprev = gy;
    gzprev = gz;
    gxglobal = gx;
    gyglobal = gy;
    gzglobal = gz;


    float fact_acc = 0.3;
    ax = (1 - fact_acc) * axprev + fact_acc * ax;
    ay = (1 - fact_acc) * ayprev + fact_acc * ay;
    az = (1 - fact_acc) * azprev + fact_acc * az;
    axprev = ax;
    ayprev = ay;
    azprev = az;


    if (print_mag == 1)
    {
        Serial.print("M     ");
        Serial.print(mx);
        Serial.print("      ,      ");
        Serial.print(my);
        Serial.print("      ,      ");
        Serial.println(mz);
    }



    float roll = atan2(ay, az);
    float pitch = atan2(-ax, sqrt(ay * ay + az * az));


    float heading;
    // Convert everything from radians to degrees:
    //heading *= 180.0 / PI;
    pitch *= 180.0 / PI;
    roll  *= 180.0 / PI;


    aroll = 0.4;
    apitch = 0.4;
    factor = 0.1;
    rollc = (1 - aroll) * (rollcprev + factor * 0.001 * PRINT_SPEED * imu.calcGyro(gx)) + aroll * roll;
    pitchc = (1 - apitch) * (pitchcprev + factor * 0.001 * PRINT_SPEED * imu.calcGyro(gy)) + apitch * pitch;

    heading  = atan2( (-my * cos(rollc * PI / 180) + mz * sin(rollc * PI / 180) ), (-mx * cos(pitchc * PI / 180) + my * sin(pitchc * PI / 180) * sin(rollc * PI / 180) + mz * sin(pitchc * PI / 180) * cos(rollc * PI / 180)) );  //the same as aboive but with -mx
    heading *= 180.0 / PI;
    heading -= DECLINATION * PI / 180;
    if (heading > PI) heading -= (2 * PI);
    else if (heading < -PI) heading += (2 * PI);

    ahead = 0.4;
    headingc = (1 - ahead) * (headingcprev + factor * 0.001 * PRINT_SPEED * imu.calcGyro(gz)) + ahead * heading;

    headingc = headingc - headingbias * 180 / PI; //headingbias is in rad, so I convert heading bias from rad to deg

    rollcprev = rollc;
    pitchcprev = pitchc;
    headingcprev = headingc;

    if (print_attitude == 1)
    {
        Serial.print("N   ");
        Serial.print(roll, 2);
        Serial.print("  ,  ");
        Serial.print(pitch, 2);
        Serial.print("  ,  ");
        Serial.println(heading, 2);

        Serial.print("C   ");
        Serial.print(-rollc, 2);
        Serial.print("  ,  ");
        Serial.print(pitch, 2);
        Serial.print("  ,  ");
        Serial.println(headingc, 2);
    }

    if (print_gyro == 1)
    {
        // imu.calcGyro(gx) gives result in deg/s
        Serial.print("G   ");
        Serial.print(imu.calcGyro(gx), 2);
        Serial.print("  ,  ");
        Serial.print(imu.calcGyro(gy), 2);
        Serial.print("  ,  ");
        Serial.println(imu.calcGyro(gz), 2);

    }


    Matrix<3, 1> rpytemp = {-rollc, pitch, headingc};

    heading = heading - headingbias * 180 / PI;
    Matrix<3, 1> rpytempbare = {-roll, pitch, heading};
    return rpytempbare;
}

Matrix<3, 1> getRPY()
{


    if ( imu.gyroAvailable() )
    {
        imu.readGyro();
    }
    if ( imu.accelAvailable())
    {
        imu.readAccel();
    }
    if ( imu.magAvailable() )
    {
        imu.readMag();
    }

    // if ((lastPrint + PRINT_SPEED) < millis())
    // {
    RPY = printAttitude(imu.ax, imu.ay, imu.az, imu.gx, imu.gy, imu.gz, imu.mx, imu.my, imu.mz );
    // lastPrint = millis(); // Update lastPrint time
    // }
    RPY = RPY * (PI / 180); //convert RPY from deg to rad
    // RPY(2)=RPY(2)-headingbias;
    return RPY;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////                                 END IMU                       ///////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
void printBLAmat300x27(Matrix<300, 27> mat)
{
    for (int i = 0; i < mat.GetColCount(); i++)
    {
        for (int j = 0; j < mat.GetRowCount(); j++)
        {
            if (mat(j, i) != 9998.98)
            {
                Serial.print(mat(j, i));
                Serial.print(" , ");
            }
        }
        Serial.println(" ");
        delay(50);
    }
}
void printBLAmat300x9(Matrix<300, 9> mat)
{
    for (int i = 0; i < mat.GetColCount(); i++)
    {
        for (int j = 0; j < mat.GetRowCount(); j++)
        {
            if (mat(j, i) != 9998.98)
            {
                Serial.print(mat(j, i));
                Serial.print(" , ");
            }
        }
        Serial.println(" ");
        delay(50);
    }
}
void printBLAmat300x3(Matrix<300, 3> mat)
{
    for (int i = 0; i < mat.GetColCount(); i++)
    {
        for (int j = 0; j < mat.GetRowCount(); j++)
        {
            if (mat(j, i) != 9998.98)
            {
                Serial.print(mat(j, i), 5);
                Serial.print(" , ");
            }
        }
        Serial.println(" ");
        delay(50);
    }
}
void printBLAmat300x4(Matrix<300, 4> mat)
{
    for (int i = 0; i < mat.GetColCount(); i++)
    {
        for (int j = 0; j < mat.GetRowCount(); j++)
        {
            if (mat(j, i) != 9998.98)
            {
                Serial.print(mat(j, i));
                Serial.print(" , ");
            }
        }
        Serial.println(" ");
        delay(50);
    }
}
void printBLAmat300x1(Matrix<300, 1> mat)
{
    for (int i = 0; i < mat.GetColCount(); i++)
    {
        for (int j = 0; j < mat.GetRowCount(); j++)
        {
            if (mat(j, i) != 9998.98)
            {
                Serial.print(mat(j, i));
                Serial.print(" , ");
            }
        }
        Serial.println(" ");
        delay(50);
    }
}

/////
void printBLAmat3x6(Matrix<3, 6> mat)
{
    for (int i = 0; i < mat.GetColCount(); i++)
    {
        for (int j = 0; j < mat.GetRowCount(); j++)
        {
            Serial.print(mat(j, i));
            Serial.print(" , ");
        }
        Serial.println(" ");
        delay(50);
    }
}
void printBLAmat3x9(Matrix<3, 9> mat)
{
    for (int i = 0; i < mat.GetColCount(); i++)
    {
        for (int j = 0; j < mat.GetRowCount(); j++)
        {
            Serial.print(mat(j, i));
            Serial.print(" , ");
        }
        Serial.println(" ");
        delay(50);
    }
}
void printBLAmat3x4(Matrix<3, 4> mat)
{
    for (int i = 0; i < mat.GetRowCount(); i++)
    {
        for (int j = 0; j < mat.GetColCount(); j++)
        {
            Serial.print(mat(i, j));
            Serial.print(' ');
        }
        Serial.println(' ');
    }
}
void printBLAmat3x1(Matrix<3, 1> mat)
{
    for (int i = 0; i < mat.GetColCount(); i++)
    {
        for (int j = 0; j < mat.GetRowCount(); j++)
        {
            Serial.print(mat(j, i), 5);
            Serial.print(" , ");
        }
        Serial.println(" ");
        delay(50);
    }
}

void printBLAmat4x1(Matrix<4, 1> mat)
{
    for (int i = 0; i < mat.GetRowCount(); i++)
    {
        for (int j = 0; j < mat.GetColCount(); j++)
        {
            Serial.print(mat(i, j), 4);
            Serial.print(' ');  //print(value, accuracy)
        }
        Serial.println(' ');
    }
}

Matrix<3, 3> skew(Matrix<3, 1> mat)
{
    Matrix<3, 3> S = {0, -mat(2), mat(1),
                      mat(2), 0, -mat(0),
                      -mat(1), mat(0), 0
                     };

    return S;
}

Matrix<4, 1> rpy2quat(float eulyaw, float eulpitch, float eulroll)
{
    // float cy = cos(yaw * 0.5);
    // float sy = sin(yaw * 0.5);
    // float cp = cos(pitch * 0.5);
    // float sp = sin(pitch * 0.5);
    // float cr = cos(roll * 0.5);
    // float sr = sin(roll * 0.5);
    // Matrix<4, 1> quaternion;
    // quaternion(0) = cr * cp * cy + sr * sp * sy;
    // quaternion(1) = sr * cp * cy - cr * sp * sy;
    // quaternion(2) = cr * sp * cy + sr * cp * sy;
    // quaternion(3) = cr * cp * sy - sr * sp * cy;
    // return quaternion;

    double c1 = cos(eulyaw * 0.5);
    double s1 = sin(eulyaw * 0.5);
    double c2 = cos(eulpitch * 0.5);
    double s2 = sin(eulpitch * 0.5);
    double c3 = cos(eulroll * 0.5);
    double s3 = sin(eulroll * 0.5);
    Matrix<4, 1> quaternion;
    quaternion(0) = c1 * c2 * c3 + s1 * s2 * s3;
    quaternion(1) = c1 * c2 * s3 - s1 * s2 * c3;
    quaternion(2) = c1 * s2 * c3 + s1 * c2 * s3;
    quaternion(3) = s1 * c2 * c3 - c1 * s2 * s3;

    return quaternion;

}

Matrix<4, 1> quatmultiply(Matrix<4, 1> q1, Matrix<4, 1> q2)
{
    float a = q1(0);
    float b = q1(1);
    float c = q1(2);
    float d = q1(3);

    float x = q2(0);
    float y = q2(1);
    float z = q2(2);
    float w = q2(3);

    Matrix<4, 1> quatres;
    quatres(0) = a * x - b * y - d * w - c * z;
    quatres(1) = a * y + b * x + c * w - d * z;
    quatres(2) = a * z - b * w + c * x + d * y;
    quatres(3) = a * w + b * z - c * y + d * x;
    return quatres;
}

float norm(Matrix<3, 1> mat)
{
    float normres = pow((pow(abs(mat(0)), 2) + pow(abs(mat(1)), 2) + pow(abs(mat(2)), 2)), 0.5);
    return normres;
}

Matrix<4, 1> quatnormalize(Matrix<4, 1> quat)
{
    Matrix<4, 1> quatres;
    float denum = pow((pow(abs(quat(0)), 2) + pow(abs(quat(1)), 2) + pow(abs(quat(2)), 2) + pow(abs(quat(3)), 2)), 0.5);
    quatres(0) = quat(0) / denum;
    quatres(1) = quat(1) / denum;
    quatres(2) = quat(2) / denum;
    quatres(3) = quat(3) / denum;
    return quatres;
}

Matrix<4, 1> quatconj(Matrix<4, 1> quat)
{
    Matrix<4, 1> quatres = {quat(0), -quat(1), -quat(2), -quat(3)};
    return quatres;
}

Matrix<3, 1> cross(Matrix<3, 1> mat1, Matrix<3, 1> mat2)
{
    float a = mat1(0);
    float b = mat1(1);
    float c = mat1(2);
    float x = mat2(0);
    float y = mat2(1);
    float z = mat2(2);
    Matrix<3, 1> crossres = {b *z - c * y, c *x - a * z, a *y - b * x};
    return crossres;
}

double det(Matrix<3, 3> m)
{
    double detres = m(0, 0) * m(1, 1) * m(2, 2) - m(0, 0) * m(1, 2) * m(2, 1) - m(0, 1) * m(1, 0) * m(2, 2) + m(0, 1) * m(1, 2) * m(2, 0) + m(0, 2) * m(1, 0) * m(2, 1) - m(0, 2) * m(1, 1) * m(2, 0);
    return detres;
}



Matrix<3, 3> J = {1, 0, 0, 0, 1, 0, 0, 0, 1};
Matrix<3, 1> w = {0, 0, 0};
Matrix<3, 1> wdot = {0, 0, 0};
Matrix<3, 1> Tc = {0, 0, 0};
Matrix<4, 1> qini = {1, 0, 0, 0};
Matrix<4, 1> q = qini;
Matrix<4, 1> qerror;
Matrix<4, 1> qprev = {1, 0, 0, 0};
Matrix<4, 1> qref = {0.7071, -0.7071, 0, 0};
Matrix<4, 1> qdot;
Matrix<4, 1> wquat;
Matrix<1, 1> slerppartA1;
Matrix<3, 1> slerppartA2;
Matrix<4, 1> slerppartA;
Matrix<3, 1> axiserror;
Matrix<3, 1> axiserrorprev;
Matrix<3, 1> Integral = {0, 0, 0};
Matrix<3, 1> u;
Matrix<3, 4> A;
Matrix<3, 1> hvector;
Matrix<3, 1> Hdot;
Matrix<4, 3> Ahash;
Matrix<3, 3> AAT;
Matrix<4, 1> ddot;
Matrix<4, 1> ddotpwm;
Matrix<300, 3> MATw;
Matrix<300, 1> MATmanip;
Matrix<300, 4> MATddot;
Matrix<300, 4> MATangles;
Matrix<3, 6> cross_sum = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
Matrix<3, 3> gravskewtot = skew({0, 0, -9.81});
Matrix<3, 1> rightterm2tot = {0, 0, 0};
Matrix<3, 1> rightterm;
Matrix<3, 9> leftterm;


Matrix<3, 3> phi;
Matrix<3, 1> dw;
Matrix<300, 9> phimat;
Matrix<300, 3> dwmat;
double fiprev;
double thetaprev ;
double psiprev  ;
double w0prev ;
double w1prev ;
double w2prev ;

float b = 0.9552187; //54.73 deg
float dt = 0.05;
float th;
float dur = 0;
float Kp = 1;
float Ki = 0;
float Kwmega = 0;
double manip;
float d1 = 0, d2 = 0, d3 = 0, d4 = 0;
Matrix<4, 1> gimbalangles;
float h01, h02, h03, h04;

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
//  MOTORS //
//GIMBAL 1
int  gimbal1_PWMpin = 29; //enable
int  gimbal1_DIRpin = 27; //phase
int  fly1_PWMpin = 3;
int  fly1_DIRpin = 1;
int  gimbal1_encA = 25;
int  gimbal1_encB = 26;
int  fly1_encA = 7;
int  fly1_encB = 6;

//GIMBAL 2
int  gimbal2_PWMpin = 35; //enable
int  gimbal2_DIRpin = 33; //phase
int  fly2_PWMpin = 30;
int  fly2_DIRpin = 28;
int  gimbal2_encA = 37;
int  gimbal2_encB = 24;
int  fly2_encA = 32;
int  fly2_encB = 31;

//GIMBAL 3
int  gimbal3_PWMpin = 22; //enable
int  gimbal3_DIRpin = 21; //phase
int  fly3_PWMpin = 36;
int  fly3_DIRpin = 34;
int  gimbal3_encA = 15;
int  gimbal3_encB = 14;
int  fly3_encA = 39;
int  fly3_encB = 38;

//GIMBAL 4
int  gimbal4_PWMpin = 2; //enable
int  gimbal4_DIRpin = 0; //phase
int  fly4_PWMpin = 20;
int  fly4_DIRpin = 23;
int  gimbal4_encA = 4;
int  gimbal4_encB = 5;
int  fly4_encA = 17;
int  fly4_encB = 16;

// ON BOARD PINS
int led = 13;
int photores = A21;
int orangeled = 9;
int greenled = 8;

//VARIABLES
float photores_threshold;
int pwmval_G1 = 70;
int pwmval_G2 = 70;
int pwmval_G3 = 70;
int pwmval_G4 = 70;
int pwmval_fly1 = 60;
int pwmval_fly2 = 60;
int pwmval_fly3 = 60;
int pwmval_fly4 = 60;
int gimbalstall = 60; //means that min pwm for gimbals is 60

//Volatile Counts
volatile int countgimbal1 = 0;
volatile int countgimbal2 = 0;
volatile int countgimbal3 = 0;
volatile int countgimbal4 = 0;
volatile int countfly1 = 0;
volatile int countfly2 = 0;
volatile int countfly3 = 0;
volatile int countfly4 = 0;
volatile int countfly1prev = 0;
volatile int countfly2prev = 0;
volatile int countfly3prev = 0;
volatile int countfly4prev = 0;

//Volatile Angles
volatile double delta1 = 0;
volatile double delta2 = 0;
volatile double delta3 = 0;
volatile double delta4 = 0;
volatile double deltafly1 = 0;
volatile double deltafly2 = 0;
volatile double deltafly3 = 0;
volatile double deltafly4 = 0;

int countprev = 0;
unsigned long previousTime = 0, currentTime = 0; //
double LOOP_FREQUENCY = 1000;
double rpms;

Matrix<4, 1> ddot2pwm()
{

    Matrix<4, 1> ddotpwmval;
    ddotpwmval(0) = 425.89 * ddot(0) - 8.3061;
    ddotpwmval(1) = 425.89 * ddot(1) - 8.3061;
    ddotpwmval(2) = 425.89 * ddot(2) - 8.3061;
    ddotpwmval(3) = 425.89 * ddot(3) - 8.3061;
    return ddotpwmval;
}

void saturateddotpwm()
{
    for (int i = 0; i < 4; ++i)
    {
        ddotpwm(i) = abs(ddotpwm(i));
        ddotpwm(i) = ddotpwm(i) > 255 ? 255 : ddotpwm(i);
        ddotpwm(i) = ddotpwm(i) < gimbalstall ? 0 : ddotpwm(i);
    }

}

int mysign(double number)
{
    if (number >= 0)
    {
        return 1;
    }
    else
    {
        return -1;
    }
}
void saturateddot()
{
    double maxddot = 0.6182;

    if (abs(ddot(0)) > maxddot)
    {
        ddot(0) = mysign(ddot(0)) * maxddot;
    }
    if (abs(ddot(1)) > maxddot)
    {
        ddot(1) = mysign(ddot(1)) * maxddot;
    }
    if (abs(ddot(2)) > maxddot)
    {
        ddot(2) = mysign(ddot(2)) * maxddot;
    }
    if (abs(ddot(3)) > maxddot)
    {
        ddot(3) = mysign(ddot(3)) * maxddot;
    }
}


double gimbalcount2theta(int count)
{
    double angleindeg = count * 360 / (12 * 1006);
    double angleinrad = angleindeg * (PI / 180);
    return angleinrad;
}
double flycount2theta(int count)
{
    double angleindeg = count * 360 / (12 * 5);
    double angleinrad = angleindeg * (PI / 180);
    return angleinrad;
}

//ISR FUNCTIONS

//GIMBAL 1
void ISR_GIMBAL1_A ()
{
    if (digitalRead(gimbal1_encA) == HIGH)
    {
        if (digitalRead(gimbal1_encB) == LOW)    //1
        {
            countgimbal1++;
        }
        else
        {
            countgimbal1--;
        }
    }
    else
    {
        if (digitalRead(gimbal1_encB) == HIGH)    //1
        {
            countgimbal1++;
        }
        else
        {
            countgimbal1--;
        }
    }
    delta1 = gimbalcount2theta(countgimbal1);
}

void ISR_GIMBAL1_B ()
{
    if (digitalRead(gimbal1_encB) == HIGH)
    {
        if (digitalRead(gimbal1_encA) == HIGH)    //1
        {
            countgimbal1++;
        }
        else
        {
            countgimbal1--;
        }
    }
    else
    {
        if (digitalRead(gimbal1_encA) == LOW)    //1
        {
            countgimbal1++;
        }
        else
        {
            countgimbal1--;
        }
    }
    delta1 = gimbalcount2theta(countgimbal1);
}

void ISR_FLY1_A ()
{
    if (digitalRead(fly1_encA) == HIGH)
    {
        if (digitalRead(fly1_encB) == LOW)    //1
        {
            countfly1++;
        }
        else
        {
            countfly1--;
        }
    }
    else
    {
        if (digitalRead(fly1_encB) == HIGH)    //1
        {
            countfly1++;
        }
        else
        {
            countfly1--;
        }
    }
    deltafly1 = flycount2theta(countfly1);
}

void ISR_FLY1_B ()
{
    if (digitalRead(fly1_encB) == HIGH)
    {
        if (digitalRead(fly1_encA) == HIGH)    //1
        {
            countfly1++;
        }
        else
        {
            countfly1--;
        }
    }
    else
    {
        if (digitalRead(fly1_encA) == LOW)    //1
        {
            countfly1++;
        }
        else
        {
            countfly1--;
        }
    }
    deltafly1 = flycount2theta(countfly1);
}


//GIMBAL 2
void ISR_GIMBAL2_A ()
{
    if (digitalRead(gimbal2_encA) == HIGH)
    {
        if (digitalRead(gimbal2_encB) == LOW)    //1
        {
            countgimbal2++;
        }
        else
        {
            countgimbal2--;
        }
    }
    else
    {
        if (digitalRead(gimbal2_encB) == HIGH)    //1
        {
            countgimbal2++;
        }
        else
        {
            countgimbal2--;
        }
    }
    delta2 = gimbalcount2theta(countgimbal2);
}

void ISR_GIMBAL2_B ()
{
    if (digitalRead(gimbal2_encB) == HIGH)
    {
        if (digitalRead(gimbal2_encA) == HIGH)    //1
        {
            countgimbal2++;
        }
        else
        {
            countgimbal2--;
        }
    }
    else
    {
        if (digitalRead(gimbal2_encA) == LOW)    //1
        {
            countgimbal2++;
        }
        else
        {
            countgimbal2--;
        }
    }
    delta2 = gimbalcount2theta(countgimbal2);
}

void ISR_FLY2_A ()
{
    if (digitalRead(fly2_encA) == HIGH)
    {
        if (digitalRead(fly2_encB) == LOW)    //2
        {
            countfly2++;
        }
        else
        {
            countfly2--;
        }
    }
    else
    {
        if (digitalRead(fly2_encB) == HIGH)    //2
        {
            countfly2++;
        }
        else
        {
            countfly2--;
        }
    }
    deltafly2 = flycount2theta(countfly2);
}

void ISR_FLY2_B ()
{
    if (digitalRead(fly2_encB) == HIGH)
    {
        if (digitalRead(fly2_encA) == HIGH)    //2
        {
            countfly2++;
        }
        else
        {
            countfly2--;
        }
    }
    else
    {
        if (digitalRead(fly2_encA) == LOW)    //2
        {
            countfly2++;
        }
        else
        {
            countfly2--;
        }
    }
    deltafly2 = flycount2theta(countfly2);
}

//GIMBAL 3
void ISR_GIMBAL3_A ()
{
    if (digitalRead(gimbal3_encA) == HIGH)
    {
        if (digitalRead(gimbal3_encB) == LOW)    //1
        {
            countgimbal3++;
        }
        else
        {
            countgimbal3--;
        }
    }
    else
    {
        if (digitalRead(gimbal3_encB) == HIGH)    //1
        {
            countgimbal3++;
        }
        else
        {
            countgimbal3--;
        }
    }
    delta3 = gimbalcount2theta(countgimbal3);
}

void ISR_GIMBAL3_B ()
{
    if (digitalRead(gimbal3_encB) == HIGH)
    {
        if (digitalRead(gimbal3_encA) == HIGH)    //1
        {
            countgimbal3++;
        }
        else
        {
            countgimbal3--;
        }
    }
    else
    {
        if (digitalRead(gimbal3_encA) == LOW)    //1
        {
            countgimbal3++;
        }
        else
        {
            countgimbal3--;
        }
    }
    delta3 = gimbalcount2theta(countgimbal3);
}

void ISR_FLY3_A ()
{
    if (digitalRead(fly3_encA) == HIGH)
    {
        if (digitalRead(fly3_encB) == LOW)    //3
        {
            countfly3++;
        }
        else
        {
            countfly3--;
        }
    }
    else
    {
        if (digitalRead(fly3_encB) == HIGH)    //3
        {
            countfly3++;
        }
        else
        {
            countfly3--;
        }
    }
    deltafly3 = flycount2theta(countfly3);
}

void ISR_FLY3_B ()
{
    if (digitalRead(fly3_encB) == HIGH)
    {
        if (digitalRead(fly3_encA) == HIGH)    //3
        {
            countfly3++;
        }
        else
        {
            countfly3--;
        }
    }
    else
    {
        if (digitalRead(fly3_encA) == LOW)    //3
        {
            countfly3++;
        }
        else
        {
            countfly3--;
        }
    }
    deltafly3 = flycount2theta(countfly3);
}

//GIMBAL 4
void ISR_GIMBAL4_A ()
{
    if (digitalRead(gimbal4_encA) == HIGH)
    {
        if (digitalRead(gimbal4_encB) == LOW)    //1
        {
            countgimbal4++;
        }
        else
        {
            countgimbal4--;
        }
    }
    else
    {
        if (digitalRead(gimbal4_encB) == HIGH)    //1
        {
            countgimbal4++;
        }
        else
        {
            countgimbal4--;
        }
    }
    delta4 = gimbalcount2theta(countgimbal4);
}

void ISR_GIMBAL4_B ()
{
    if (digitalRead(gimbal4_encB) == HIGH)
    {
        if (digitalRead(gimbal4_encA) == HIGH)    //1
        {
            countgimbal4++;
        }
        else
        {
            countgimbal4--;
        }
    }
    else
    {
        if (digitalRead(gimbal4_encA) == LOW)    //1
        {
            countgimbal4++;
        }
        else
        {
            countgimbal4--;
        }
    }
    delta4 = gimbalcount2theta(countgimbal4);
}

void ISR_FLY4_A ()
{
    if (digitalRead(fly4_encA) == HIGH)
    {
        if (digitalRead(fly4_encB) == LOW)    //4
        {
            countfly4++;
        }
        else
        {
            countfly4--;
        }
    }
    else
    {
        if (digitalRead(fly4_encB) == HIGH)    //4
        {
            countfly4++;
        }
        else
        {
            countfly4--;
        }
    }
    deltafly4 = flycount2theta(countfly4);
}

void ISR_FLY4_B ()
{
    if (digitalRead(fly4_encB) == HIGH)
    {
        if (digitalRead(fly4_encA) == HIGH)    //4
        {
            countfly4++;
        }
        else
        {
            countfly4--;
        }
    }
    else
    {
        if (digitalRead(fly4_encA) == LOW)    //4
        {
            countfly4++;
        }
        else
        {
            countfly4--;
        }
    }
    deltafly4 = flycount2theta(countfly4);
}

void setgimbaldirections()
{
    if (ddotpwm(0) > 0)
    {
        digitalWrite(gimbal1_DIRpin, LOW);
    }
    else
    {
        digitalWrite(gimbal1_DIRpin, HIGH);
    }
    if (ddotpwm(1) > 0)
    {
        digitalWrite(gimbal2_DIRpin, LOW);
    }
    else
    {
        digitalWrite(gimbal2_DIRpin, HIGH);
    }
    if (ddotpwm(2) > 0)
    {
        digitalWrite(gimbal3_DIRpin, LOW);
    }
    else
    {
        digitalWrite(gimbal3_DIRpin, HIGH);
    }
    if (ddotpwm(3) > 0)
    {
        digitalWrite(gimbal4_DIRpin, LOW);
    }
    else
    {
        digitalWrite(gimbal4_DIRpin, HIGH);
    }
}

/////////////// VOID - PRINT FUNCTIONS ////////////////////////
void printddot()
{
    Serial.print("D   ");
    Serial.print(ddot(0), 3);
    Serial.print("  ,  ");
    Serial.print(ddot(1), 3);
    Serial.print("  ,  ");
    Serial.print(ddot(2), 3);
    Serial.print("  ,  ");
    Serial.println(ddot(3), 3);
}

void printmanipulability()
{
    Serial.print("M   ");
    Serial.println(manip, 3);

}
void printgimbalangles()
{
    Serial.print("G   ");
    Serial.print(gimbalangles(0), 3);
    Serial.print("  ,  ");
    Serial.print(gimbalangles(1), 3);
    Serial.print("  ,  ");
    Serial.print(gimbalangles(2), 3);
    Serial.print("  ,  ");
    Serial.println(gimbalangles(3), 3);
}
void printaxiserror()
{
    Serial.print("E   ");
    Serial.print(axiserror(0), 3);
    Serial.print("  ,  ");
    Serial.print(axiserror(1), 3);
    Serial.print("  ,  ");
    Serial.println(axiserror(2), 3);
}
void printquaternion(Matrix<4, 1> p)
{
    Serial.print("Q   ");
    Serial.print(p(0), 3);
    Serial.print("  ,  ");
    Serial.print(p(1), 3);
    Serial.print("  ,  ");
    Serial.print(p(2), 3);
    Serial.print("  ,  ");
    Serial.println(p(3), 3);
    // Serial.println("  ;  ");
}
void printMATs()
{
    Serial.println("W");
    printBLAmat300x3(MATw);
    Serial.println("M");
    printBLAmat300x1(MATmanip);
    Serial.println("D");
    printBLAmat300x4(MATddot);
    Serial.println("A");
    printBLAmat300x4(MATangles);
    Serial.println("END");
}

void printbalanceMATs()
{
    Serial.println("P");
    printBLAmat300x9(phimat);
    Serial.println("D");
    printBLAmat300x3(dwmat);
    Serial.println("END");
}
//////////////// END ISR FUNCTIONS /////////////////////////////
void myblink()
{
    for (int i = 0; i < 4; ++i)
    {
        digitalWrite(orangeled, HIGH);
        digitalWrite(greenled, LOW);
        delay(500);
        digitalWrite(orangeled, LOW);
        digitalWrite(greenled, HIGH);
        delay(500);
    }
    digitalWrite(orangeled, HIGH);
    digitalWrite(greenled, LOW);

}
void resetcounts()
{
    countgimbal1 = 0;
    countgimbal2 = 0;
    countgimbal3 = 0;
    countgimbal4 = 0;
    countfly1 = 0;
    countfly2 = 0;
    countfly3 = 0;
    countfly4 = 0;
}
void returntoinitialconfiguration()
{
    double acc = 0.3;
    double refangle = 0;
    int speed = 200;
    double temp1 = fmod(abs(delta1 - refangle), (2 * PI));
    double temp2 = fmod(abs(delta2 - refangle), (2 * PI));
    double temp3 = fmod(abs(delta3 - refangle), (2 * PI));
    double temp4 = fmod(abs(delta4 - refangle), (2 * PI));
    while ( temp1 > acc)
    {
        digitalWrite(gimbal1_DIRpin, LOW);
        analogWrite(gimbal1_PWMpin, speed * temp1);
        temp1 = fmod(abs(delta1 - refangle), (2 * PI));

    }
    analogWrite(gimbal1_PWMpin, 0);

    while ( temp2 > acc)
    {
        digitalWrite(gimbal2_DIRpin, LOW);
        analogWrite(gimbal2_PWMpin, speed * temp2);
        temp2 = fmod(abs(delta2 - refangle), (2 * PI));

    }
    analogWrite(gimbal2_PWMpin, 0);

    while ( temp3 > acc)
    {
        digitalWrite(gimbal3_DIRpin, LOW);
        analogWrite(gimbal3_PWMpin, speed * temp3);
        temp3 = fmod(abs(delta3 - refangle), (2 * PI));

    }
    analogWrite(gimbal3_PWMpin, 0);

    while ( temp4 > acc)
    {
        digitalWrite(gimbal4_DIRpin, LOW);
        analogWrite(gimbal4_PWMpin, speed * temp4);
        temp4 = fmod(abs(delta4 - refangle), (2 * PI));

    }
    analogWrite(gimbal4_PWMpin, 0);
}

void spinupflywheels()
{
    digitalWrite(fly1_DIRpin, HIGH);
    digitalWrite(fly2_DIRpin, HIGH);
    digitalWrite(fly3_DIRpin, HIGH);
    digitalWrite(fly4_DIRpin, HIGH);

    int flystall = 60;
    for (int i = flystall; i < 255; i = i + 20)
    {
        analogWrite(fly1_PWMpin, i);
        analogWrite(fly2_PWMpin, i);
        analogWrite(fly3_PWMpin, i);
        analogWrite(fly4_PWMpin, i);
        delay(100);
    }
}
void initializegimbalangles()
{
    digitalWrite(greenled, HIGH);
    digitalWrite(orangeled, LOW);

    do
    {

        //INITIALIZE CMG 1
        if ((countfly1 - countfly1prev) > 0)
        {
            digitalWrite(gimbal1_DIRpin, HIGH);
            analogWrite(gimbal1_PWMpin, 70 * (countfly1 - countfly1prev));
        }
        else if ((countfly1 - countfly1prev) < 0)
        {
            digitalWrite(gimbal1_DIRpin, LOW);
            analogWrite(gimbal1_PWMpin, 70 * abs(countfly1 - countfly1prev));
        }
        else
        {
            analogWrite(gimbal1_PWMpin, 0);
        }
        countfly1prev = countfly1;

        //INITIALIZE CMG 2
        if ((countfly2 - countfly2prev) > 0)
        {
            digitalWrite(gimbal2_DIRpin, HIGH);
            analogWrite(gimbal2_PWMpin, 70 * (countfly2 - countfly2prev));
        }
        else if ((countfly2 - countfly2prev) < 0)
        {
            digitalWrite(gimbal2_DIRpin, LOW);
            analogWrite(gimbal2_PWMpin, 70 * abs(countfly2 - countfly2prev));
        }
        else
        {
            analogWrite(gimbal2_PWMpin, 0);
        }
        countfly2prev = countfly2;

        //INITIALIZE CMG 3
        if ((countfly3 - countfly3prev) > 0)
        {
            digitalWrite(gimbal3_DIRpin, HIGH);
            analogWrite(gimbal3_PWMpin, 70 * (countfly3 - countfly3prev));
        }
        else if ((countfly3 - countfly3prev) < 0)
        {
            digitalWrite(gimbal3_DIRpin, LOW);
            analogWrite(gimbal3_PWMpin, 70 * abs(countfly3 - countfly3prev));
        }
        else
        {
            analogWrite(gimbal3_PWMpin, 0);
        }
        countfly3prev = countfly3;


        //INITIALIZE CMG 4
        if ((countfly4 - countfly4prev) > 0)
        {
            digitalWrite(gimbal4_DIRpin, HIGH);
            analogWrite(gimbal4_PWMpin, 70 * (countfly4 - countfly4prev));
        }
        else if ((countfly4 - countfly4prev) < 0)
        {
            digitalWrite(gimbal4_DIRpin, LOW);
            analogWrite(gimbal4_PWMpin, 70 * abs(countfly4 - countfly4prev));
        }
        else
        {
            analogWrite(gimbal4_PWMpin, 0);
        }
        countfly4prev = countfly4;



        delay(50);

    }
    while (photoresistortouched() == 0);
    digitalWrite(greenled, HIGH);
    digitalWrite(orangeled, HIGH);
    delay(2000);


}

void calibrate_photoresistor()
{
    int val = 0;
    int samples = 200;
    for (int i = 0; i < samples; ++i)
    {
        val = val + analogRead(photores);
    }
    val = val / samples;
    photores_threshold = val / 3;
}


int photoresistortouched()
{
    // Serial.print(analogRead(photores)); Serial.print("   ");
    // Serial.println(photores_threshold);
    if (analogRead(photores) < photores_threshold)
    {
        digitalWrite(orangeled, LOW);
        digitalWrite(greenled, HIGH);
        return 1;
    }
    else
    {
        return 0;
    }
}

void stopeverything()
{
    digitalWrite(orangeled, HIGH);
    digitalWrite(greenled, LOW);
    // STOP 1st CMG
    do
    {
        ddotpwm(0) = ddotpwm(0) - 10;
        analogWrite(gimbal1_PWMpin, ddotpwm(0));
        delay(25);
    }
    while (ddotpwm(0) > 60);
    ddotpwm(0) = 0;
    analogWrite(gimbal1_PWMpin, ddotpwm(0));

    do
    {
        pwmval_fly1 = pwmval_fly1 - 10;
        analogWrite(fly1_PWMpin, pwmval_fly1);
        delay(50);
    }
    while (pwmval_fly1 > 60);
    pwmval_fly1 = 0;
    analogWrite(fly1_PWMpin, pwmval_fly1);

    // STOP 2nd CMG
    do
    {
        ddotpwm(1) = ddotpwm(1) - 10;
        analogWrite(gimbal2_PWMpin, ddotpwm(1));
        delay(25);
    }
    while (ddotpwm(1) > 60);
    ddotpwm(1) = 0;
    analogWrite(gimbal2_PWMpin, ddotpwm(1));

    do
    {
        pwmval_fly2 = pwmval_fly2 - 10;
        analogWrite(fly2_PWMpin, pwmval_fly2);
        delay(50);
    }
    while (pwmval_fly2 > 60);
    pwmval_fly2 = 0;
    analogWrite(fly2_PWMpin, pwmval_fly2);

    // STOP 3rd CMG
    do
    {
        ddotpwm(2) = ddotpwm(2) - 10;
        analogWrite(gimbal3_PWMpin, ddotpwm(2));
        delay(25);
    }
    while (ddotpwm(2) > 60);
    ddotpwm(2) = 0;
    analogWrite(gimbal3_PWMpin, ddotpwm(2));

    do
    {
        pwmval_fly3 = pwmval_fly3 - 10;
        analogWrite(fly3_PWMpin, pwmval_fly3);
        delay(50);
    }
    while (pwmval_fly3 > 60);
    pwmval_fly3 = 0;
    analogWrite(fly3_PWMpin, pwmval_fly3);

    // STOP 4th CMG
    do
    {
        ddotpwm(3) = ddotpwm(3) - 10;
        analogWrite(gimbal4_PWMpin, ddotpwm(3));
        delay(25);
    }
    while (ddotpwm(3) > 60);
    ddotpwm(3) = 0;
    analogWrite(gimbal4_PWMpin, ddotpwm(3));

    do
    {
        pwmval_fly4 = pwmval_fly4 - 10;
        analogWrite(fly4_PWMpin, pwmval_fly4);
        delay(50);
    }
    while (pwmval_fly4 > 60);
    pwmval_fly4 = 0;
    analogWrite(fly4_PWMpin, pwmval_fly4);

    // returntoinitialconfiguration();

}
void savetobalanceMATs()
{
    for (int ii = 0; ii < 3; ++ii)
    {
        for (int jj = 0; jj < 3; ++jj)
        {
            phimat(dur, 3 * ii + jj) = phi(ii, jj);
        }

    }

    dwmat(dur, 0) = dw(0);
    dwmat(dur, 1) = dw(1);
    dwmat(dur, 2) = dw(2);
}
void savetoMATs()
{
    MATw(dur, 0) = w(0);
    MATw(dur, 1) = w(1);
    MATw(dur, 2) = w(2);
    MATmanip(dur) = manip;
    MATddot(dur, 0) = ddot(0);
    MATddot(dur, 1) = ddot(1);
    MATddot(dur, 2) = ddot(2);
    MATddot(dur, 3) = ddot(3);
    MATangles(dur, 0) = gimbalangles(0);
    MATangles(dur, 1) = gimbalangles(1);
    MATangles(dur, 2) = gimbalangles(2);
    MATangles(dur, 3) = gimbalangles(3);
}

void setup()
{

    pinMode(ledpin, OUTPUT);
    digitalWrite(ledpin, HIGH);
    Wire.begin();
    if (imu.begin() == false) // with no arguments, this uses default addresses (AG:0x6B, M:0x1E) and i2c port (Wire).
    {
        Serial.println("Failed to communicate with LSM9DS1.");
        Serial.println("Double-check wiring.");
        Serial.println("Default settings in this sketch will " \
                       "work for an out of the box LSM9DS1 " \
                       "Breakout, but may need to be modified " \
                       "if the board jumpers are.");
        while (1);
    }
    Serial.println("================================");
    Serial.println("accel/gyro Calibration: Started");
    delay(1000);
    imu.calibrate(1);
    Serial.println("accel/gyro Calibration: Fininshed");
    Serial.println("================================");
    delay(1000);


    Serial.begin(115200); //Start Serial at 115200bps
    pinMode(led,  OUTPUT);
    pinMode(gimbal1_DIRpin, OUTPUT);
    pinMode(gimbal1_PWMpin, OUTPUT);
    pinMode(gimbal2_DIRpin, OUTPUT);
    pinMode(gimbal2_PWMpin, OUTPUT);
    pinMode(gimbal3_DIRpin, OUTPUT);
    pinMode(gimbal3_PWMpin, OUTPUT);
    pinMode(gimbal4_DIRpin, OUTPUT);
    pinMode(gimbal4_PWMpin, OUTPUT);
    pinMode(fly1_PWMpin, OUTPUT);
    pinMode(fly1_DIRpin, OUTPUT);
    pinMode(fly2_PWMpin, OUTPUT);
    pinMode(fly2_DIRpin, OUTPUT);
    pinMode(fly3_PWMpin, OUTPUT);
    pinMode(fly3_DIRpin, OUTPUT);
    pinMode(fly4_PWMpin, OUTPUT);
    pinMode(fly4_DIRpin, OUTPUT);

    pinMode(photores, INPUT);
    pinMode(orangeled,  OUTPUT);
    pinMode(greenled,  OUTPUT);

    //interrupt-counting_pins
    pinMode(gimbal1_encA, INPUT);
    pinMode(gimbal1_encB, INPUT);
    pinMode(gimbal2_encA, INPUT);
    pinMode(gimbal2_encB, INPUT);
    pinMode(gimbal3_encA, INPUT);
    pinMode(gimbal3_encB, INPUT);
    pinMode(gimbal4_encA, INPUT);
    pinMode(gimbal4_encB, INPUT);
    pinMode(fly1_encA, INPUT);
    pinMode(fly1_encB, INPUT);
    pinMode(fly2_encA, INPUT);
    pinMode(fly2_encB, INPUT);
    pinMode(fly3_encA, INPUT);
    pinMode(fly3_encB, INPUT);
    pinMode(fly4_encA, INPUT);
    pinMode(fly4_encB, INPUT);

    // SOS !!!  : Το digitalPinToInterrupt δε το θελει  στο τεενσυ
    attachInterrupt(gimbal1_encA, ISR_GIMBAL1_A, CHANGE);
    attachInterrupt(gimbal1_encB, ISR_GIMBAL1_B, CHANGE);
    attachInterrupt(fly1_encA, ISR_FLY1_A, CHANGE);
    attachInterrupt(fly1_encB, ISR_FLY1_B, CHANGE);

    attachInterrupt(gimbal2_encA, ISR_GIMBAL2_A, CHANGE);
    attachInterrupt(gimbal2_encB, ISR_GIMBAL2_B, CHANGE);
    attachInterrupt(fly2_encA, ISR_FLY2_A, CHANGE);
    attachInterrupt(fly2_encB, ISR_FLY2_B, CHANGE);

    attachInterrupt(gimbal3_encA, ISR_GIMBAL3_A, CHANGE);
    attachInterrupt(gimbal3_encB, ISR_GIMBAL3_B, CHANGE);
    attachInterrupt(fly3_encA, ISR_FLY3_A, CHANGE);
    attachInterrupt(fly3_encB, ISR_FLY3_B, CHANGE);

    attachInterrupt(gimbal4_encA, ISR_GIMBAL4_A, CHANGE);
    attachInterrupt(gimbal4_encB, ISR_GIMBAL4_B, CHANGE);
    attachInterrupt(fly4_encA, ISR_FLY4_A, CHANGE);
    attachInterrupt(fly4_encB, ISR_FLY4_B, CHANGE);
    //variable initialization


    digitalWrite(led, HIGH);
    calibrate_photoresistor();
    resetcounts();
    initializegimbalangles();
    resetcounts();
    delay(500);

    Serial.println("================================");
    RPY = getRPY();
    headingbias = RPY(2);
    Serial.print("headbias= ");
    Serial.println(headingbias);
    Serial.println("================================");

    // spinupflywheels();

    MATw.Fill(9998.98);
    MATddot.Fill(9998.98);
    MATangles.Fill(9998.98);
    MATmanip.Fill(9998.98);

    phimat.Fill(9998.98);
    dwmat.Fill(9998.98);

    resetcounts();
}
int cled = 0;
int kk = 0;
void loop()
{



    if (photoresistortouched() == 1 )
    {
        stopeverything();
        delay(1000);

        while(1)
        {
            if (cled < 1)
            {
                myblink();
            }
            else
            {
                digitalWrite(orangeled, HIGH);
                digitalWrite(greenled, LOW);
            }
            if (photoresistortouched() == 1 )
            {
                cled = cled + 1;
                digitalWrite(orangeled, HIGH);
                digitalWrite(greenled, LOW);

                Serial.println("======= Printing Results =========");
                printbalanceMATs();
            }
        }
    }
    else
    {

        if ((lastPrint + PRINT_SPEED) < millis()) //accurate time-keeping
        {
            // Serial.println(millis()-lastPrint);



            RPY = getRPY();
            // w(0) = (RPY(0) - RPYprev(0)) / dt;
            // w(1) = (RPY(1) - RPYprev(1)) / dt;
            // w(2) = (RPY(2) - RPYprev(2)) / dt;
            w(0) = imu.calcGyro(gxglobal) * (PI / 180); //convert from deg 2 rad
            w(1) = -imu.calcGyro(gyglobal) * (PI / 180); //convert from deg 2 rad
            w(2) = imu.calcGyro(gzglobal) * (PI / 180); //convert from deg 2 rad
            // double w1 = w(0);
            // double w2 = w(1);
            // double w3 = w(2);
            // Matrix<3, 6> WMEGA = {w1, 0, 0, w2, w3, 0,
            //                       0, w2, 0, w1, 0, w3,
            //                       0, 0, w3, 0, w1, w2
            //                      };
            // Matrix<3, 1> grav = {0, 0, -9.81};
            // Matrix<3, 3> wskew = skew(w);
            // Matrix<3, 3> gravskew = skew(grav);

            // cross_sum = cross_sum + (wskew * WMEGA) * dt;

            // Matrix<3, 6> leftterm1 = WMEGA + cross_sum;

            // gravskew = skew({9.81 * sin(RPY(1)), -9.81 * cos(RPY(1))*sin(RPY(0)), -9.81 * cos(RPY(1))*cos(RPY(0))}); //from paper 29 in air bearing folder

            // gravskewtot = gravskewtot + gravskew * dt;

            // leftterm = leftterm1 || gravskewtot;



            // h01 =  0.001041332; //J1*flywheelspeed1
            // h02 =  0.001041332; //J2*flywheelspeed2
            // h03 =  0.001041332; //J3*flywheelspeed3
            // h04 =  0.001041332; //J4*flywheelspeed4

            // d1 = PI / 2;
            // d2 = PI / 2;
            // d3 = PI / 2;
            // d4 = PI / 2;

            // // d1 = 0;
            // // d2 = 0;
            // // d3 = 0;
            // // d4 = 0;

            // hvector = {-cos(b) *sin(d1) *h01 - cos(d2) *h02 + cos(b) *sin(d3) *h03 + cos(d4) *h04,
            //            cos(d1) *h01 - cos(b) *sin(d2) *h02 - cos(d3) *h03 + cos(b) *sin(d4) *h04,
            //            sin(b) *sin(d1) *h01 + sin(b) *sin(d2) *h02 + sin(b) *sin(d3) *h03 + sin(b) *sin(d4) *h04
            //           };

            // Matrix<3, 1> rightterm2now = cross(w, hvector) ;
            // rightterm2tot = rightterm2tot + rightterm2now * dt;
            // rightterm = - hvector - rightterm2tot;

            // // Serial.println("L");
            // // printBLAmat3x9(leftterm);
            // // Serial.println("R");
            // // printBLAmat3x1(rightterm);

            double M = 0.678;
            double g = 9.81;
            double Ixx = 0.0046;
            double Iyy = 0.0047;
            double Izz = 0.0034;
            double fi = RPY(0);
            double theta = RPY(1);
            double psi = RPY(2);
            if (dur == 0)
            {
                fiprev = fi;
                thetaprev = theta;
                psiprev = psi;
                w0prev = w(0);
                w1prev = w(1);
                w2prev = w(2);
            }
            double ph12 = - M * g * dt * ( cos(fi) * cos(theta) + cos(fiprev) * cos(thetaprev) ) / (2 * Ixx);

            double ph13 = M * g * dt * ( sin(fi) * cos(theta) + sin(fiprev) * cos(thetaprev) ) / (2 * Ixx);

            double ph21 = M * g * dt * ( cos(fi) * cos(theta) + cos(fiprev) * cos(thetaprev) ) / (2 * Iyy);

            double ph23 = M * g * dt * ( sin(theta) + sin(thetaprev) ) / (2 * Iyy);

            double ph31 = - M * g * dt * ( sin(fi) * cos(theta) + sin(fiprev) * cos(thetaprev) ) / (2 * Izz);

            double ph32 = - M * g * dt * ( sin(theta) + sin(thetaprev) ) / (2 * Izz);


            double dwx = w(0) - w0prev;
            double dwy = w(1) - w1prev;
            double dwz = w(2) - w2prev;

            w0prev = w(0);
            w1prev = w(1);
            w2prev = w(2);
            fiprev = fi;
            thetaprev = theta;
            psiprev = psi;

            phi = {0, ph12, ph13, ph21, 0, ph23, ph31, ph32, 0};
            dw = {dwx, dwy, dwz};


            Matrix < 900, 3 > phiv;
            Matrix < 900, 1 > dwv;
            for (int i = kk; i < kk + 3; ++i)
            {
                for (int j = 0; j < 3; ++j)
                {
                    phiv(i, j) = phi(i - kk, j);
                }

                dwv(i) = dw(i - kk);
            }
            kk = kk + 3;

            Invert(phiv);

            Matrix < 3, 1 >rxyz = phiv * dwv;

            savetobalanceMATs();
            // savetoMATs();





            dur = dur + 1;
            lastPrint = millis(); // Update lastPrint time

        }

    }





}
