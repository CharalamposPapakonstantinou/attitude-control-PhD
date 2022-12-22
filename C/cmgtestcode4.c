////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////           ---- NOTES ----                  //////////////////////////////
// -- Make returntoinitialconfiguration function to return gimbals to initial-zero configuration
// -- FIX printing the matrices/values after finishing the manoeuvre / maybe insert delay in prints
// -- fix controlflywheelspeed to reach steady state before starting
// -- change h0 in h matrix to J*w
// -- warning: 'errorfly1' may be used uninitialized in this function  errorfly1prev = errorfly1;
// check the axiserror sign. propably is correct but check it anyway
// change b
// CALCULATING W FROM GYRO ONLY I DONT USE THE IMU ROLL PITCH YAW AT ALL!!!
// Do initializegimbals first and roll pitch yaw calibration after.
////////////////////////////////////////////////////////////////////////////////////////////////////////

// CODE COMMENT
// cmgtestcode2.ino includes the MATs with savetoMATs and printMATS functions along the SD card

////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////                                 START IMU                     ///////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
#include <Wire.h>
#include <SPI.h>
#include <SparkFunLSM9DS1.h>
#include <BasicLinearAlgebra.h>
#include <SD.h>
using namespace BLA;
LSM9DS1 imu;
#define PRINT_RAW
#define PRINT_SPEED 50 //50 // 250 ms between prints
double PRINT_SPEEDCFS = 100;
double dtCFS = PRINT_SPEEDCFS * 0.001;
static unsigned long lastPrint = 0; // Keep track of print time
static unsigned long lastPrinthead = 0; // Keep track of print time
static unsigned long lastPrintCFS = 0;
#define DECLINATION 4.53
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
double mxprev = 0;
double myprev = 0;
double mzprev = 0;
int dur = 0;
double dt = (PRINT_SPEED) * 0.001 ;
////////////////////////// OPTIONS ////////////////////////////////////////////////////////////////
int ONLYGYRO = 0;   //for 1 I use the usegyroinloop , for 0 I use getRPY
int writetoSD = 1;
int fhywheelsON = 1;
double gimbalsON = 1;
int print_attitude = 0;
double IMUcalibrate = 0;
int print_gyro = 0;
int print_mag = 0; //set to 1 to get data for calibration, normally must be 0 to get the correct mag vlues
int print_rpms = 0;
int printgimbals = 0;
////////////////////////////////////////////////////////////////////////////////////////////////
Matrix<3, 1> RPY;
Matrix<3, 1> RPYprev;
double rpy0 = 0;
double rpy1 = 0;
double rpy2 = 0;
double noise_gx = 0;
double noise_gy = 0;
double noise_gz = 0;
double noise_ax = 0;
double noise_ay = 0;
double noise_az = 0;
float noise_counterg = 0; //wste na einai 1 otan fhywheelsON=0 kai den kaleitai h spinup
float noise_countera = 0;



double AX = 0;
double AY = 0;
double AZ = 0;

double mxwrite;
double mywrite;
double mzwrite;


Matrix<3, 1> printAttitude(float ax, float ay, float az, float gx, float gy, float gz, float mx, float my, float mz)
{
    if (print_mag == 0) //dont use this when i do calibration
    {

        //WITH cable
        // Matrix<3, 1> mtemp = {1532.287 , 741.8238 , 385.5341};
        // mx = mx - (mtemp(0));
        // my = my - (mtemp(1));
        // mz = mz - (mtemp(2));
        // Matrix<3, 3> Mc = {1 , 0.041334 , -0.0096088 , 0.041334 , 0.98117 , 0.0087635 , -0.0096088 , 0.0087635 , 1.0528};

        //WITHOUT cable
        Matrix<3, 1> mtemp = {1880.0329, -7966.7361, -1971.0795};
        mx = mx - (mtemp(0));
        my = my - (mtemp(1));
        mz = mz - (mtemp(2));
        Matrix<3, 3> Mc = {1, 0.039117, -0.0032794, 0.039117, 0.97401, 0.0096616, -0.0032794, 0.0096616, 1.032};


        mx =  Mc(0, 0) * mx +  Mc(1, 0) * my +  Mc(2, 0) * mz;
        my =  Mc(0, 1) * mx +  Mc(1, 1) * my +  Mc(2, 1) * mz;
        mz =  Mc(0, 2) * mx +  Mc(1, 2) * my +  Mc(2, 2) * mz;
    }
    double magnorm = sqrt(mx * mx + my * my + mz * mz);
    mx = mx / magnorm;
    my = my / magnorm;
    mz = mz / magnorm;

    gx = gx - (1 - IMUcalibrate) * imu.gBiasRaw[0];
    gy = gy - (1 - IMUcalibrate) * imu.gBiasRaw[1];
    gz = gz - (1 - IMUcalibrate) * imu.gBiasRaw[2];


    // noise is calculated inside getgyroaccelnoise
    if (noise_counterg > 0)
    {
        gx = gx - noise_gx / noise_counterg;
        gy = gy - noise_gy / noise_counterg;
        gz = gz - noise_gz / noise_counterg;
    }

    float fact_gyro = 0.25; //0.25 //0.15 for PRINT_SPEED= 10ms
    float fact_acc = 0;  //0        # i dont use it
    float fact_mag = 0;  //0.03   # i dont use it
    gx = (1 - fact_gyro) * gxprev + fact_gyro * gx;
    gy = (1 - fact_gyro) * gyprev + fact_gyro * gy;
    gz = (1 - fact_gyro) * gzprev + fact_gyro * gz;

    gxprev = gx;
    gyprev = gy;
    gzprev = gz;
    gxglobal = gx;
    gyglobal = gy;
    gzglobal = gz;

    // noise is calculated inside getgyroaccelnoise
    if (noise_countera > 0)
    {
        ax = ax - noise_ax / noise_countera;
        ay = ay - noise_ay / noise_countera;
        az = az - noise_az / noise_countera;
    }

    ax = (1 - fact_acc) * axprev + fact_acc * ax;
    ay = (1 - fact_acc) * ayprev + fact_acc * ay;
    az = (1 - fact_acc) * azprev + fact_acc * az;
    axprev = ax;
    ayprev = ay;
    azprev = az;


    mx = (1 - fact_mag) * mxprev + fact_mag * mx;
    my = (1 - fact_mag) * myprev + fact_mag * my;
    mz = (1 - fact_mag) * mzprev + fact_mag * mz;
    mxprev = mx;
    myprev = my;
    mzprev = mz;


    if (print_mag == 1)
    {
        // Serial.print("M     ");
        Serial.print(mx);
        Serial.print("      ,      ");
        Serial.print(my);
        Serial.print("      ,      ");
        Serial.print(mz);
        Serial.println("      ;      ");
    }


    float roll = atan2(ay, az);
    float pitch = atan2(-ax, sqrt(ay * ay + az * az));


    pitch = 0;
    roll = 0;
    // these are the correct equations
    float magx = (-mx * cos(pitch) + my * sin(pitch) * sin(roll) + mz * sin(pitch) * cos(roll));
    float magy = (-my * cos(roll) + mz * sin(roll));
    float heading  = atan2(magy, magx); //the same as aboive but with -mx


    heading -= DECLINATION * PI / 180;
    if (heading > PI) heading -= (2 * PI);
    else if (heading < -PI) heading += (2 * PI);

    aroll = 0.05; //0.05    # i dont use it
    apitch = 0.05; //0.05   # i dont use it
    ahead = 0; //0.05


    rollc = (1 - aroll) * (rollcprev + dt * imu.calcGyro(gx) * PI / 180) + aroll * roll;
    pitchc = (1 - apitch) * (pitchcprev +  dt * imu.calcGyro(gy) * PI / 180) + apitch * pitch;
    headingc = (1 - ahead) * (headingcprev + dt * (-imu.calcGyro(gz)) * PI / 180) + ahead * heading;


    rollcprev = rollc;
    pitchcprev = pitchc;
    headingcprev = headingc;

    // heading *= 180.0 / PI;
    pitchc *= 180.0 / PI;
    rollc  *= 180.0 / PI;
    headingc *= 180.0 / PI;

    double finalroll = -rollc * 0;
    double finalpitch = pitchc * 0;
    double finalheading = -headingc - (headingbias) * 180 / PI;

    if (print_attitude == 1 )
    {
        // Serial.print("N   ");
        // Serial.print(roll, 2);
        // Serial.print("  ,  ");
        // Serial.print(pitch, 2);
        // Serial.print("  ,  ");
        // Serial.println(heading, 2);

        Serial.print("C   ");
        Serial.print(finalroll, 2);
        Serial.print("  ,  ");
        Serial.print(finalpitch, 2);
        Serial.print("  ,  ");
        Serial.println(finalheading, 2);
    }

    if (print_gyro == 1)
    {
        // imu.calcGyro(gx) gives result in deg/s
        Serial.print("G   ");
        Serial.print(imu.calcGyro(gx), 2);
        Serial.print("  ,  ");
        Serial.print(-imu.calcGyro(gy), 2);
        Serial.print("  ,  ");
        Serial.println(imu.calcGyro(gz), 2);

    }


    Matrix<3, 1> rpytemp = {finalroll, finalpitch, finalheading};
    return rpytemp;
}

Matrix<3, 1> getRPY()
{


    if ( imu.gyroAvailable() )
    {
        imu.readGyro();
    }
    else
    {
        Serial.println("gyro not ready");
    }
    if ( imu.accelAvailable())
    {
        imu.readAccel();
    }
    else
    {
        Serial.println("accel not ready");
    }
    if ( imu.magAvailable() )
    {
        imu.readMag();
    }
    else
    {
        // Serial.println("magn not ready");
    }


    RPY = printAttitude(imu.ax, imu.ay, imu.az, imu.gx, imu.gy, imu.gz, imu.mx, imu.my, imu.mz );

    RPY = RPY * (PI / 180); //convert RPY from deg to rad

    return RPY;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////                                 END IMU                       ///////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
void printBLAmat300x3(Matrix<1200, 3> mat)
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
        delay(70);
    }
}
void printBLAmat300x4(Matrix<1200, 4> mat)
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
        delay(70);
    }
}
void printBLAmat300x1(Matrix<1200, 1> mat)
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
        delay(70);
    }
}

/////

void printBLAmat3x3(Matrix<3, 3> mat)
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
void printBLAmat3x4(Matrix<3, 4> mat)
{
    for (int i = 0; i < mat.GetRowCount(); i++)
    {
        for (int j = 0; j < mat.GetColCount(); j++)
        {
            Serial.print(mat(i, j), 8);
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
            Serial.print(mat(i, j));
            Serial.print(' ');
        }
        Serial.println(' ');
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

Matrix<3, 3> abs_skew(Matrix<3, 1> mat)
{
    Matrix<3, 3> S = {0, abs(mat(2)), abs(mat(1)),
                      abs(mat(2)), 0, abs(mat(0)),
                      abs(mat(1)), abs(mat(0)), 0
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
Matrix<3, 1> quat2eul(Matrix<4, 1> quat) // in ZYX sequence
{

    double qw = quat(0);
    double qx = quat(1);
    double qy = quat(2);
    double qz = quat(3);

    double aSinInput = -2 * (qx * qz - qw * qy);
    if  (aSinInput > 1)
    {
        aSinInput = 1;
    }

    Matrix<3, 1> res = { atan2( 2 * (qx * qy + qw * qz), pow(qw, 2) + pow(qx, 2) - pow(qy, 2) - pow(qz, 2) ),
                         asin( aSinInput ),
                         atan2( 2 * (qy * qz + qw * qx), pow(qw, 2) - pow(qx, 2) - pow(qy, 2) + pow(qz, 2) )
                       };
    return res;
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
Matrix<4, 1> qref; //{0.7071, 0, 0, 0.7071};
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
Matrix <4, 1> ddotpwmprev = {0, 0, 0, 0};
Matrix<4, 1> drefprev = {0, 0, 0, 0};
Matrix<4, 1> dref = {0, 0, 0, 0};
Matrix<3, 1> eulererror;
Matrix<3, 1> atangle;
Matrix<1200, 3> MATw; // 1 value in 0.01s ---> 300 values in 15 sec --> Ara 0.05*300=15s einai h maximum diarkeia
Matrix<1200, 1> MATmanip;
Matrix<1200, 4> MATddot;
Matrix<1200, 4> MATangles;
Matrix<1200, 3> MATeulererror;
Matrix<1200, 3> MATmagnetic;
int overflowlimit = 30 / (PRINT_SPEED * 0.001);
int overflow = 0;
double errorfly1sum = 0;
double errorfly2sum = 0;
double errorfly3sum = 0;
double errorfly4sum = 0;
double errorfly1prev = 0;
double errorfly2prev = 0;
double errorfly3prev = 0;
double errorfly4prev = 0;
double rpmsfly1;
double rpmsfly2;
double rpmsfly3;
double rpmsfly4;

float b = 0.9552187 ; //for b=54.73deg --> b=0.9552187;
float th;
float Kp = 6; //50
float Ki = 0.001; //0.01
float Kwmega = 6; //70
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
int pwmval_fly1 = 00;
int pwmval_fly2 = 00;
int pwmval_fly3 = 00;
int pwmval_fly4 = 00;
int gimbalstall = 40; //means that min pwm for gimbals is 40
int flystall = 40;
double flywhellcontrollerpwm1;
double flywhellcontrollerpwm2;
double flywhellcontrollerpwm3;
double flywhellcontrollerpwm4;

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
double rpms;
int ctimer = 0;
double ctimerlim;
int cled = 0;

double steady_statecounter = 0;
double steady_statelim = 3 / (PRINT_SPEED * 0.001); // 1sec/0.01 --> 100 cycles
double manoeuvreaccuracy = 1.5 * (PI / 180);


// SD CARD
const int chipSelect = BUILTIN_SDCARD;
File dataFile;
//




// ==================== GRI ========================
int GRI = 0;

Matrix<4, 3> AhashGRI(Matrix<3, 4> A)
{
    double l_gain = 100;
    double mi = -100;
    double e0 = 1.5e-1;
    double freq = 0.4;
    double phi1 = PI / 3;
    double phi2 = PI / 2;
    double phi3 = PI;
    double lamda;
    double e1, e2, e3;
    Matrix<3, 3> V;
    Matrix<3, 3> temp;
    Matrix<4, 4> W;
    Matrix<4, 3> resAhash;

    lamda = l_gain * exp(mi * manip);
    // Serial.println(lamda);
    e1 = e0 * sin(freq * dur + phi1);
    e2 = e0 * sin(freq * dur + phi2);
    e3 = e0 * sin(freq * dur + phi3);
    V =  abs_skew({e1, e2, e3})*lamda;
    W = {1.1, lamda, lamda, lamda,
         lamda, 1.2, lamda, lamda,
         lamda, lamda, 1.3, lamda,
         lamda, lamda, lamda, 1.4
        };
    temp = (A * W * (~A) + V);
    Invert(temp);
    resAhash = W * (~A) * temp;
    return resAhash;
}
// ==================== END GRI ====================

Matrix<4, 1> ddot2pwm()
{

    Matrix<4, 1> ddotpwmval;

    if (ddotpwmval(0) >= 0)
    {
        ddotpwmval(0) = 210.7857 * ddot(0) + 9.6999;
    }
    else
    {
        ddotpwmval(0) = -(210.7857 * (-ddot(0)) + 9.6999);
    }
    //////////
    if (ddotpwmval(1) >= 0)
    {
        ddotpwmval(1) = 206.8227 * ddot(1) + 15.3940;
    }
    else
    {
        ddotpwmval(1) = -(206.8227 * (-ddot(1)) + 15.3940);
    }
    //////////
    if (ddotpwmval(2) >= 0)
    {
        ddotpwmval(2) = 215.7739 * ddot(2) + 12.0957;
    }
    else
    {
        ddotpwmval(2) = -(215.7739 * (-ddot(2)) + 12.0957);
    }
    //////////
    if (ddotpwmval(3) >= 0)
    {
        ddotpwmval(3) = 207.5022 * ddot(2) + 16.0053;
    }
    else
    {
        ddotpwmval(3) = -(207.5022 * (-ddot(2)) + 16.0053);
    }

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

double myabsmax(Matrix<4, 1> mat)
{
    double maxval = -100000;
    for (int i = 0; i < 4; ++i)
    {
        mat(i) = abs(mat(i));
        if (mat(i) > maxval)
        {
            maxval = mat(i);
        }
    }

    return maxval;
}

void saturateddot()
{
    double ddotlim = 1.1257; // it is the min of the max velocties, see code "rpm2pwm.m"
    double ddotmax = myabsmax(ddot);
    double ddgamma = ddotmax / ddotlim;
    if (ddgamma > 1)
    {
        for (int i = 0; i < 4; ++i)
        {
            ddot(i) = ddot(i) / ddgamma;
        }
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
    Serial.print(gimbalangles(0) * 180 / PI, 3);
    Serial.print("  ,  ");
    Serial.print(gimbalangles(1) * 180 / PI, 3);
    Serial.print("  ,  ");
    Serial.print(gimbalangles(2) * 180 / PI, 3);
    Serial.print("  ,  ");
    Serial.println(gimbalangles(3) * 180 / PI, 3);
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
    // Serial.print("Q   ");
    Serial.print(p(0), 3);
    Serial.print("  ,  ");
    Serial.print(p(1), 3);
    Serial.print("  ,  ");
    Serial.print(p(2), 3);
    Serial.print("  ,  ");
    Serial.print(p(3), 3);
    Serial.println("  ;  ");
}
void printMATs()
{
    Serial.flush();
    Serial.println("W");
    printBLAmat300x3(MATw);
    Serial.println("M");
    printBLAmat300x1(MATmanip);
    Serial.println("D");
    printBLAmat300x4(MATddot);
    Serial.println("A");
    printBLAmat300x4(MATangles);
    Serial.println("R");
    printBLAmat300x3(MATeulererror);
    if (print_mag == 1)
    {
        Serial.println("G"); //G for Geomagnetic field
        printBLAmat300x3(MATmagnetic);
    }
    Serial.println("END");
}
/////////////// END PRINT FUNCTIONS /////////////////////////////
static unsigned long lastPrintgimb = 0;
double PRINT_SPEEDgimb = 3;

void gimbalsPIDcontrol()
{
    double KG_Pgain = 5, KG_Dgain = 15, KG_Igain = 0.00;
    int gi = 0;
    double errorgimbal1, errorgimbal2, errorgimbal3, errorgimbal4;
    double errorgimbal1prev = 0, errorgimbal2prev = 0, errorgimbal3prev = 0, errorgimbal4prev = 0, errorgimbal1sum = 0, errorgimbal2sum = 0, errorgimbal3sum = 0, errorgimbal4sum = 0;

    Matrix<4, 1> gimbalcontrollerpwm;
    while(gi < 5)
    {
        if ((lastPrintgimb + PRINT_SPEEDgimb) <= millis())
        {
            lastPrintgimb = millis();

            errorgimbal1 = (dref(0) - delta1);
            errorgimbal2 = (dref(1) - delta2);
            errorgimbal3 = (dref(2) - delta3);
            errorgimbal4 = (dref(3) - delta4);
            errorgimbal1sum = errorgimbal1sum + errorgimbal1;
            gimbalcontrollerpwm(0) = ddotpwm(0) + KG_Pgain * errorgimbal1 + KG_Dgain * (errorgimbal1 - errorgimbal1prev) + KG_Igain * errorgimbal1sum;
            gimbalcontrollerpwm(1) = ddotpwm(1) + KG_Pgain * errorgimbal2 + KG_Dgain * (errorgimbal2 - errorgimbal2prev) + KG_Igain * errorgimbal2sum;
            gimbalcontrollerpwm(2) = ddotpwm(2) + KG_Pgain * errorgimbal3 + KG_Dgain * (errorgimbal3 - errorgimbal3prev) + KG_Igain * errorgimbal3sum;
            gimbalcontrollerpwm(3) = ddotpwm(3) + KG_Pgain * errorgimbal4 + KG_Dgain * (errorgimbal4 - errorgimbal4prev) + KG_Igain * errorgimbal4sum;

            // Serial.print(dref(0) *180/PI); Serial.print("   ");
            // Serial.println(delta1 *180/PI); Serial.print("   ");
            // Serial.print(errorgimbal1 * 180 / PI);
            // Serial.print("   ");
            // Serial.print(errorgimbal2 * 180 / PI);
            // Serial.print("   ");
            // Serial.print(errorgimbal3 * 180 / PI);
            // Serial.print("   ");
            // Serial.println(errorgimbal4 * 180 / PI);

            if (gimbalcontrollerpwm(0) > 0)
            {
                digitalWrite(gimbal1_DIRpin, LOW);
            }
            else
            {
                digitalWrite(gimbal1_DIRpin, HIGH);
            }
            if (gimbalcontrollerpwm(1) > 0)
            {
                digitalWrite(gimbal2_DIRpin, LOW);
            }
            else
            {
                digitalWrite(gimbal2_DIRpin, HIGH);
            }
            if (gimbalcontrollerpwm(2) > 0)
            {
                digitalWrite(gimbal3_DIRpin, LOW);
            }
            else
            {
                digitalWrite(gimbal3_DIRpin, HIGH);
            }
            if (gimbalcontrollerpwm(3) > 0)
            {
                digitalWrite(gimbal4_DIRpin, LOW);
            }
            else
            {
                digitalWrite(gimbal4_DIRpin, HIGH);
            }


            for (int i = 0; i < 4; ++i)
            {
                gimbalcontrollerpwm(i) = abs(gimbalcontrollerpwm(i));
                gimbalcontrollerpwm(i) = gimbalcontrollerpwm(i) > 255 ? 255 : gimbalcontrollerpwm(i);
                gimbalcontrollerpwm(i) = gimbalcontrollerpwm(i) < gimbalstall ? 0 : gimbalcontrollerpwm(i);
            }


            analogWrite(gimbal1_PWMpin, gimbalcontrollerpwm(0));
            analogWrite(gimbal2_PWMpin, gimbalcontrollerpwm(1));
            analogWrite(gimbal3_PWMpin, gimbalcontrollerpwm(2));
            analogWrite(gimbal4_PWMpin, gimbalcontrollerpwm(3));


            errorgimbal1prev = errorgimbal1;

        }

        gi = gi + 1;
    }

}

void low_pass_gimbal_rates()
{
    float fact_ddotpwm = 0.05;
    for (int ddi = 0; ddi < 4; ++ddi)
    {
        ddotpwm(ddi) = (1 - fact_ddotpwm) * ddotpwmprev(ddi) + fact_ddotpwm * ddotpwm(ddi);
        ddotpwmprev(ddi) = ddotpwm(ddi);
    }
}
void getgyroaccelnoise()
{

    if ( imu.gyroAvailable() )
    {
        imu.readGyro();
        noise_gx = noise_gx + (imu.gx - (1 - IMUcalibrate) * imu.gBiasRaw[0]);
        noise_gy = noise_gy + (imu.gy - (1 - IMUcalibrate) * imu.gBiasRaw[1]);
        noise_gz = noise_gz + (imu.gz - (1 - IMUcalibrate) * imu.gBiasRaw[2]);
        noise_counterg = noise_counterg + 1;
    }

    if ( imu.accelAvailable() )
    {
        imu.readAccel();
        noise_ax = noise_ax + (AX - imu.ax);
        noise_ay = noise_ay + (AY - imu.ay);
        noise_az = noise_az + (AZ - imu.az);
        noise_countera = noise_countera + 1;
    }


}
void controlflywheelspeed()
{

    double errorfly1, errorfly2, errorfly3, errorfly4;



    if (ctimer == 1)
    {
        errorfly1 = 0;
        errorfly2 = 0;
        errorfly3 = 0;
        errorfly4 = 0;

        errorfly1prev = errorfly1;
        errorfly2prev = errorfly2;
        errorfly3prev = errorfly3;
        errorfly4prev = errorfly4;

    }

    if ((lastPrintCFS + PRINT_SPEEDCFS) <= millis())
    {

        // Serial.print(millis() - lastPrintCFS);
        // Serial.print("  ");

        // Serial.println(millis() - lastPrintCFS);
        // if (ctimer > 0)
        // {
        //     dtCFS = 0.001 * (millis() - lastPrintCFS);
        // }

        lastPrintCFS = millis();


        double rpmsref = 4000;
        double Kgain = 0.05  ; //0.05
        double Kigain = 0.01 ; //0.02
        double Kdgain = 0.1 * ctimer / ctimerlim; //0.1
        if (ctimer >= ctimerlim)
        {
            Kdgain = 0.1;
        }
        ctimer = ctimer + 1;


        rpmsfly1 = abs(60 * (6 * (countfly1 - countfly1prev) / dtCFS) / 360);
        countfly1prev = countfly1;
        errorfly1 = (rpmsref - rpmsfly1);
        errorfly1sum = errorfly1sum + errorfly1;
        flywhellcontrollerpwm1 = 200 + Kgain * errorfly1 + Kigain * errorfly1sum + Kdgain * (errorfly1 - errorfly1prev);
        analogWrite(fly1_PWMpin, flywhellcontrollerpwm1);
        errorfly1prev = errorfly1;


        rpmsfly2 = abs(60 * (6 * (countfly2 - countfly2prev) / dtCFS) / 360);
        countfly2prev = countfly2;
        errorfly2 = (rpmsref - rpmsfly2);
        errorfly2sum = errorfly2sum + errorfly2;
        flywhellcontrollerpwm2 = 200 + Kgain * errorfly2 + Kigain * errorfly2sum + Kdgain * (errorfly2 - errorfly2prev);
        analogWrite(fly2_PWMpin, flywhellcontrollerpwm2);
        errorfly2prev = errorfly2;

        // Serial.print(Kgain * errorfly2 + Kigain * errorfly2sum + Kdgain * (errorfly2 - errorfly2prev));
        // Serial.print("  ");

        rpmsfly3 = abs(60 * (6 * (countfly3 - countfly3prev) / dtCFS) / 360);
        countfly3prev = countfly3;
        errorfly3 = (rpmsref - rpmsfly3);
        errorfly3sum = errorfly3sum + errorfly3;
        flywhellcontrollerpwm3 = 200 + Kgain * errorfly3 + Kigain * errorfly3sum + Kdgain * (errorfly3 - errorfly3prev);
        analogWrite(fly3_PWMpin, flywhellcontrollerpwm3 );
        errorfly3prev = errorfly3;


        rpmsfly4 = abs(60 * (6 * (countfly4 - countfly4prev) / dtCFS) / 360);
        countfly4prev = countfly4;
        errorfly4 = (rpmsref - rpmsfly4);
        errorfly4sum = errorfly4sum + errorfly4;
        flywhellcontrollerpwm4 = 200 + Kgain * errorfly4 + Kigain * errorfly4sum + Kdgain * (errorfly4 - errorfly4prev);
        analogWrite(fly4_PWMpin, flywhellcontrollerpwm4);
        errorfly4prev = errorfly4;

        if (print_rpms == 1)
        {
            Serial.print(rpmsfly1);
            Serial.print("  ");
            Serial.print(rpmsfly2);
            Serial.print("  ");
            Serial.print(rpmsfly3);
            Serial.print("  ");
            Serial.print(rpmsfly4);
            Serial.println(" ; ");
        }

    }

}
void myblink()
{
    digitalWrite(orangeled, HIGH);
    digitalWrite(greenled, LOW);
    delay(300);
    digitalWrite(orangeled, LOW);
    digitalWrite(greenled, HIGH);
    delay(300);
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
    countfly1prev = 0;
    countfly2prev = 0;
    countfly3prev = 0;
    countfly4prev = 0;

    delta1 = 0;
    delta2 = 0;
    delta3 = 0;
    delta4 = 0;
}
void returntoinitialconfiguration()
{
    double acc = 0.5;
    double refangle = 0;
    int speed = 180;
    double temp1 = fmod(abs(delta1 - refangle), (2 * PI));
    double temp2 = fmod(abs(delta2 - refangle), (2 * PI));
    double temp3 = fmod(abs(delta3 - refangle), (2 * PI));
    double temp4 = fmod(abs(delta4 - refangle), (2 * PI));
    int returned1, returned2, returned3, returned4, allreturned;
    do
    {
        if ( temp1 > acc)
        {
            if (mysign(delta1)*temp1 > PI)
            {
                digitalWrite(gimbal1_DIRpin, HIGH);
            }
            else
            {
                digitalWrite(gimbal1_DIRpin, LOW);
            }
            analogWrite(gimbal1_PWMpin, speed );
            temp1 = fmod(abs(delta1 - refangle), (2 * PI));
            returned1 = 0;
        }
        else
        {
            analogWrite(gimbal1_PWMpin, 0);
            returned1 = 1;
        }

        if ( temp2 > acc)
        {
            if (mysign(delta2)*temp2 > PI)
            {
                digitalWrite(gimbal2_DIRpin, HIGH);
            }
            else
            {
                digitalWrite(gimbal2_DIRpin, LOW);
            }
            analogWrite(gimbal2_PWMpin, speed );
            temp2 = fmod(abs(delta2 - refangle), (2 * PI));
            returned2 = 0;

        }
        else
        {
            analogWrite(gimbal2_PWMpin, 0);
            returned2 = 1;
        }

        if ( temp3 > acc)
        {
            if (mysign(delta3)*temp3 > PI)
            {
                digitalWrite(gimbal3_DIRpin, HIGH);
            }
            else
            {
                digitalWrite(gimbal3_DIRpin, LOW);
            }
            analogWrite(gimbal3_PWMpin, speed );
            temp3 = fmod(abs(delta3 - refangle), (2 * PI));
            returned3 = 0;

        }
        else
        {
            analogWrite(gimbal3_PWMpin, 0);
            returned3 = 1;
        }

        if ( temp4 > acc)
        {
            if (mysign(delta4)*temp4 > PI)
            {
                digitalWrite(gimbal4_DIRpin, HIGH);
            }
            else
            {
                digitalWrite(gimbal4_DIRpin, LOW);
            }
            analogWrite(gimbal4_PWMpin, speed );
            temp4 = fmod(abs(delta4 - refangle), (2 * PI));
            returned4 = 0;

        }
        else
        {
            analogWrite(gimbal4_PWMpin, 0);
            returned4 = 1;
        }

        allreturned = returned1 * returned2 * returned3 * returned4;
    }
    while(allreturned == 0);

}

void spinupflywheels()
{

    digitalWrite(orangeled, HIGH);
    digitalWrite(greenled, LOW);

    digitalWrite(fly1_DIRpin, HIGH);
    digitalWrite(fly2_DIRpin, HIGH);
    digitalWrite(fly3_DIRpin, HIGH);
    digitalWrite(fly4_DIRpin, HIGH);


    for (int i = flystall; i < 200; i = i + 10) //220 here and 200 in controlflywheelspeed
    {
        analogWrite(fly1_PWMpin, i);
        analogWrite(fly2_PWMpin, i);
        analogWrite(fly3_PWMpin, i);
        analogWrite(fly4_PWMpin, i);
        delay(400); //300
    }

    ctimerlim = 3 / dtCFS;
    resetcounts();
    do
    {
        controlflywheelspeed();
        getgyroaccelnoise(); //an vgalw ayto tote den kanoyn spike ta rpms twn flywheel
    }
    while(ctimer < ctimerlim);
    controlflywheelspeed();
    // delay(PRINT_SPEEDCFS-10); // = 51(ms)

    digitalWrite(orangeled, HIGH);
    digitalWrite(greenled, HIGH);

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
    float photores_denum = 2.5; //the lower the more sensitive
    photores_threshold = val / photores_denum;
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
    int flystopdelay = 400;
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


    int flywheelsstoped1;
    int flywheelsstoped2;
    int flywheelsstoped3;
    int flywheelsstoped4;
    int allflywheelsstoped;

    do
    {

        if(flywhellcontrollerpwm1 > flystall)
        {
            flywhellcontrollerpwm1 = flywhellcontrollerpwm1 - 10;
            analogWrite(fly1_PWMpin, flywhellcontrollerpwm1);
            flywheelsstoped1 = 0;
        }
        else
        {
            flywhellcontrollerpwm1 = 0;
            analogWrite(fly1_PWMpin, flywhellcontrollerpwm1);
            flywheelsstoped1 = 1;
        }


        if (flywhellcontrollerpwm2 > flystall)
        {
            flywhellcontrollerpwm2 = (flywhellcontrollerpwm2 - 10);
            analogWrite(fly2_PWMpin, flywhellcontrollerpwm2);
            flywheelsstoped2 = 0;
        }
        else
        {
            flywhellcontrollerpwm2 = 0;
            analogWrite(fly2_PWMpin, flywhellcontrollerpwm2);
            flywheelsstoped2 = 1;
        }

        if  (flywhellcontrollerpwm3 > flystall)
        {
            flywhellcontrollerpwm3 = (flywhellcontrollerpwm3 - 10);
            analogWrite(fly3_PWMpin, flywhellcontrollerpwm3);
            flywheelsstoped3 = 0;
        }
        else
        {
            flywhellcontrollerpwm3 = 0;
            analogWrite(fly3_PWMpin, flywhellcontrollerpwm3);
            flywheelsstoped3 = 1;
        }

        if (flywhellcontrollerpwm4 > flystall)
        {
            flywhellcontrollerpwm4 = (flywhellcontrollerpwm4 - 10);
            analogWrite(fly4_PWMpin, flywhellcontrollerpwm4);
            flywheelsstoped4 = 0;
        }
        else
        {
            flywhellcontrollerpwm4 = 0;
            analogWrite(fly4_PWMpin, flywhellcontrollerpwm4);
            flywheelsstoped4 = 1;
        }

        allflywheelsstoped = flywheelsstoped1 * flywheelsstoped2 * flywheelsstoped3 * flywheelsstoped4;
        delay(flystopdelay);

    }
    while (allflywheelsstoped == 0);

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
    MATeulererror(dur, 0) = eulererror(0);
    MATeulererror(dur, 1) = eulererror(1);
    MATeulererror(dur, 2) = eulererror(2);

    if (print_mag == 1)
    {
        MATmagnetic(dur, 0) = imu.mx;
        MATmagnetic(dur, 1) = imu.my;
        MATmagnetic(dur, 2) = imu.mz;
    }
}

void calcheadingbias()
{
    Serial.println("================================");
    double headingsamples = 100; //1000
    double headingbiastemp = 0;
    double garbagesamples = 50; //500
    int iter = 1;
    digitalWrite(orangeled, HIGH);
    digitalWrite(greenled, LOW);
    while (iter <= headingsamples)
    {
        if ((lastPrinthead + PRINT_SPEED) < millis()) //accurate time-keeping
        {
            lastPrinthead = millis();
            RPY = getRPY();
            if (iter > garbagesamples)
            {
                headingbiastemp = headingbiastemp + RPY(2);
            }
            iter = iter + 1;
        }
    }
    digitalWrite(orangeled, HIGH);
    digitalWrite(greenled, HIGH);
    headingbias = headingbiastemp / (headingsamples - garbagesamples);
    Serial.print("headbias= ");
    Serial.println(headingbias * 180 / PI);
    Serial.println("================================");
}
void prepareSD()
{
    Serial.print("Initializing SD card...");
    if (!SD.begin(chipSelect))
    {
        Serial.println("Card failed, or not present");
        while (1)
        {
            digitalWrite(ledpin, HIGH);
            delay(600);
            digitalWrite(ledpin, LOW);
            delay(600);
        }
    }
    Serial.println("card initialized.");
    if (SD.exists("CMGDATA.txt"))
    {
        SD.remove("CMGDATA.txt");
    }
    dataFile = SD.open("CMGDATA.txt", FILE_WRITE);


}
void writedatatoSD()
{
    if (dataFile)
    {
        int digits = 4; //3
        dataFile.print(w(0), digits);
        dataFile.print(" , ");
        dataFile.print(w(1), digits);
        dataFile.print(" , ");
        dataFile.print(w(2), digits);
        dataFile.print(" , ");
        dataFile.print(manip, digits);
        dataFile.print(" , ");
        dataFile.print(ddot(0), digits);
        dataFile.print(" , ");
        dataFile.print(ddot(1), digits);
        dataFile.print(" , ");
        dataFile.print(ddot(2), digits);
        dataFile.print(" , ");
        dataFile.print(ddot(3), digits);
        dataFile.print(" , ");
        dataFile.print(gimbalangles(0), digits);
        dataFile.print(" , ");
        dataFile.print(gimbalangles(1), digits);
        dataFile.print(" , ");
        dataFile.print(gimbalangles(2), digits);
        dataFile.print(" , ");
        dataFile.print(gimbalangles(3), digits);
        dataFile.print(" , ");
        dataFile.print(eulererror(0), digits);
        dataFile.print(" , ");
        dataFile.print(eulererror(1), digits);
        dataFile.print(" , ");
        dataFile.println(eulererror(2), digits);
    }
    else
    {
        Serial.println("error opening the .txt file");
    }
}

void usegyroinloop()
{
    if ( imu.gyroAvailable() )
    {
        imu.readGyro();
    }
    else
    {
        Serial.println("gyro not available");
    }

    double gx = imu.gx - (1 - IMUcalibrate) * imu.gBiasRaw[0]  - noise_gx / noise_counterg;
    double gy = imu.gy - (1 - IMUcalibrate) * imu.gBiasRaw[1]  - noise_gy / noise_counterg;
    double gz = imu.gz - (1 - IMUcalibrate) * imu.gBiasRaw[2]  - noise_gz / noise_counterg;

    double fact_gyro = 0.15; // 0.15
    gx = (1 - fact_gyro) * gxprev + fact_gyro * gx;
    gy = (1 - fact_gyro) * gyprev + fact_gyro * gy;
    gz = (1 - fact_gyro) * gzprev + fact_gyro * gz;

    gxprev = gx;
    gyprev = gy;
    gzprev = gz;

    gxglobal = gx;
    gyglobal = gy;
    gzglobal = gz;
    if (print_gyro == 1)
    {
        // imu.calcGyro(gx) gives result in deg/s
        Serial.print("G   ");
        Serial.print(imu.calcGyro(gxglobal), 2);
        Serial.print("  ,  ");
        Serial.print(-imu.calcGyro(gyglobal), 2);
        Serial.print("  ,  ");
        Serial.println(imu.calcGyro(gzglobal), 2);

    }
    if (print_attitude == 1)
    {
        Serial.print("rpy    ");
        Serial.print(rpy0 * 180 / PI, 2);
        Serial.print("  ,  ");
        Serial.print(rpy1 * 180 / PI, 2);
        Serial.print("  ,  ");
        Serial.println(rpy2 * 180 / PI, 2);
    }

}


// these are the help variable for the calculatepwm2rpm function
int pwmval = gimbalstall;
int samples = 1;
double c = 0;
double rpms1, rpms2, rpms3, rpms4;
int countgimbal1prev = 0;
int countgimbal2prev = 0;
int countgimbal3prev = 0;
int countgimbal4prev = 0;
double sumrpms1 = 0;
double sumrpms2 = 0;
double sumrpms3 = 0;
double sumrpms4 = 0;
void calculatepwm2rpm()
{

    if (samples <= 25 && pwmval <= 255 && photoresistortouched() == 0)
    {


        digitalWrite(gimbal1_DIRpin, LOW);
        digitalWrite(gimbal2_DIRpin, LOW);
        digitalWrite(gimbal3_DIRpin, LOW);
        digitalWrite(gimbal4_DIRpin, LOW);
        analogWrite(gimbal1_PWMpin, pwmval);
        analogWrite(gimbal2_PWMpin, pwmval);
        analogWrite(gimbal3_PWMpin, pwmval);
        analogWrite(gimbal4_PWMpin, pwmval);

        // detachInterrupt(inter_pin_1);
        rpms1 = 60 * (0.029821 * (countgimbal1 - countgimbal1prev) / dt) / 360; // the number 6 is derived by 360/(12*shaft_ratio) (=360/(12*5)), for gimbal is 360/(12*1006)=0.029821
        countgimbal1prev = countgimbal1;
        rpms2 = 60 * (0.029821 * (countgimbal2 - countgimbal2prev) / dt) / 360; // the number 6 is derived by 360/(12*shaft_ratio) (=360/(12*5)), for gimbal is 360/(12*1006)=0.029821
        countgimbal2prev = countgimbal2;
        rpms3 = 60 * (0.029821 * (countgimbal3 - countgimbal3prev) / dt) / 360; // the number 6 is derived by 360/(12*shaft_ratio) (=360/(12*5)), for gimbal is 360/(12*1006)=0.029821
        countgimbal3prev = countgimbal3;
        rpms4 = 60 * (0.029821 * (countgimbal4 - countgimbal4prev) / dt) / 360; // the number 6 is derived by 360/(12*shaft_ratio) (=360/(12*5)), for gimbal is 360/(12*1006)=0.029821
        countgimbal4prev = countgimbal4;
        // attachInterrupt( inter_pin_1, ISR_1, CHANGE);

        if (samples > 5)
        {
            // Serial.println(rpms);  //in rpm
            sumrpms1 = sumrpms1 + rpms1;
            sumrpms2 = sumrpms2 + rpms2;
            sumrpms3 = sumrpms3 + rpms3;
            sumrpms4 = sumrpms4 + rpms4;
            c = c + 1;
        }


        samples = samples + 1;

    }
    else
    {
        Serial.print("PWM= ");
        Serial.print(pwmval);
        Serial.print(" , mean rpms: ");
        Serial.print(sumrpms1 / c);
        Serial.print(" , ");
        Serial.print(sumrpms2 / c);
        Serial.print(" , ");
        Serial.print(sumrpms3 / c);
        Serial.print(" , ");
        Serial.println(sumrpms4 / c);
        pwmval = pwmval + 10;
        c = 0;
        samples = 1;
        sumrpms1 = 0;
        sumrpms2 = 0;
        sumrpms3 = 0;
        sumrpms4 = 0;
    }



}


void setup()
{

    pinMode(ledpin, OUTPUT);
    digitalWrite(ledpin, HIGH);
    Serial.begin(115200); //Start Serial at 115200bps

    if (writetoSD == 1)
    {
        prepareSD();  //initialize sd card and creat txt file
    }

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
    imu.calibrate(IMUcalibrate);
    Serial.println("accel/gyro Calibration: Fininshed");
    Serial.println("================================");
    delay(1000);

    imu.setMagODR(7); // 80Hz
    imu.setGyroODR(6);
    imu.setAccelODR(6);


    while (imu.accelAvailable())  //mia fora mpainei mono
    {
        imu.readAccel();
        AX = imu.ax;
        AY = imu.ay;
        AZ = imu.az;
    }

    if (fhywheelsON == 0)  //an fhywheelsON=0 pairnei apo edw ta noises, alliws h getgyroaccelnoise kaleitai mesa sthn spinup
    {
        for (int i = 0; i < 100; ++i) //i<100
        {
            getgyroaccelnoise();
            delay(10);
        }
    }

    calcheadingbias();
    RPY = getRPY();
    q = rpy2quat(RPY(2), RPY(1), RPY(0));
    qprev = q;
    atangle = RPY;
    Serial.println("---q_ini---");
    printquaternion(q);
    Serial.println("------");

    double y = (180) * PI / 180;
    double p = 0;
    double r = 0;
    qref = rpy2quat(y, p, r); //yaw,pitch,roll
    Serial.println("---q_ref---");
    printquaternion(qref);
    Serial.println("------");
    // MATw.Fill(9998.98);
    // MATddot.Fill(9998.98);
    // MATangles.Fill(9998.98);
    // MATmanip.Fill(9998.98);
    // MATeulererror.Fill(9998.98);
    // MATmagnetic.Fill(9998.98);

    if (fhywheelsON == 1)
    {
        spinupflywheels();
    }

}

void loop()
{



    if (photoresistortouched() == 1 || overflow == 1 || steady_statecounter >= steady_statelim)
    {
        dataFile.close();
        stopeverything();
        // returntoinitialconfiguration();
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
                // printMATs();
                delay(100000000);
            }
        }
    }
    else
    {

        if ((lastPrint + PRINT_SPEED) <= millis()) //accurate time-keeping
        {
            // Serial.println(millis() - lastPrint);
            if (dur > 0)
            {
                dt = 0.001 * (millis() - lastPrint);
            }
            lastPrint = millis();
            if (fhywheelsON == 1)
            {
                controlflywheelspeed();
            }


            if (ONLYGYRO == 1)
            {

                usegyroinloop();



                // SOS----- THIS WAY I DONT USE THE IMU ROLL PITCH YAW AT ALL!!!
                w(0) =  imu.calcGyro(gxglobal) * (PI / 180); //convert from deg 2 rad
                w(1) = -imu.calcGyro(gyglobal) * (PI / 180); //convert from deg 2 rad
                w(2) =  imu.calcGyro(gzglobal) * (PI / 180); //convert from deg 2 rad
                // SOS----- THIS WAY I DONT USE THE IMU ROLL PITCH YAW AT ALL!!!

                wquat = {0, w(0), w(1), w(2)};
                // printquaternion(q);
                qdot = quatmultiply(q, wquat) * 0.5;
                if (norm(w) != 0)
                {
                    slerppartA1 = cos(norm(w) * dt / 2);
                    slerppartA2 = (w / norm(w)) * sin(norm(w) * dt / 2);
                    slerppartA =  slerppartA1 && slerppartA2;
                    q = quatmultiply(qprev, slerppartA);
                }
                RPYprev = RPY;
                qprev = q;

                atangle = quat2eul(q);
                rpy0 = atangle(2);
                rpy1 = atangle(1);
                rpy2 = atangle(0);
            }
            else if (ONLYGYRO == 0)
            {
                //den xreiazetai kapoy to w stoys ypologismoys alla to vazw gia na to enhmerwnw kai na to grafw sthn SD
                w(0) =  imu.calcGyro(gxglobal) * (PI / 180); //convert from deg 2 rad
                w(1) = -imu.calcGyro(gyglobal) * (PI / 180); //convert from deg 2 rad
                w(2) =  imu.calcGyro(gzglobal) * (PI / 180); //convert from deg 2 rad
                RPY = getRPY();
                q = rpy2quat(RPY(2), RPY(1), RPY(0));
            }


            qerror = quatmultiply(quatconj(quatnormalize(qref)), quatnormalize((q)));

            eulererror = quat2eul(qerror);
            if (abs(eulererror(0)) < manoeuvreaccuracy)
            {
                steady_statecounter = steady_statecounter + 1;
            }
            else
            {
                steady_statecounter = 0;
            }

            th = 2 * acos(qerror(0)); //maybe there is a faster function than this(acos)
            axiserror = {qerror(1), qerror(2), qerror(3)};


            if (qerror(0) < 0)
            {
                axiserror = -axiserror;
            }

            Integral = Integral + axiserror;

            if (dur == 0)
            {
                axiserrorprev = axiserror;
            }

            u = -axiserror * Kp - Integral * Ki - w * Kwmega;

            A = {-cos(b) *cos(d1),  sin(d2),     cos(b) *cos(d3),  -sin(d4),
                 -sin(d1),    -cos(b) *cos(d2), sin(d3),     cos(b) *cos(d4),
                 sin(b) *cos(d1),  sin(b) *cos(d2),    sin(b) *cos(d3),  sin(b) *cos(d4)
                };
            // Kanonika prepei na pollaplasiasw to A me ta h01,h02...



            h01 = 2.04779008e-6 * rpmsfly1 * (2 * PI / 60); //J1*flywheelspeed1
            h02 = 2.04779008e-6  * rpmsfly2 * (2 * PI / 60); //J2*flywheelspeed2
            h03 = 2.04779008e-6  * rpmsfly3 * (2 * PI / 60); //J3*flywheelspeed3
            h04 = 2.04779008e-6  * rpmsfly4 * (2 * PI / 60); //J4*flywheelspeed4
            hvector = {-cos(b) *sin(d1) *h01 - cos(d2) *h02 + cos(b) *sin(d3) *h03 + cos(d4) *h04,
                       cos(d1) *h01 - cos(b) *sin(d2) *h02 - cos(d3) *h03 + cos(b) *sin(d4) *h04,
                       sin(b) *sin(d1) *h01 + sin(b) *sin(d2) *h02 + sin(b) *sin(d3) *h03 + sin(b) *sin(d4) *h04
                      };

            Hdot = -u - cross(w, hvector);

            AAT = A * (~A);
            manip = sqrt(det(AAT));


            if (GRI == 0)
            {
                // Matrix<3,3> AATtemp = AAT.Inverse();
                Invert(AAT); //  inverts a matrix and sets the result to the matrix itself,
                Ahash = (~A) * (AAT);
            }
            else
            {
                Ahash = AhashGRI(A);
            }
            // Serial.println("-----");
            // printBLAmat3x4(~((~A) * (AAT)));
            // printBLAmat3x4(~Ahash);
            // Serial.println("-----");
            ddot = Ahash * Hdot * gimbalsON; //ddot in rad/sec
            saturateddot();
            // printddot();


            dref =  ddot * dt + drefprev;
            drefprev = dref;



            ddotpwm = ddot2pwm();
            // low_pass_gimbal_rates();
            gimbalsPIDcontrol();



            // setgimbaldirections();
            // saturateddotpwm();
            // analogWrite(gimbal1_PWMpin, ddotpwm(0));
            // analogWrite(gimbal2_PWMpin, ddotpwm(1));
            // analogWrite(gimbal3_PWMpin, ddotpwm(2));
            // analogWrite(gimbal4_PWMpin, ddotpwm(3));

            d1 = delta1;
            d2 = delta2;
            d3 = delta3;
            d4 = delta4;

            gimbalangles = {d1, d2, d3, d4};
            if (printgimbals == 1)
            {
                printgimbalangles();
            }


            dur = dur + 1;


            // calculatepwm2rpm(); //mono gia otan thelw na dw pws syndeontai ta pwm me ta rpm.

            if (dur >= overflowlimit)
            {
                overflow = 1;
            }


            if (writetoSD == 1)
            {
                writedatatoSD();
            }



        }

    }





}
