/*****************************************************************
  LSM9DS1_Basic_I2C.ino
  SFE_LSM9DS1 Library Simple Example Code - I2C Interface
  Jim Lindblom @ SparkFun Electronics
  Original Creation Date: April 30, 2015
  https://github.com/sparkfun/LSM9DS1_Breakout

  The LSM9DS1 is a versatile 9DOF sensor. It has a built-in
  accelerometer, gyroscope, and magnetometer. Very cool! Plus it
  functions over either SPI or I2C.

  This Arduino sketch is a demo of the simple side of the
  SFE_LSM9DS1 library. It'll demo the following:
  How to create a LSM9DS1 object, using a constructor (global
  variables section).
  How to use the begin() function of the LSM9DS1 class.
  How to read the gyroscope, accelerometer, and magnetometer
  using the readGryo(), readAccel(), readMag() functions and
  the gx, gy, gz, ax, ay, az, mx, my, and mz variables.
  How to calculate actual acceleration, rotation speed,
  magnetic field strength using the calcAccel(), calcGyro()
  and calcMag() functions.
  How to use the data from the LSM9DS1 to calculate
  orientation and heading.

  Hardware setup: This library supports communicating with the
  LSM9DS1 over either I2C or SPI. This example demonstrates how
  to use I2C. The pin-out is as follows:
  LSM9DS1 --------- Arduino
   SCL ---------- SCL (A5 on older 'Duinos')
   SDA ---------- SDA (A4 on older 'Duinos')
   VDD ------------- 3.3V
   GND ------------- GND
  (CSG, CSXM, SDOG, and SDOXM should all be pulled high.
  Jumpers on the breakout board will do this for you.)

  The LSM9DS1 has a maximum voltage of 3.6V. Make sure you power it
  off the 3.3V rail! I2C pins are open-drain, so you'll be
  (mostly) safe connecting the LSM9DS1's SCL and SDA pins
  directly to the Arduino.

  Development environment specifics:
  IDE: Arduino 1.6.3
  Hardware Platform: SparkFun Redboard
  LSM9DS1 Breakout Version: 1.0

  This code is beerware. If you see me (or any other SparkFun
  employee) at the local, and you've found our code helpful,
  please buy us a round!

  Distributed as-is; no warranty is given.
*****************************************************************/
// The SFE_LSM9DS1 library requires both Wire and SPI be
// included BEFORE including the 9DS1 library.
#include <Wire.h>
#include <SPI.h>
#include <SparkFunLSM9DS1.h>
#include <BasicLinearAlgebra.h>
using namespace BLA;

//////////////////////////
// LSM9DS1 Library Init //
//////////////////////////
// Use the LSM9DS1 class to create an object. [imu] can be
// named anything, we'll refer to that throught the sketch.
LSM9DS1 imu;

///////////////////////
// Example I2C Setup //
///////////////////////
// SDO_XM and SDO_G are both pulled high, so our addresses are:
// #define LSM9DS1_M  0x1E // Would be 0x1C if SDO_M is LOW
// #define LSM9DS1_AG 0x6B // Would be 0x6A if SDO_AG is LOW

////////////////////////////
// Sketch Output Settings //
////////////////////////////
//#define PRINT_CALCULATED
#define PRINT_RAW
#define PRINT_SPEED 50 //50 // 250 ms between prints
static unsigned long lastPrint = 0; // Keep track of print time

// Earth's magnetic field varies by location. Add or subtract
// a declination to get a more accurate heading. Calculate
// your's here:
// http://www.ngdc.noaa.gov/geomag-web/#declination
#define DECLINATION 4.56 // Declination (degrees) in Boulder, CO.

//Function definitions
void printGyro();
void printAccel();
void printMag();
void printAttitude(float ax, float ay, float az, float gx, float gy, float gz, float mx, float my, float mz);
void findbiases();

int ledpin = 13;
double rollc;
double pitchc;
double headingc;
double rollcprev;
double pitchcprev;
double headingcprev;
double a, factor;
double headingcfprev = 0;

double gxmean = 0;
double gymean = 0;
double gzmean = 0;
double axmean = 0;
double aymean = 0;
double azmean = 0;
double gxprev = 0;
double gyprev = 0;
double gzprev = 0;
double axprev = 0;
double ayprev = 0;
double azprev = 0;


struct biases
{
    double ax;
    double ay;
    double az;
    double gx;
    double gy;
    double gz;
    double mx;
    double my;
    double mz;
};

typedef struct biases Struct;

struct biases bias;

void setup()
{

    rollcprev = 0;
    pitchcprev = 0;
    headingcprev = 0;
    pinMode(ledpin, OUTPUT);
    digitalWrite(ledpin, HIGH);
    Serial.begin(115200);

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

    //  Serial.println("compass calibration started");
    //  delay(500);
    //  imu.calibrateMag(1);
    //  Serial.println("compass calibration finished");
    //  Serial.println(imu.mBiasRaw[0]);
    //  Serial.println(imu.mBiasRaw[1]);
    //  Serial.println(imu.mBiasRaw[2]);



    Serial.println("accel/gyro calibration started");
    delay(1000);
    imu.calibrate(1);
    //  find_ac_gyro_biases();
    Serial.println("accel/gyro Calibration fininshed");
    Serial.println("================================");
    delay(1000);


}

void loop()
{
    // Update the sensor values whenever new data is available
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




    if ((lastPrint + PRINT_SPEED) < millis())
    {

        // printGyro();  // Print "G: gx, gy, gz"
        // printAccel(); // Print "A: ax, ay, az"
        //    printMag();   // Print "M: mx, my, mz"
        // Print the heading and orientation for fun!
        // Call print attitude. The LSM9DS1's mag x and y
        // axes are opposite to the accelerometer, so my, mx are
        // substituted for each other.

        printAttitude(imu.ax, imu.ay, imu.az, imu.gx, imu.gy, imu.gz,
                      imu.mx, imu.my, imu.mz );


        //    Serial.print("M     ");
        //    Serial.print(imu.mx - mxav);
        //    Serial.print("      ,      ");
        //    Serial.print(imu.my - myav);
        //    Serial.print("      ,      ");
        //    Serial.println(imu.mz - mzav);


        lastPrint = millis(); // Update lastPrint time
    }
}

void printGyro()
{
    // Now we can use the gx, gy, and gz variables as we please.
    // Either print them as raw ADC values, or calculated in DPS.
    Serial.print("G: ");
#ifdef PRINT_CALCULATED
    // If you want to print calculated values, you can use the
    // calcGyro helper function to convert a raw ADC value to
    // DPS. Give the function the value that you want to convert.
    Serial.print(imu.calcGyro(imu.gx), 2);
    Serial.print(", ");
    Serial.print(imu.calcGyro(imu.gy), 2);
    Serial.print(", ");
    Serial.print(imu.calcGyro(imu.gz), 2);
    Serial.println(" deg/s");
#elif defined PRINT_RAW
    Serial.print(imu.gx);
    Serial.print(", ");
    Serial.print(imu.gy);
    Serial.print(", ");
    Serial.println(imu.gz);
#endif
}

void printAccel()
{
    // Now we can use the ax, ay, and az variables as we please.
    // Either print them as raw ADC values, or calculated in g's.
    Serial.print("A: ");
#ifdef PRINT_CALCULATED
    // If you want to print calculated values, you can use the
    // calcAccel helper function to convert a raw ADC value to
    // g's. Give the function the value that you want to convert.
    Serial.print(imu.calcAccel(imu.ax), 2);
    Serial.print(", ");
    Serial.print(imu.calcAccel(imu.ay), 2);
    Serial.print(", ");
    Serial.print(imu.calcAccel(imu.az), 2);
    Serial.println(" g");
#elif defined PRINT_RAW
    Serial.print(imu.ax);
    Serial.print(", ");
    Serial.print(imu.ay);
    Serial.print(", ");
    Serial.println(imu.az);
#endif

}

void printMag()
{
    // Now we can use the mx, my, and mz variables as we please.
    // Either print them as raw ADC values, or calculated in Gauss.
    Serial.print("M: ");
#ifdef PRINT_CALCULATED
    // If you want to print calculated values, you can use the
    // calcMag helper function to convert a raw ADC value to
    // Gauss. Give the function the value that you want to convert.
    Serial.print(imu.calcMag(imu.mx), 2);
    Serial.print(", ");
    Serial.print(imu.calcMag(imu.my), 2);
    Serial.print(", ");
    Serial.print(imu.calcMag(imu.mz), 2);
    Serial.println(" gauss");
#elif defined PRINT_RAW
    Serial.print(imu.mx);
    Serial.print(", ");
    Serial.print(imu.my);
    Serial.print(", ");
    Serial.println(imu.mz);
#endif
}


void find_ac_gyro_biases()
{

    int sample = 1000;

    for (int i = 0; i < sample; ++i)
    {
        if ( imu.gyroAvailable() )
        {
            imu.readGyro();
            gxmean = gxmean + imu.gx;
            gymean = gymean + imu.gy;
            gzmean = gzmean + imu.gz;
        }
    }

    gxmean = gxmean / sample;
    gymean = gymean / sample;
    gzmean = gzmean / sample;

}

// Calculate pitch, roll, and heading.
// Pitch/roll calculations take from this app note:
// http://cache.freescale.com/files/sensors/doc/app_note/AN3461.pdf?fpsp=1
// Heading calculations taken from this app note:
// http://www51.honeywell.com/aero/common/documents/myaerospacecatalog-documents/Defense_Brochures-documents/Magnetic__Literature_Application_notes-documents/AN203_Compass_Heading_Using_Magnetometers.pdf
void printAttitude(float ax, float ay, float az, float gx, float gy, float gz, float mx, float my, float mz)
{


    // gx = gx - gxmean;
    // gy = gy - gymean;
    // gz = gz - gzmean;


    // mx = mx -  1357.7;
    // my = my -   635.0;
    // mz = mz - (-14.9);
    // float M11 =  1.000;
    // float M12 =  0.0291;
    // float M13 = -0.0002;
    // float M21 =  0.0291;
    // float M22 =  0.9665;
    // float M23 =  0.0088;
    // float M31 = -0.0002;
    // float M32 =  0.0088;
    // float M33 =  1.0516;
    mx = mx - (1.4462 * 1000);
    my = my - (0.3519 * 1000);
    mz = mz - (0.7036 * 1000);
    Matrix<3, 3> Mc = {  1.0000,   0.0331,  -0.0086,   0.0331,  0.9803,   0.0223, -0.0086,  0.0223, 1.0126 };
    mx =  Mc(0, 0) * mx +  Mc(1, 0) * my +  Mc(2, 0) * mz;
    my =  Mc(0, 1) * mx +  Mc(1, 1) * my +  Mc(2, 1) * mz;
    mz =  Mc(0, 2) * mx +  Mc(1, 2) * my +  Mc(2, 2) * mz;



    float  fact_gyro = 0.1;
    gx = (1 - fact_gyro) * gxprev + fact_gyro * gx;
    gy = (1 - fact_gyro) * gyprev + fact_gyro * gy;
    gz = (1 - fact_gyro) * gzprev + fact_gyro * gz;
    gxprev = gx;
    gyprev = gy;
    gzprev = gz;

    float fact_acc = 0.6;
    ax = (1 - fact_acc) * axprev + fact_acc * ax;
    ay = (1 - fact_acc) * ayprev + fact_acc * ay;
    az = (1 - fact_acc) * azprev + fact_acc * az;
    axprev = ax;
    ayprev = ay;
    azprev = az;



    // Serial.print("M     ");
    // Serial.print(mx);
    // Serial.print("      ,      ");
    // Serial.print(my);
    // Serial.print("      ,      ");
    // Serial.println(mz);


    float roll = atan2(ay, az);
    float pitch = atan2(-ax, sqrt(ay * ay + az * az));


    float heading;
    // Convert everything from radians to degrees:
    //heading *= 180.0 / PI;
    pitch *= 180.0 / PI;
    roll  *= 180.0 / PI;


    a = 0.3;
    factor = 0.1;
    rollc = (1 - a) * (rollcprev + factor * 0.001 * PRINT_SPEED * imu.calcGyro(gx)) + a * roll;
    pitchc = (1 - a) * (pitchcprev + factor * 0.001 * PRINT_SPEED * imu.calcGyro(gy)) + a * pitch;

    heading  = atan2( (-my * cos(rollc * PI / 180) + mz * sin(rollc * PI / 180) ), (-mx * cos(pitchc * PI / 180) + my * sin(pitchc * PI / 180) * sin(rollc * PI / 180) + mz * sin(pitchc * PI / 180) * cos(rollc * PI / 180)) );  //the same as aboive but with -mx
    heading *= 180.0 / PI;
    heading -= DECLINATION * PI / 180;
    if (heading > PI) heading -= (2 * PI);
    else if (heading < -PI) heading += (2 * PI);

    headingc = (1 - a) * (headingcprev + factor * 0.001 * PRINT_SPEED * imu.calcGyro(gz)) + a * heading;


    rollcprev = rollc;
    pitchcprev = pitchc;
    headingcprev = headingc;

    Serial.print("N   ");
    Serial.print(roll, 2);
    Serial.print("  ,  ");
    Serial.print(-pitch, 2);
    Serial.print("  ,  ");
    Serial.println(-heading, 2);

    Serial.print("C   ");
    Serial.print(rollc, 2);
    Serial.print("  ,  ");
    Serial.print(-pitch, 2);
    Serial.print("  ,  ");
    Serial.println(-headingc, 2);

}
