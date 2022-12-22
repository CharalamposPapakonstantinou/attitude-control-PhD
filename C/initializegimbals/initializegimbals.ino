

//INTERRUPT COUNTERS
int inter_pin_1 = 7;
int inter_pin_2 = 6;
int  PWMA = 2; //enable
int  GPIOA = 0; //phase
int  PWMB = 3; //enable
int  GPIOB = 1; //phase
int  photores = A21;
int orangeled = 9;
int greenled = 8;
float photores_threshold;

int led = 13;

volatile int count1 = 0; //left wheel

int countprev = 0;
unsigned long previousTime = 0, currentTime = 0; //
double LOOP_FREQUENCY = 1000;
double rpms;
double dt;


//left wheel int 0  SERVICE ROUTINE
void ISR_1 ()      //also avail waveform on pin 1
{
    if (digitalRead(inter_pin_1) == HIGH)
    {
        if (digitalRead(inter_pin_2) == LOW)    //1
        {
            count1++;
        }
        else
        {
            count1--;
        }
    }
    else
    {
        if (digitalRead(inter_pin_2) == HIGH)    //1
        {
            count1++;
        }
        else
        {
            count1--;
        }
    }
}

void ISR_2 ()      //also avail waveform on pin 1
{
    if (digitalRead(inter_pin_2) == HIGH)
    {
        if (digitalRead(inter_pin_1) == HIGH)    //1
        {
            count1++;
        }
        else
        {
            count1--;
        }
    }
    else
    {
        if (digitalRead(inter_pin_1) == LOW)    //1
        {
            count1++;
        }
        else
        {
            count1--;
        }
    }
}

float calibrate_photoresistor()
{
    int val = 0;
    int samples = 200;
    for (int i = 0; i < samples; ++i)
    {
        val = val + analogRead(photores);
    }
    val = val / samples;
    photores_threshold = val / 2;
}

void initializegimbalangles()
{
    digitalWrite(orangeled, HIGH);
    if ((count1 - countprev) > 0)
    {
        digitalWrite(GPIOA, LOW);
        analogWrite(PWMA, 70 * (count1 - countprev));
    }
    else if ((count1 - countprev) < 0)
    {
        digitalWrite(GPIOA, HIGH);
        analogWrite(PWMA, 70 * abs(count1 - countprev));
    }
    else
    {
        analogWrite(PWMA, 0);
    }
    countprev = count1;
    delay(50);
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

//////////////////////////////////////////////////////    SETUP
void setup()
{
    Serial.begin(115200); //Start Serial at 115200bps
    pinMode(led,  OUTPUT);
    pinMode(PWMA, OUTPUT);
    pinMode(GPIOA, OUTPUT);
    pinMode(PWMB, OUTPUT);
    pinMode(GPIOB, OUTPUT);
    pinMode(photores, INPUT);
    pinMode(orangeled,  OUTPUT);
    pinMode(greenled,  OUTPUT);

    //interrupt-counting_pins

    pinMode(inter_pin_1, INPUT);       //SOS diff from ARDUINO , NEEDs to be declared as INPUT before attachInterrupt !!!
    pinMode(inter_pin_2, INPUT);

    attachInterrupt( digitalPinToInterrupt(inter_pin_1), ISR_1, CHANGE);     //left wheel
    attachInterrupt( digitalPinToInterrupt(inter_pin_2), ISR_2, CHANGE);   //right wheel

    //variable initialization
    count1 = 0;

    digitalWrite(led, HIGH) ;

    dt = LOOP_FREQUENCY * 0.001;

    //////////////////////////////////////////////////////////////
    calibrate_photoresistor();

    do
    {
        initializegimbalangles();
        // Serial.println(photoresistortouched());
    }
    while (photoresistortouched() == 0);
}



int samples = 1;
float sumrpms = 0;
int pwmval = 50;
int c = 1;

void loop()
{


}

















