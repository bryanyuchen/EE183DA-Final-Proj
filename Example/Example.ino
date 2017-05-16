#include <Kalmf.h>
#include <LevMar.h>
#include <MatrixMath.h>
LevMar levmar(20.0f);
Kalmf kalmf;

//measurement array
#define NUM_MEASUREMENTS (20)
float measurements[NUM_MEASUREMENTS][3] =
    {
        {-0.557700,   0.028600,  -0.543400},
        {-0.557700,  -0.028600,  -0.543400},
        {-0.572000,   0.028600,  -0.557700},
        {-0.572000,   0.028600,  -0.543400},
        {-0.529100,   0.028600,  -0.529100},
        {-0.557700,   0.028600,  -0.529100},
        {-0.614900,   0.028600,  -0.529100},
        {-0.543400,  -0.028600,  -0.572000},
        {-0.557700,  -0.042900,  -0.572000},
        {-0.557700,  -0.028600,  -0.557700},
        {-0.529100,  -0.042900,  -0.543400},
        {-0.572000,  -0.057200,  -0.572000},
        {-0.543400,  -0.042900,  -0.557700},
        {-0.543400,   0.028600,  -0.529100},
        {-0.529100,  -0.042900,  -0.557700},
        {-0.514800,  -0.042900,  -0.572000},
        {-0.543400,  -0.042900,  -0.557700},
        {-0.529100,  -0.057200,  -0.572000},
        {-0.543400,  -0.028600,  -0.557700},
        {-0.529100,   0.028600,  -0.529100}
    };
    
void setup() {
 Serial.begin(9600); 
}

void loop() {
//float measurements[3] = {-0.577883,0.117166,-0.460717}; // ideal measurement

for (int i = 3; i < NUM_MEASUREMENTS; i++) {
//solves nonlinear LS problem
float* state_est = levmar.Run(measurements[i]);
Matrix.Print((float*)state_est,2,1,"LevMar State Estimation");

//estimates state evolution based on levmar solution and state memory
float* Kalmf_return_state = kalmf.Run(state_est);
Matrix.Print((float*)Kalmf_return_state,2,1,"Kalmf State Estimation"); 
}
delay(100000);
}

