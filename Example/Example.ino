#include <LevMar.h>
#include <MatrixMath.h>
LevMar levmar(20.0f);

void setup() {
Serial.begin(9600); 
}

void loop() {
float measurements[3] = {-0.577883,0.117166,-0.460717};
float* state_est = levmar.Run(measurements);
Matrix.Print((float*)state_est,2,1,"LevMar State Estimation");
delay(100000);
}
