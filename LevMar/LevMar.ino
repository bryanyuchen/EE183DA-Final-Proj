#include <MatrixMath.h>

/********************************************************************/
/*
 * Global Variables
 */

// 2D Microphone Array
#define N  (2)
float c = 34.3;
uint8_t mdist = 20;
float m1[N]; // vector of mic 1
float m2[N]; // vector of mic 2
float m3[N]; // vector of mic 3

// LevMar Parameters
uint8_t n_iters = 10; // # of iterations for LM
float lambda= 0.01; // initial damping factor 
uint8_t updateJ=1; //update?

//LevMar Variables
float state_est[N] = {1,1}; //state estimate (maybe can be optimized)
float state_lm[N]; //Lev_Marstate estimate
float dp[N]; //step
float y_est [3]; //y estimate, derived from state_est
float y_est_lm [3]; //y estimate, derived from LevMar
float d[3]; //Residual
float d_lm[3]; //Residual after LevMar step
float e; //Error (residual^2)
float e_lm; //Error (residual^2) after LevMar step
float J [3][2]; //Jacobian
float J_transpose [2][3]; //Transpose of Jacobian
float temp[2]; //temp matrix to help J'*d
float H[2][2]; //Hessian
float H_lm[2][2];
float res1[N]; //result vector for Jacobian math
float res2[N]; //result vector for Jacobian math
float res3[N]; //result vector for Jacobian math
float eye[N][N];

/********************************************************************/
/*
 * Setup
 */
void setup() {
Serial.begin(9600);      // open the serial port at 9600 bps:  

// 2D Microphone Array 
m1[0] = 0;
m1[1] = 0;
m2[0] = mdist*-0.5f;
m2[1] = mdist*-0.8660254f;
m3[0] = mdist*0.5f;
m3[1] = mdist*-0.8660254f;

// identity matrix init
eye[0][0] = 1;
eye[0][1] = 0;
eye[1][0] = 0;
eye[1][1] = 1;
}

/********************************************************************/
/*
 * Norm Function
 */
float norm(float* vector, uint8_t vector_size) {
  float total = 0;
  for (int i = 0; i < vector_size; i++) {
    total += (vector[i]*vector[i]);
  }
  return sqrt(total);
}

/********************************************************************/
/*
 * Levenberg-Marquardt Algorithm
 */
void LevMar(float* measurements) {
for (uint8_t it = 1; it <= n_iters; it++) {
  Serial.println("starting levmar");
  // if updated, calculate new parameters
  if (updateJ == 1) {
    // Evaluate the Jacobian matrix at the current parameters 
    Matrix.Subtract((float*)state_est,(float*)m1,N,1,(float*)res1);
    Matrix.Subtract((float*)state_est,(float*)m2,N,1,(float*)res2);
    Matrix.Subtract((float*)state_est,(float*)m3,N,1,(float*)res3);

    J[0][0] = (1/c) * ((state_est[0]-m1[0])/norm(res1,N) - (state_est[0]-m2[0])/norm(res2,N));
    J[0][1] = (1/c) * ((state_est[1]-m1[1])/norm(res1,N) - (state_est[1]-m2[1])/norm(res2,N));
    J[1][0] = (1/c) * ((state_est[0]-m2[0])/norm(res2,N) - (state_est[0]-m3[0])/norm(res3,N));
    J[1][1] = (1/c) * ((state_est[1]-m2[1])/norm(res2,N) - (state_est[1]-m3[1])/norm(res3,N));
    J[2][0] = (1/c) * ((state_est[0]-m1[0])/norm(res1,N) - (state_est[0]-m3[0])/norm(res3,N));
    J[2][1] = (1/c) * ((state_est[1]-m1[1])/norm(res1,N) - (state_est[1]-m3[1])/norm(res3,N));

    // Evaluate the distance error at the current parameters
    y_est[0] = (1/c) * (norm(res1,N) - norm(res2,N));
    y_est[1] = (1/c) * (norm(res2,N) - norm(res3,N));
    y_est[2] = (1/c) * (norm(res1,N) - norm(res3,N));

    // Evaluate residual matrix
    Matrix.Subtract(measurements,(float*)y_est,3,1,(float*)d);
  
    // compute the approximated Hessian matrix, J’ is the transpose of J
    Matrix.Transpose((float*)J,3,2,(float*)J_transpose);
    Matrix.Multiply((float*)J_transpose,(float*)J,2,3,2,(float*)H);

    //if the first iteration, compute the total error    
    if (it==1) { 
      //slight modification to error calculation using abs() to save memory and improve precision
      e = (d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
      Serial.print("e = ");
      Serial.println(e,7);
    }
  }

  // Apply the damping factor to the Hessian matrix
  eye[0][0] = lambda;
  eye[1][1] = lambda;
  Matrix.Add((float*)H,(float*)eye,N,N,(float*)H_lm);
  
  // Compute the updated parameters
  Matrix.Invert((float*)H_lm,2);
  Matrix.Multiply((float*)J_transpose,(float*)d,2,3,1,(float*)temp);
  Matrix.Multiply((float*)H_lm,(float*)temp,2,2,1,(float*)dp);

  Matrix.Add((float*)state_est,(float*)dp,2,1,(float*)state_lm);
  
  //going towards direction of wrong minima
  if (state_lm[1] < 1.0f) {
      state_lm[1] = state_est[1]-(10*dp[1]);
  }
  
  //Evaluate the total distance error at the updated parameters
  Matrix.Subtract((float*)state_lm,(float*)m1,N,1,(float*)res1);
  Matrix.Subtract((float*)state_lm,(float*)m2,N,1,(float*)res2);
  Matrix.Subtract((float*)state_lm,(float*)m3,N,1,(float*)res3);
      
  y_est_lm[0] = (1/c) * (norm(res1,N) - norm(res2,N));
  y_est_lm[1] = (1/c) * (norm(res2,N) - norm(res3,N));
  y_est_lm[2] = (1/c) * (norm(res1,N) - norm(res3,N));
  
  Matrix.Subtract((float*)measurements,(float*)y_est_lm,3,1,d_lm);
  e_lm = (d_lm[0]*d_lm[0] + d_lm[1]*d_lm[1] + d_lm[2]*d_lm[2]);
  Serial.print("e_lm = ");
  Serial.println(e_lm,7);
  Matrix.Print((float*)state_est,2,1,"state_est");
  
  //check total distance error
  if (e_lm < e) { 
   Serial.println("update step");
   lambda /= 10;
   state_est[0] = state_lm[0];
   state_est[1] = state_lm[1];
   e = e_lm;
   Serial.print("error: ");
   Serial.println(e,7);
   updateJ = 1;
  }
  else {
    Serial.println("increase damping");
    updateJ = 0;
    lambda *= 10;
  }
}
return;
}
/********************************************************************/

void loop() {

float measurements[3] = {-0.577883,0.117166,-0.460717};

LevMar(measurements);
Matrix.Print((float*)state_est,2,1,"LevMar State Estimation");
  delay(100000);
}
