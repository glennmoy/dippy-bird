/*********************************************************
This is a Dippy Bird Simulator similar to my Physics FYP

*********************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>



//Function to calculate dw/dt
double update_omega(double t,double theta, double omega, double alpha, double beta){

	double k;
	k=-sin(theta)-beta*omega+alpha*t;
	return k;
}



void main(int argc, char* argv[]){

	//argv are going to be the initial conditions
	//and the parameters alpha, beta, gamma

	
	double t=0.0	//absolute time
	double t0=0.0;	//time since last dip
	double dt=0.01;	 //time-step
	
	//input parameters
	double alpha=0.01; //relative strength of driving force
	double beta=0.08;  //relative strength of damping
	double t_max=2000; //max time

	//malloc theta, omega arrays
	double* theta=malloc(sizeof(double)*(double)N/dt);
	double* omega=malloc(sizeof(double)*(double)N/dt);

	//initial conditions
	theta[0]=M_PI/4;
	omega[0]=0;


	printf("%g\t%.3f\t%.3f\t%.3f\n",t,theta[0],omega[0],0.0);

	int i=1;		//array index
	double k1,k2,k3,k4;	//runge kutta parameters
	
	while(t<t_max){
		
		//update runge kutta values
		k1=update_omega(t0,theta[i-1],omega[i-1],alpha,beta);
		k2=update_omega(t0+0.5*dt,theta[i-1],omega[i-1]+0.5*dt*k1,alpha,beta);
		k2=update_omega(t0+0.5*dt,theta[i-1],omega[i-1]+0.5*dt*k2,alpha,beta);
		k3=update_omega(t0+0.5*dt,theta[i-1],omega[i-1]+0.5*dt*k3,alpha,beta);

		//solve eoms
		omega[i]=omega[i-1]+0.1666*dt*(k1+2*k2+2*k3+k4);
		theta[i]=theta[i-1]+omega[i]*dt;		
	
		printf("%g\t%.3f\t%.3f\t%.3f\n",t,theta[i],omega[i],t0);

		//update time and index
		t+=dt; t0+=dt;
		i+=1;
		
		//check if theta has gone overboard
		if(theta[i]>max_angle){
			t0=0;	
		}
	}

	free(theta);
	free(omega);
} 
