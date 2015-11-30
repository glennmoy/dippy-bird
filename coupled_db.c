/*********************************************************
This is a Dippy Bird Simulator similar to my Physics FYP

*********************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>


//Function to calculate dw/dt
double update_omega(double t,double theta, double omega, double alpha, double beta, double gamma, double theta1){

	double k;
	double d,x,y,x1,y1,phi;
	
	//assuming l1=l2=1
	x=sin(theta);
	y=cos(theta);
	x1=2+sin(theta1); 
	y1=sin(theta1);

	d=sqrt(pow(x-x1,2)+pow(y-y1,2));
	phi=atan((y1-y)/(x1-x));

	k=-sin(theta)-beta*omega+alpha*t+gamma*d*cos(theta-phi);
	return k;
}



int main(int argc, char* argv[]){

	//argv are going to be the initial conditions
	//and the parameters alpha, beta, gamma
	
	//Also assuming for now parameters are equal
	//for both birds and l1=l2=1

	
	double t=0.0;	//absolute time
	double t0=0.0;	//time since last dip
	double t1=0.0;
	double dt=0.01;	 //time-step
	
	//input parameters
	double alpha=0.01; //relative strength of driving force
	double beta=0.09;  //relative strength of damping
	double gamma=0.05;
	double t_max=200.0; //max time
	double max_angle=M_PI/3;

	//malloc theta, omega arrays
	double* theta=malloc(sizeof(double)*t_max/dt);
	double* omega=malloc(sizeof(double)*t_max/dt);
	double* theta1=malloc(sizeof(double)*t_max/dt);
	double* omega1=malloc(sizeof(double)*t_max/dt);

	printf("%g\n\n\n",t_max/dt);

	//initial conditions
	theta[0]=0.4;
	omega[0]=-1;
	theta1[0]=1.2;
	omega1[0]=0;

	printf("%g\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n",t,theta[0],omega[0],theta1[0],omega1[0],t0,t1);


	int i=1;		//array index
	double k1,k2,k3,k4;	//runge kutta parameters
	
	while(t<t_max){
		
		//update runge kutta values
		k1=update_omega(t0,theta[i-1],omega[i-1],alpha,beta,gamma,theta1[i-1]);
		k2=update_omega(t0+0.5*dt,theta[i-1],omega[i-1]+0.5*dt*k1,alpha,beta,gamma,theta1[i-1]);
		k2=update_omega(t0+0.5*dt,theta[i-1],omega[i-1]+0.5*dt*k2,alpha,beta,gamma,theta1[i-1]);
		k3=update_omega(t0+0.5*dt,theta[i-1],omega[i-1]+0.5*dt*k3,alpha,beta,gamma,theta1[i-1]);

		//solve eoms
		omega[i]=omega[i-1]+0.1666*dt*(k1+2*k2+2*k3+k4);
		theta[i]=theta[i-1]+omega[i]*dt;		

		
		k1=update_omega(t1,theta1[i-1],omega1[i-1],alpha,beta,gamma,theta[i-1]);
		k2=update_omega(t1+0.5*dt,theta1[i-1],omega1[i-1]+0.5*dt*k1,alpha,beta,gamma,theta[i-1]);
		k2=update_omega(t1+0.5*dt,theta1[i-1],omega1[i-1]+0.5*dt*k2,alpha,beta,gamma,theta[i-1]);
		k3=update_omega(t1+0.5*dt,theta1[i-1],omega1[i-1]+0.5*dt*k3,alpha,beta,gamma,theta[i-1]);

		omega1[i]=omega1[i-1]+0.1666*dt*(k1+2*k2+2*k3+k4);
		theta1[i]=theta1[i-1]+omega1[i]*dt;		

	
		printf("%g\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n",t,theta[i],omega[i],theta1[i],omega1[i],t0,t1);

		//update time and index
		t+=dt;	t0+=dt;  t1+=dt;
		//check if theta has gone overboard
		if(theta[i]>max_angle){
			t0=0;		//if not update t0
		}
		if(theta1[i]>max_angle){
			t1=0;		//else set to 0
		}
		
		i+=1;
	}

	free(theta); free(theta1);
	free(omega); free(omega1);
} 
