#include "math.h"
#include "mex.h"
#include <conio.h>
#include "stdlib.h "
const double PI = 3.1416;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	//*******Variable declarations****//
	double *T;
	int Tm;
	int Tn;
	double z;
	double *p;
	double gamma;
	double beta;
	double Dt;
	mxArray *tData;
	double *zData;
	mxArray *pData;
	double *gammaData;
	double *betaData;
	double *DtData;
	double *outSd;
    double *outd;
    double *outv;
    double *outa;
    
    double dAbs;
	double m;
	double *w;
	double *c;
	double *k;
	int lengthT;
	int lengthP;
	int pm,pn;
	int i,j;
	double *d;
	double *v;
	double *a;
	double *fs;
    double *Sd;
    
    double maxD;
	double *A;
	double *B;
	int index;
    int next;
	double DPi;
	double ki;
	double Ki;
	double Ddi;

	double Df;
	double DR;

	double Dvi;
	///////////////////////////////////

	//*******Get data from MATLAB***************//
    
	tData = prhs[0];
	zData = mxGetPr(prhs[1]);
	pData = prhs[2];
	gammaData = mxGetPr(prhs[3]);
	betaData = mxGetPr(prhs[4]);
	DtData = mxGetPr(prhs[5]);

	T = mxGetPr(tData);
	Tm = mxGetM(tData); // dimensions of T
	Tn = mxGetN(tData);
	lengthT = Tm>=Tn?Tm:Tn;
	z = zData[0];
	p = mxGetPr(pData);
	pm = mxGetM(pData); // dimensions of p
	pn = mxGetN(pData);
	lengthP = pm>=pn?pm:pn; //P can be row or column vector
	gamma = gammaData[0];
	beta = betaData[0];
	Dt = DtData[0];
	////////////////////////////////
	
	//********Main Program********//

	//setup
	m = 1;
    // Memory allocations
	w = (double*)malloc(lengthT*sizeof(double));
	c = (double*)malloc(lengthT*sizeof(double));
	k = (double*)malloc(lengthT*sizeof(double));
	A = (double*)malloc(lengthT*sizeof(double));
	B = (double*)malloc(lengthT*sizeof(double));

	for(i = 0 ; i < lengthT ; i++){
		w[i] = 2*PI/T[i];
		c[i] = z*2*m*w[i];
		k[i] = (w[i]*w[i])*m;
		A[i] = 1/(beta*Dt)*m+gamma/beta*c[i];
		B[i] = 1/(2*beta)*m + Dt*(gamma/(2*beta)-1)*c[i];
	}

	d = (double*)malloc(lengthP*sizeof(double));
	v = (double*)malloc(lengthP*sizeof(double));
	a = (double*)malloc(lengthP*sizeof(double));
	fs = (double*)malloc(lengthP*sizeof(double));
    Sd = (double*)malloc(lengthT*sizeof(double));
    
	for(i = 0; i<lengthT; i++){
		a[i] = (p[0])/m; // a[0] = (p-c*v[0]-fs[0])/m
		v[i] = 0;
		d[i] = 0;
		fs[i] = 0;
	}

	//*************Time Stepping ************//

	for(j = 0; j<lengthT; j++){
        maxD = 0.0;
		for(i = 0; i<lengthP-1; i++){
			index = i;
            next =(i+1);
			DPi = p[i+1]-p[i]+A[j]*v[index]+B[j]*a[index];
			ki = k[j];
			Ki = ki + A[j]/Dt; // = ki + GAMMA/(BETA*Dt)*c + 1/(BETA*Dt^2)*m
            Ddi = DPi/Ki;
            fs[next] = fs[index] + ki*Ddi;
            
			d[next] = d[index] + Ddi;
            
			Dvi = gamma/(beta*Dt)*Ddi - gamma/beta*v[index] + Dt*(1-gamma/(2*beta))*a[index];
			v[next] = v[index] + Dvi;
			//a[next] = (p[i+1] - c[j]*v[next] - fs[next])/m;
            a[next] = a[index] + 1/(beta*Dt*Dt)*Ddi - 1/(beta*Dt)*v[index] - 1/(2*beta)*a[index];
            //printf("d = %f maxD = %f \n",abs(d[next]),maxD);
            if(d[next]<0){
                dAbs = -1.0*d[next];
            }else{
                dAbs = d[next];
            }
            if(dAbs > maxD){
                maxD = dAbs;
            }
		}
        Sd[j] = maxD;
	}


	///////////////////////////////

	//SendData
	plhs[0] = mxCreateDoubleMatrix(lengthT,1,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(lengthP,1,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(lengthP,1,mxREAL);
    plhs[3] = mxCreateDoubleMatrix(lengthP,1,mxREAL);
    
	outSd = mxGetPr(plhs[0]);
    outd = mxGetPr(plhs[1]);
    outv = mxGetPr(plhs[2]);
    outa = mxGetPr(plhs[3]);
    
    
	memcpy(outSd,Sd,lengthT*sizeof(double));
    memcpy(outd,d,lengthP*sizeof(double));
    memcpy(outv,v,lengthP*sizeof(double));
    memcpy(outa,a,lengthP*sizeof(double));
	
	free(w);
	free(c);
	free(k);
	free(A);
	free(B);
	free(d);
	free(v);
	free(a);
	free(fs);
}