/*
Author Mike van der Naald, July 18, 2016:
* This program implements algorithms to classify ellipses.  Completely 
*taken from the paper "Calculating ellipse Overlap Areas."  By 
* Gary Hughes and Mohcine Chraibi.
*/

#include <iostream>
#include <string>
#include <math.h>
#include<stdlib.h>
#include <vector>
#include <complex>
#include <boost/numeric/ublas/matrix.hpp>



/* These two functions will be needed later.  They're just linear algebra functions I didn't want to 
 * import another library. */
std::complex<double> twoByTwoDet( boost::numeric::ublas::matrix< std::complex<double> > & matrix )
{
	return matrix(0,0)*matrix(1,1)-matrix(0,1)*matrix(1,0);
}

std::complex<double> threeByThreeDet( boost::numeric::ublas::matrix< std::complex<double> > & matrix )
{
	return matrix(0,0)*matrix(1,1)*matrix(2,2)+matrix(0,1)*matrix(1,2)*matrix(2,0)+
	matrix(0,2)*matrix(1,0)*matrix(2,1)-matrix(0,2)*matrix(1,1)*matrix(2,0)-
	matrix(0,1)*matrix(1,0)*matrix(2,2)-matrix(0,0)*matrix(1,2)*matrix(2,1);
}




/*This calculates coefficients that will be needed in the case finder below */
std::vector<std::complex<double>> coefficientFinder(double phi,double A,double B,double h,double k)
{
	std::vector<std::complex<double>> coefficientHolder (6) ;
	coefficientHolder[0] = (cos(phi)*cos(phi))/(A*A)+(sin(phi)*sin(phi))/(B*B);
	coefficientHolder[1] = 2*sin(phi)*cos(phi)/(A*A)-2*sin(phi)*cos(phi)/(B*B);
	coefficientHolder[2] = (sin(phi)*sin(phi))/(A*A)+(cos(phi)*cos(phi))/(B*B);
	coefficientHolder[3] = -2*cos(phi)*(h*cos(phi)+k*sin(phi))/(A*A)+2*sin(phi)*(k*cos(phi)-h*sin(phi))/(B*B);
	coefficientHolder[4] = -2*sin(phi)*(h*cos(phi)+k*sin(phi))/(A*A)+2*cos(phi)*(h*sin(phi)-k*cos(phi))/(B*B);
	coefficientHolder[5] = (h*cos(phi)+k*cos(phi))*(h*cos(phi)+k*cos(phi))/(A*A) + (h*sin(phi)-k*sin(phi))*(h*sin(phi)-k*sin(phi))/(B*B) - 1;
	return coefficientHolder;
}


int ellipseOverlapCaseFinder(std::complex<double> phi0,std::complex<double> A0,std::complex<double> B0,std::complex<double> h0,std::complex<double> k0,std::complex<double> phi1,std::complex<double> A1,std::complex<double> B1,std::complex<double> h1,std::complex<double> k1)
{
	/* These are the coefficients AA BB CC DD EE FF for each ellipse. */
	std::vector<std::complex<double>> ellipseOne = coefficientFinder(phi0,A0,B0,h0,k0);
	std::vector<std::complex<double>> ellipseTwo = coefficientFinder(phi1,A1,B1,h1,k1);
	
	/* NOTE: ellipseOne[0]=AA,ellipseOne[1]=BB,ellipseOne[2]=CC,ellipseOne[3]=DD, ellipseOne[4]=EE, */
	
	std::complex<double> d = ellipseOne[0]*(ellipseOne[2]*ellipseOne[5]-ellipseOne[4]*ellipseOne[4])-
	(ellipseOne[2]*ellipseOne[3]*ellipseOne[3]
	-2*ellipseOne[1]*ellipseOne[3]*ellipseOne[4]+ellipseOne[5]*ellipseOne[1]*ellipseOne[1]);
	
	/* division is expensive and 1/d is used for a,b and c so I just keep it in memory. */
	std::complex<double> oneOverd = 1/d;
	 
	std::complex<double> a = oneOverd*(ellipseOne[0]*(ellipseOne[2]*ellipseTwo[5]-2*ellipseOne[4]*ellipseTwo[4]+ellipseOne[5]*ellipseTwo[2]) +
	2*ellipseOne[1]*(ellipseOne[4]*ellipseTwo[3]-ellipseOne[5]*ellipseTwo[1]+ellipseOne[3]*ellipseTwo[5])+2*ellipseOne[3]*(
	ellipseOne[4]*ellipseTwo[1]-ellipseOne[2]*ellipseTwo[3])-(ellipseOne[1]*ellipseOne[1]*ellipseTwo[5]+ellipseOne[3]*ellipseOne[3]*ellipseTwo[2]+
	ellipseOne[5]*ellipseOne[5]*ellipseTwo[0])+(ellipseOne[2]*ellipseOne[5]*ellipseTwo[0]));
	
	std::complex<double> b = oneOverd*(ellipseOne[0]*(ellipseTwo[2]*ellipseTwo[5]-ellipseTwo[4]*ellipseTwo[4])+2*ellipseOne[1]*
	(ellipseTwo[4]*ellipseTwo[3]-ellipseTwo[5]*ellipseTwo[1])+2*ellipseOne[3]*(ellipseTwo[4]*ellipseTwo[1]-ellipseTwo[2]*
	ellipseTwo[3])+ellipseOne[2]*(ellipseTwo[0]*ellipseTwo[5]-ellipseTwo[3]*ellipseTwo[3])+2*ellipseOne[4]*(ellipseTwo[1]*
	ellipseTwo[3]-ellipseTwo[0]*ellipseTwo[4])+ellipseOne[5]*(ellipseTwo[0]*ellipseTwo[2]-ellipseTwo[1]*ellipseTwo[1]));
	
	std::complex<double> c = 	oneOverd*(ellipseTwo[0]*(ellipseTwo[2]*ellipseTwo[5]-ellipseTwo[4]*ellipseTwo[4])-(ellipseTwo[1]*ellipseTwo[1]
	*ellipseTwo[5]-2*ellipseTwo[1]*ellipseTwo[3]*ellipseTwo[4]+ellipseTwo[3]*ellipseTwo[3]*ellipseTwo[2]));
	
	std::complex<double> s4 = -27*pow(c,3)+18*c*a*b+a*a*b*b-4*pow(a,3)*c-4*pow(b,3);
	
	
	/*We're almost done with the necessary constants we just have to craft matrices A and B 
	 * For these I'm going to use boosts linear algebra library.*/
	boost::numeric::ublas::matrix<complex<double>> aMatrix (3, 3);
	boost::numeric::ublas::matrix<complex<double>> bMatrix (3, 3);
	
	aMatrix (0,0) = ellipseOne[0];
	aMatrix (0,1) = .5*ellipseOne[1];
	aMatrix (1,0) = .5*ellipseOne[1];
	aMatrix (1,1) = ellipseOne[2];
	aMatrix (0,2) = .5*ellipseOne[3];
	aMatrix (2,0) = .5*ellipseOne[3];
	aMatrix (1,2) = .5*ellipseOne[4];
	aMatrix (2,1) = .5*ellipseOne[4];
	aMatrix (2,2) = ellipseOne[5];
	
	bMatrix (0,0) = ellipseTwo[0];
	bMatrix (0,1) = .5*ellipseTwo[1];
	bMatrix (1,0) = .5*ellipseTwo[1];
	bMatrix (1,1) = ellipseTwo[2];
	bMatrix (0,2) = .5*ellipseTwo[3];
	bMatrix (2,0) = .5*ellipseTwo[3];
	bMatrix (1,2) = .5*ellipseTwo[4];
	bMatrix (2,1) = .5*ellipseTwo[4];
	bMatrix (2,2) = ellipseTwo[5];
	
	
	
	
	
	
	
	
	int returner;
	
	if (s4>0)
		returner = 2;
	else if (s4<0)
		
		double s1 = a;
		double s2 = a*a-3*b;
		double s3 = 3*a*c+b*a*a-4*b*b;
		if ((s1>0)&&(s2>0)&&(s3>0))
			std::complex<double> u = -a-pow(s2,(1/2))/3;
			std::complex<double> v = -a+pow(s2,(1/2))/3;
			boost::numeric::ublas::matrix< std::complex<double> > mMatrix (3, 3);
			boost::numeric::ublas::matrix< std::complex<double> > nMatrix (3, 3);
			mMatrix = u*A+B;
			nMatrix = v*A+B;
			boost::numeric::ublas::matrix< std::complex<double> > mMatrix (3, 3);
			
			boost::numeric::ublas::matrix< std::complex<double> > m11minor (2, 2);
			boost::numeric::ublas::matrix< std::complex<double> > n11minor (2, 2);
			
			m11minor(0,0)=mMatrix(1,1);
			m11minor(1,0)=mMatrix(2,1);
			m11minor(0,1)=mMatrix(1,2);
			m11minor(1,1)=mMatrix(2,2);
			
			n11minor(0,0)=nMatrix(1,1);
			n11minor(1,0)=nMatrix(2,1);
			n11minor(0,1)=nMatrix(1,2);
			n11minor(1,1)=nMatrix(2,2);
			
			
			if ((mMatrix(1,1)*threeByThreeDet(mMatrix)>0)&&(twoByTwoDet(mTopMinor)>0))||((nMatrix(1,1)*threeByThreeDet(nMatrix)>0)&&(twoByTwoDet(nTopMinor)>0))
				returner = 4;
			else
				returner = 1;
		else
			returner = 3;
	else if (s4==0)
		double s1 = a;
		double s2 = a*a-3*b;
		double s3 = 3*a*c+b*a*a-4*b*b;
		
		if ((s1>0)&&(s2>0)&&(s3<0))
			returner =6;
		else if ((s1>0)&&(s2>0)&&(s3>0))
			std::complex<double> beta = (9*c-a*b)/(2*s2);
			std::complex<double> alpha = (4*a*b-a*a*a-9*c)/s2;
			boost::numeric::ublas::matrix< std::complex<double> > mMatrix (3, 3);
			boost::numeric::ublas::matrix< std::complex<double> > nMatrix (3, 3);
			
			M = beta*A+B;
			N = alpha*A+b;
			
			boost::numeric::ublas::matrix< std::complex<double> > N33Minor(2,2);
			boost::numeric::ublas::matrix< std::complex<double> > M33Minor(2,2);
			boost::numeric::ublas::matrix< std::complex<double> > M11Minor(2,2);
			boost::numeric::ublas::matrix< std::complex<double> > M22Minor(2,2); 
			
			N33Minor(0,0) = nMatrix(0,0);
			N33Minor(1,0) = nMatrix(1,0);
			N33Minor(0,1) = nMatrix(0,1);
			N33Minor(1,1) = nMatrix(1,1);
			
			M33Minor(0,0) = mMatrix(0,0);
			M33Minor(1,0) = mMatrix(1,0);
			M33Minor(0,1) = mMatrix(0,1);
			M33Minor(1,1) = mMatrix(1,1);
			
			M11Minor(0,0) = mMatrix(0,0);
			M11Minor(0,1) = mMatrix(0,1);
			M11Minor(1,0) = mMatrix(1,0);
			M11Minor(1,1) = mMatrix(1,1); 
			
			M22Minor(0,0) = mMatrix(0,0);
			M22Minor(1,0) = mMatrix(2,0);
			M22Minor(0,1) = mMatrix(0,2);
			M22Minor(1,1) = mMatrix(2,2);
			
			
			
			if (twoByTwoDet(M33Minor)>0)
				if (twoByTwoDet(N33Minor)>0)
					returner = 5;
				else
					returner = 7;
			else if ((twoByTwoDet(M11Minor)+twoByTwoDet(M22Minor))>0)
				returner = 4;
			else if ((twoByTwoDet(M11Minor)+twoByTwoDet(M22Minor))<0)
				returner = 8;
			
			else if (twoByTwoDet(M33Minor)<0)
				double newAlpha = -a/3;
				
			
				
			  
	
	
	
	/*Most of the computation is done now, now we v */
}

	

int main () {
	
	
	
	int RR=3;
	
	int CC=3;
	std::vector< std::vector< double > > matrix(3);
	for ( int i = 0 ; i < RR ; i++ )
		matrix[i].resize(CC);

	double b = threeByThreeDet(matrix);
	std::cout << b << std::endl ;
	return 0;
}
