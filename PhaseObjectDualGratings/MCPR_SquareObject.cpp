#define _USE_MATH_DEFINES
#define C 3.0e+08
#define h 6.62607004e-34

#include <cmath>
#include <vector>
#include <complex>
#include <algorithm>  

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <omp.h>

using namespace std;

vector<double> avgFilter(vector<double>orignalData, vector<double>filterSet);
void printVector(vector<double>toPrint, string vectorName);

class Grating {
	double pitchWidth;
	double gSR;
	double gLength;
	double gDutyCycle;

	vector<double> gratingGenerator(double pitchLength, double samplingRate, double gratingLength, double dutyCycle);
	void createGrating();

public:
	Grating();
	Grating(double gratingPitch, double samplingRate, double gratingLength, double gratingDutyCycle);

	vector<double> gGrid;
	vector<double> gAmp;
};

Grating::Grating() {
	pitchWidth = 1;
	gSR = 0.001;
	gLength = 100;
	gDutyCycle = 0.5;

	createGrating();
}

Grating::Grating(double gratingPitch, double samplingRate, double gratingLength, double gratingDutyCycle) {
	pitchWidth = gratingPitch;
	gSR = samplingRate;
	gLength = gratingLength;
	gDutyCycle = gratingDutyCycle;

	createGrating();
}

vector<double> Grating::gratingGenerator(double pitchLength, double SR, double gratingLength, double dutyCycle) {
	int onepitchPointNum = (int)floor(pitchLength / SR);
	if (onepitchPointNum % 2 != 0) {
		onepitchPointNum += 1;
	}

	int pitchNum = (int)floor(gratingLength / pitchLength);
	if (pitchNum % 2 == 0) {
		pitchNum += 1;
	}

	int gPointNum = pitchNum * onepitchPointNum + floor(onepitchPointNum * (1 - dutyCycle));
	vector<double> grating(gPointNum, 0);

	for (int j = 0; j < pitchNum; j++) {
		for (int k = 0; k < onepitchPointNum; k++) {
			if (k >= floor(onepitchPointNum * (1 - dutyCycle))) {
				grating[k + j * onepitchPointNum] = 1;
			}
		}
	}

	return grating;
}

void Grating::createGrating() {
	gAmp = gratingGenerator(pitchWidth, gSR, gLength, gDutyCycle);
	gGrid.resize(gAmp.size(), 0);

	for (int i = 0; i < gAmp.size(); i++) {
		gGrid[i] = (-floor(gAmp.size() / 2) + i) * gSR;
	}
}

class TestObject {
	int objectOption;
	double objectDim;
	double objectPos;
	int ptN;
	double oSR;

	void objectGenerator(double angle);

public:
	TestObject(int gratingN, double angle);
	TestObject(int gratingN, int objectIndex, double objectLength, double objectCenter, double samplingRate, double angle);

	vector<double> objectSample;
};

TestObject::TestObject(int gratingN, double angle) {
	objectOption = 1;
	objectDim = 100;
	objectPos = 0;
	ptN = gratingN;
	oSR = 0.001;

	objectGenerator(angle);
}

TestObject::TestObject(int gratingN, int objectIndex, double objectLength, double objectCenter, double samplingRate, double angle) {
	objectOption = objectIndex;
	objectDim = objectLength;
	objectPos = objectCenter;
	ptN = gratingN;
	oSR = samplingRate;

	objectGenerator(angle);
}

void TestObject::objectGenerator(double angle) {
	objectSample.resize(ptN, 0);

	int centerOptIndex = (int) (floor(ptN / 2) + 1) + ceil(objectPos / oSR);
	int objPtNum = ceil(objectDim / oSR);
	int startIndex = (int) centerOptIndex - floor(objPtNum / 2) + 1;

	if (objectOption == 1) {  // Ramp
		double adder = tan(angle*M_PI/180);
		for (int i = startIndex; i < startIndex + objPtNum; i++) {
			objectSample[i] = (double)(i - startIndex + adder) / objPtNum;
		}
	}
	else if (objectOption == 2) {  // Sphere
		for (int i = startIndex; i < 2 * centerOptIndex - startIndex + 1; i++) {
			objectSample[i] = 2 * sqrt(pow(centerOptIndex - startIndex, 2) - pow(i - centerOptIndex, 2)) / (2 * (centerOptIndex - startIndex));
		}
	}
	else {  // Thin Slab
		for (int i = startIndex; i < startIndex + objPtNum; i++) {
			objectSample[i] = 1;
		}
	}
}

int main(int argc, char** argv)
{

/*	
        if (argc < 4) {
	  cout << argv[0] << " D  lambda dir" << "\n";
	  cout << " Notes: D = distance between two gratings mm "  << "\n";
	  cout << "        lambda =  wavelength in nm"  << "\n";
	  cout << "        dir =  directory to save files"  << "\n";
	  exit(0);
	}
*/

	double D1 = atof(argv[1]); // Grating-to-Grating distance in mm
	double D = D1*1000; // Grating-to-Grating Distance in micron 

	double lambda1  = atof(argv[2]);
	double lambda = lambda1/1000; // in microns

	//char dir[1000];
	//dir = strcpy(argv[3]);   
  //  string dir = argv[3];
  // Phase Grating G1 Setup
	double G1Pitch = 2.4; 
	double G1SR = 0.002; 
	//double G1SR = 0.2; 
	double G1Length = 2460;   
	double G1DC = 0.5;

	Grating G1(G1Pitch, G1SR, G1Length, G1DC);
	int N1 = G1.gGrid.size();
  int N = N1;

	double G2Pitch = 2.4;
	double G2SR = 0.002;
	//double G2SR = 0.2;
	double G2Length = 2460; 
	double G2DC = 0.5;

	Grating G2(G2Pitch, G2SR, G2Length, G2DC);
	int N2 = G2.gGrid.size(); 


  // Phase Object Setup
  double scanItem1 = 3;
 	double objectSize1 = 1600; 
	double objectWeight1 = 1;
	double objectPl1 =  0;
	double oShift = M_PI;
	double angle = 60;
	double scanItem2 = 3;
	double objectSize2 = 500;
	double objectWeight2 = 1;
	double objectPl2 = 0;




	// Detector Setup
	vector<double> dGrid;
	double dLength = 4000;
	double dSR = 0.01;
	int dPtNum = floor(dLength / dSR);
	if (dPtNum % 2 == 0) {
		dPtNum += 1;
	}

	dGrid.resize(dPtNum, 0);
	for (int i = 0; i < dPtNum; i++) {
		dGrid[i] = (-floor(dPtNum / 2) + i) * dSR;
	}
	int M = dGrid.size();  


	// Object for scan
	int objPtNum = floor(4100 / dSR);
	TestObject sampleTest1(objPtNum, scanItem1, objectSize1, objectPl1, dSR, angle);
	TestObject sampleTest2(objPtNum, scanItem2, objectSize2, objectPl2, dSR, angle);
 
  vector<double> objectAll(objPtNum, 0);
  //vector<double> object1(N, 0);
  //vector<double> side(N, 0);
  
  for (int om = 0; om < objPtNum; om++) {
   // objectAll[om] = 0.78 * sampleTest.objectSample[om]; 
    //objectAll[n] = 0.587 * sampleTest.objectSample[n]; //+  0.159 * sampleTest2.objectSample[n] -0.4138 *  sampleTest3.objectSample[n];
    objectAll[om] = sampleTest1.objectSample[om];  
  }
  
  int N3 = objectAll.size();

/*
	// Pixel Filter
	vector<double> pxlFilter(M, 0);
	double pixelSize = 100;
	int cX = floor(M / 2);
	int halfPixelPointNum = ceil(pixelSize / 2 / dSR);
	for (int fi = cX - halfPixelPointNum; fi <= cX + halfPixelPointNum; fi++) {
		pxlFilter[fi] = (double)1 / (halfPixelPointNum * 2);
	}
*/


//Monochromatic set-up

  double pShiftConstant= lambda/0.44e-3;
	//double lambda = 0.80e-3;  // lambda is an input modifed dparameter
	double k = 2 * M_PI / lambda; 
	double pShift = -1 * 0.135 * pShiftConstant* lambda ; 

	double L1 = 120e+04;; // Sorce-to-Grating Distance
	//	double D = 120e+02; // Grating-to-Grating Distance // D is being read as an input variable 
	double L2 = 299e+04 - L1 - D; // Grating-to-Detector Distance

	double ZOD = 450e+03; // Object-to-Detector Distance ; 500 is 50 cm
  double ZGO = L2 - ZOD ; // Grating-to-Object Distance 


//Polychromatic set-up
/*
  double pShiftConstant= lambda/0.50e-3;
	//double lambda = 0.80e-3;  // lambda is an input modifed dparameter
	double k = 2 * M_PI / lambda; 
	double pShift = -1 * 0.10 * pShiftConstant* lambda ; 
 
	double L1 = 465e+04;; // Sorce-to-Grating Distance
	//	double D = 120e+02; // Grating-to-Grating Distance // D is being read as an input variable 
	double L2 = 880e+04 - L1 - D; // Grating-to-Detector Distance  
*/

	complex<double> unitImag(0, 1);
	vector<complex<double> > sg1(N1, complex<double>(0, 0));
	vector<complex<double> > g1g2(N2, complex<double>(0, 0));
//	vector<complex<double> > dA(M, complex<double>(0, 0)); // Amplitude on detector
	vector<complex<double> > go(N3, complex<double>(0, 0));
  vector<complex<double> > phaseShiftObj(objPtNum, complex<double>(0, 0));
 	vector<complex<double> > dAblank(M, complex<double>(0, 0));
	vector<complex<double> > dAobj(M, complex<double>(0, 0));
	vector<double> dIblank(M, 0);
	vector<double> dIobj(M, 0);

 

#pragma omp parallel shared(N1,N2,N3,N,objPtNum,M,lambda,k,pShift,oShift,L1,D,L2,ZOD,ZGO,unitImag,sg1,g1g2,phaseShiftObj,dAblank,dAobj,dIblank,dIobj) 
  {
	double rSG1 = 0;
  #pragma omp for
	for (int n1 = 0; n1 < N1; n1++) {
		rSG1 = sqrt(pow(L1, 2) + pow(G1.gGrid[n1], 2));
		sg1[n1] = exp(unitImag * k * rSG1) / rSG1;
	}

  #pragma omp for 
	for (int n1 = 0; n1 < N1; n1++) {
		sg1[n1] = sg1[n1] * exp(unitImag * k * pShift * G1.gAmp[n1]);
	}

	double rG1G2 = 0;
  #pragma omp for 
	for (int n2 = 0; n2 < N2; n2++) {
		for (int n1 = 0; n1 < N1; n1++) {
			rG1G2 = sqrt(pow(D, 2) + pow(G2.gGrid[n2] - G1.gGrid[n1], 2));
      if ((D / rG1G2) < 0.9998  ) {                 // theta=1 degrees (1 degrees acceptance) 
            continue;
        }
			g1g2[n2] = g1g2[n2] + sg1[n1] * exp(unitImag * k * rG1G2) / pow(rG1G2, 2);  
		}
	}

	#pragma omp for 
  for (int n2 = 0; n2 < N2; n2++) {
		g1g2[n2] = g1g2[n2] * exp(unitImag * k * pShift * G2.gAmp[n2]) * D / unitImag / lambda;
	}

	double rG2D = 0;
 	double pOD = 0;
	int ODindex = 0;
 
  #pragma omp for 
 	for (int om = 0; om < objPtNum; om++) {
		phaseShiftObj[om] = exp(unitImag * oShift * objectAll[om]);
	}
 
  #pragma omp for 
 	for (int m = 0; m < M; m++) {
		for (int n = 0; n < N; n++) {
			rG2D = sqrt(pow(L2, 2) + pow(dGrid[m] - G2.gGrid[n], 2));
			pOD = G2.gGrid[n] + (L2 - ZOD) * (dGrid[m] - G2.gGrid[n]) / L2;
			ODindex = floor(pOD / dSR) + (floor(objPtNum / 2) + 1);

			dAblank[m] = dAblank[m] + g1g2[n] * exp(unitImag * k * rG2D) / pow(rG2D, 2);
			dAobj[m] = dAobj[m] + phaseShiftObj[ODindex] * g1g2[n] * exp(unitImag * k * rG2D) / pow(rG2D, 2);
		}
	}
 
 #pragma omp for
	for (int m = 0; m < M; m++) {
		dAblank[m] = dAblank[m] * L2 / lambda / unitImag;
		dAobj[m] = dAobj[m] * L2 / lambda / unitImag;
	}


 #pragma omp for
	for (int m = 0; m < M; m++) {
		dIblank[m] = pow(abs(dAblank[m]), 2);
		dIobj[m] = pow(abs(dAobj[m]), 2);
	}
 
}
//  vector<double> dIfiltered = avgFilter(dI, pxlFilter); // Intensity filtered by pixel


	// Output file
	printVector(dIblank,"dIblank_D12mm_S4");
	printVector(dIobj,"dIobj_D12mm_S4");
	printVector(dGrid,"dGrid_D12mm_S4");
	printVector(objectAll,"objectAll_D12mm_S4");
	printVector(G1.gAmp,"G1_gAmp_D12mm_S4");
	printVector(G2.gAmp,"G2_gAmp_D12mm_S4");
	printVector(G1.gGrid,"G1_gGrid_D12mm_S4");
	printVector(G2.gGrid,"G2_gGrid_D12mm_S4");



 int lambdaint= static_cast<int>(round(lambda1*1000)); 
 int distanceint= static_cast<int>(round(D1));
 

	return 0;
}
 


void printVector(vector<double> toPrint, string vectorName)
{
	vectorName.append(".bin");
	ofstream testFileOutput2(vectorName.c_str(), iostream::out | ofstream::binary);
	testFileOutput2.write(reinterpret_cast<char *>(&toPrint[0]), toPrint.size() * sizeof(double));
	testFileOutput2.close();

}

/*
vector<double> avgFilter(vector<double>orignalData, vector<double>filterSet) {
	int m = orignalData.size();
	int n = filterSet.size();

	int filteredlength = m + n - 1;
	vector<double> filteredData(filteredlength, 0);

	int j1 = 0;
	int j2 = 0;

	for (int k = 1; k <= filteredlength; k++) {
		j1 = max(1, k + 1 - m);
		j2 = min(k, n);

		for (int j = j1; j <= j2; j++) {
			filteredData[k - 1] = filteredData[k - 1] + filterSet[j - 1] * orignalData[k - j];
		}
	}

	vector<double> filteredDataCut(m, 0);
	int cutOffindex = (n - 1) / 2;
	for (int i = 0; i < m; i++) {
		filteredDataCut[i] = filteredData[i + cutOffindex];
	}

	return filteredDataCut;

}
*/
 
