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
	double unitWidth;
	int unitNum;

	vector<double> gratingGenerator(double pitchLength, double samplingRate, double unitLength, int uN);
	void createGrating();

public:
	Grating();
	Grating(double gratingPitch, double samplingRate, double unitLength, int uN);

	vector<double> gGrid;
	vector<double> gAmp;
};

Grating::Grating() {
	pitchWidth = 1;
	gSR = 0.001;
	unitWidth = 60;
	unitNum = 3;

	createGrating();
}

Grating::Grating(double gratingPitch, double samplingRate, double unitLength, int uN) {
	pitchWidth = gratingPitch;
	gSR = samplingRate;
	unitWidth = unitLength;
	unitNum = uN;

	createGrating();
}

vector<double> Grating::gratingGenerator(double pitchLength, double SR, double unitLength, int uN) {
  int pitchPtNum = ceil(pitchLength / SR);
  int pitchNum = ceil(unitLength / pitchLength);
	int gPtNum = pitchPtNum * pitchNum * (uN * 2 + 1) + floor(pitchPtNum / 2);

	vector<double> grating(gPtNum, 0);
 
  for (int i = 0; i < uN * 2 + 1; i++) {
    for (int j = 0; j < pitchNum; j++) {
      for(int k = 0; k < pitchPtNum; k++) {
        
        if (k >= floor(pitchPtNum / 2)) {
          if (i % 2 == 0) {
            grating[i * pitchNum * pitchPtNum + j * pitchPtNum + k] = 1;
          }
          else {
            grating[i * pitchNum * pitchPtNum + j * pitchPtNum + k] = 4;
          }      
        } 
        
      }        
    } 
  } 
  
	return grating;
}

void Grating::createGrating() {
	gAmp = gratingGenerator(pitchWidth, gSR, unitWidth, unitNum);
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

	double L2I = atof(argv[1]); // Grating-to-Grating distance in mm
	double L2 = L2I * 10000; // Grating-to-Grating Distance in micron 

	double lambda1  = atof(argv[2]);
	double lambda = lambda1/1000; // in microns

	double gPitch = 2.1;
	double gSampleRate = 0.001;
	double uWidth = 35;
	double uNum = 65;

	Grating phasorG(gPitch, gSampleRate, uWidth, uNum);
	int N = phasorG.gGrid.size();



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
    objectAll[om] = 0.25 * sampleTest1.objectSample[om] + 0.50 * sampleTest2.objectSample[om]; 
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
	double pShift = -1 * pShiftConstant * lambda / 8 ; 

	double L1 = 100e+04; // Sorce-to-Grating Distance
	//	double D = 120e+02; // Grating-to-Grating Distance // D is being read as an input variable 
  //	double L2 = 350e+04; // Grating-to-Detector Distance

	double ZOG = 450e+03; // Object-to-Detector Distance ; 500 is 50 cm
  double ZSO = L1 - ZOG ; // Grating-to-Object Distance 


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
	vector<complex<double> > sAblank(N, complex<double>(0, 0));
 	vector<complex<double> > sAobj(N, complex<double>(0, 0));
	vector<complex<double> > go(N3, complex<double>(0, 0));
  vector<complex<double> > phaseShiftObj(objPtNum, complex<double>(0, 0));
 	vector<complex<double> > dAblank(M, complex<double>(0, 0));
	vector<complex<double> > dAobj(M, complex<double>(0, 0));
	vector<double> dIblank(M, 0);
	vector<double> dIobj(M, 0);

 

#pragma omp parallel shared(N,objPtNum,M,lambda,k,pShift,oShift,L1,L2,ZOG,ZSO,unitImag,sAblank,sAobj,phaseShiftObj,dAblank,dAobj,dIblank,dIobj) 
  {
  
	double rSG = 0;
	double pOD = 0;
	int ODindex = 0;
 
  #pragma omp for 
 	for (int om = 0; om < objPtNum; om++) {
		phaseShiftObj[om] = exp(unitImag * oShift * objectAll[om]);
	}
  #pragma omp for
  for (int n = 0; n < N; n++) {
		rSG = sqrt(pow(L1, 2) + pow(phasorG.gGrid[n], 2));
		pOD = (L1 - ZOG) * (phasorG.gGrid[n]) / L1;
		ODindex = floor(pOD / dSR) + (floor(objPtNum / 2) + 1);
		sAblank[n] = exp(unitImag * k * rSG) / rSG;
    sAobj[n] = phaseShiftObj[ODindex] * exp(unitImag * k * rSG) / rSG;
	}

  #pragma omp for
	for (int n = 0; n < N; n++) {
		sAblank[n] = sAblank[n] * exp(unitImag * k * pShift * phasorG.gAmp[n]);
    sAobj[n] = sAobj[n] * exp(unitImag * k * pShift * phasorG.gAmp[n]);
	}

	double rG2D = 0;
 
  #pragma omp for 
 	for (int m = 0; m < M; m++) {
		for (int n = 0; n < N; n++) {
			rG2D = sqrt(pow(L2, 2) + pow(dGrid[m] - phasorG.gGrid[n], 2));
			dAblank[m] = dAblank[m] + sAblank[n] * exp(unitImag * k * rG2D) / pow(rG2D, 2);
			dAobj[m] = dAobj[m] +  sAobj[n] * exp(unitImag * k * rG2D) / pow(rG2D, 2);
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
	printVector(dIblank,"dIblank_ObjMPG_S1");
	printVector(dIobj,"dIobj_ObjMPG_S1");
	printVector(dGrid,"dGrid_ObjMPG_S1");
	printVector(objectAll,"objectAll_ObjMPG_S1");
	printVector(phasorG.gAmp,"gAmp_ObjMPG_S1");
	printVector(phasorG.gGrid,"gGrid_ObjMPG_S1");



 int lambdaint= static_cast<int>(round(lambda1*1000)); 
 int distanceint= static_cast<int>(round(L2I));
 

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
 
