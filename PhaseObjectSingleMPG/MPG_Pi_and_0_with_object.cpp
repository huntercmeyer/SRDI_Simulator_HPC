#define _USE_MATH_DEFINES
#define C 3.0e+08
#define h 6.62607004e-34

#include <cmath>
#include <vector>
#include <complex>
#include <algorithm>  

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <omp.h>

using namespace std;

vector<double> avgFilter(vector<double>orignalData, vector<double>filterSet);
void printVector(vector<double>toPrint, string vectorName, string parameters);

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
            grating[i * pitchNum * pitchPtNum + j * pitchPtNum + k] = 0.5;
          }
          else {
            grating[i * pitchNum * pitchPtNum + j * pitchPtNum + k] = 2;
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
  int startIndexPyramid = (int) centerOptIndex - objPtNum  + 1;
  int startIndexTrapezoid = (int) centerOptIndex - 3 * objPtNum / 2  + 1; 
  int startIndexParabolaRamp = (int) centerOptIndex - 3 * objPtNum / 2  + 1; 

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
 
  else if (objectOption == 3) {  // Thin Slab
    for (int i = startIndex; i < startIndex + objPtNum; i++) {
      objectSample[i] = 1;
    }
  }
  else if (objectOption == 4) {  // Pyramid
    double adder = tan(angle*M_PI/180);
    for (int i = startIndexPyramid; i < centerOptIndex + objPtNum ; i++) {
      if (i <= centerOptIndex ) {
	objectSample[i] = (double)(i - startIndexPyramid + adder) / objPtNum;
      } else {
	objectSample[i] = (double)(objPtNum - i + centerOptIndex ) / objPtNum;
      }
    }
  }
 
  else if (objectOption == 5) {  // Parabola (negative a) vertex form
    for (int i = startIndex; i < 2 * centerOptIndex - startIndex + 1; i++) {
      objectSample[i] =  -1 /pow(centerOptIndex - startIndex, 2) * pow(i - centerOptIndex, 2) + 1; 
    }
  }
 
  else if (objectOption == 6) {  // Parabola (positive a) vertex form
    for (int i = startIndex; i < 2 * centerOptIndex - startIndex + 1; i++) {
      objectSample[i] =  1 /pow(centerOptIndex - startIndex, 2) * pow(i - centerOptIndex, 2) + 1; 
    }
  }

  else if (objectOption == 7) {  // Parabola (positive a) vertex form with ramps on sides 
    double adder = tan(angle*M_PI/180);
    for (int i = startIndexParabolaRamp; i < centerOptIndex + 3 * objPtNum /2; i++) {
      if (i <= centerOptIndex - objPtNum / 2) {
	objectSample[i] = (double)(i - startIndexParabolaRamp + adder) / objPtNum;
      } 
      else if (i <= centerOptIndex + objPtNum / 2) {
	objectSample[i] =  (1 /pow(centerOptIndex - startIndex, 2) * pow(i - centerOptIndex, 2) + 1.0) / 2; 
      }
      else if (i > centerOptIndex + objPtNum / 2){
	objectSample[i] = (double)(objPtNum - i + centerOptIndex + objPtNum / 2 )  / objPtNum;
      }   
    }
  }
 
 
  else {  // Trapezoid 
    double adder = tan(angle*M_PI/180);
    for (int i = startIndexTrapezoid; i < centerOptIndex + 3 * objPtNum /2; i++) {
      if (i <= centerOptIndex - objPtNum / 2) {
	objectSample[i] = (double)(i - startIndexTrapezoid + adder) / objPtNum;
      } 
      else if (i <= centerOptIndex + objPtNum / 2) {
	objectSample[i] = (double)(objPtNum + adder) / objPtNum;
      }
      else if (i > centerOptIndex + objPtNum / 2){
	objectSample[i] = (double)(objPtNum - i + centerOptIndex + objPtNum / 2 ) / objPtNum;
      }
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

  double L2I = atof(argv[1]); // Grating-to-Grating distance in mm (Hunter: I believe this is actually cm)
  double L2 = L2I * 10000; // Grating-to-Grating Distance in micron 

  double lambda1  = atof(argv[2]); // Wavelength in nm
  double lambda = lambda1/1000; // in microns

  double gPitch = atof(argv[3]); // pitch in microns
  double gSampleRate = 0.001;
  double W = atof(argv[4]); // grating period in microns
  double uWidth = W/2; // W/2
  double uNum = 50;

  Grating phasorG(gPitch, gSampleRate, uWidth, uNum);
  int N = phasorG.gGrid.size();

  int pitchPtNum = ceil(gPitch / gSampleRate);
  int pitchNum = ceil(uWidth / gPitch);
  int gPtNum = pitchPtNum * pitchNum * (uNum * 2 + 1) + floor(pitchPtNum / 2);

  // Phase Object Setup
  /*
    double scanItem1 = 3;
    double objectSize1 = 1400; 
    double objectWeight1 = 1;
    double objectPl1 =  0;
    double oShift = M_PI;
    double angle = 60;
    double scanItem2 = 3;
    double objectSize2 = 700;
    double objectWeight2 = 1;
    double objectPl2 = 0;
  */



  double scanItem1 = atof(argv[5]); // Which object to use

  double objectSize1 = 1600.;
  // If the scan item is a pyramid (triangle = 4), then the object size should be halved.
  // TODO: Add this functionality for the trapezoid.  The trapezoid should be 1/3 of what you want.
  if(scanItem1 == 4) { // Triangle
    objectSize1 /= 2.;
  } else if(scanItem1 == 8) { // Trapezoid
    objectSize1 /= 3.;
  }

  // double objectSize1 = 800;  
  double objectWeight1 = 1;
  double objectPl1 =  0;
  double oShift = M_PI;
  double angle = 60;



  // Detector Setup
  vector<double> dGrid;
  double dLength = 13000;
  double dSR = 0.1;
  int dPtNum = floor(dLength / dSR);
  if (dPtNum % 2 == 0) {
    dPtNum += 1;
  }

  dGrid.resize(dPtNum, 0);
  for (int i = 0; i < dPtNum; i++) {
    dGrid[i] = (-floor(dPtNum / 2) + i) * dSR;
  }
  int M = dGrid.size();  

  double ObjSR = 0.1;
  
  // Object for scan
  int objPtNum = floor(dLength / ObjSR);
  TestObject sampleTest1(objPtNum, scanItem1, objectSize1, objectPl1, ObjSR, angle);
  //	TestObject sampleTest2(objPtNum, scanItem2, objectSize2, objectPl2, dSR, angle); 
 
  vector<double> objectAll(objPtNum, 0);
  //vector<double> object1(N, 0);
  //vector<double> side(N, 0);
  
  for (int om = 0; om < objPtNum; om++) {

    objectAll[om] = 16.0 * sampleTest1.objectSample[om] ; 
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
  // double pShiftConstant= lambda / 0.22e-3;
  //double lambda = 0.80e-3;  // lambda is an input modifed dparameter
  double k = 2 * M_PI / lambda; 
  double pShift = -1 * pShiftConstant * lambda / 8 ; 

  double L1 = atof(argv[6])*1e+04; // Source to grating distance in microns.  argv[6] must be in cm.
  // double L1 = 100e+04; // Sorce-to-Grating Distance ; Hunter: 100e+04 is in 1e6 microns = 1 meter
  //	double D = 120e+02; // Grating-to-Grating Distance // D is being read as an input variable 
  //	double L2 = 350e+04; // Grating-to-Detector Distance

  // double ZOD = 250e+04; // Object-to-Detector Distance ; 500 is 50 cm ; Hunter: I think this is wrong.  I think 500 is 500 cm
  //  double ZOD = 350e+04; // Object-to-Detector Distance ; 500 is 50 cm
  // double ZGO = L2 - ZOD ; // Grating-to-Object Distance 

  double ZGO = 13.5e+04; // Grating-to-Object Distance in micron
  double ZOD = L2 - ZGO; // Object-to-Detector Distance in micron;

  complex<double> unitImag(0, 1);
  vector<complex<double> > sA(N, complex<double>(0, 0));
  vector<complex<double> > go(N3, complex<double>(0, 0));
  vector<complex<double> > phaseShiftObj(objPtNum, complex<double>(0, 0));
  vector<complex<double> > dAblank(M, complex<double>(0, 0));
  vector<complex<double> > dAobj(M, complex<double>(0, 0));
  vector<double> dIblank(M, 0);
  vector<double> dIobj(M, 0);

 

#pragma omp parallel shared(N,objPtNum,M,lambda,k,pShift,oShift,L1,L2,ZOD,ZGO,unitImag,sA,phaseShiftObj,dAblank,dAobj,dIblank,dIobj) 
  {
    double rSG = 0;
#pragma omp for
    for (int n = 0; n < N; n++) {
      rSG = sqrt(pow(L1, 2) + pow(phasorG.gGrid[n], 2));
      sA[n] = exp(unitImag * k * rSG) / rSG;
    }

#pragma omp for
    for (int n = 0; n < N; n++) {
      sA[n] = sA[n] * exp(unitImag * k * pShift * phasorG.gAmp[n]);
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
	rG2D = sqrt(pow(L2, 2) + pow(dGrid[m] - phasorG.gGrid[n], 2));
	pOD = phasorG.gGrid[n] + (L2 - ZOD) * (dGrid[m] - phasorG.gGrid[n]) / L2;
	//			ODindex = floor(pOD / ObjSR) + (floor(objPtNum / 2) + 1);
	ODindex = round((pOD / ObjSR) + (objPtNum / 2)) ;
	dAblank[m] = dAblank[m] + sA[n] * exp(unitImag * k * rG2D) / pow(rG2D, 2);
	if (ODindex >= 0 && ODindex < objPtNum ) {
	  dAobj[m] = dAobj[m] + phaseShiftObj[ODindex] * sA[n] * exp(unitImag * k * rG2D) / pow(rG2D, 2);
	}
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

  stringstream parameters_ss;
  parameters_ss << fixed << setprecision(0) << "W_" << W << "_p_" << gPitch << "_DSG_" << L1/10000 << "_DGO_" << ZGO/10000 << "_DSD_" << (L1 + L2)/10000 << "_SHAPE_" << scanItem1 << "_SR100_gSR1_16pi_ObjSize1600New_MPG3000_Det14000_PIand0PS_1"; // fixed prevents it from using scientific notation
  string parameters = parameters_ss.str();
  // string parameters = "W_120_p_2_DSG100_DGO30_DSD300_TriangleSR100_gSR1_16pi_ObjSize1600New_MPG3000_Det14000_PIand0PS_1";
 
  // Output file
  printVector(dIblank,"dIblank_MPG_", parameters);
  printVector(dIobj,"dIobj_MPG_", parameters);
  printVector(dGrid,"dGrid_MPG_", parameters);
  printVector(objectAll,"objectAll_MPG_", parameters);
  printVector(phasorG.gAmp,"gAmp_MPG_", parameters);
  printVector(phasorG.gGrid,"gGrid_MPG_", parameters);

  string fileName3 = "InputVars_";
  fileName3.append(parameters);
  fileName3.append(".txt"); 
  

  //Print map of vars to file
  ofstream varFile;
  varFile.open(fileName3.c_str());
  varFile << "G1 Length:\t\t" << uWidth * uNum  << '\n';

  varFile << "G1 SR:\t\t" << gSampleRate<< '\n';
	
  varFile << "G1 Pitch:\t\t" << gPitch<< '\n';

  varFile << "dLength:\t\t" << dLength << '\n';
  varFile << "dSR:\t\t" << dSR << '\n';

  varFile << "ScanItem1:\t\t" << scanItem1<< '\n';
  varFile << "objectSize1:\t\t" << objectSize1<< '\n';
  varFile << " ObjSR:\t\t" << ObjSR<< '\n';


  varFile << "L1:\t\t" << L1 << '\n';
  varFile << "L2:\t\t" << L2 << '\n';
  varFile << "zOD:\t\t" << ZOD << '\n';
  varFile.close();


  int lambdaint= static_cast<int>(round(lambda1*1000)); 
  int distanceint= static_cast<int>(round(L2I));
 


  return 0;
} // main
 
 

void printVector(vector<double> toPrint, string vectorName, string parameters)
{
  vectorName.append("parameters_");
  vectorName.append(parameters);
  vectorName.append(".bin");
  // vectorName.append("parameters_W_120_DSG1_DSD4_TriangleSR100_gSR1_16pi_ObjSize1600New_MPG4000_Det13000_PIand0PS_1.bin");
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
 
