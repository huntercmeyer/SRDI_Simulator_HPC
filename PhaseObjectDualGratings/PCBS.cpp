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

  double lambda1  = atof(argv[2]); // Wavelength in nm
  double lambda = lambda1/1000; // in microns

  //char dir[1000];
  //dir = strcpy(argv[3]);   
  //  string dir = argv[3];
  // Phase Grating G1 Setup
  double G1Pitch = 1.0; 
  double G1SR = 0.002; 
  //double G1SR = 0.2; 
  double G1Length = 3500;
  double G1DC = 0.5;
  // TODO: consider decreasing grating size, change both d length and grating length to 1000 to test
  Grating G1(G1Pitch, G1SR, G1Length, G1DC);
  int N1 = G1.gGrid.size();

  double G2Pitch = 1.0;
  double G2SR = 0.002;
  //double G2SR = 0.2;
  double G2Length = 3500; 
  double G2DC = 0.5;
  double phaseStep = 0.1; // 0.1 micron
  double phaseStepIndex = atof(argv[4]); // HM phase step converted to micron
  double G2phase = phaseStepIndex*phaseStep;

  Grating G2(G2Pitch, G2SR, G2Length, G2DC);
  int N2 = G2.gGrid.size(); 

  // HM: Acceptance Angle is used to determine how large of an angle will be accepted in the G1->G2 propagation
  double acceptanceAngle = atof(argv[5]); // degrees

  // Detector Setup
  vector<double> dGrid;
  double dLength = 2000;
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


  // Pixel Filter
  vector<double> pxlFilter(M, 0);
  double pixelSize = 100;
  int cX = floor(M / 2);
  int halfPixelPointNum = ceil(pixelSize / 2 / dSR);
  for (int fi = cX - halfPixelPointNum; fi <= cX + halfPixelPointNum; fi++) {
    pxlFilter[fi] = (double)1 / (halfPixelPointNum * 2);
  }

  //Monochromatic set-up

  //double pShiftConstant= lambda/0.44e-3; // in case we switch back to neutrons
  //double lambda = 0.80e-3;  // lambda is an input modifed dparameter
  double k = 2 * M_PI / lambda; 
  double pShift = -1 * 0.5 * lambda; // Ultimately corresponds to a pi-shift grating
  // For now, I'll hard-code the Silicon height used in Oragnista paper
  double height = 28.17; // micron
  double density = 2.32; // g / cm^3
  double muOverRho = 3.191; // cm^2 / g
  double mu = density*muOverRho; // cm^-1
  mu *= 1e-4; // micron^-1
  mu = 0; // testing for now
  // Eventually, it will something like this:
  // height = pi / (delta*k);
  // mu = k*beta;

  // double L1 = 5e+04; // Sorce-to-Grating Distance
  double L1 = atof(argv[3])*1e+04; // Source-to-Grating converted to micron (from cm)
  //	double D = 120e+02; // Grating-to-Grating Distance // D is being read as an input variable 
  // double L2 = 299e+04 - L1 - D; // Grating-to-Detector Distance  
  //double L2 = 32e+04;
  // double L2 = atof(argv[4])*1e+04; // Grating to Detector Distance converted to micron (from cm)
  double L2 = 100e+04 - L1 - D; // Force a 100 cm Source-to-Detector Distance (Organista paper had this)

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
  vector<complex<double> > sg1(N1, complex<double>(0, 0)); // Ampltiude after g1
  vector<complex<double> > g1g2(N2, complex<double>(0, 0)); // Amplitude after g2
  vector<complex<double> > dA(M, complex<double>(0, 0)); // Amplitude on detector
  vector<double> dI(M, 0); // Intensity on detector

#pragma omp parallel shared(N1,N2,M,lambda,k,pShift,L1,D,L2,unitImag,sg1,g1g2,dA,dI) 
  {
    double rSG1 = 0;
#pragma omp for
    // Propagate from source to G1
    for (int n1 = 0; n1 < N1; n1++) {
      rSG1 = sqrt(pow(L1, 2) + pow(G1.gGrid[n1], 2));
      sg1[n1] = exp(unitImag * k * rSG1) / rSG1;
    }

#pragma omp for 
    // Pass through G1 (considers both amplitude and phase)
    for (int n1 = 0; n1 < N1; n1++) {
      sg1[n1] = sg1[n1] * exp(unitImag * k * pShift * G1.gAmp[n1] - mu*height*G1.gAmp[n1]);
    }

    double rG1G2 = 0;
#pragma omp for 
    // Propagate from G1 to G2
    for (int n2 = 0; n2 < N2; n2++) {
      for (int n1 = 0; n1 < N1; n1++) {
	rG1G2 = sqrt(pow(D, 2) + pow((G2.gGrid[n2] - G2phase) - G1.gGrid[n1], 2));
	//if ((D / rG1G2) < cos(acceptanceAngle*M_PI/180.0) ) { // acceptance angle is used to avoid numerical errors when gratings are close together, we're not sure why this is such an issue, but we believe it's a numerical problem, not a physical one.
	//  continue;
	// }
	g1g2[n2] = g1g2[n2] + sg1[n1] * exp(unitImag * k * rG1G2) / pow(rG1G2, 2);
      }
    }

#pragma omp for 
    // Pass through G2 (considers both amplitude and phase) and account for scaling missed in previous propagation (Eq. 4-9 in Goodman)
    for (int n2 = 0; n2 < N2; n2++) {
      g1g2[n2] = g1g2[n2] * exp(unitImag * k * pShift * G2.gAmp[n2] - mu*height*G2.gAmp[n2]) * D / unitImag / lambda;
    }

    double rG2D = 0;
#pragma omp for 
    // Propagate from G2 to Detector
    for (int m = 0; m < M; m++) {
      for (int n2 = 0; n2 < N2; n2++) {
	rG2D = sqrt(pow(L2, 2) + pow(dGrid[m] - (G2.gGrid[n2] - G2phase), 2));
	dA[m] = dA[m] + g1g2[n2] * exp(unitImag * k * rG2D) / pow(rG2D, 2);
      }
    }

#pragma omp for 
    // Account for scaling missed in previous propagation (Eq. 4-9 in Goodman)
    for (int m = 0; m < M; m++) {
      dA[m] = dA[m] * L2 / unitImag / lambda;
    }

#pragma omp for 
    for (int m = 0; m < M; m++) {
      dI[m] = pow(abs(dA[m]), 2);
    }
  }
  //  vector<double> dIfiltered = avgFilter(dI, pxlFilter); // Intensity filtered by pixel

  stringstream sl;
  sl << M;
  string sIlength = sl.str();

  int lambdaint= static_cast<int>(round(lambda1*1000)); 
  int distanceint= static_cast<int>(round(D1));
  int L1int = static_cast<int>(round(L1));

  stringstream parameters;
  parameters << fixed << setprecision(0) << "phaseStepIndex_" << phaseStepIndex << "_L1_" << L1*1e-3 << "mm_D_" << D*1e-3 << "mm_L2_" << L2*1e-3 << "mm_lambda_" << lambda*1e6 << "pm_dLength_" << dLength << "um_gLength_" << G1Length << "um_" << acceptanceAngle << "degrees";
  string fileName1 = "Intensity_";
  fileName1.append(parameters.str());
  fileName1.append(".bin");

  string fileName2 = "dGrid_";
  fileName2.append(parameters.str());
  fileName2.append(".bin");
 
  //	string fileName1 = dir;
  //  fileName1.append("/");
  // char fn[100];
  // sprintf(fn,"MCBS_L1_$d_D_%d_Lambda_%d_", L1int, distanceint, lambdaint); // use the input params in mm and pm
  // string fileName1 = fn;
  //	fileName1.append(fn);
  // fileName1.append(sIlength);
  // fileName1.append(".bin");
  // cout << "save " << fileName1 << "\n";;
  ofstream testFileOutput1(fileName1.c_str(), iostream::out | ofstream::binary);
  testFileOutput1.write(reinterpret_cast<char *>(&dI[0]), dI.size() * sizeof(double));
  testFileOutput1.close();

  // string fileName2 = "dGrid_"; 
  // fileName2.append(fn);
  // fileName2.append(sIlength);
  // fileName2.append(".bin");
  // cout << "save " << fileName2 << "\n";;
  ofstream testFileOutput2(fileName2.c_str(), iostream::out | ofstream::binary);
  testFileOutput2.write(reinterpret_cast<char*>(&dGrid[0]), dGrid.size() * sizeof(double)); 
  testFileOutput2.close();

  return 0;
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
 
