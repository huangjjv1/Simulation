// Translate this file with
//
// g++ -O3 --std=c++11 assignment-2019.c -o assignment
//
// Run it with
//
// ./assignment
//
// There should be a result.pvd file that you can open with Paraview.
// Sometimes, Paraview requires to select the representation "Point Gaussian"
// to see something meaningful.
//
// (C) 2019 Tobias Weinzierl
#include <omp.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <cmath>
#include <limits>
#include <iomanip>
#include <algorithm>

double t          = 0;
double tFinal     = 0;
double tPlot      = 0;
double tPlotDelta = 0;

int NumberOfBodies = 0;

/**
 * Pointer to pointers. Each pointer in turn points to three coordinates, i.e.
 * each pointer represents one molecule/particle/body.
 */
double** x;

/**
 * Equivalent to x storing the velocities.
 */
double** v;

/**
 * Global time step size used.
 */
double   timeStepSize = 0.0001;

/**
 * Maximum velocity of all particles.
 */
double   maxV;

/**
 * Minimum distance between two elements.
 */
double   minDx;


/**
 * Set up scenario from the command line.
 *
 * This operation is not to be changed in the assignment.
 */
void setUp(int argc, char** argv) {
  NumberOfBodies = (argc-2) / 6;

  x    = new double*[NumberOfBodies];
  v    = new double*[NumberOfBodies];

  int readArgument = 1;

  tPlotDelta  = std::stof(argv[readArgument]); readArgument++;//1
  tFinal      = std::stof(argv[readArgument]); readArgument++;//2

  for (int i=0; i<NumberOfBodies; i++) {
    x[i] = new double[3];
    v[i] = new double[3];

    x[i][0] = std::stof(argv[readArgument]); readArgument++;//3
    x[i][1] = std::stof(argv[readArgument]); readArgument++;//4
    x[i][2] = std::stof(argv[readArgument]); readArgument++;//5

    v[i][0] = std::stof(argv[readArgument]); readArgument++;//6
    v[i][1] = std::stof(argv[readArgument]); readArgument++;//7
    v[i][2] = std::stof(argv[readArgument]); readArgument++;//8
  }

  std::cout << "created setup with " << NumberOfBodies << " bodies" << std::endl;
  
  if (tPlotDelta<=0.0) {
    std::cout << "plotting switched off" << std::endl;
    tPlot = tFinal + 1.0;
  }
  else {
    std::cout << "plot initial setup plus every " << tPlotDelta << " time units" << std::endl;
    tPlot = 0.0;
  }
}


std::ofstream videoFile;


/**
 * This operation is not to be changed in the assignment.
 */
void openParaviewVideoFile() {
  videoFile.open( "result.pvd" );
  videoFile << "<?xml version=\"1.0\"?>" << std::endl
            << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl
            << "<Collection>";
}





/**
 * This operation is not to be changed in the assignment.
 */
void closeParaviewVideoFile() {
  videoFile << "</Collection>"
            << "</VTKFile>" << std::endl;
}


/**
 * The file format is documented at http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
 *
 * This operation is not to be changed in the assignment.
 */
void printParaviewSnapshot() {
  static int counter = -1;
  counter++;
  std::stringstream filename;
  filename << "result-" << counter <<  ".vtp";
  std::ofstream out( filename.str().c_str() );
  out << "<VTKFile type=\"PolyData\" >" << std::endl
      << "<PolyData>" << std::endl
      << " <Piece NumberOfPoints=\"" << NumberOfBodies << "\">" << std::endl
      << "  <Points>" << std::endl
      << "   <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">";

  for (int i=0; i<NumberOfBodies; i++) {
    out << x[i][0]
        << " "
        << x[i][1]
        << " "
        << x[i][2]
        << " ";
  }

  out << "   </DataArray>" << std::endl
      << "  </Points>" << std::endl
      << " </Piece>" << std::endl
      << "</PolyData>" << std::endl
      << "</VTKFile>"  << std::endl;

  videoFile << "<DataSet timestep=\"" << counter << "\" group=\"\" part=\"0\" file=\"" << filename.str() << "\"/>" << std::endl;
}



/**
 * This is the only operation you are allowed to change in the assignment.
 */
void updateBody() {
  maxV   = 0.0;
  minDx  = std::numeric_limits<double>::max();

  const double sigma   = 3.4e-10;
  const double epsilon = 1.64e-21;
  const double mass    = 39.948;

  double* force0 = new double[NumberOfBodies];
  double* force1 = new double[NumberOfBodies];
  double* force2 = new double[NumberOfBodies];

#pragma omp parallel
    {
  for (int i=0; i<NumberOfBodies; i++) {//每次计算两个粒子之间的距离，和作用力。 投影到三个方向上，force 0，1，2
    force0[i] = 0.0; //force x
    force1[i] = 0.0; //force y
    force2[i] = 0.0; //force z

#pragma omp for
    for (int j=0; j<NumberOfBodies; j++) {
      if (i!=j) {
        const double distance = sqrt(
          (x[j][0]-x[i][0]) * (x[j][0]-x[i][0]) +
          (x[j][1]-x[i][1]) * (x[j][1]-x[i][1]) +
          (x[j][2]-x[i][2]) * (x[j][2]-x[i][2])
        );

        minDx = std::min( minDx,distance );
        double quantity = 4.0  * epsilon * ( -12.0 * std::pow(sigma,12.0) / std::pow(distance,12.0) + 6.0 * std::pow(sigma,6.0) / std::pow(distance,6.0) ) / distance;

          
        force0[i] += (x[j][0]-x[i][0]) * quantity / distance ;
        force1[i] += (x[j][1]-x[i][1]) * quantity / distance ;
        force2[i] += (x[j][2]-x[i][2]) * quantity / distance ;
          
      }
    }
  }//计算结束
    }
    
    
#pragma omp parallel
    {
#pragma omp for
  for (int i=0; i<NumberOfBodies; i++) { //更新位置
    x[i][0] = x[i][0] + timeStepSize * v[i][0];
    x[i][1] = x[i][1] + timeStepSize * v[i][1];
    x[i][2] = x[i][2] + timeStepSize * v[i][2];
  }
    }
#pragma omp parallel
    {
#pragma omp for
  for (int i=0; i<NumberOfBodies; i++) { //更新速度
    v[i][0] = v[i][0] + timeStepSize * force0[i] / mass;
    v[i][1] = v[i][1] + timeStepSize * force1[i] / mass;
    v[i][2] = v[i][2] + timeStepSize * force2[i] / mass;

    double thisV = std::sqrt( v[i][0]*v[0][0] + v[i][1]*v[0][1] + v[i][2]*v[0][2] );

    maxV = std::max(maxV,thisV);

  }
    }
    
  delete[] force0;
  delete[] force1;
  delete[] force2;

  t += timeStepSize;
    
}


/**
 * Main routine.
 *
 * Not to be changed in assignment.
 */
int main(int argc, char** argv) {
  if (argc==1) {
    std::cerr << "usage: " + std::string(argv[0]) + " snapshot final-time objects" << std::endl
              << "  snapshot        interval after how many time units to plot. Use 0 to switch off plotting" << std::endl
              << "  final-time      simulated time (greater 0)" << std::endl
              << std::endl
              << "Examples:" << std::endl
              << "10.0  10000.0  0 0 0 0 0 0  2e-9 0 0 0 0 0  0.9e-9 1e-9 0 0 0 0  \t Three body setup" << std::endl
              << std::endl;

    return -1;
  }
  else if ( (argc-3)%6!=0 ) {
    std::cerr << "error in arguments: each planet is given by six entries (position, velocity)" << std::endl;
    return -2;
  }

  setUp(argc,argv);

  openParaviewVideoFile();

  int snapshotCounter = 0;
  if (t > tPlot) {
    printParaviewSnapshot();
    std::cout << "plotted initial setup" << std::endl;
    tPlot = tPlotDelta;
  }

  int timeStepCounter = 0;
  while (t<=tFinal) {
    updateBody();
    timeStepCounter++;
    if (t >= tPlot) {
      printParaviewSnapshot();
      std::cout << "plot next snapshot"
                << ",\t time step=" << timeStepCounter
                << ",\t t="         << t
                << ",\t dt="        << timeStepSize
                << ",\t v_max="     << maxV
                << ",\t dx_min="    << minDx
                << std::endl;

      tPlot += tPlotDelta;
    }
  }

  closeParaviewVideoFile();

  return 0;
}

