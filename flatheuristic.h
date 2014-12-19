#ifndef FLATHEURISTIC_H
#define FLATHEURISTIC_H

#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include "gdl.h"

#define Default_Flat_Distance_Threshold 2
#define Default_Sample_Dimensions 10
#define Default_Disc_Radius 10

class FlatHeuristic
{
    public:
        FlatHeuristic(vtkSmartPointer<vtkPolyData> newMesh, int dim=Default_Sample_Dimensions, double disc=Default_Disc_Radius, double dist=Default_Flat_Distance_Threshold);
        void GenerateFlatSet();
        bool IsPointFlat(double *point, int numNeighbors = 10);
        vtkSmartPointer<vtkPoints> ArePointsFlat (vtkSmartPointer<vtkPoints>, int numNeighbors = 10);
        void SetSampleDimensions(int dim);
        void SetDiscRadius(double disc);
        void SetDistanceThreshold(double dist);
        int GetSampleDimensions();
        double GetDiscRadius();
        double GetDistanceThreshold();
        void VisualizeFlatPoints();
        void VisualizeSampledPoints(int increment = 100);
    private:
        vtkSmartPointer<vtkPolyData> originalMesh;
        vtkSmartPointer<vtkKdTreePointLocator> kdTree;
        int sampleDimensions;
        double discRadius;
        double distanceThreshold;
        std::vector<vtkIdType> flatSet;
        std::vector<vtkIdType> FlatSet();
};

FlatHeuristic::FlatHeuristic(vtkSmartPointer<vtkPolyData> newMesh, int dim, double disc, double dist)
{
    originalMesh = vtkSmartPointer<vtkPolyData>::New();
    originalMesh->DeepCopy(newMesh);
    distanceThreshold = dist;
    discRadius = disc;
    sampleDimensions = dim;
    kdTree = vtkSmartPointer<vtkKdTreePointLocator>::New();
    kdTree->SetDataSet(originalMesh);
    kdTree->BuildLocator();
    kdTree->Update();
}

void FlatHeuristic::GenerateFlatSet()
{
    flatSet = FlatSet();
}

// Render mesh and concave points
void FlatHeuristic::VisualizeFlatPoints()
{
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for (int i=0;i<flatSet.size();i++)
    {
        double *point = new double[3];
        originalMesh->GetPoint(flatSet[i],point);
        points->InsertNextPoint(point);
    }
    VisualizePointsOnMesh(originalMesh, points);
}

void FlatHeuristic::VisualizeSampledPoints(int increment)
{
    // Output
    vtkSmartPointer<vtkPoints> pointsToVisualize = vtkSmartPointer<vtkPoints>::New();

    // Get point normals for mesh
    bool hasPointNormals = GetPointNormals(originalMesh);
    // If normals don't exist - generate them
    if(!hasPointNormals)
    {
        std::cout << "No point normals were found. Computing normals..." << std::endl;

        vtkSmartPointer<vtkPolyDataNormals> normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();
    #if VTK_MAJOR_VERSION <= 5
        normalGenerator->SetInput(originalMesh);
    #else
        normalGenerator->SetInputData(originalMesh);
    #endif
        normalGenerator->ComputePointNormalsOn();
        normalGenerator->ComputeCellNormalsOff();
        normalGenerator->Update();

        originalMesh = normalGenerator->GetOutput();
    }

    // Double normals in an array
    vtkFloatArray* normalDataFloat = vtkFloatArray::SafeDownCast(originalMesh->GetPointData()->GetArray("Normals"));

    // Array creation success
    if(normalDataFloat)
    {

     // create point set of circle sample points
     double *samplePoints [sampleDimensions * sampleDimensions];
     double sampleStep = 2 * discRadius/sampleDimensions;
     for (int i=0;i<sampleDimensions;i++)
     {
         for (int j=0;j<sampleDimensions;j++)
         {
             double *temp = new double[3];
             temp[0] = i*sampleStep - discRadius;
             temp[1] = j*sampleStep - discRadius;
             temp[2] = 0;
             if (temp[0] * temp[0] + temp[1] * temp[1] <= discRadius * discRadius)
             {
                samplePoints[(i*sampleDimensions) + j] = temp;
             }
             else
             {
                 double *zero = new double[3];
                 zero[0] = 0;
                 zero[1] = 0;
                 zero[2] = 0;
                 samplePoints[(i*sampleDimensions) + j] = zero;
             }
         }
     }

     // Iterate through normals, transpose sample points, and add them to visualization
     int numPoints = originalMesh->GetNumberOfPoints();
     vtkSmartPointer<vtkImplicitPolyDataDistance> implicitPolyDataDistance = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
     implicitPolyDataDistance->SetInput(originalMesh);
     for (vtkIdType i=0;i<numPoints;i+=increment)
     {

         // assume 3-D normal vector
         float *normal = new float[3];
         normalDataFloat->GetTupleValue(i,normal);

         // make the normal a unit vector
         double length = sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
         normal[0] = normal[0] / length;
         normal[1] = normal[1] / length;
         normal[2] = normal[2] / length;

         // a = [0,0,1]
         double v[3] = {-normal[1], normal[0], 0}; // v = a x normal
         double s = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]); // sine = len(v)
         double c = normal[2]; // cos = a dot normal

         double vConst = (1-c) / (s * s);
         double *R = new double[9]; // 3 x 3 matrix, row by row

         // create rotation matrices
         // if vConst is defined (ie, either normal 0 or normal 1 nonzero)
         if (s) {
             R[0] = 1 + vConst * (-v[1] * v[1] - v[2] * v[2]);
             R[1] = -v[2] + vConst * (v[0] * v[1]);
             R[2] = v[1] + vConst * (v[0] * v[2]);
             R[3] = v[2] + vConst * (v[0] * v[1]);
             R[4] = 1 + vConst * (-v[0] * v[0] - v[2] * v[2]);
             R[5] = -v[0] + vConst * (v[1] * v[2]);
             R[6] = -v[1] + vConst * (v[0] * v[2]);
             R[7] = v[0] + vConst * (v[1] * v[2]);
             R[8] = 1 + vConst * (-v[0] * v[0] - v[1] * v[1]);
         }
         // axis is aligned with [0,0,1], so no rotation necessary. Identity matrix
         else{
             R[0] = 1;
             R[1] = 0;
             R[2] = 0;
             R[3] = 0;
             R[4] = 1;
             R[5] = 0;
             R[6] = 0;
             R[7] = 0;
             R[8] = 1;
         }
         double *currentPoint = new double[3];
         int id = i;
         originalMesh->GetPoint(id, currentPoint);
         int addPoint = 1;

         // iterate over sample points, transposing to frame of point in question and then, if all sample points are below some distance
         // from the mesh, add the current point to the set of flat points
         for (int j=0;j<(sampleDimensions*sampleDimensions);j++)
         {
             double *temp = new double[3];
             temp[0] = currentPoint[0] + samplePoints[j][0] * R[0] + samplePoints[j][1] * R[1] + samplePoints[j][2] * R[2];
             temp[1] = currentPoint[1] + samplePoints[j][0] * R[3] + samplePoints[j][1] * R[4] + samplePoints[j][2] * R[5];
             temp[2] = currentPoint[2] + samplePoints[j][0] * R[6] + samplePoints[j][1] * R[7] + samplePoints[j][2] * R[8];
             pointsToVisualize->InsertNextPoint(temp);
         }
     }
   }

    VisualizePointsOnMesh(originalMesh, pointsToVisualize);
}

// check if point is flat by checking closest numNeighbors neighbors
bool FlatHeuristic::IsPointFlat(double *point, int numNeighbors)
{

    vtkSmartPointer<vtkIdList> result =
         vtkSmartPointer<vtkIdList>::New();
    result->Reset();
    kdTree->FindClosestNPoints(numNeighbors,point,result);

    int numHit = 0;
    for(vtkIdType j = 0; j < numNeighbors; j++)
    {
        vtkIdType point_ind = result->GetId(j);
        double p[3];
        originalMesh->GetPoint(point_ind, p);
        if (std::find(flatSet.begin(), flatSet.end(), point_ind) != flatSet.end())
        {
            numHit++;
        }
    }

    return numHit > numNeighbors - numHit;
}

// check if vtkPoints are flat
vtkSmartPointer<vtkPoints> FlatHeuristic::ArePointsFlat (vtkSmartPointer<vtkPoints> pointsToCheck, int numNeighbors)
{
    return ArePointsInSet(originalMesh, flatSet, pointsToCheck, numNeighbors);
}

void FlatHeuristic::SetSampleDimensions(int dim)
{
    sampleDimensions = dim;
}

void FlatHeuristic::SetDiscRadius(double disc)
{
    discRadius = disc;
}

void FlatHeuristic::SetDistanceThreshold(double dist)
{
    distanceThreshold = dist;
}

int FlatHeuristic::GetSampleDimensions()
{
    return sampleDimensions;
}

double FlatHeuristic::GetDiscRadius()
{
    return discRadius;
}

double FlatHeuristic::GetDistanceThreshold()
{
    return distanceThreshold;
}

// Create the set of points that are flat on the mesh
std::vector<vtkIdType> FlatHeuristic::FlatSet ()
{

        // Output
        std::vector<vtkIdType> flatIdVector;

        // Get point normals for mesh
        bool hasPointNormals = GetPointNormals(originalMesh);
        // If normals don't exist - generate them
        if(!hasPointNormals)
        {
            std::cout << "No point normals were found. Computing normals..." << std::endl;

            vtkSmartPointer<vtkPolyDataNormals> normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();
        #if VTK_MAJOR_VERSION <= 5
            normalGenerator->SetInput(originalMesh);
        #else
            normalGenerator->SetInputData(originalMesh);
        #endif
            normalGenerator->ComputePointNormalsOn();
            normalGenerator->ComputeCellNormalsOff();
            normalGenerator->Update();

            originalMesh = normalGenerator->GetOutput();
        }

        // Double normals in an array
        vtkFloatArray* normalDataFloat = vtkFloatArray::SafeDownCast(originalMesh->GetPointData()->GetArray("Normals"));

        // Array creation success
        if(normalDataFloat)
        {

         // create point set of circle sample points
         double *samplePoints [sampleDimensions * sampleDimensions];
         double sampleStep = 2 * discRadius/sampleDimensions;
         for (int i=0;i<sampleDimensions;i++)
         {
             for (int j=0;j<sampleDimensions;j++)
             {
                 double *temp = new double[3];
                 temp[0] = i*sampleStep - discRadius;
                 temp[1] = j*sampleStep - discRadius;
                 temp[2] = 0;
                 if (temp[0] * temp[0] + temp[1] * temp[1] <= discRadius * discRadius)
                 {
                    samplePoints[(i*sampleDimensions) + j] = temp;
                 }
                 else
                 {
                     double *zero = new double[3];
                     zero[0] = 0;
                     zero[1] = 0;
                     zero[2] = 0;
                     samplePoints[(i*sampleDimensions) + j] = zero;
                 }
             }
         }

         // Iterate through normals, transpose sample points, test sum mesh distance against threshold, add to id vector if flat
         int numPoints = originalMesh->GetNumberOfPoints();
         vtkSmartPointer<vtkImplicitPolyDataDistance> implicitPolyDataDistance = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
         implicitPolyDataDistance->SetInput(originalMesh);
         for (vtkIdType i=0;i<numPoints;i++)
         {

             // assume 3-D normal vector
             float *normal = new float[3];
             normalDataFloat->GetTupleValue(i,normal);

             // make the normal a unit vector
             double length = sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
             normal[0] = normal[0] / length;
             normal[1] = normal[1] / length;
             normal[2] = normal[2] / length;

             // a = [0,0,1]
             double v[3] = {-normal[1], normal[0], 0}; // v = a x normal
             double s = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]); // sine = len(v)
             double c = normal[2]; // cos = a dot normal

             double vConst = (1-c) / (s * s);
             double *R = new double[9]; // 3 x 3 matrix, row by row

             // create rotation matrices
             // if vConst is defined (ie, either normal 0 or normal 1 nonzero)
             if (s) {
                 R[0] = 1 + vConst * (-v[1] * v[1] - v[2] * v[2]);
                 R[1] = -v[2] + vConst * (v[0] * v[1]);
                 R[2] = v[1] + vConst * (v[0] * v[2]);
                 R[3] = v[2] + vConst * (v[0] * v[1]);
                 R[4] = 1 + vConst * (-v[0] * v[0] - v[2] * v[2]);
                 R[5] = -v[0] + vConst * (v[1] * v[2]);
                 R[6] = -v[1] + vConst * (v[0] * v[2]);
                 R[7] = v[0] + vConst * (v[1] * v[2]);
                 R[8] = 1 + vConst * (-v[0] * v[0] - v[1] * v[1]);
             }
             // axis is aligned with [0,0,1], so no rotation necessary. Identity matrix
             else{
                 R[0] = 1;
                 R[1] = 0;
                 R[2] = 0;
                 R[3] = 0;
                 R[4] = 1;
                 R[5] = 0;
                 R[6] = 0;
                 R[7] = 0;
                 R[8] = 1;
             }
             double *currentPoint = new double[3];
             int id = i;
             originalMesh->GetPoint(id, currentPoint);
             int addPoint = 1;

             // iterate over sample points, transposing to frame of point in question and then, if all sample points are below some distance
             // from the mesh, add the current point to the set of flat points
             for (int j=0;j<(sampleDimensions*sampleDimensions);j++)
             {
                 double *temp = new double[3];
                 temp[0] = currentPoint[0] + samplePoints[j][0] * R[0] + samplePoints[j][1] * R[1] + samplePoints[j][2] * R[2];
                 temp[1] = currentPoint[1] + samplePoints[j][0] * R[3] + samplePoints[j][1] * R[4] + samplePoints[j][2] * R[5];
                 temp[2] = currentPoint[2] + samplePoints[j][0] * R[6] + samplePoints[j][1] * R[7] + samplePoints[j][2] * R[8];
                 // this line right here is veryyyyyy slow
                 if (fabs(implicitPolyDataDistance->EvaluateFunction(temp)) > distanceThreshold)
                 {
                    addPoint = 0;
                    break;
                 }
             }
             if (addPoint){
                 flatIdVector.push_back(i);
             }
         }
       }

     return flatIdVector;
 }

#endif // FLATHEURISTIC_H
