#ifndef CONCAVEHEURISTIC_H
#define CONCAVEHEURISTIC_H

#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkDelaunay3D.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkPoints.h>
#include "gdl.h"

#define Default_Concave_Distance_Threshold 50

class ConcaveHeuristic
{ 
    public:
        ConcaveHeuristic(vtkSmartPointer<vtkPolyData> newMesh, double dist=Default_Concave_Distance_Threshold);
        void GenerateConcaveSet();
        double GetDistThresh();
        double SetDistThresh(double dist);
        std::vector<std::vector<vtkIdType>*> GetConcaveSets();
        bool IsPointConcave(double *point, int numNeighbors = 10);
        vtkSmartPointer<vtkPoints> ArePointsConcave (vtkSmartPointer<vtkPoints>, int numNeighbors = 10);
        void VisualizeConcavePoints();
    private:
        vtkSmartPointer<vtkPolyData> originalMesh;
        vtkSmartPointer<vtkKdTreePointLocator> kdTree;
        double distThresh;
        std::vector<std::vector<vtkIdType>*> concaveSets;
        std::vector<std::vector<vtkIdType>*> ConcaveSets();
        std::vector<vtkIdType> convergedConcaveSet;
};

ConcaveHeuristic::ConcaveHeuristic(vtkSmartPointer<vtkPolyData> newMesh, double dist)
{
    originalMesh = vtkSmartPointer<vtkPolyData>::New();
    originalMesh->DeepCopy(newMesh);
    distThresh = dist;
//    concaveSets = new std::vector<std::vector<vtkIdType>*>();
//    convergedConcaveSet = new std::vector<vtkIdType>();
    kdTree = vtkSmartPointer<vtkKdTreePointLocator>::New();
    kdTree->SetDataSet(originalMesh);
    kdTree->BuildLocator();
    kdTree->Update();
}

// Render mesh and concave points
void ConcaveHeuristic::VisualizeConcavePoints()
{
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for (int i=0;i<convergedConcaveSet.size();i++)
    {
        double *point = new double[3];
        originalMesh->GetPoint(convergedConcaveSet[i],point);
        points->InsertNextPoint(point);
    }
    VisualizePointsOnMesh(originalMesh, points);
}

// Generate originalMesh's concave sets
void ConcaveHeuristic::GenerateConcaveSet()
{
    concaveSets = ConcaveSets();
    convergedConcaveSet = ConvergeSets(concaveSets);
}

// check if point is concave by checking closest numNeighbors neighbors
bool ConcaveHeuristic::IsPointConcave(double *point, int numNeighbors)
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
        if (std::find(convergedConcaveSet.begin(), convergedConcaveSet.end(), point_ind) != convergedConcaveSet.end())
        {
            numHit++;
        }
    }

    return numHit > numNeighbors - numHit;
}

vtkSmartPointer<vtkPoints> ConcaveHeuristic::ArePointsConcave (vtkSmartPointer<vtkPoints> pointsToCheck, int numNeighbors)
{
    return ArePointsInSet(originalMesh, convergedConcaveSet, pointsToCheck, numNeighbors);
}

double ConcaveHeuristic::GetDistThresh()
{
    return distThresh;
}

double ConcaveHeuristic::SetDistThresh(double dist)
{
    distThresh = dist;
}

std::vector<std::vector<vtkIdType>*> ConcaveHeuristic::GetConcaveSets()
{
    return concaveSets;
}

// Finds the sets of concave points on the originalMesh. This is done by testing a point to see if it is
// a certain distance (distThresh) away from the convex hull of the mesh.
std::vector<std::vector<vtkIdType>*> ConcaveHeuristic::ConcaveSets()
{
     // Generate a tetrahedral mesh from the input points. By
     // default, the generated volume is the convex hull of the points.
     vtkSmartPointer<vtkDelaunay3D> delaunay = vtkSmartPointer<vtkDelaunay3D>::New();
     delaunay->SetInputDataObject(originalMesh);
     delaunay->Update();
     vtkSmartPointer<vtkDataSetSurfaceFilter> surfaceFilter = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
     surfaceFilter->SetInputConnection(delaunay->GetOutputPort());
     surfaceFilter->Update();
     vtkSmartPointer<vtkPolyData> convexMesh = vtkSmartPointer<vtkPolyData>::New();
     convexMesh = surfaceFilter->GetOutput();

     // Find difference between original mesh and convex mesh to compute concave point set
     std::vector<double*> concavePoints;
     for (int i=0;i<originalMesh->GetNumberOfPoints();i++){

         // Get closest point in convex mesh
         vtkIdType closestID = convexMesh->FindPoint(originalMesh->GetPoint(i));
         double *p0 = new double[3];
         double *p1 = new double[3];
         convexMesh->GetPoint(closestID,p0);
         originalMesh->GetPoint(i,p1);

         // Find the squared distance between the points.
         double squaredDistance = vtkMath::Distance2BetweenPoints(p0, p1);

         // Threshold distance
         if (squaredDistance > distThresh){
              concavePoints.push_back(p1);
         }
     }

    // Iterate through concave points, expanding and building sets
    std::vector<std::vector<vtkIdType>*> connectionsVector;
    while(concavePoints.size() > 0){

        // Creating a new connection set
        std::vector<vtkIdType> *connectionSet = new std::vector<vtkIdType>();

        // get new start point by popping from convex set and get connections
        double *firstPoint = concavePoints.back();
        vtkIdType id = originalMesh->FindPoint(firstPoint);
        concavePoints.pop_back();
        vtkSmartPointer<vtkIdList> connections = GetConnectedVertices(originalMesh, id);
        connectionSet->push_back(id);

        // add each connection still in connection set to vector to be checked and pop from connection set
        std::vector<vtkIdType> checkList;
        for(int i=0; i<connections->GetNumberOfIds(); i++){

            double *originalMeshPoint = new double[3];
            originalMesh->GetPoint(connections->GetId(i), originalMeshPoint);

            // check to see if point in convex set
            int isConvex = 0;
            int concaveIndex = 0;
            for (int j=0;j<concavePoints.size();j++)
            {
              if (SamePoint(originalMeshPoint, concavePoints[j]))
              {
                  isConvex = 1;
                  concaveIndex = j;
                  break;
              }
            }
            if (isConvex)
            {
                  connectionSet->push_back(connections->GetId(i));
                  checkList.push_back(connections->GetId(i));
                  concavePoints.erase(concavePoints.begin() + concaveIndex);
            }
        }

        // while points to be checked vector size > 0 (expanding connections)
        // find point's connections, add them to connection set to be checked, and pop them from convex points to be checked vector (same as above)
        while (checkList.size()>0){
            connections = GetConnectedVertices(originalMesh,checkList.back());
            checkList.pop_back();

            for(int i=0; i<connections->GetNumberOfIds(); i++){
                double *originalMeshPoint = new double[3];
                originalMesh->GetPoint(connections->GetId(i), originalMeshPoint);

                // check to see if point in convex point set
                int isConvex = 0;
                int concaveIndex = 0;
                for (int j=0;j<concavePoints.size();j++)
                {
                  if (SamePoint(originalMeshPoint, concavePoints[j]))
                  {
                      isConvex = 1;
                      concaveIndex = j;
                      break;
                  }
                }
                if (isConvex)
                {
                      connectionSet->push_back(connections->GetId(i));
                      checkList.push_back(connections->GetId(i));
                      concavePoints.erase(concavePoints.begin() + concaveIndex);
                }
            }
        }

        // add connectionSet to connectionsVector
        connectionsVector.push_back(connectionSet);
    }
    return connectionsVector;
 }

#endif // CONCAVEHEURISTIC_H
