#ifndef GDL_H
#define GDL_H

#include <vtkVersion.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkLinearSubdivisionFilter.h>
#include <vtkCleanPolyData.h>
#include <vtkKdTreePointLocator.h>
#include <algorithm>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPointSet.h>

bool GetPointNormals(vtkPolyData* polydata);
bool IsPointInSet(vtkSmartPointer<vtkPolyData> originalMesh, std::vector<vtkIdType> flatPoints, double *pointToCheck, int numNeighbors);
vtkSmartPointer<vtkPoints> ArePointsInSet(vtkSmartPointer<vtkPolyData> originalMesh, std::vector<vtkIdType> flatPoints, vtkSmartPointer<vtkPoints> pointsToCheck, int numNeighbors);
vtkSmartPointer<vtkIdList> GetConnectedVertices(vtkSmartPointer<vtkPolyData> mesh, int id);
std::vector<vtkIdType> ConvergeSets(std::vector<std::vector<vtkIdType>*> setOfSets);
vtkSmartPointer<vtkPolyData> DiscretizeMesh (vtkSmartPointer<vtkPolyData> originalMesh, int numberOfSubdivisions);
int SamePoint (const double* p0, const double* p1);
void VisualizePointsOnMesh(vtkSmartPointer<vtkPolyData> mesh,  vtkSmartPointer<vtkPoints> points);

//  Looks for multiple types of normals in a polydata and generates them if they do not exist.
//  Taken from vtk.org
bool GetPointNormals(vtkPolyData* polydata)
{
     std::cout << "In GetPointNormals: " << polydata->GetNumberOfPoints() << std::endl;
     std::cout << "Looking for point normals..." << std::endl;

     // Count points
     vtkIdType numPoints = polydata->GetNumberOfPoints();
     std::cout << "There are " << numPoints << " points." << std::endl;

     // Count triangles
     vtkIdType numPolys = polydata->GetNumberOfPolys();
     std::cout << "There are " << numPolys << " polys." << std::endl;

     ////////////////////////////////////////////////////////////////
     // Double normals in an array
     vtkDoubleArray* normalDataDouble =
     vtkDoubleArray::SafeDownCast(polydata->GetPointData()->GetArray("Normals"));

     if(normalDataDouble)
       {
       int nc = normalDataDouble->GetNumberOfTuples();
       std::cout << "There are " << nc
               << " components in normalDataDouble" << std::endl;
       return true;
       }

     ////////////////////////////////////////////////////////////////
     // Double normals in an array
     vtkFloatArray* normalDataFloat = vtkFloatArray::SafeDownCast(polydata->GetPointData()->GetArray("Normals"));

     if(normalDataFloat)
       {
       int nc = normalDataFloat->GetNumberOfTuples();
       std::cout << "There are " << nc
               << " components in normalDataFloat" << std::endl;
       return true;
       }

     ////////////////////////////////////////////////////////////////
     // Point normals
     vtkDoubleArray* normalsDouble =
       vtkDoubleArray::SafeDownCast(polydata->GetPointData()->GetNormals());

     if(normalsDouble)
       {
       std::cout << "There are " << normalsDouble->GetNumberOfComponents()
                 << " components in normalsDouble" << std::endl;
       return true;
       }

     ////////////////////////////////////////////////////////////////
     // Point normals
     vtkFloatArray* normalsFloat =
       vtkFloatArray::SafeDownCast(polydata->GetPointData()->GetNormals());

     if(normalsFloat)
       {
       std::cout << "There are " << normalsFloat->GetNumberOfComponents()
                 << " components in normalsFloat" << std::endl;
       return true;
       }

     /////////////////////////////////////////////////////////////////////
     // Generic type point normals
     vtkDataArray* normalsGeneric = polydata->GetPointData()->GetNormals(); //works
     if(normalsGeneric)
       {
       std::cout << "There are " << normalsGeneric->GetNumberOfTuples()
                 << " normals in normalsGeneric" << std::endl;

       double testDouble[3];
       normalsGeneric->GetTuple(0, testDouble);

       std::cout << "Double: " << testDouble[0] << " "
                 << testDouble[1] << " " << testDouble[2] << std::endl;

       // Can't do this:
       /*
       float testFloat[3];
       normalsGeneric->GetTuple(0, testFloat);

       std::cout << "Float: " << testFloat[0] << " "
                 << testFloat[1] << " " << testFloat[2] << std::endl;
       */
       return true;
       }

     // If the function has not yet quit, there were none of these types of normals
     std::cout << "Normals not found!" << std::endl;
     return false;
}

// Takes a mesh and returns a discretized version. numberOfSubdivisions is used as input for the VTK subdivison function
vtkSmartPointer<vtkPolyData> DiscretizeMesh (vtkSmartPointer<vtkPolyData> originalMesh, int numberOfSubdivisions){

    vtkSmartPointer<vtkPolyDataAlgorithm> subdivisionFilter;
    subdivisionFilter = vtkSmartPointer<vtkLinearSubdivisionFilter>::New();
    dynamic_cast<vtkLinearSubdivisionFilter *> (subdivisionFilter.GetPointer())->SetNumberOfSubdivisions(numberOfSubdivisions);
    subdivisionFilter->SetInputDataObject(originalMesh);
    subdivisionFilter->Update();
    vtkSmartPointer<vtkPolyData> discreteMesh = vtkSmartPointer<vtkPolyData>::New();
    discreteMesh = subdivisionFilter->GetOutput();

    // Clean new mesh
    vtkCleanPolyData *Clean = vtkCleanPolyData::New();
    Clean->SetInputDataObject(discreteMesh);
    Clean->SetTolerance(0.0);
    Clean->PointMergingOn();
    Clean->Update();
    vtkSmartPointer<vtkPolyData> outputMesh = vtkSmartPointer<vtkPolyData>::New();
    outputMesh = subdivisionFilter->GetOutput();

    return outputMesh;
}

// Check if two points are exactly the same
int SamePoint (const double* p0, const double* p1)
{
   double squaredDistance = vtkMath::Distance2BetweenPoints(p0, p1);
   return (squaredDistance == 0);
}

// Find connected vertices in mesh
// Taken from vtk.org
vtkSmartPointer<vtkIdList> GetConnectedVertices(vtkSmartPointer<vtkPolyData> mesh, int id)
{
 vtkSmartPointer<vtkIdList> connectedVertices =
     vtkSmartPointer<vtkIdList>::New();

 //get all cells that vertex 'id' is a part of
 vtkSmartPointer<vtkIdList> cellIdList =
     vtkSmartPointer<vtkIdList>::New();
 mesh->GetPointCells(id, cellIdList);

 for(vtkIdType i = 0; i < cellIdList->GetNumberOfIds(); i++)
   {

   vtkSmartPointer<vtkIdList> pointIdList =
     vtkSmartPointer<vtkIdList>::New();
   mesh->GetCellPoints(cellIdList->GetId(i), pointIdList);

   if(pointIdList->GetId(0) != id)
     {
     connectedVertices->InsertNextId(pointIdList->GetId(0));
     }
   else
     {
     connectedVertices->InsertNextId(pointIdList->GetId(1));
     }
   }

 return connectedVertices;
}

// Takes a vector of vectors, each which contains vtkIds, and moves all of them to one vector.
// This is used to converge sets calculated by the concave heuristic class
std::vector<vtkIdType> ConvergeSets(std::vector<std::vector<vtkIdType>*> setOfSets)
{
    std::vector<vtkIdType> output;
    for (int i=0;i<setOfSets.size();i++)
    {
        for (int j=0;j<(setOfSets[i])->size();j++)
        {
            output.push_back((*setOfSets[i])[j]);
        }
    }
    return output;
}

// Check if a point should be included in a class of points, defined by a vector of points. This is done by checking the point's
// nearest neighbors in that class and then
bool IsPointInSet(vtkSmartPointer<vtkPolyData> originalMesh, std::vector<vtkIdType> pointSet, double *pointToCheck, int numNeighbors)
{

    // Kd tree to find nearest neighbors
    vtkSmartPointer<vtkKdTreePointLocator> Tree = vtkSmartPointer<vtkKdTreePointLocator>::New();
    Tree->SetDataSet(originalMesh);
    Tree->BuildLocator();

    vtkSmartPointer<vtkIdList> result =
         vtkSmartPointer<vtkIdList>::New();
    result->Reset();
    Tree->Update();
    Tree->FindClosestNPoints(numNeighbors,pointToCheck,result);

    int numHit = 0;
    for(vtkIdType j = 0; j < numNeighbors; j++)
    {
        vtkIdType point_ind = result->GetId(j);
        double p[3];
        originalMesh->GetPoint(point_ind, p);
        if (std::find(pointSet.begin(), pointSet.end(), point_ind) != pointSet.end())
        {
            numHit++;
        }
    }

    return numHit > numNeighbors - numHit;
}

vtkSmartPointer<vtkPoints> ArePointsInSet(vtkSmartPointer<vtkPolyData> originalMesh, std::vector<vtkIdType> pointSet, vtkSmartPointer<vtkPoints> pointsToCheck, int numNeighbors)
{
    vtkSmartPointer<vtkPoints> output = vtkSmartPointer<vtkPoints>::New();
    // Kd tree to find nearest neighbors
    vtkSmartPointer<vtkKdTreePointLocator> Tree = vtkSmartPointer<vtkKdTreePointLocator>::New();
    Tree->SetDataSet(originalMesh);
    Tree->BuildLocator();

    vtkSmartPointer<vtkIdList> result =
         vtkSmartPointer<vtkIdList>::New();
    Tree->Update();
    for (int i=0;i<pointsToCheck->GetNumberOfPoints();i++)
    {
        double *point = new double[3];
        pointsToCheck->GetPoint(i,point);
        result->Reset(); // move back if slows
        Tree->FindClosestNPoints(numNeighbors,point,result);

        int numHit = 0;
        for(vtkIdType j = 0; j < numNeighbors; j++)
        {
            vtkIdType point_ind = result->GetId(j);
            double p[3];
            originalMesh->GetPoint(point_ind, p);
            if (std::find(pointSet.begin(), pointSet.end(), point_ind) != pointSet.end())
            {
                numHit++;
            }
        }
        if (numHit > numNeighbors - numHit)
        {
            output->InsertNextPoint(point);
        }
    }
    return output;
}

// Render points of surface
void VisualizePointsOnMesh(vtkSmartPointer<vtkPolyData> mesh, vtkSmartPointer<vtkPoints> points)
{
         vtkSmartPointer<vtkRenderer> renderer =
           vtkSmartPointer<vtkRenderer>::New();
         vtkSmartPointer<vtkRenderWindow> renderWindow =
           vtkSmartPointer<vtkRenderWindow>::New();
         renderWindow->AddRenderer(renderer);
         vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
           vtkSmartPointer<vtkRenderWindowInteractor>::New();
         renderWindowInteractor->SetRenderWindow(renderWindow);

         vtkSmartPointer<vtkPolyData> pointsPolydata =
           vtkSmartPointer<vtkPolyData>::New();

         pointsPolydata->SetPoints(points);

               vtkSmartPointer<vtkVertexGlyphFilter> vertexFilter =
                  vtkSmartPointer<vtkVertexGlyphFilter>::New();

           #if VTK_MAJOR_VERSION <= 5
             vertexFilter->SetInputConnection(pointsPolydata->GetProducerPort());
           #else
             vertexFilter->SetInputData(pointsPolydata);
           #endif
             vertexFilter->Update();

             vtkSmartPointer<vtkPolyData> polydata =
               vtkSmartPointer<vtkPolyData>::New();
             polydata->ShallowCopy(vertexFilter->GetOutput());

             vtkSmartPointer<vtkPolyDataMapper> mapper2 =
                 vtkSmartPointer<vtkPolyDataMapper>::New();
             #if VTK_MAJOR_VERSION <= 5
               mapper2->SetInputConnection(polydata->GetProducerPort());
             #else
               mapper2->SetInputData(polydata);
             #endif

               vtkSmartPointer<vtkActor> actor2 =
                 vtkSmartPointer<vtkActor>::New();
               actor2->SetMapper(mapper2);
               actor2->GetProperty()->SetPointSize(10);
               renderer->AddActor(actor2);

           // Visualize
           vtkSmartPointer<vtkPolyDataMapper> mapper =
               vtkSmartPointer<vtkPolyDataMapper>::New();
           mapper->SetInputDataObject(mesh);

           vtkSmartPointer<vtkActor> actor =
               vtkSmartPointer<vtkActor>::New();
             actor->SetMapper(mapper);

             renderer->AddActor(actor);

             renderWindow->Render();
             renderWindowInteractor->Start();
}

#endif // GDL_H
