#include <math.h>
#include <stdlib.h>
#include <vtkCellArray.h>
#include <vtkProperty.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolygon.h>
#include <vtkSmartPointer.h>
#include <vtkPointSet.h>
#include <vtkMath.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkCleanPolyData.h>
#include <vtkDelaunay3D.h>
#include <vtkPolyDataReader.h>
#include <vtkTriangle.h>
#include <vtkTriangleFilter.h>
#include <vtkPLYReader.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkPLYWriter.h>
#include <vtkLinearSubdivisionFilter.h>
#include <vtkPolyDataNormals.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkImplicitPolyDataDistance.h>
#include <vtkTuple.h>

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkTriangle.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkVertexGlyphFilter.h>

#include "concaveheuristic.h"
#include "flatheuristic.h"
#include "gdl.h"
#include <time.h>

// Usage: first argument is input mesh file,
int main (int argc, char** argv)
{

    // Read in meshes //

    //Read the file
   vtkSmartPointer<vtkPLYReader> reader = vtkSmartPointer<vtkPLYReader>::New();
   reader->SetFileName(argv[1]);

   vtkSmartPointer<vtkPolyData> originalMesh = vtkSmartPointer<vtkPolyData>::New();
   vtkSmartPointer<vtkPolyData> testMesh = vtkSmartPointer<vtkPolyData>::New();

   // Subdivision filters only work on triangles
   vtkSmartPointer<vtkTriangleFilter> triangles = vtkSmartPointer<vtkTriangleFilter>::New();
   triangles->SetInputConnection(reader->GetOutputPort());
   triangles->Update();
   originalMesh = triangles->GetOutput();
   testMesh->DeepCopy(originalMesh);
   testMesh = DiscretizeMesh(testMesh,1); // Discretize test mesh so it has more (unseen) points



   // Code example //

   clock_t t1,t2;




   // Concave Class //

   // Initialize with new mesh
   ConcaveHeuristic *ch = new ConcaveHeuristic(originalMesh);

   // Generate concave points (timed)
   t1=clock();
   ch->GenerateConcaveSet(); // HERE
   t2=clock();
   float diff = ((float)t2-(float)t1);
   float seconds = diff / CLOCKS_PER_SEC;
   cout << "Generating concave points: " << seconds << " seconds" << endl;

   // Display concave points
   cout << "Now showing points labeled concave on original mesh" << endl;
   ch->VisualizeConcavePoints();

   // Test points in test mesh (timed)
   vtkSmartPointer<vtkPoints> concavePoints = vtkSmartPointer<vtkPoints>::New();
   seconds = 0;
   for (int i=0;i<testMesh->GetNumberOfPoints();i++)
   {
       double *point = new double[3];
       testMesh->GetPoint(i,point);
       t1=clock();
       bool temp = ch->IsPointConcave(point); // HERE is the check
       t2=clock();
       diff = ((float)t2-(float)t1);
       seconds += diff / CLOCKS_PER_SEC;
       if(temp)
       {
           concavePoints->InsertNextPoint(point[0],point[1],point[2]);
       }
   }
   cout << "Testing " << testMesh->GetNumberOfPoints() << " points for concavity:" << seconds << " seconds" << endl;

   // Visualize Points labeled concave
   cout << "Now showing points labeled concave on test mesh" << endl;
   VisualizePointsOnMesh(originalMesh, concavePoints);






   // Flat Class //

   // Initialize with new mesh
   FlatHeuristic *fh = new FlatHeuristic(DiscretizeMesh(originalMesh,1));
   testMesh = DiscretizeMesh(testMesh,1); // Discretize test mesh so it has more (unseen) points

   // Generate flat points (timed)
   t1=clock();
   fh->GenerateFlatSet(); // HERE
   t2=clock();
   diff = ((float)t2-(float)t1);
   seconds = diff / CLOCKS_PER_SEC;
   cout << "Generating flat points: " << seconds << " seconds (AFTER DISCRETIZING: IE MORE POINTS)" << endl;

   // Display flat points
   cout << "Now showing sampling points for flat point generation on original mesh" << endl;
   fh->VisualizeSampledPoints();

   // Display flat points
   cout << "Now showing points labeled flat on original mesh" << endl;
   fh->VisualizeFlatPoints();

   // Test points in test mesh (timed)
   vtkSmartPointer<vtkPoints> flatPoints = vtkSmartPointer<vtkPoints>::New();
   seconds = 0;
   for (int i=0;i<testMesh->GetNumberOfPoints();i++)
   {
       double *point = new double[3];
       testMesh->GetPoint(i,point);
       t1=clock();
       bool temp = fh->IsPointFlat(point); // HERE is the check
       t2=clock();
       diff = ((float)t2-(float)t1);
       seconds += diff / CLOCKS_PER_SEC;
       if(temp)
       {
           flatPoints->InsertNextPoint(point[0],point[1],point[2]);
       }
   }
   cout << "Testing " << testMesh->GetNumberOfPoints() << " points for flatness:" << seconds << " seconds" << endl;

   // Visualize Points labeled concave
   cout << "Now showing points labeled flat on test mesh" << endl;
   VisualizePointsOnMesh(originalMesh, flatPoints);

    return (0);
}

// NOTES: what is slowing down the flat set calculation is the distance from point to mesh VTK function.
