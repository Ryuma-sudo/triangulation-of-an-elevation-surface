#include "filter_qhull.h"
#include "qhull_tools.h"
#include <Eigen/Dense>
#include <vcg/complex/algorithms/convex_hull.h>

using namespace std;
using namespace vcg;

QhullPlugin::QhullPlugin()
{
	typeList = {FP_TRIANGULATION};

	for (ActionIDType tt : types())
		actionList.push_back(new QAction(filterName(tt), this));
}

QhullPlugin::~QhullPlugin()
{
}

QString QhullPlugin::pluginName() const
{
	return "FilterTriangulation";
}

QString QhullPlugin::filterName(ActionIDType f) const
{
	switch (f) {
	case FP_TRIANGULATION:
		return QString("Triangulation: Elevation surface structured by a sparse 3D grid");
	default: assert(0); return QString();
	}
}

QString QhullPlugin::pythonFilterName(ActionIDType filterId) const
{
	switch (filterId) {
	case FP_TRIANGULATION:
		return QString(
			"Method to triangulate an elevation surface structured in a sparse 3D grid by using a "
			"fast search algorithm based on the 2D Delaunay triangulation");
	default: assert(0); return QString();
	}
}

QString QhullPlugin::filterInfo(ActionIDType filterId) const
{
	switch (filterId) {
	case FP_TRIANGULATION:
		return QString(
			"Method to triangulate an elevation surface structured in a sparse 3D grid by using a "
			"fast search algorithm based on the 2D Delaunay triangulation");
	default: assert(0);
	}
	return QString("Error: Unknown Filter");
}

QhullPlugin::FilterClass QhullPlugin::getClass(const QAction* a) const
{
	switch (ID(a)) {
	case FP_TRIANGULATION: return FilterClass(FilterPlugin::Remeshing);
	default: assert(0);
	}
	return FilterClass(0);
}

RichParameterList QhullPlugin::initParameterList(const QAction* action, const MeshModel& m)
{
	RichParameterList parlst;
	switch (ID(action)) {
	case FP_TRIANGULATION:
	default: break; // do not add any parameter for the other filters
	}
	return parlst;
}

/*
Flat the 2d array of coordinate to 1d array by, one point has three element respectively to its axis
Thus, we need to allocate numpoints * dimension * sizeof(coordT) to save the coordinate of all point
*/
coordT* readpointsFromMesh(int* numpoints, int* dimension, MeshModel& m)
{
	coordT *points, *coords;

	coords = points = (coordT*) malloc((*numpoints) * (*dimension) * sizeof(coordT));

	int                    cnt = 0;
	CMeshO::VertexIterator vi;
	for (vi = m.cm.vert.begin(); vi != m.cm.vert.end(); ++vi)
		if (!(*vi).IsD()) {
			for (int ii = 0; ii < *dimension; ++ii)
				*(coords++) = (*vi).P()[ii];
			++cnt;
		}
	assert(cnt == m.cm.vn);

	return (points);
}

const int    dX[] = {0, 0, 1, -1};
const int    dY[] = {1, -1, 0, 0};
const coordT EPS  = 1e-9;

struct MyPoint
{
	coordT x, y;
};

struct MyEdge
{
	MyPoint a, b;
	// adjacent vertex;
	MyPoint adjLeft, adjRight = {-1.0, -1.0};

	// Define equality operator for MyEdge struct
	bool operator==(const MyEdge& other) const
	{
		return (a.x == other.a.x && a.y == other.a.y && b.x == other.b.x && b.y == other.b.y) ||
			   (a.x == other.b.x && a.y == other.b.y && b.x == other.a.x && b.y == other.a.y);
	}
};

bool isInCircumcircle(MyPoint P1, MyPoint P2, MyPoint P3, MyPoint P)
{
	coordT centerX, centerY;
	coordT a1, b1, d1, a2, b2, d2;
	coordT radius, checkDistance;

	// Make the below calculation easier
	a1 = -2.0 * P1.x + 2.0 * P2.x;
	b1 = -2.0 * P1.y + 2.0 * P2.y;
	d1 = -P1.x * P1.x - P1.y * P1.y + P2.x * P2.x + P2.y * P2.y;
	a2 = -2.0 * P1.x + 2.0 * P3.x;
	b2 = -2.0 * P1.y + 2.0 * P3.y;
	d2 = -P1.x * P1.x - P1.y * P1.y + P3.x * P3.x + P3.y * P3.y;

	// Using determinant to calculate the coordinate x, and y of center point
	centerX = (d1 * b2 - d2 * b1) / (a1 * b2 - a2 * b1);
	centerY = (d2 * a1 - d1 * a2) / (a1 * b2 - a2 * b1);
	radius  = sqrt((centerX - P1.x) * (centerX - P1.x) + (centerY - P1.y) * (centerY - P1.y));

	// Distance between check point and center point
	checkDistance = sqrt((centerX - P.x) * (centerX - P.x) + (centerY - P.y) * (centerY - P.y));
	return checkDistance < radius;
}

bool isLeft(MyPoint A, MyPoint B, MyPoint P)
{
	return (B.x - A.x) * (P.y - A.y) - (B.y - A.y) * (P.x - A.x) >= 0;
}

MyPoint centerOf(MyEdge edge)
{
	MyPoint M = {(edge.a.x + edge.b.x) / 2, (edge.a.y + edge.b.y) / 2};

	return{{M.x - (edge.a.y - edge.b.y) / 2}, {M.y + (edge.a.x - edge.b.x) / 2}};
}

std::vector<MyPoint> inSquare(MyEdge edge, MyPoint center, coordT maxX, coordT maxY)
{
	vector<MyPoint> list, finalList;
	coordT          radius = sqrt((edge.a.x - center.x) * (edge.a.x - center.x) +
					(edge.a.y - center.y) * (edge.a.y - center.y));

	// find all point inside the circle

	for (int x = (int) (- radius + center.x); x <= (int) (radius + center.x) + 1; x++) {
		if (x < 0 || x > (int) maxX) continue;
		for (int y = (int) (- radius + center.y); y <= (int) (radius + center.y) + 1; y++) {
			if (y < 0 || y > (int) maxY) continue;
			if ((x - center.x) * (x - center.x) + (y - center.y) * (y - center.y) <= radius * radius) {
				list.push_back({(coordT) x ,(coordT) y});
			}
		}
	}

	// remove all point outside the square
	
	// coordinates of the square

	MyPoint A = edge.b;
	MyPoint B = edge.a;
	MyPoint C = {2 * center.x - A.x, 2 * center.y - A.y};
	MyPoint D = {2 * center.x - B.x, 2 * center.y - B.y};

	// check all point

	for (int i = 0; i < list.size(); i++) {
		MyPoint p = list[i];
		// remove if it's outside the square or the same as one of the edge's vertex 
		if (!isLeft(A, B, p) || !isLeft(B, C, p) || !isLeft(C, D, p) || !isLeft(D, A, p) ||
			((p.x == edge.a.x && p.y == edge.a.y) || (p.x == edge.b.x && p.y == edge.b.y))) {
			continue;
		}
		finalList.push_back(list[i]);
	}
	return finalList;
}

bool specialEdge(MyEdge edge, coordT maxX,coordT maxY)
{
	MyEdge sEdge1 = {{0, 1}, {1, 0}};
	MyEdge sEdge2 = {{maxX - 1, 0}, {maxX, 1}};
	MyEdge sEdge3 = {{0, maxY - 1}, {1, maxY}};
	MyEdge sEdge4 = {{maxX - 1, maxY}, {maxX, maxY - 1}};
	if (edge == sEdge1 || edge == sEdge2 || edge == sEdge3 || edge == sEdge4) {
		return true;
	}
	return false;
}

std::map<std::string, QVariant> QhullPlugin::applyFilter(
	const QAction*           filter,
	const RichParameterList& par,
	MeshDocument&            md,
	unsigned int& /*postConditionMask*/,
	vcg::CallBackPos* /* cb*/)
{
	qhT  qh_qh = {};
	qhT* qh    = &qh_qh;

	switch (ID(filter)) {
	case FP_TRIANGULATION: {
		MeshModel& m         = *md.mm();
		MeshModel& nm        = *md.addNewMesh("", "Triangulation 2d");
		MeshModel& fm        = *md.addNewMesh("", "Final Triangulation");
		int        dim       = 3;
		int        numpoints = m.cm.vn;
		coordT*    points;
		points = readpointsFromMesh(&numpoints, &dim, m);

		// vector<MyPoint> list1 = inSquare({{1, 0}, {1, 1}}, centerOf({{1, 0}, {1, 1}}));

		// for (MyPoint p : list1) {
		//	log("%f, %f\n", p.x, p.y);
		// }

		// vector<MyPoint> list2 = inSquare({{1, 1}, {0, 1}}, centerOf({{1, 1}, {0, 1}}));

		// for (MyPoint p : list2) {
		//	log("%f, %f\n", p.x, p.y);
		// }

		coordT maxX = -1.0, maxY = -1.0;
		for (int i = 0; i < numpoints; i++) {
			maxX = max(maxX, points[i * dim]);
			maxY = max(maxY, points[i * dim + 1]);
		}
		int row = (int) maxX;
		int col = (int) maxY;
		row++, col++;
		vector<vector<int>>    lookup, boundary;
		vector<vector<coordT>> zAxis;
		lookup.resize(row);
		boundary.resize(row);
		zAxis.resize(row);
		for (int i = 0; i < row; i++) {
			lookup[i].resize(col, 0);
			boundary[i].resize(col, 0);
			zAxis[i].resize(col, 0.0);
		}

		for (int i = 0; i < numpoints; i++) {
			int curX           = (int) points[i * dim];
			int curY           = (int) points[i * dim + 1];
			lookup[curX][curY] = 1;
			if (curX == 0 || curX == row - 1 || curY == 0 || curY == col - 1)
				boundary[curX][curY] = 1;
		}

		 for (int i = 0; i < numpoints; i++) {
			int curX = (int) points[i * dim];
			int curY = (int) points[i * dim + 1];
			if (boundary[curX][curY] == 0) {
				for (int j = 0; j < 4; j++) {
					int nxtX = curX + dX[j];
					int nxtY = curY + dY[j];
					if (nxtX >= 0 && nxtX <= row - 1 && nxtY >= 0 && nxtY <= col - 1)
						if (!lookup[nxtX][nxtY]) {
							boundary[curX][curY] = 2;
							break;
						}
				}
			}
		 }

		// Convert 3d surface to 2d grid by eliminating the z coordinate

		for (int i = 0; i < numpoints; i++) {
			Point3d curVertex = {points[i * dim], points[i * dim + 1], 0};
			coordT  curX      = points[i * dim];
			coordT  curY      = points[i * dim + 1];
			// Use sparse matrix as a lookup table
			zAxis[(int) curX][(int) curY] = points[i * dim + 2];
			tri::Allocator<CMeshO>::AddVertex(nm.cm, curVertex);
		}

		/* Create seed triangle */
		// (0, 0) (0, 1) and (1, 0)

		Point3d p0 = {0, 0, 0};
		Point3d p1 = {0, 1, 0};
		Point3d p2 = {1, 0, 0};

		// Edges of this seed triangle

		MyEdge seedE3 = {{1, 0}, {0, 1}, {0, 0}};
		tri::Allocator<CMeshO>::AddFace(nm.cm, p0, p1, p2);

		/* Triangulating a surface */

		queue<MyEdge>  EdgePool;
		vector<MyEdge> visited;
		// Add edge of seed triangle into edges pool

		EdgePool.push(seedE3);
		int cnt = 0;
		while (!EdgePool.empty()) {
			 //if (cnt == 8)
				//break;
			 //cnt++;
			//  Take the front edge
			MyEdge curEdge = EdgePool.front();
			EdgePool.pop();

			bool isVisited = false;
			for (int i = 0; i < visited.size(); i++) {
				if (curEdge.a.x == visited[i].b.x && curEdge.a.y == visited[i].b.y &&
					curEdge.b.x == visited[i].a.x && curEdge.b.y == visited[i].a.y) {
					isVisited = true;
					break;
				}
			}
			if (isVisited) {
				continue;
			}
			log("%f, %f, %d **", curEdge.a.x, curEdge.a.y, boundary[curEdge.a.x][curEdge.a.y]);
			log("%f, %f, %d **", curEdge.b.x, curEdge.b.y, boundary[curEdge.b.x][curEdge.b.y]);

			// only consider if the edge is not the boundary edge
			if (!(boundary[curEdge.a.x][curEdge.a.y] == 1 && boundary[curEdge.b.x][curEdge.b.y] == 1) ||
				!(boundary[curEdge.a.x][curEdge.a.y] == 2 && boundary[curEdge.b.x][curEdge.b.y] == 2) ||
				specialEdge(curEdge, maxX, maxY)) {
				// find the right adjacent vertex of curEdge

				vector<MyPoint> neighborPoints = inSquare(curEdge, centerOf(curEdge), maxX, maxY);

				// find the best point to form new triangle

				coordT bestCompact = -1e18;
				coordT bestPX = -1.0, bestPY = -1.0;
				for (int i = 0; i < neighborPoints.size(); i++) {
					MyPoint curPoint = neighborPoints[i];
					// Check if this point achieve delaunay criteration
					if (!lookup[(int) neighborPoints[i].x][(int) neighborPoints[i].y] &&
						boundary[neighborPoints[i].x][neighborPoints[i].y] != 2) {
						continue;
					}
					bool isFalse = false;

					for (int ii = 0; ii < neighborPoints.size(); ii++) {
						if (!lookup[(int) neighborPoints[ii].x][(int) neighborPoints[ii].y] &&
							boundary[neighborPoints[i].x][neighborPoints[i].y] != 2) {
							continue;
						}
						log("%f, %f *", curEdge.a.x, curEdge.a.y);
						log("%f, %f *", curEdge.b.x, curEdge.b.y);
						log("%f, %f *", curPoint.x, curPoint.y);
						log("%f, %f *", neighborPoints[ii].x, neighborPoints[ii].y);
						log("");
						if (isInCircumcircle(curEdge.a, curEdge.b, curPoint, neighborPoints[ii])) {
							isFalse = true;
							break;
						}
					}

					if (!isFalse) {
						coordT l0 = sqrt(
							(curEdge.a.x - curEdge.b.x) * (curEdge.a.x - curEdge.b.x) +
							(curEdge.a.y - curEdge.b.y) * (curEdge.a.y - curEdge.b.y));
						coordT l1 = sqrt(
							(curEdge.a.x - curPoint.x) * (curEdge.a.x - curPoint.x) +
							(curEdge.a.y - curPoint.y) * (curEdge.a.y - curPoint.y));
						coordT l2 = sqrt(
							(curEdge.b.x - curPoint.x) * (curEdge.b.x - curPoint.x) +
							(curEdge.b.y - curPoint.y) * (curEdge.b.y - curPoint.y));
						coordT p          = (l0 + l1 + l2) / 2.0;
						coordT A          = sqrt(p * (p - l0) * (p - l1) * (p - l2));
						coordT curCompact = (4.0 * sqrt(3.0) * A) / (l0 * l0 + l1 * l1 + l2 * l2);
						log("compact: %f", curCompact);
						if (curCompact > bestCompact) {
							bestCompact = curCompact;
							bestPX      = curPoint.x;
							bestPY      = curPoint.y;
						}
					}
					else
						log("false 1");
				}
				neighborPoints.clear();

				// Add new triangle
				if (bestPX != -1.0 && bestPY != -1.0) {
					// New edge
					MyEdge newEdge1;
					newEdge1.a       = curEdge.a;
					newEdge1.b       = {bestPX, bestPY};
					newEdge1.adjLeft = curEdge.b;
					MyEdge newEdge2;
					newEdge2.a       = {bestPX, bestPY};
					newEdge2.b       = curEdge.b;
					newEdge2.adjLeft = curEdge.a;
					bool isVisited   = false;

					EdgePool.push(newEdge1);
					EdgePool.push(newEdge2);
					visited.push_back(newEdge1);
					visited.push_back(newEdge2);
					Point3d p0 = {curEdge.a.x, curEdge.a.y, 0.0};
					Point3d p1 = {curEdge.b.x, curEdge.b.y, 0.0};
					Point3d p2 = {bestPX, bestPY, 0.0};
					log("%f, %f", bestPX, bestPY);

					tri::Allocator<CMeshO>::AddFace(nm.cm, p0, p1, p2);

				}
				else
					log("false 2");
			}
			else
				log("false 3");
		}

		// Delete duplicate vertices and faces
		tri::Clean<CMeshO>::RemoveDuplicateVertex(nm.cm);
		tri::Clean<CMeshO>::RemoveDuplicateFace(nm.cm);

		// Restore z
		CMeshO::FaceIterator fi;
		for (fi = nm.cm.face.begin(); fi != nm.cm.face.end(); fi++) {
			if (!(*fi).IsD()) {
				coordT cur1X, cur1Y, cur2X, cur2Y, cur3X, cur3Y, cur1Z, cur2Z, cur3Z;
				cur1X = (*fi).P(0)[0];
				cur1Y = (*fi).P(0)[1];
				cur1Z = zAxis[cur1X][cur1Y];
				cur2X = (*fi).P(1)[0];
				cur2Y = (*fi).P(1)[1];
				cur2Z = zAxis[cur2X][cur2Y];
				cur3X = (*fi).P(2)[0];
				cur3Y = (*fi).P(2)[1];
				cur3Z = zAxis[cur3X][cur3Y];

				Point3d p0 = {cur1X, cur1Y, cur1Z};
				Point3d p1 = {cur2X, cur2Y, cur2Z};
				Point3d p2 = {cur3X, cur3Y, cur3Z};

				tri::Allocator<CMeshO>::AddFace(fm.cm, p0, p1, p2);
			}
		}

		// Delete duplicate vertices and faces
		tri::Clean<CMeshO>::RemoveDuplicateVertex(fm.cm);
		tri::Clean<CMeshO>::RemoveDuplicateFace(fm.cm);

	} break;
	default: wrongActionCalled(filter);
	}
	return std::map<std::string, QVariant>();
}
MESHLAB_PLUGIN_NAME_EXPORTER(QhullPlugin)
