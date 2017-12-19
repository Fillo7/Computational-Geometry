// Author: Filip Petrovic (422334)

// Voronoi diagram variables
boolean showVoronoi;
ArrayList<Edge> voronoiEdges;
ArrayList<Point> voronoiPoints;

Edge getOutgoingVoronoiEdge(Edge triangleEdge, Triangle triangle)
{
    Point intersection = triangleEdge.getNearestPoint(triangle.circleCenter);
    PVector vector = new PVector(intersection.x - triangle.circleCenter.x, intersection.y - triangle.circleCenter.y);
    vector.mult(100.0f);
    
    Edge test = new Edge(new Point(triangle.circleCenter.x, triangle.circleCenter.y), new Point(vector.x, vector.y));
    if (!triangle.pointInside(triangle.circleCenter))
    {
        Edge closestIntersecting = test.getClosestIntersectingEdge(delaunayEdges);
        if (closestIntersecting != null && closestIntersecting.sharesPoints(triangleEdge))
        {
            vector.mult(-1.0f);
        }
    }
    
    return new Edge(new Point(triangle.circleCenter.x, triangle.circleCenter.y), new Point(vector.x, vector.y));
}

void calculateVoronoi(ArrayList<Triangle> triangles)
{
    for (Triangle triangle : triangles)
    {
        for (Triangle other : triangles)
        {
            if (triangle.equals(other))
            {
                continue;
            }
            
            Edge shared = triangle.getSharedEdge(other);
            if (shared == null)
            {
                continue;
            }
            
            shared.onConvexHull = false;
            voronoiPoints.add(new Point(triangle.circleCenter.x, triangle.circleCenter.y));
            voronoiPoints.add(new Point(other.circleCenter.x, other.circleCenter.y));
            voronoiEdges.add(new Edge(new Point(triangle.circleCenter.x, triangle.circleCenter.y), new Point(other.circleCenter.x, other.circleCenter.y)));
        }
    }
    
    for (Triangle triangle : triangles)
    {
        if (!triangle.circleCenter.onScreen())
        {
            continue;
        }
      
        if (triangle.first.onConvexHull)
        {
            Edge outgoing = getOutgoingVoronoiEdge(triangle.first, triangle);
            voronoiEdges.add(outgoing);
        }
        
        if (triangle.second.onConvexHull)
        {
            Edge outgoing = getOutgoingVoronoiEdge(triangle.second, triangle);
            voronoiEdges.add(outgoing);
        }
        
        if (triangle.third.onConvexHull)
        {
            Edge outgoing = getOutgoingVoronoiEdge(triangle.third, triangle);
            voronoiEdges.add(outgoing);
        }
    }
}

void drawVoronoi(PVector voronoiColor)
{
    stroke(voronoiColor.x, voronoiColor.y, voronoiColor.z);
    fill(voronoiColor.x, voronoiColor.y, voronoiColor.z);
    
    for (Point point : voronoiPoints)
    {
        ellipse(point.x, point.y, pointRadius * 1.5f, pointRadius * 1.5f);
    }
    
    for (Edge edge: voronoiEdges)
    {
        line(edge.start.x, edge.start.y, edge.end.x, edge.end.y);
    }
    
    stroke(pointColor);
    fill(pointColor);
}