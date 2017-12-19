// Author: Filip Petrovic (422334)

// Delaunay triangulation variables
boolean showDelaunay;
ArrayList<Edge> delaunayEdges;
ArrayList<Triangle> delaunayTriangles;

Circle getCircumscribedCircle(Point first, Point second, Point third)
{
    float cross = crossProduct(first, second, third);
    
    if (cross == 0.0f)
    {
        return null; 
    }
    
    float firstSq = first.x * first.x + first.y * first.y;
    float secondSq = second.x * second.x + second.y * second.y;
    float thirdSq = third.x * third.x + third.y * third.y;
    
    float helperX = firstSq * (second.y - third.y) + secondSq * (third.y - first.y) + thirdSq * (first.y - second.y);
    float centerX = helperX / (2.0f * cross);
    
    float helperY = firstSq * (third.x - second.x) + secondSq * (first.x - third.x) + thirdSq * (second.x - first.x);
    float centerY = helperY / (2.0f * cross);
    
    Point center = new Point(centerX, centerY);
    float radius = getDistance(first, center);
    
    return new Circle(center, radius);
}

Point getClosestPoint(Point target, ArrayList<Point> points)
{
    Point result = null;
    float distance = Float.MAX_VALUE;
    
    for (Point point : points)
    {
        float currentDistance = getDistance(target, point);
        if (!point.equals(target) && currentDistance < distance)
        {
            result = point;
            distance = currentDistance;
        }
    }
    
    return result;
}

Point getClosestDelaunayDistancePoint(Edge target, ArrayList<Point> points)
{
    Point result = null;
    float bestDistance = Float.MAX_VALUE;
    
    for (Point point : points)
    {
        if (getOrientation(target.start, target.end, point) >= 0.0f)
        {
            continue;
        }
        
        Circle circle = getCircumscribedCircle(target.start, target.end, point);
        float distance = getDistance(circle.center, point);
        
        if (getOrientation(target.start, target.end, circle.center) > 0.0f)
        {
            distance = -distance;
        }
        
        if (distance < bestDistance)
        {
            result = point;
            bestDistance = distance;
            target.circleCenter = circle.center;
        }
    }
    
    return result;
}

void triangulateDelaunay(ArrayList<Point> polygon)
{
    Queue<Edge> activeEdgeList = new LinkedList<Edge>();
    
    Point first = polygon.get(0);
    Point second = getClosestPoint(first, polygon);
    Edge firstEdge = new Edge(first, second);
    
    Point third = getClosestDelaunayDistancePoint(firstEdge, polygon);
    
    if (third == null)
    {
        firstEdge.swapOrientation();
        third = getClosestDelaunayDistancePoint(firstEdge, polygon);
    }
    
    Edge secondEdge = new Edge(firstEdge.end, third);
    Edge thirdEdge = new Edge(third, firstEdge.start);
    
    activeEdgeList.add(firstEdge);
    activeEdgeList.add(secondEdge);
    activeEdgeList.add(thirdEdge);
    delaunayTriangles.add(new Triangle(firstEdge, secondEdge, thirdEdge, firstEdge.circleCenter));
    
    while (!activeEdgeList.isEmpty())
    {
        Edge firstNew = activeEdgeList.poll();
        firstNew.swapOrientation();
        
        Point point = getClosestDelaunayDistancePoint(firstNew, polygon);
        if (point != null)
        {
            Edge secondNew = new Edge(firstNew.end, point);
            Edge thirdNew = new Edge(point, firstNew.start);
            delaunayTriangles.add(new Triangle(firstNew, secondNew, thirdNew, firstNew.circleCenter));
            
            Edge secondReverse = new Edge(point, firstNew.end);
            if (!activeEdgeList.contains(secondNew) && !activeEdgeList.contains(secondReverse) && !delaunayEdges.contains(secondNew) && !delaunayEdges.contains(secondReverse))
            {
                activeEdgeList.add(secondNew);
            }
            
            Edge thirdReverse = new Edge(firstNew.start, point);
            if (!activeEdgeList.contains(thirdNew) && !activeEdgeList.contains(thirdReverse) && !delaunayEdges.contains(thirdNew) && !delaunayEdges.contains(thirdReverse))
            {
                activeEdgeList.add(thirdNew);
            }
        }
        
        delaunayEdges.add(firstNew);
    }
}