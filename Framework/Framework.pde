// Author: Filip Petrovic (422334)

import java.awt.geom.Line2D;
import java.awt.Polygon;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.Objects;
import java.util.Stack;
import java.util.Queue;

enum Path
{
    Left,
    Right
};

public class Point
{
    public Point(float x, float y)
    {
        this.x = x;
        this.y = y;
        this.radius = 1.0f;
        this.grahamAngle = Float.MAX_VALUE;
        edgeStart = null;
        edgeEnd = null;
        path = Path.Left;
    }
  
    public Point(float x, float y, float radius)
    {
        this.x = x;
        this.y = y;
        this.radius = radius;
        this.grahamAngle = Float.MAX_VALUE;
        edgeStart = null;
        edgeEnd = null;
        path = Path.Left;
    }
    
    public void move(float x, float y)
    {
        this.x = x;
        this.y = y;
    }
    
    public boolean sharesEdge(Point other)
    {
        if (edgeStart == other.edgeStart || edgeStart == other.edgeEnd || edgeEnd == other.edgeStart || edgeEnd == other.edgeEnd)
        {
            return true;
        }
        
        return false;
    }
    
    public boolean sharesPath(Point other)
    {
        return path == other.path;
    }
    
    public boolean onScreen()
    {
        return x >= 0.0f && x <= windowWidth && y >= 0.0f && y <= windowHeight;
    }
    
    @Override
    public int hashCode()
    {
        return Objects.hash(x + y);
    }

    @Override
    public boolean equals(final Object other)
    {
        if (!(other instanceof Point))
        {
            return false;
        }
        
        float differenceX = Math.abs(x - ((Point)other).x);
        float differenceY = Math.abs(y - ((Point)other).y);
        return differenceX < 0.001 && differenceY < 0.001;
    }
    
    float x;
    float y;
    float radius;
    float grahamAngle;
    Edge edgeStart;
    Edge edgeEnd;
    Path path;
};

public class PositionComparator implements Comparator<Point>
{
    @Override
    public int compare(Point first, Point second)
    {
        if (first.y < second.y || first.y == second.y && first.x > second.x)
        {
            return -1;
        }
        return 1;
    }
};

public class AngleComparator implements Comparator<Point>
{
    @Override
    public int compare(Point first, Point second)
    {
        if (first.grahamAngle < second.grahamAngle)
        {
            return -1;
        }
        return 1;
    }
};

public class XComparator implements Comparator<Point>
{
    @Override
    public int compare(Point first, Point second)
    {
        if (first.x < second.x)
        {
            return -1;
        }
        return 1;
    }
};

public class YComparator implements Comparator<Point>
{
    @Override
    public int compare(Point first, Point second)
    {
        if (first.y < second.y)
        {
            return -1;
        }
        return 1;
    }
};

public class Edge
{
    public Edge(Point start, Point end)
    {
        this.start = start;
        this.end = end;
        this.circleCenter = null;
        onConvexHull = true;
    }
    
    public boolean contains(Point point)
    {
        return point == start || point == end;
    }
    
    public void swapOrientation()
    {
        Point swap = start;
        start = end;
        end = swap;
    }
    
    public boolean sharesPoints(Edge other)
    {
        return start.equals(other.start) && end.equals(other.end) || start.equals(other.end) && end.equals(other.start);
    }
    
    public Point getNearestPoint(Point other)
    {
        double apx = other.x - start.x;
        double apy = other.y - start.y;
        double abx = end.x - start.x;
        double aby = end.y - start.y;
    
        double ab2 = abx * abx + aby * aby;
        double apab = apx * abx + apy * aby;
        double t = apab / ab2;

        Point result = new Point((float)(start.x + abx * t), (float)(start.y + aby * t));
        return result;
    }
    
    Point getIntersection(Edge other)
    {
        double a1 = end.y - start.y;
        double b1 = start.x - end.x;
        double c1 = a1 * start.x + b1 * start.y;
      
        double a2 = other.end.y - other.start.y;
        double b2 = other.start.x - other.end.x;
        double c2 = a2 * other.start.x + b2 * other.start.y;
      
        double determinant = a1 * b2 - a2 * b1;
      
        if (determinant == 0)
        {
            return null;
        }
        else
        {
            double x = (b2 * c1 - b1 * c2) / determinant;
            double y = (a1 * c2 - a2 * c1) / determinant;
            return new Point((float)x, (float)y);
        }
    }
    
    Edge getClosestIntersectingEdge(ArrayList<Edge> edges)
    {
        Edge result = null;
        float shortestDistance = Float.MAX_VALUE;
        Line2D current = new Line2D.Float(start.x, start.y, end.x, end.y);
        
        for (Edge other : edges)
        {
            Line2D otherLine = new Line2D.Float(other.start.x, other.start.y, other.end.x, other.end.y);
            if (current.intersectsLine(otherLine))
            {
                Point intersection = getIntersection(other);
                float distance = getLength(getVector(start, intersection));
                if (distance < shortestDistance)
                {
                    result = other;
                    shortestDistance = distance;
                }
            }
        }
        
        return result;
    }
    
    @Override
    public int hashCode()
    {
        return Objects.hash(start.hashCode() + end.hashCode());
    }

    @Override
    public boolean equals(final Object other)
    {
        if (!(other instanceof Edge))
        {
            return false;
        }
        
        return start.equals(((Edge)other).start) && end.equals(((Edge)other).end);
    }
    
    Point start;
    Point end;
    Point circleCenter;
    boolean onConvexHull;
};

public class Circle
{
    public Circle()
    {
       center = new Point(0.0f, 0.0f);
       radius = 1.0f;
    }
    
    public Circle(Point center, float radius)
    {
       this.center = center;
       this.radius = radius;
    }
    
    public boolean inside(Point point)
    {
        return getDistance(point, center) < radius;
    }
    
    Point center;
    float radius;
};

public class Triangle
{
    public Triangle()
    {
        this.first = null;
        this.second = null;
        this.third = null;
        this.circleCenter = null;
    }
  
    public Triangle(Edge first, Edge second, Edge third)
    {
        this.first = first;
        this.second = second;
        this.third = third;
        this.circleCenter = null;
    }
    
    public Triangle(Edge first, Edge second, Edge third, Point circleCenter)
    {
        this.first = first;
        this.second = second;
        this.third = third;
        this.circleCenter = circleCenter;
    }
    
    public Edge getSharedEdge(Triangle other)
    {
        if (first.sharesPoints(other.first) || first.sharesPoints(other.second) || first.sharesPoints(other.third))
        {
            return first;
        }
        else if (second.sharesPoints(other.first) || second.sharesPoints(other.second) || second.sharesPoints(other.third))
        {
            return second;
        }
        else if (third.sharesPoints(other.first) || third.sharesPoints(other.second) || third.sharesPoints(other.third))
        {
            return third;
        }
        return null;
    }
    
    boolean pointInside(Point point)
    {
        Polygon area = new Polygon();
        area.addPoint((int)first.start.x, (int)first.start.y);
        area.addPoint((int)first.end.x, (int)first.end.y);
        
        if (second.start.equals(first.start) || second.start.equals(first.end))
        {
            area.addPoint((int)second.end.x, (int)second.end.y);
        }
        else
        {
            area.addPoint((int)second.start.x, (int)second.start.y);
        }
        
        if (area.contains(point.x, point.y))
        {
            return true;
        }
        
        return false;
    }
    
    @Override
    public int hashCode()
    {
        return Objects.hash(first.hashCode() + second.hashCode() + third.hashCode() + circleCenter.hashCode());
    }

    @Override
    public boolean equals(final Object other)
    {
        if (!(other instanceof Triangle))
        {
            return false;
        }
        
        return first.equals(((Triangle)other).first) && second.equals(((Triangle)other).second) && third.equals(((Triangle)other).third) && circleCenter.equals(((Triangle)other).circleCenter);
    }
    
    Edge first;
    Edge second;
    Edge third;
    Point circleCenter;
};

public class KdNode
{
    public KdNode()
    {
        depth = 0;
        id = null;
        parent = null;
        left = null;
        right = null;
    }
    
    public KdNode(int depth, Point id, KdNode parent, KdNode left, KdNode right)
    {
        this.depth = depth;
        this.id = id;
        this.parent = parent;
        this.left = left;
        this.right = right;
    }
    
    int depth;
    Point id;
    KdNode parent;
    KdNode left;
    KdNode right;
};

public class KdTree
{
    public KdTree()
    {
        horizontalEdges = null;
        verticalEdges = null;
        root = null;
    }
    
    public KdTree(ArrayList<Point> points)
    {
        horizontalEdges = new ArrayList<Edge>();
        verticalEdges = new ArrayList<Edge>();
        root = buildTreeRecursive(points, 0);
    }
    
    private KdNode buildTreeRecursive(ArrayList<Point> points, int depth)
    {
        if (points.size() < 1)
        {
            return null;
        }
        
        if (points.size() == 1)
        {
            KdNode newNode = new KdNode(depth, points.get(0), null, null, null);
            return newNode;
        }
        
        int medianIndex;
        if (points.size() % 2 == 0)
        {
            medianIndex = points.size() / 2 - 1;
        }
        else
        {
            medianIndex = points.size() / 2;
        }
        
        ArrayList<Point> firstPart = new ArrayList<Point>();
        ArrayList<Point> secondPart = new ArrayList<Point>();
        
        if (depth % 2 == 0)
        {
            Collections.sort(points, new XComparator());
        }
        else
        {
            Collections.sort(points, new YComparator());
        }
        
        for (int i = 0; i <= medianIndex; i++)
        {
            firstPart.add(points.get(i));
        }
        for (int i = medianIndex + 1; i < points.size(); i++)
        {
            secondPart.add(points.get(i));
        }
        
        Point id = points.get(medianIndex);
        if (depth % 2 == 0)
        {
            float startX = id.x;
            float endX = id.x;
            Point intersection = findIntersectionVertical(id);
            float startY = intersection.x;
            float endY = intersection.y;
            verticalEdges.add(new Edge(new Point(startX, startY), new Point(endX, endY)));
        }
        else
        {
            float startY = id.y;
            float endY = id.y;
            Point intersection = findIntersectionHorizontal(id);
            float startX = intersection.x;
            float endX = intersection.y;
            horizontalEdges.add(new Edge(new Point(startX, startY), new Point(endX, endY)));
        }
        
        KdNode left = buildTreeRecursive(firstPart, depth + 1);
        KdNode right = buildTreeRecursive(secondPart, depth + 1);
        
        KdNode newNode = new KdNode(depth, id, null, left, right);
        
        if (left != null)
        {
            left.parent = newNode;
        }
        if (right != null)
        {
            right.parent = newNode;
        }
        return newNode;
    }
    
    public void drawEdges()
    {
        stroke(255, 255, 0);
        strokeWeight(3);
        for (int i = 0; i < horizontalEdges.size(); i++)
        {
            Edge current = horizontalEdges.get(i);
            line(current.start.x, current.start.y, current.end.x, current.end.y);
        }
        
        stroke(0, 255, 255);
        for (int i = 0; i < verticalEdges.size(); i++)
        {
            Edge current = verticalEdges.get(i);
            line(current.start.x, current.start.y, current.end.x, current.end.y);
        }
        
        stroke(pointColor);
        strokeWeight(1);
    }
    
    private Point findIntersectionHorizontal(Point originPoint)
    {
        float start = 0.0f;
        float end = windowWidth;
        
        for (int i = 0; i < verticalEdges.size(); i++)
        {
            Edge currentEdge = verticalEdges.get(i);
            if (currentEdge.start.y <= originPoint.y && currentEdge.end.y >= originPoint.y && currentEdge.start.x < originPoint.x && currentEdge.start.x > start)
            {
                start = currentEdge.start.x;
            }
            
            if (currentEdge.start.y <= originPoint.y && currentEdge.end.y >= originPoint.y && currentEdge.start.x >= originPoint.x && currentEdge.start.x < end)
            {
                end = currentEdge.start.x;
            }
        }
        
        return new Point(start, end);
    }
    
    private Point findIntersectionVertical(Point originPoint)
    {
        float start = 0.0f;
        float end = windowHeight;
        
        for (int i = 0; i < horizontalEdges.size(); i++)
        {
            Edge currentEdge = horizontalEdges.get(i);
            if (currentEdge.start.x <= originPoint.x && currentEdge.end.x >= originPoint.x && currentEdge.start.y < originPoint.y && currentEdge.start.y > start)
            {
                start = currentEdge.start.y;
            }
            
            if (currentEdge.start.x <= originPoint.x && currentEdge.end.x >= originPoint.x && currentEdge.start.y >= originPoint.y && currentEdge.start.y < end)
            {
                end = currentEdge.start.y;
            }
        }
        
        return new Point(start, end);
    }
    
    KdNode root;
    ArrayList<Edge> horizontalEdges;
    ArrayList<Edge> verticalEdges;
};

// Basic variables
int windowWidth; // window resolution in size() method needs to be modified separately
int windowHeight;
int backgroundColor;
int pointColor;

// Point variables
boolean showPoints;
boolean pointSelected;
boolean pointDragged;
float pointRadius;
Point selectedPoint;
ArrayList<Point> points;

// Hull variables
boolean showHull;
PVector hullColor;
ArrayList<Point> hullPoints;

// Triangulation variables
boolean showTriangulation;
ArrayList<Edge> triangulationEdges;

// Delaunay triangulation variables
boolean showDelaunay;
ArrayList<Edge> delaunayEdges;

// Voronoi diagram variables
boolean showVoronoi;
ArrayList<Triangle> delaunayTriangles;
ArrayList<Edge> voronoiEdges;
ArrayList<Point> voronoiPoints;

// k-D tree variables
boolean showTree;
KdTree kdTree;

void setup()
{
    windowWidth = 1280;
    windowHeight = 720;
    backgroundColor = 255;
    pointColor = 0;
    
    showPoints = true;
    pointSelected = false;
    pointDragged = false;
    pointRadius = 7.0f;
    selectedPoint = null;
    points = new ArrayList<Point>();
    
    showHull = false;
    hullColor = new PVector(0, 0, 255);
    hullPoints = new ArrayList<Point>();
    
    showTriangulation = false;
    triangulationEdges = new ArrayList<Edge>();
    
    showDelaunay = false;
    delaunayEdges = new ArrayList<Edge>();
    
    showVoronoi = false;
    delaunayTriangles = new ArrayList<Triangle>();
    voronoiEdges = new ArrayList<Edge>();
    voronoiPoints = new ArrayList<Point>();
    
    showTree = false;
    kdTree = new KdTree();
    
    size(1280, 720);
    background(backgroundColor);
    fill(pointColor);
    textSize(14);
}

void draw()
{
    background(backgroundColor);
    
    if (showPoints)
    {
        for (Point point : points)
        {
            ellipse(point.x, point.y, point.radius * 2.0f, point.radius * 2.0f);
        }
        
        if (!pointDragged)
        {
            pointSelected = checkPointSelection();
        }
        
        if (pointSelected && selectedPoint != null)
        {
            fill(128);
            ellipse(selectedPoint.x, selectedPoint.y, selectedPoint.radius * 2.0f, selectedPoint.radius * 2.0f);
            fill(pointColor);
        }
    }
    
    if (showTriangulation)
    {
        drawEdges(triangulationEdges, new PVector(255, 0, 0));
        
        if (!showHull)
        {
            drawLines(points, new PVector(0, 0, 255));
        }
    }
    
    if (showDelaunay)
    {
        drawEdges(delaunayEdges, new PVector(255, 0, 0));
        
        if (showVoronoi)
        {
            drawVoronoi(new PVector(0, 150, 0));
        }
    }
    
    if (showHull)
    {
        drawLines(hullPoints, hullColor);
    }
    
    if (showTree)
    {
        kdTree.drawEdges();
    }
}

void mousePressed()
{
    if (pointSelected && selectedPoint != null)
    {
        pointDragged = true;
        clearStructures();
    }
}

void mouseDragged()
{
    if (pointDragged)
    {
        selectedPoint.move(min(max(pointRadius, mouseX), windowWidth - pointRadius), min(max(pointRadius, mouseY), windowHeight - pointRadius));
    }
}

void mouseReleased()
{
    pointDragged = false;
}

void keyPressed()
{
    // Add point
    if (key == 'a')
    {
        points.add(new Point(mouseX, mouseY, pointRadius));
        clearStructures();
    }
    
    // Delete point
    if (key == 'd' && pointSelected && selectedPoint != null)
    {
        points.remove(selectedPoint);
        selectedPoint = null;
        clearStructures();
    }
    
    // Reload scene
    if (key == 'r')
    {
        points.clear();
        clearStructures();
    }
    
    // Generate 5 points
    if (key == 'g')
    {
        generateRandomPoints(5);
        clearStructures();
    }
    
    // Toggle points
    if (key == 'p')
    {
        showPoints = !showPoints;
    }
    
    // Compute convex hull (Gift wrapping)
    if (key == 'c' && points.size() > 2)
    {
        clearStructures();
        hullSimple();
        hullColor = new PVector(0, 0, 255);
        showHull = true;
    }
    
    // Compute convex hull (Graham Scan)
    if (key == 'v' && points.size() > 2)
    {
        clearStructures();
        hullGraham();
        hullColor = new PVector(0, 255, 0);
        showHull = true;
    }
    
    // Perform triangulation
    if (key == 't' && points.size() > 2)
    {
        if (showHull)
        {
            triangulate(hullPoints);
        }
        else
        {
            triangulate(points);
        }
        showTriangulation = true;
    }
    
    // Perform Delaunay triangulation
    if (key == 'y' && points.size() > 2)
    {
        if (showHull)
        {
            triangulateDelaunay(hullPoints);
        }
        else
        {
            triangulateDelaunay(points);
        }
        showDelaunay = true;
    }
    
    // Draw Voronoi diagram
    if (key == 'u' && showDelaunay)
    {
        calculateVoronoi(delaunayTriangles);
        showVoronoi = true;
    }
    
    // Compute k-D tree
    if (key == 'k' && points.size() > 0)
    {
        clearStructures();
        buildTree(points);
        showTree = true;
    }
}

boolean floatEquals(float a, float b)
{
    float difference = Math.abs(a - b);
    return difference <= Math.ulp(1.0f);
}

void generateRandomPoints(int count)
{
    for (int i = 0; i < count; i++)
    {
        points.add(new Point(random(pointRadius, (float)windowWidth - pointRadius), random(pointRadius, (float)windowHeight - pointRadius), pointRadius));
    }
}

boolean checkPointSelection()
{
    for (Point point : points)
    {
        if (getDistance(point, new Point(mouseX, mouseY)) < point.radius)
        {
            selectedPoint = point;
            return true;
        }
    }
    return false;
}

void clearStructures()
{
    hullPoints.clear();
    showHull = false;
    
    triangulationEdges.clear();
    showTriangulation = false;
    
    delaunayEdges.clear();
    showDelaunay = false;
    
    delaunayTriangles.clear();
    voronoiEdges.clear();
    voronoiPoints.clear();
    showVoronoi = false;
    
    kdTree = new KdTree();
    showTree = false;
}

ArrayList<Edge> pointsToEdges(ArrayList<Point> points)
{
    ArrayList<Edge> edges = new ArrayList<Edge>();
    
    for (int i = 0; i < points.size(); i++)
    {
        Point start = points.get(i);
        Point end = points.get(0);
        
        if (i + 1 != points.size())
        {
            end = points.get(i + 1);
        }
        
        Edge edge = new Edge(start, end);
        start.edgeStart = edge;
        end.edgeEnd = edge;
        edges.add(edge);
    }
    
    return edges;
}

void drawLines(ArrayList<Point> points, PVector lineColor)
{
    ArrayList<Edge> edges = pointsToEdges(points);
    drawEdges(edges, lineColor);
}

void drawEdges(ArrayList<Edge> edges, PVector lineColor)
{
    stroke(lineColor.x, lineColor.y, lineColor.z);
    
    for (Edge edge : edges)
    {
        line(edge.start.x, edge.start.y, edge.end.x, edge.end.y);
    }
    
    stroke(pointColor);
}

float getDistance(Point first, Point second)
{
    float dx = second.x - first.x;
    float dy = second.y - first.y;
    
    return sqrt(dx * dx + dy * dy);
}

float getLength(PVector vector)
{
    return sqrt(vector.x * vector.x + vector.y * vector.y);
}

PVector getVector(Point first, Point second)
{
    return new PVector(second.x - first.x, second.y - first.y);
}

float getCosine(PVector first, PVector second)
{
    float dot = first.dot(second);
    float firstLength = getLength(first);
    float secondLength = getLength(second);
    
    float cosine = dot / (firstLength * secondLength);
    return cosine;
}

float getAngle(PVector first, PVector second)
{
   return acos(getCosine(first, second)); 
}

/* Possible results:
 * = 0 - points lie on one line
 * > 0 - points are oriented clockwise
 * < 0 - points are oriented counter-clockwise
 */
float getOrientation(Point first, Point second, Point third)
{
    PVector firstToSecond = getVector(first, second);
    PVector firstToThird = getVector(first, third);
    
    return firstToSecond.cross(firstToThird).z;
}

void hullSimple()
{
    Point initial = points.get(0);
    
    for (Point point : points)
    {
        if (point.x > initial.x)
        {
            initial = point;
        }
    }
    
    Point previous = new Point(initial.x, initial.y - 1.0f);
    Point current = initial;
    
    do
    {
        hullPoints.add(current);
        
        PVector v1 = getVector(previous, current);
        Point next = null;
        float cosineBest = 0.0f;
        
        for (Point point : points)
        {
            if (point == current || point == previous)
            {
                continue;
            }
            
            PVector v2 = getVector(current, point);
            float cosineCurrent = getCosine(v1, v2);
            
            if (next == null || cosineCurrent > cosineBest)
            {
                next = point;
                cosineBest = cosineCurrent;
            }
        }
        
        previous = current;
        current = next;
    }
    while (current != initial);
}

void hullGraham()
{
    Point initial = points.get(0);
    
    for (Point point : points)
    {
        if (point.y < initial.y)
        {
            initial = point;
        }
    }
    
    PVector xAxis = new PVector(1.0f, 0.0f);
    initial.grahamAngle = 0.0f;
    
    for (Point point : points)
    {
        if (point == initial)
        {
            continue;
        }
        
        PVector pointConnection = getVector(initial, point);
        point.grahamAngle = getAngle(xAxis, pointConnection);
    }
    
    Collections.sort(points, new AngleComparator());
    
    Stack<Point> stack = new Stack<Point>();
    stack.push(points.get(0));
    stack.push(points.get(1));
    stack.push(points.get(2));
    
    for (int i = 3; i <= points.size(); i++)
    {
        Point third = stack.pop();
        Point second = stack.pop();
        Point first = stack.pop();
        
        while (getOrientation(first, second, third) <= 0.0f && stack.size() > 0) // we keep removing middle point until orientation is correct or stack becomes empty
        {
            second = first;
            first = stack.pop();
        }
        
        stack.push(first);
        if (getOrientation(first, second, third) > 0.0f)
        {
            stack.push(second);
        }
        stack.push(third);
        
        if (i < points.size())
        {
            stack.push(points.get(i));
        }
    }
    
    hullPoints.addAll(stack);
}

void triangulate(ArrayList<Point> polygon)
{
    ArrayList<Point> points = new ArrayList<Point>();
    points.addAll(polygon);
    pointsToEdges(points); // initialize edge variables inside points
    
    Collections.sort(points, new PositionComparator());
    
    // Get 2 points which are connected to the top point
    Point nextOne = points.get(0).edgeStart.end;
    Point nextTwo = points.get(0).edgeEnd.start;
    
    // We only modify right path, because all points have path set to left by default
    if (nextOne.x < nextTwo.x)
    {
        while (nextTwo != points.get(points.size() - 1))
        {
            nextTwo.path = Path.Right;
            nextTwo = nextTwo.edgeEnd.start;
        }
    }
    else
    {
        while (nextOne != points.get(points.size() - 1))
        {
            nextOne.path = Path.Right;
            nextOne = nextOne.edgeStart.end;
        }
    }
    
    Stack<Point> stack = new Stack<Point>();
    stack.push(points.get(0));
    stack.push(points.get(1));
    
    for (int i = 2; i < points.size(); i++)
    {
        if (stack.peek().sharesPath(points.get(i)))
        {
            Point top = stack.pop();
            
            if (top.path == Path.Left)
            {
                while (stack.size() > 0 && getOrientation(points.get(i), top, stack.peek()) > 0.0f)
                {
                    triangulationEdges.add(new Edge(points.get(i), stack.peek()));
                    top = stack.pop();
                }
            }
            else
            {
                while (stack.size() > 0 && getOrientation(points.get(i), top, stack.peek()) <= 0.0f)
                {
                    triangulationEdges.add(new Edge(points.get(i), stack.peek()));
                    top = stack.pop();
                }
            }
            stack.push(top);
        }
        else
        {
            Point top = stack.pop();
            triangulationEdges.add(new Edge(points.get(i), top));
            
            while (!stack.empty())
            {
                triangulationEdges.add(new Edge(points.get(i), stack.pop()));
            }
            
            stack.push(top);
        }
        stack.push(points.get(i));
    }
}

float crossProduct(Point first, Point second, Point third)
{
    float u1 = second.x - first.x;
    float v1 = second.y - first.y;
    float u2 = third.x - first.x;
    float v2 = third.y - first.y;
    return u1 * v2 - v1 * u2;
}

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
            Point intersection = triangle.first.getNearestPoint(triangle.circleCenter);
            PVector vector = new PVector(intersection.x - triangle.circleCenter.x, intersection.y - triangle.circleCenter.y);
            vector.mult(100.0f);
            
            Edge test = new Edge(new Point(triangle.circleCenter.x, triangle.circleCenter.y), new Point(vector.x, vector.y));
            if (!triangle.pointInside(triangle.circleCenter))
            {
                Edge closestIntersecting = test.getClosestIntersectingEdge(delaunayEdges);
                if (closestIntersecting != null && closestIntersecting.sharesPoints(triangle.first))
                {
                    vector.mult(-1.0f);
                }
            }
            
            voronoiEdges.add(new Edge(new Point(triangle.circleCenter.x, triangle.circleCenter.y), new Point(vector.x, vector.y)));
        }
        
        if (triangle.second.onConvexHull)
        {
            Point intersection = triangle.second.getNearestPoint(triangle.circleCenter);
            PVector vector = new PVector(intersection.x - triangle.circleCenter.x, intersection.y - triangle.circleCenter.y);
            vector.mult(100.0f);
            
            Edge test = new Edge(new Point(triangle.circleCenter.x, triangle.circleCenter.y), new Point(vector.x, vector.y));
            if (!triangle.pointInside(triangle.circleCenter))
            {
                Edge closestIntersecting = test.getClosestIntersectingEdge(delaunayEdges);
                if (closestIntersecting != null && closestIntersecting.sharesPoints(triangle.second))
                {
                    vector.mult(-1.0f);
                }
            }
            
            voronoiEdges.add(new Edge(new Point(triangle.circleCenter.x, triangle.circleCenter.y), new Point(vector.x, vector.y)));
        }
        
        if (triangle.third.onConvexHull)
        {
            Point intersection = triangle.third.getNearestPoint(triangle.circleCenter);
            PVector vector = new PVector(intersection.x - triangle.circleCenter.x, intersection.y - triangle.circleCenter.y);
            vector.mult(100.0f);
            
            Edge test = new Edge(new Point(triangle.circleCenter.x, triangle.circleCenter.y), new Point(vector.x, vector.y));
            if (!triangle.pointInside(triangle.circleCenter))
            {
                Edge closestIntersecting = test.getClosestIntersectingEdge(delaunayEdges);
                if (closestIntersecting != null && closestIntersecting.sharesPoints(triangle.third))
                {
                    vector.mult(-1.0f);
                }
            }
            
            voronoiEdges.add(new Edge(new Point(triangle.circleCenter.x, triangle.circleCenter.y), new Point(vector.x, vector.y)));
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

void buildTree(ArrayList<Point> points)
{
    ArrayList<Point> copy = new ArrayList<Point>();
    copy.addAll(points);
    kdTree = new KdTree(copy);
}