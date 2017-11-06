// Author: Filip Petrovic (422334)

import java.util.Collections;
import java.util.Comparator;
import java.util.Objects;
import java.util.Stack;

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
    }
    
    public boolean contains(Point point)
    {
        return point == start || point == end;
    }
    
    Point start;
    Point end;
};

enum Path
{
    Left,
    Right
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
        root = null;
    }
    
    public KdTree(ArrayList<Point> points)
    {
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
        
        for (int i = 0; i < medianIndex; i++)
        {
            firstPart.add(points.get(i));
        }
        for (int i = medianIndex + 1; i < points.size(); i++)
        {
            secondPart.add(points.get(i));
        }
        
        KdNode left = buildTreeRecursive(firstPart, depth + 1);
        KdNode right = buildTreeRecursive(secondPart, depth + 1);
        
        KdNode newNode = new KdNode(depth, points.get(medianIndex), null, left, right);
        
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
    
    KdNode root;
};

// Basic variables
int windowWidth; // window resolution in size() method needs to be modified separately
int windowHeight;
int backgroundColor;
int pointColor;
boolean showHint;

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

// k-D tree variables
boolean showTree;
KdTree kdTree;

void setup()
{
    windowWidth = 1280;
    windowHeight = 720;
    backgroundColor = 255;
    pointColor = 0;
    showHint = true;
    
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
    
    if (showHint)
    {
        text("Controls:", 10, 20); 
        text("A - Add new point under mouse cursor", 10, 40);
        text("D - Delete existing point under mouse cursor", 10, 60);
        text("Mouse click - Drag existing point to another place", 10, 80);
        text("G - Generate 5 random points", 10, 100);
        text("C - Convex hull (Gift wrapping)", 10, 120);
        text("V - Convex hull (Graham Scan)", 10, 140);
        text("T - Triangulation of a polygon defined by all points in order of their addition, the last added point will be automatically connected to the first, the polygon has to be y-monotone", 10, 160);
        text("T (after C/V) - Triangulation of convex hull", 10, 180);
        text("K - k-D tree", 10, 200);
        text("R - Reload", 10, 220);
        text("H - Toggle control hint", 10, 240);
        text("P - Toggle points", 10, 260);
    }
    
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
    
    if (showHull)
    {
        drawLines(hullPoints, hullColor); //<>//
    }
    
    if (showTree)
    {
        drawTree(kdTree);
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
    
    // Toggle hint
    if (key == 'h')
    {
        showHint = !showHint;
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
    
    // Compute k-D tree
    if (key == 'k' && points.size() > 0)
    {
        buildTree(points);
        showTree = true; //<>//
    }
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
    Point current = initial; //<>//
    
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

void hullGraham() //<>//
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

void buildTree(ArrayList<Point> points)
{
    ArrayList<Point> copy = new ArrayList<Point>();
    copy.addAll(points);
    kdTree = new KdTree(copy);
}

void drawLineHorizontal(Point point)
{
    stroke(255, 255, 0);
    line(0, point.y, windowWidth, point.y);
    stroke(pointColor);
}

void drawLineVertical(Point point)
{
    stroke(0, 255, 255);
    line(point.x, 0, point.x, windowHeight);
    stroke(pointColor);
}

void drawTree(KdTree kdTree)
{
    // to do
}