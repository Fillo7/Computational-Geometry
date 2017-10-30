// Author: Filip Petrovic (422334)
// Version: 30.10.2017 (Bug fixes in triangulation function)

import java.util.Collections;
import java.util.Comparator;
import java.util.Objects;
import java.util.Stack;

public class Point
{
    Point(float x, float y)
    {
        this.x = x;
        this.y = y;
        this.radius = 1.0f;
        this.grahamAngle = Float.MAX_VALUE;
    }
  
    Point(float x, float y, float radius)
    {
        this.x = x;
        this.y = y;
        this.radius = radius;
        this.grahamAngle = Float.MAX_VALUE;
    }
    
    public void move(float x, float y)
    {
        this.x = x;
        this.y = y;
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

public class Edge
{
    Edge(Point start, Point end)
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
PVector triangulationColor;
ArrayList<Edge> triangulationEdges;

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
    triangulationColor = new PVector(255, 0, 0);
    triangulationEdges = new ArrayList<Edge>();
    
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
        text("T - Triangulation of a body defined by all points, the last added point will be automatically connected to the first", 10, 160);
        text("T (after C/V) - Triangulation of convex hull", 10, 180);
        text("R - Reload", 10, 200);
        text("H - Toggle control hint", 10, 220);
        text("P - Toggle points", 10, 240);
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
     
    if (showHull)
    {
        drawLines(hullPoints, hullColor); //<>//
    }
    
    if (showTriangulation)
    {
        drawEdges(triangulationEdges, triangulationColor);
        
        if (!showHull)
        {
            drawLines(points, hullColor);
        }
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
        
        edges.add(new Edge(start, end));
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

boolean sharesEdge(Point first, Point second, ArrayList<Edge> edges)
{
    for (Edge edge : edges)
    {
        if (edge.contains(first) && edge.contains(second))
        {
            return true;
        }
    }
    
    return false;
}

void triangulate(ArrayList<Point> body)
{
    ArrayList<Point> points = new ArrayList<Point>();
    points.addAll(body);
    ArrayList<Edge> edges = pointsToEdges(points);
    
    Collections.sort(points, new PositionComparator());
    
    Stack<Point> stack = new Stack<Point>();
    stack.push(points.get(0));
    stack.push(points.get(1));
    
    for (int i = 2; i < points.size(); i++)
    {
        if (sharesEdge(stack.peek(), points.get(i), edges))
        {
            Point top = stack.pop();
            
            while (stack.size() > 0 && getOrientation(points.get(i), top, stack.peek()) > 0.0f)
            {
                triangulationEdges.add(new Edge(points.get(i), stack.peek()));
                top = stack.pop();
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