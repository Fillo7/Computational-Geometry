// Author: Filip Petrovic (422334)

import java.awt.geom.Line2D;
import java.awt.Polygon;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.Objects;
import java.util.Stack;
import java.util.Queue;

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

float crossProduct(Point first, Point second, Point third)
{
    float u1 = second.x - first.x;
    float v1 = second.y - first.y;
    float u2 = third.x - first.x;
    float v2 = third.y - first.y;
    return u1 * v2 - v1 * u2;
}