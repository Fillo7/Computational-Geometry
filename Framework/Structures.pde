// Author: Filip Petrovic (422334)

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