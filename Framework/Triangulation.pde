// Author: Filip Petrovic (422334)

// Triangulation variables
boolean showTriangulation;
ArrayList<Edge> triangulationEdges;

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