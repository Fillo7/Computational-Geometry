// Author: Filip Petrovic (422334)

// k-D tree variables
boolean showTree;
KdTree kdTree;

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

void buildTree(ArrayList<Point> points)
{
    ArrayList<Point> copy = new ArrayList<Point>();
    copy.addAll(points);
    kdTree = new KdTree(copy);
}