// Author: Filip Petrovic (422334)

// Hull variables
boolean showHull;
PVector hullColor;
ArrayList<Point> hullPoints;

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