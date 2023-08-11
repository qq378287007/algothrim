#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

class Point
{
public:
    double x, y;
};

class LineSegment
{
public:
    Point p1, p2;
};

// p点为起点，向x轴负向做射线，求与线段l的交点
//-1，无交点
// 0，p点位于线段上
// 1, 交点位于线段中部
// 2，交点位于线段端点
// 3，射线与线段重叠
int insertType(const Point &p, const LineSegment &l)
{
    Point p1 = l.p1;
    Point p2 = l.p2;
    if (p1.y > p2.y)
    {
        p1 = l.p2;
        p2 = l.p1;
    }

    // p点位于l矩形框的下侧或上侧
    if (p.y < p1.y || p.y > p2.y)
        return -1;

    double x1 = p2.x - p1.x;
    double y1 = p2.y - p1.y;

    double x2 = p.x - p1.x;
    double y2 = p.y - p1.y;

    double cross = x1 * y2 - x2 * y1;
    if (cross > 0.0)
    {
        return -1;
    }
    else if (cross == 0.0)
    {
        if (p1.y != p2.y)
            return 0;

        if (p.x < p1.x && p.x < p2.x)
            return -1;
        else if (p.x > p1.x && p.x > p2.x)
            return 3;
        else
            return 0;
    }
    else
    {
        if (p1.y == p2.y)
            return -1;
        else if (p.y == p1.y || p.y == p2.y)
            return 2;
        else
            return 1;
    }
}

// 判断点p是否在item区域内
// 0，位于边上
// 1, 位于内部
//-1, 位于外部
int inItem(const Point &p, const vector<LineSegment> &item)
{
    int type;
    int n1 = 0;
    int n2 = 0;
    for (auto iter = item.cbegin(); iter != item.cend(); iter++)
    {
        type = insertType(p, *iter);
        switch (type)
        {
        case 0:
            return 0;
        case 1:
        case 3:
            n1++;
            break;
        case 2:
            n2++;
            break;
        }
    }

    n1 += n2 / 2;
    if (n1 % 2 == 1)
        return 1;

    return -1;
}

int main()
{
    Point p;
    p.x = 2.5;
    p.y = 2.0;

    vector<LineSegment> item{
        {{-1.0, -1.0}, {1.0, -1.0}},
        {{1.0, -1.0}, {1.0, 1.0}},
        {{1.0, 1.0}, {-1.0, 1.0}},
        {{-1.0, 1.0}, {-1.0, -1.0}},
        {{-2.0, -2.0}, {2.0, -2.0}},
        {{2.0, -2.0}, {2.0, 2.0}},
        {{2.0, 2.0}, {-2.0, 2.0}},
        {{-2.0, 2.0}, {-2.0, -2.0}},
    };

    int flag = inItem(p, item);
    cout << "flag: " << flag << endl;

    return 0;
}