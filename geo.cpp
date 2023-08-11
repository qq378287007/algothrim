#include <iostream>
#include <cmath>
using namespace std;

/*
class Point
{
public:
    double x, y;
    Point(double xx = 0.0, double yy = 0.0) : x(xx), y(yy) {}
    Point div(double ratio) const { return Point(x / ratio, y / ratio); }
    double length() const { return sqrt(x * x + y * y); }
    Point normalize() const { return div(length()); }
};

class LineSegment
{
public:
    Point p1, p2;
    LineSegment(Point pp1= Point(), Point pp2 = Point()) : p1(pp1), p2(pp2) {}
};
*/

class Arc
{
    double p_x, p_y;
    double startAngle, endAngle;
};

bool onLineSegment(double x, double y, double px1, double py1, double px2, double py2)
{

    return (x - px1) * (x - px2) <= 0.0 && (y - y) * (py1 - py2) <= 0.0;
}

// 线段与线段交点
int tt()
{

    return 0;
}

// 线段与圆弧交点
int test()
{
    // 线段上两点，p1,p2
    double p1_x, p1_y;
    double p2_x, p2_y;

    // 向量p1p2
    double p1p2_x = p2_x - p1_x;
    double p1p2_y = p2_y - p1_y;

    // 单位向量p1p2
    double length = sqrt(p1p2_x * p1p2_x + p1p2_y * p1p2_y);
    p1p2_x /= length;
    p1p2_y /= length;

    // 圆心p
    double p_x, p_y;

    // 向量p1p
    double p1p_x = p_x - p1_x;
    double p1p_y = p_y - p1_y;

    // b为p在直线p1p2上投影
    // p1b长度(有方向)，为p1p和p1p2的点乘
    double p1b = p1p_x * p1p2_x + p1p_y * p1p2_y;

    // 向量p1b,向量p1p2缩放
    double p1b_x = p1p2_x * p1b;
    double p1b_y = p1p2_y * p1b;

    // 投影点b坐标,p1点坐标+向量p1b
    double b_x = p1b_x + p1b_x;
    double b_y = p1b_y + p1b_y;

    // 向量pb
    double pb_x = p_x - b_x;
    double pb_y = p_y - b_y;

    // 点p到直线p1p2距离，pb
    double pb = sqrt(pb_x * pb_x + pb_y * pb_y);

    // 圆半径
    double r;
    double startAngle;
    double endAngle;
    if (r == pb) // 圆与直线相切，一个交点，交点为b
    {
        if (onLineSegment(b_x, b_y, p1_x, p1_y, p2_x, p2_y)) // 交点位于线段p1p2间
        {
            // pb与x轴的夹角
            double angle = atan2(pb_y, pb_x);
            if (angle < M_PI)
                angle += 2 * M_PI; // 0~2pi

            if ((angle >= startAngle && angle <= endAngle) || (angle <= endAngle - 360.0)) // 交点位于圆弧之间
                return 1;
        }
    }
    else if (r > pb) // 圆与直线相交，两个交点
    {
        int n = 0;
        // b与交点的距离
        double b_d = sqrt(r * r - pb * pb);

        // 交点bl坐标
        double bl_x = b_x + p1p2_x * b_d;
        double bl_y = b_y + p1p2_y * b_d;
        if (onLineSegment(bl_x, b_y, bl_y, p1_y, p2_x, p2_y)) // 交点位于线段p1p2间
        {
            // pbl与x轴的夹角
            double angle = atan2(bl_y - p_y, bl_x - p_x);
            if (angle < M_PI)
                angle += 2 * M_PI; // 0~2pi

            if ((angle >= startAngle && angle <= endAngle) || (angle <= endAngle - 360.0)) // 交点位于圆弧之间
                n++;
        }

        // 交点br坐标
        double br_x = b_x - p1p2_x * b_d;
        double br_y = b_y - p1p2_y * b_d;
        if (onLineSegment(br_x, b_y, br_y, p1_y, p2_x, p2_y)) // 交点位于线段p1p2间
        {
            // pbr与x轴的夹角
            double angle = atan2(br_y - p_y, br_x - p_x);
            if (angle < M_PI)
                angle += 2 * M_PI; // 0~2pi

            if ((angle >= startAngle && angle <= endAngle) || (angle <= endAngle - 360.0)) // 交点位于圆弧之间
                n++;
        }

        return n;
    }

    return 0;
}

// 圆弧与圆弧的交点
int tt2()
{
    return 0;
}

int main()
{
    return 0;
}