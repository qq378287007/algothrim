#include <algorithm>
#include <cmath>
#include <cstring>
#include <list>
#include <utility>
#include <vector>

#include <iostream>
#include <ctime>
#include <random>
#include <fstream>
using namespace std;

static const double EPS = 1.0e-8;
static int cmp(double v) { return fabs(v) > EPS ? (v > 0.0 ? 1 : -1) : 0; }

struct Point
{
    double x, y;
    int id;
    Point(double a = 0.0, double b = 0.0, int c = -1) : x(a), y(b), id(c) {}
    bool operator<(const Point &a) const { return x < a.x || (cmp(x - a.x) == 0 && y < a.y); }
    bool operator==(const Point &a) const { return cmp(x - a.x) == 0 && cmp(y - a.y) == 0; }
    double dist2(const Point &b) { return (x - b.x) * (x - b.x) + (y - b.y) * (y - b.y); }
};

struct Point3D
{
    double x, y, z;
    Point3D(double a = 0.0, double b = 0.0, double c = 0.0) : x(a), y(b), z(c) {}
    Point3D(const Point &p) { x = p.x, y = p.y, z = p.x * p.x + p.y * p.y; }
    Point3D operator-(const Point3D &a) const { return Point3D(x - a.x, y - a.y, z - a.z); }
    double dot(const Point3D &a) { return x * a.x + y * a.y + z * a.z; }
};

struct Edge
{
    list<Edge>::iterator c; // 用于删除交叉边（LR与LL|RR）
    int id;
    Edge(int d = -1) : id(d) {}
};

//>0逆时针，<0顺时针，=0共线
double cross(const Point &o, const Point &a, const Point &b) { return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x); }

// 平面法向量
Point3D cross(const Point3D &a, const Point3D &b) { return Point3D(a.y * b.z - a.z * b.y, -a.x * b.z + a.z * b.x, a.x * b.y - a.y * b.x); }

//>0圆外，<0圆内，=0圆上
int inCircle(const Point &a, Point b, Point c, const Point &p)
{
    if (cross(a, b, c) < 0)
        swap(b, c);
    Point3D a3(a), b3(b), c3(c), p3(p);
    b3 = b3 - a3, c3 = c3 - a3, p3 = p3 - a3;
    Point3D f = cross(b3, c3);
    return cmp(p3.dot(f));
}

// 线段ab和线段cd是否相交（不含公用一个端点）
static bool intersection(const Point &a, const Point &b, const Point &c, const Point &d)
{
    return cmp(cross(a, c, b)) * cmp(cross(a, b, d)) > 0 && cmp(cross(c, a, d)) * cmp(cross(c, d, b)) > 0;
}

const int MAXV = 10000;
class Delaunay
{
public:
    Point *p;
    int n;
    list<Edge> head[MAXV]; // graph

    void init(int n, Point p[])
    {
        sort(p, p + n);
        for (int i = 0; i < n; i++)
            p[i].id = i;
        this->p = p;
        this->n = n;

        divide(0, n - 1);
    }

    void addEdge(int u, int v)
    {
        head[u].emplace_front(Edge(v));
        head[v].emplace_front(Edge(u));
        head[u].begin()->c = head[v].begin();
        head[v].begin()->c = head[u].begin();
    }

    void divide(int l, int r)
    {
        if (r - l <= 2) // 小于等于3个点
        {
            for (int i = l; i <= r; i++)
                for (int j = i + 1; j <= r; j++)
                    addEdge(i, j);
            return;
        }

        int mid = l + (r - l) / 2;
        divide(l, mid);
        divide(mid + 1, r);

        int nowl = l;
        int nowr = r;

        bool update = true;
        while (update) // 求左凸和右凸的下公切线
        {
            update = false;
            Point ptL = p[nowl], ptR = p[nowr];
            for (list<Edge>::const_iterator it = head[nowl].cbegin(); it != head[nowl].cend(); it++)
            {
                Point t = p[it->id];
                double v = cross(ptR, ptL, t);
                if (cmp(v) > 0 || (cmp(v) == 0 && ptR.dist2(t) < ptR.dist2(ptL)))
                {
                    nowl = it->id, update = true; // 找到左边端点的最下部和最右部的点
                    break;
                }
            }

            if (update)
                continue;

            for (list<Edge>::const_iterator it = head[nowr].cbegin(); it != head[nowr].cend(); it++)
            {
                Point t = p[it->id];
                double v = cross(ptL, ptR, t);
                if (cmp(v) < 0 || (cmp(v) == 0 && ptL.dist2(t) < ptL.dist2(ptR)))
                {
                    nowr = it->id, update = true; // 找到右边端点的最下部和最左部的点
                    break;
                }
            }
        }

        addEdge(nowl, nowr); // 添加下公切线

        while (true)
        {
            Point ptL = p[nowl], ptR = p[nowr];
            int ch = -1, side = 0;
            for (list<Edge>::const_iterator it = head[nowl].cbegin(); it != head[nowl].cend(); it++)
            {
                if (cmp(cross(ptL, ptR, p[it->id])) > 0 && (ch == -1 || inCircle(ptL, ptR, p[ch], p[it->id]) < 0))
                {
                    ch = it->id, side = -1; // 筛选左候选点
                }
            }

            for (list<Edge>::const_iterator it = head[nowr].cbegin(); it != head[nowr].cend(); it++)
            {
                if (cmp(cross(ptR, p[it->id], ptL)) > 0 && (ch == -1 || inCircle(ptL, ptR, p[ch], p[it->id]) < 0))
                {
                    ch = it->id, side = 1; // 筛选右候选点
                }
            }

            if (ch == -1)
                break; // 上公切线

            if (side == -1)
            {
                for (list<Edge>::const_iterator it = head[nowl].cbegin(); it != head[nowl].cend();)
                {
                    if (intersection(ptL, p[it->id], ptR, p[ch])) // LL与LR相交，移除LL线段
                    {
                        head[it->id].erase(it->c);
                        head[nowl].erase(it++);
                    }
                    else
                    {
                        it++;
                    }
                }
                nowl = ch;
            }
            else
            {
                for (list<Edge>::const_iterator it = head[nowr].cbegin(); it != head[nowr].cend();)
                {
                    if (intersection(ptR, p[it->id], ptL, p[ch])) // RR与LR相交，移除RR线段
                    {
                        head[it->id].erase(it->c);
                        head[nowr].erase(it++);
                    }
                    else
                    {
                        it++;
                    }
                }
                nowr = ch;
            }
            addEdge(nowl, nowr); // 添加LR线
        }
    }

    vector<pair<int, int>> getEdge()
    {
        vector<pair<int, int>> ret;
        ret.reserve(n);

        for (int i = 0; i < n; i++) // 遍历点
        {
            const list<Edge> &lEdge = head[i]; // 点相连的所有边
            for (list<Edge>::const_iterator it = lEdge.cbegin(); it != lEdge.cend(); it++)
                if (it->id > i) // 去重
                    ret.push_back(make_pair(p[i].id, p[it->id].id));
        }
        return ret;
    }
};

int main()
{
    default_random_engine e;
    e.seed(time(0));
    normal_distribution<double> u(5.0, 5.0); // 均值为5.0，标准差为5.0

    static const int n = 500;
    Point p[n];
    for (int i = 0; i < n; i++)
    {
        p[i].x = u(e);
        p[i].y = u(e);
    }

    Delaunay d;
    d.init(n, p);

    vector<pair<int, int>> edges = d.getEdge();

    ofstream ofs;
    ofs.open("out2.txt", ios::out);

    ofs << n << endl;
    for (int i = 0; i < n; i++)
        ofs << d.p[i].x << " " << d.p[i].y << endl;

    int m = edges.size();
    ofs << m << endl;
    for (int i = 0; i < m; i++)
        ofs << d.p[edges[i].first].x << " " << d.p[edges[i].first].y << " "
            << d.p[edges[i].second].x << " " << d.p[edges[i].second].y << endl;

    ofs.close();

    return 0;
}
