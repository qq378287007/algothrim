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

const double EPS = 1e-8;
const int MAXV = 10000;

int cmp(double v) { return fabs(v) > EPS ? (v > 0 ? 1 : -1) : 0; }

struct Point
{
    double x, y;
    int id;
    Point(double a = 0.0, double b = 0.0, int c = -1) : x(a), y(b), id(c) {}
    bool operator<(const Point &a) const { return x < a.x || (fabs(x - a.x) < EPS && y < a.y); }
    bool operator==(const Point &a) const { return fabs(x - a.x) < EPS && fabs(y - a.y) < EPS; }
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
    list<Edge>::iterator c;
    int id;
    Edge(int d = 0) : id(d) {}
};

double cross(const Point &o, const Point &a, const Point &b) { return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x); }

Point3D cross(const Point3D &a, const Point3D &b) { return Point3D(a.y * b.z - a.z * b.y, -a.x * b.z + a.z * b.x, a.x * b.y - a.y * b.x); }

int inCircle(const Point &a, Point b, Point c, const Point &p)
{
    if (cross(a, b, c) < 0)
        swap(b, c);
    Point3D a3(a), b3(b), c3(c), p3(p);
    b3 = b3 - a3, c3 = c3 - a3, p3 = p3 - a3;
    Point3D f = cross(b3, c3);
    return cmp(p3.dot(f)); // check same direction, in: < 0, on: = 0, out: > 0
}

// seg(a, b) and seg(c, d)
int intersection(const Point &a, const Point &b, const Point &c, const Point &d)
{
    return cmp(cross(a, c, b)) * cmp(cross(a, b, d)) > 0 && cmp(cross(c, a, d)) * cmp(cross(c, d, b)) > 0;
}

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
        head[u].push_front(Edge(v));
        head[v].push_front(Edge(u));
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

        list<Edge>::iterator it;
        int nowl = l;
        int nowr = r;

        for (int update = 1; update;)
        {
            // find left and right convex, lower common tangent
            update = 0;
            Point ptL = p[nowl], ptR = p[nowr];
            for (it = head[nowl].begin(); it != head[nowl].end(); it++)
            {
                Point t = p[it->id];
                double v = cross(ptR, ptL, t);
                if (cmp(v) > 0 || (cmp(v) == 0 && ptR.dist2(t) < ptR.dist2(ptL)))
                {
                    nowl = it->id, update = 1;
                    break;
                }
            }

            if (update)
                continue;
                
            for (it = head[nowr].begin(); it != head[nowr].end(); it++)
            {
                Point t = p[it->id];
                double v = cross(ptL, ptR, t);
                if (cmp(v) < 0 || (cmp(v) == 0 && ptL.dist2(t) < ptL.dist2(ptR)))
                {
                    nowr = it->id, update = 1;
                    break;
                }
            }
        }

        addEdge(nowl, nowr); // add tangent

        for (int update = 1; true;)
        {
            update = 0;
            Point ptL = p[nowl], ptR = p[nowr];
            int ch = -1, side = 0;
            for (it = head[nowl].begin(); it != head[nowl].end(); it++)
            {
                if (cmp(cross(ptL, ptR, p[it->id])) > 0 && (ch == -1 || inCircle(ptL, ptR, p[ch], p[it->id]) < 0))
                {
                    ch = it->id, side = -1;
                }
            }
            for (it = head[nowr].begin(); it != head[nowr].end(); it++)
            {
                if (cmp(cross(ptR, p[it->id], ptL)) > 0 && (ch == -1 || inCircle(ptL, ptR, p[ch], p[it->id]) < 0))
                {
                    ch = it->id, side = 1;
                }
            }
            if (ch == -1)
                break; // upper common tangent
            if (side == -1)
            {
                for (it = head[nowl].begin(); it != head[nowl].end();)
                {
                    if (intersection(ptL, p[it->id], ptR, p[ch]))
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
                addEdge(nowl, nowr);
            }
            else
            {
                for (it = head[nowr].begin(); it != head[nowr].end();)
                {
                    if (intersection(ptR, p[it->id], ptL, p[ch]))
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
                addEdge(nowl, nowr);
            }
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
            {
                if (it->id > i) // 去重
                    ret.push_back(make_pair(p[i].id, p[it->id].id));
            }
        }
        return ret;
    }
};

int main()
{
    std::default_random_engine e;
    std::normal_distribution<double> u(5.0, 5.0); // 均值为5.0，标准差为5.0
    e.seed(time(0));

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
