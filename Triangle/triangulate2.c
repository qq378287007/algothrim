#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#define EPSILON 1.0e-14
#define TRUE 1
#define FALSE 0

typedef struct
{
    int p1, p2, p3;
} ITRIANGLE; // 三角形顶点序号

typedef struct
{
    int p1, p2;
} IEDGE; // 边端点序号

typedef struct
{
    double x, y, z;
} XYZ; // 点坐标，z未使用

int XYZCompare(const void *v1, const void *v2) // 用于qsort
{
    XYZ *p1 = (XYZ *)v1;
    XYZ *p2 = (XYZ *)v2;
    if (p1->x < p2->x)
        return -1;
    else if (p1->x > p2->x)
        return 1;
    else
        return 0;
}

// 点是否在三角形外接圆内（含圆上）
int CircumCircle(double xp, double yp,                                             // 输入点
                 double x1, double y1, double x2, double y2, double x3, double y3, // 三角形顶点
                 double *xc, double *yc,                                           // 外心
                 double *rsqr)                                                     // 半径平方
{
    double fabsy1y2 = fabs(y1 - y2);
    double fabsy2y3 = fabs(y2 - y3);

    if (fabsy1y2 < EPSILON && fabsy2y3 < EPSILON) // 三点y值相同
        return FALSE;

    double m1, m2, mx1, mx2, my1, my2;
    if (fabsy1y2 < EPSILON) // y1和y2相同
    {
        m2 = -(x3 - x2) / (y3 - y2);
        mx2 = (x2 + x3) / 2.0;        // x2和x3中间
        my2 = (y2 + y3) / 2.0;        // y2和y3中间
        *xc = (x2 + x1) / 2.0;        // x1和x2中间
        *yc = m2 * (*xc - mx2) + my2; // 线段(mx2, my2)(xc, yc)，与线段(x2, y2)(x3, y3)垂直，矢量点乘为0
    }
    else if (fabsy2y3 < EPSILON) // y2和y3相同
    {
        m1 = -(x2 - x1) / (y2 - y1);
        mx1 = (x1 + x2) / 2.0;
        my1 = (y1 + y2) / 2.0;
        *xc = (x3 + x2) / 2.0;
        *yc = m1 * (*xc - mx1) + my1;
    }
    else
    {
        m1 = -(x2 - x1) / (y2 - y1);
        m2 = -(x3 - x2) / (y3 - y2);
        mx1 = (x1 + x2) / 2.0;
        mx2 = (x2 + x3) / 2.0;
        my1 = (y1 + y2) / 2.0;
        my2 = (y2 + y3) / 2.0;
        *xc = (m1 * mx1 - m2 * mx2 + my2 - my1) / (m1 - m2);
        if (fabsy1y2 > fabsy2y3)
            *yc = m1 * (*xc - mx1) + my1;
        else
            *yc = m2 * (*xc - mx2) + my2;
    }

    double dx, dy, drsqr;
    dx = x2 - *xc;
    dy = y2 - *yc;
    *rsqr = dx * dx + dy * dy;

    dx = xp - *xc;
    dy = yp - *yc;
    drsqr = dx * dx + dy * dy;

    return ((drsqr - *rsqr) <= EPSILON ? TRUE : FALSE);
}

// 输入nv个点pxyz（实际大小nv+3，按x升序排列）
// 生成ntri个三角形v（按3*nv大小分配，v中三角形以一致的顺时针顺序排列）
int Triangulate(int nv, XYZ *pxyz, ITRIANGLE *v, int *ntri)
{
    int trimax = 4 * nv;

    int status = 0;
    int *complete = NULL;
    IEDGE *edges = NULL;

    // Allocate memory for the completeness list, flag for each triangle
    if ((complete = malloc(trimax * sizeof(int))) == NULL)
    {
        status = 1;
        goto skip;
    }

    int emax = 200;
    // Allocate memory for the edge list
    if ((edges = malloc(emax * (long)sizeof(IEDGE))) == NULL)
    {
        status = 2;
        goto skip;
    }

    //  Find the maximum and minimum vertex bounds. This is to allow calculation of the bounding triangle
    double xmin, xmax, ymin, ymax, xmid, ymid;
    double dx, dy, dmax;
    xmin = pxyz[0].x;
    ymin = pxyz[0].y;
    xmax = xmin;
    ymax = ymin;
    for (int i = 1; i < nv; i++)
    {
        if (pxyz[i].x < xmin)
            xmin = pxyz[i].x;
        if (pxyz[i].x > xmax)
            xmax = pxyz[i].x;
        if (pxyz[i].y < ymin)
            ymin = pxyz[i].y;
        if (pxyz[i].y > ymax)
            ymax = pxyz[i].y;
    }
    dx = xmax - xmin;
    dy = ymax - ymin;
    dmax = (dx > dy) ? dx : dy;
    xmid = (xmax + xmin) / 2.0;
    ymid = (ymax + ymin) / 2.0;

    // Set up the supertriangle
    // This is a triangle which encompasses all the sample points.
    // The supertriangle coordinates are added to the end of the vertex list.
    // The supertriangle is the first triangle in the triangle list.
    pxyz[nv + 0].x = xmid - 20 * dmax;
    pxyz[nv + 0].y = ymid - dmax;
    pxyz[nv + 0].z = 0.0;
    pxyz[nv + 1].x = xmid;
    pxyz[nv + 1].y = ymid + 20 * dmax;
    pxyz[nv + 1].z = 0.0;
    pxyz[nv + 2].x = xmid + 20 * dmax;
    pxyz[nv + 2].y = ymid - dmax;
    pxyz[nv + 2].z = 0.0;
    v[0].p1 = nv;
    v[0].p2 = nv + 1;
    v[0].p3 = nv + 2;
    complete[0] = FALSE;
    *ntri = 1;

    //  Include each point one at a time into the existing mesh
    for (int i = 0; i < nv; i++)
    {
        double xp = pxyz[i].x;
        double yp = pxyz[i].y;

        int nedge = 0;

        // Set up the edge buffer.
        // If the point (xp, yp) lies inside the circumcircle then
        // the three edges of that triangle are added to the edge buffer
        // and that triangle is removed.
        for (int j = 0; j < (*ntri); j++)
        {
            if (complete[j])
                continue;
            double x1 = pxyz[v[j].p1].x;
            double y1 = pxyz[v[j].p1].y;
            double x2 = pxyz[v[j].p2].x;
            double y2 = pxyz[v[j].p2].y;
            double x3 = pxyz[v[j].p3].x;
            double y3 = pxyz[v[j].p3].y;
            double xc, yc, r;
            int inside = CircumCircle(xp, yp, x1, y1, x2, y2, x3, y3, &xc, &yc, &r);
            if (xc < xp && ((xp - xc) * (xp - xc)) > r)
                complete[j] = TRUE;
            if (inside)
            {
                // Check that we haven't exceeded the edge list size
                if (nedge + 3 >= emax)
                {
                    emax += 100;
                    if ((edges = realloc(edges, emax * (long)sizeof(IEDGE))) == NULL)
                    {
                        status = 3;
                        goto skip;
                    }
                }
                edges[nedge + 0].p1 = v[j].p1;
                edges[nedge + 0].p2 = v[j].p2;
                edges[nedge + 1].p1 = v[j].p2;
                edges[nedge + 1].p2 = v[j].p3;
                edges[nedge + 2].p1 = v[j].p3;
                edges[nedge + 2].p2 = v[j].p1;
                nedge += 3;
                v[j] = v[(*ntri) - 1];
                complete[j] = complete[(*ntri) - 1];
                (*ntri)--;
                j--;
            }
        }

        // Tag multiple edges
        // Note: if all triangles are specified anticlockwise then
        // all interior edges are opposite pointing in direction.
        for (int j = 0; j < nedge - 1; j++)
        {
            for (int k = j + 1; k < nedge; k++)
            {
                if ((edges[j].p1 == edges[k].p2) && (edges[j].p2 == edges[k].p1))
                {
                    edges[j].p1 = -1;
                    edges[j].p2 = -1;
                    edges[k].p1 = -1;
                    edges[k].p2 = -1;
                }
                // Shouldn't need the following, see note above
                if ((edges[j].p1 == edges[k].p1) && (edges[j].p2 == edges[k].p2))
                {
                    edges[j].p1 = -1;
                    edges[j].p2 = -1;
                    edges[k].p1 = -1;
                    edges[k].p2 = -1;
                }
            }
        }

        // Form new triangles for the current point Skipping over any tagged edges.
        // All edges are arranged in clockwise order.
        for (int j = 0; j < nedge; j++)
        {
            if (edges[j].p1 < 0 || edges[j].p2 < 0)
                continue;
            if ((*ntri) >= trimax)
            {
                status = 4;
                goto skip;
            }
            v[*ntri].p1 = edges[j].p1;
            v[*ntri].p2 = edges[j].p2;
            v[*ntri].p3 = i;
            complete[*ntri] = FALSE;
            (*ntri)++;
        }
    }

    // Remove triangles with supertriangle vertices
    // These are triangles which have a vertex number greater than nv
    for (int i = 0; i < (*ntri); i++)
    {
        if (v[i].p1 >= nv || v[i].p2 >= nv || v[i].p3 >= nv)
        {
            v[i] = v[(*ntri) - 1];
            (*ntri)--;
            i--;
        }
    }

skip:
    free(edges);
    free(complete);
    return status;
}

int main(int argc, char **argv)
{
    srand((unsigned)time(NULL));

    static const int nv = 500;

    XYZ *p = malloc(sizeof(XYZ) * (nv + 3));
    for (int i = 0; i < nv; i++)
    {
        p[i].x = 10.0 * (rand() % 1000) / (1000 - 1);
        p[i].y = 10.0 * (rand() % 1000) / (1000 - 1);
        p[i].z = 0.0;
    }
    fprintf(stderr, "Read %d points\n", nv);

    qsort(p, nv, sizeof(XYZ), XYZCompare);

    ITRIANGLE *v = malloc(3 * nv * sizeof(ITRIANGLE));

    int ntri = 0;
    Triangulate(nv, p, v, &ntri);
    fprintf(stderr, "Formed %d triangles\n", ntri);

    FILE *fptr = fopen("out.txt", "w");
    fprintf(fptr, "%d\n", nv);
    for (int i = 0; i < nv; i++)
        fprintf(fptr, "%g %g\n", p[i].x, p[i].y);
    fprintf(fptr, "%d\n", ntri);
    for (int i = 0; i < ntri; i++)
        fprintf(fptr, "%g %g %g %g %g %g\n",
                p[v[i].p1].x, p[v[i].p1].y,
                p[v[i].p2].x, p[v[i].p2].y,
                p[v[i].p3].x, p[v[i].p3].y);
    fclose(fptr);

    free(v);
    free(p);
    return 0;
}
