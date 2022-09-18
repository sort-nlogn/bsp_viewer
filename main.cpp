#include <stdio.h>
#include <stdlib.h>
#include <graphics.h>
#include <vector> 
#include <time.h>
#include <iostream>
#include <math.h>

using namespace std;

#define EPS 0.000000001

#define near_dist 1.0

#define width 1080
#define height 720
#define ratio (float)height / width

#define norm3(v) sqrt(v.x*v.x + v.y*v.y + v.z*v.z)

bool selected_prism = false;
double centroid_prism[3] = {0.0, 0.0, -2.26666667};
double centroid_cube[3] = {-1.5, 0.5, -2.5};

double ground = -1.5;
double light_dir[3] = {0.0, -1.0, 0.0};

struct bsp_tree{
    struct point3d{double x, y, z;};
    struct plane{double a, b, c, d;};

    struct polygon{vector<point3d> points; int color;};

    struct node{
        node *l, *r;
        vector<polygon> polygons;
        node(): l(nullptr), r(nullptr), polygons(vector<polygon>()){}
    };

    enum p_sign{
        ON = 0,
        NEG = 1, 
        POS = 2,
        HALF = 3
    };

    node *root;

    point3d cross(point3d u, point3d v){
        point3d res;
        res.x = u.y * v.z - v.y * u.z;
        res.y = v.x * u.z - u.x * v.z;
        res.z = u.x * v.y - v.x * u.y;
        return res;
    }

    double dot(point3d p1, point3d p2){
        return p1.x * p2.x + p1.y * p2.y + p1.z * p2.z;
    }

    plane construct_plane(point3d p1, point3d p2, point3d p3){
        point3d u = {p2.x - p1.x, p2.y - p1.y, p2.z - p1.z};
        point3d v = {p3.x - p1.x, p3.y - p1.y, p3.z - p1.z};
        point3d n = cross(u, v);
        plane p = {n.x, n.y, n.z, -dot(p1, n)};
        return p;
    }

    bool is_visible(polygon &p, point3d view){
        point3d p1 = p.points[0], p2 = p.points[1], p3 = p.points[2];
        point3d u = {p2.x - p1.x, p2.y - p1.y, p2.z - p1.z};
        point3d v = {p3.x - p1.x, p3.y - p1.y, p3.z - p1.z};
        point3d n = cross(u, v);
        point3d to_view = {view.x - p1.x, view.y - p1.y, view.z - p1.z};
        return dot(n, to_view) > EPS;
    }

    vector<polygon> backface_culling(vector<polygon> polygons, point3d view){
        vector<polygon> visible;

        for(int i = 0; i < polygons.size(); i++){
            if(is_visible(polygons[i], view))
                visible.push_back(polygons[i]);
        }

        return visible;
    }

    void construct_tree(vector<polygon> &polygons, point3d view){
        root = nullptr;
        if (polygons.size() == 0)
            return;

        root = new node();
        vector<polygon> v = backface_culling(polygons, view);
        split_recursively(v, root);
    }

    double get_sign(point3d p, plane pl){
        return pl.a*p.x + pl.b*p.y + pl.c*p.z + pl.d;
    }

    p_sign poly_type(polygon &p, plane pl){
        bool has_pos = false;
        bool has_neg = false;

        for(int i = 0; i < p.points.size(); i++){
            double sgn = get_sign(p.points[i], pl);
            if(sgn > EPS)  has_pos = true;
            if(sgn < -EPS) has_neg = true;
        }

        if(has_pos && has_neg) return HALF;
        if(has_pos) return POS;
        if(has_neg) return NEG;
        return ON;
    }

    point3d plane_to_line_intersection(point3d p1, point3d p2, plane pl){
        point3d rd = {p2.x - p1.x, p2.y - p1.y, p2.z - p1.z};
        point3d n = {pl.a, pl.b, pl.c}; 
        double t = -(pl.d + dot(p1, n)) / dot(rd, n);
        point3d res = {p1.x + t*rd.x, p1.y + t*rd.y, p1.z + t*rd.z};
        return res;
    }

    pair<polygon, polygon> split_polygon(polygon &p, plane pl){
        polygon a, b; a.color = p.color; b.color = p.color;

        for(int i = 0; i < p.points.size(); i++){
            double sgn1 = get_sign(p.points[i], pl), sgn2 = get_sign(p.points[(i + 1) % p.points.size()], pl);
            if(sgn1 > EPS)
                a.points.push_back(p.points[i]);
            else if(sgn1 < -EPS)
                b.points.push_back(p.points[i]);
            else
                continue;
            if(abs(sgn2) < EPS){
                a.points.push_back(p.points[(i + 1) % p.points.size()]);
                b.points.push_back(p.points[(i + 1) % p.points.size()]);
            }else if(sgn1 > EPS && sgn2 < -EPS || sgn1 < -EPS && sgn2 > EPS){
                point3d i_point = plane_to_line_intersection(p.points[i], p.points[(i + 1) % p.points.size()], pl);
                a.points.push_back(i_point);
                b.points.push_back(i_point);
            }
        }
        return pair<polygon, polygon>(a, b);
    }

    void split_recursively(vector<polygon> &polygons, node *start_n){
        int div = polygons.size() / 2;
        start_n->polygons.push_back(polygons[div]);
        plane pl = construct_plane(polygons[div].points[0], polygons[div].points[1], polygons[div].points[2]);
        vector<polygon> neg, pos;

        for(int i = 0; i < polygons.size(); i++){
            if (i == div)
                continue;
            p_sign sgn = poly_type(polygons[i], pl); 
            if(sgn == ON){
                start_n->polygons.push_back(polygons[i]);
            }else if(sgn == POS){
                pos.push_back(polygons[i]);
            }else if(sgn == NEG){
                neg.push_back(polygons[i]);
            }else if(sgn == HALF){
                pair<polygon, polygon> chunks= split_polygon(polygons[i], pl);
                neg.push_back(chunks.second);
                pos.push_back(chunks.first);
            }
        }

        if(neg.size() > 0){
            start_n->l = new node();
            split_recursively(neg, start_n->l);
        }
        if(pos.size() > 0){
            start_n->r = new node();
            split_recursively(pos, start_n->r);
        }
    }

    void dealloc_nodes(node *start_n){
        if(start_n->l)
            dealloc_nodes(start_n->l);
        if(start_n->r)
            dealloc_nodes(start_n->r);
        delete start_n;
    }

    ~bsp_tree(){
        dealloc_nodes(root);
        root = nullptr;
    }
};

void render_shape(bsp_tree::polygon &poly){
    int *shape = (int *)malloc((poly.points.size() + 1) * 2 * sizeof(int));
    for(int i = 0; i < poly.points.size(); i++){
        double x = poly.points[i].x * ratio, y = poly.points[i].y, z = poly.points[i].z;
        x = near_dist * x / (-z*0.8);
        y = near_dist * y / (-z*0.8);
        shape[2 * i] = (int)((x * 0.5 + 0.5) * width);
        shape[2 * i + 1] = (int)((-y * 0.5 + 0.5) * height);
    }

    shape[2 * poly.points.size()] = shape[0];
    shape[2 * poly.points.size() + 1] = shape[1];

    // tone mapping
    int icol = poly.color;
    // float fcol[3] = {(float)GetRValue(icol) / 255, (float)GetGValue(icol) / 255, (float)GetBValue(icol) / 255};
    // fcol[0] /= (fcol[0] + 1.0); fcol[0] = pow(fcol[0], 1.0 / 2.2);
    // fcol[1] /= (fcol[1] + 1.0); fcol[1] = pow(fcol[1], 1.0 / 2.2);
    // fcol[2] /= (fcol[2] + 1.0); fcol[2] = pow(fcol[2], 1.0 / 2.2);
    // icol = RGB((int)(fcol[0] * 255), (int)(fcol[1] * 255), (int)(fcol[2] * 255));

    setcolor(icol);
    setfillstyle(SOLID_FILL, icol);
    fillpoly(poly.points.size() + 1, shape);
    free(shape);
}

void painters_algorithm(bsp_tree *t, bsp_tree::node *v, bsp_tree::point3d view_p){
    if(v == nullptr)
        return;
    if (v->l == nullptr && v->r == nullptr){
        for(int i = 0; i < v->polygons.size(); i++)
            render_shape(v->polygons[i]);
        return;
    }
    auto pl = t->construct_plane(v->polygons[0].points[0], v->polygons[0].points[1], v->polygons[0].points[2]);
    double view_sgn = t->get_sign(view_p, pl);
    if(view_sgn > EPS){
        painters_algorithm(t, v->l, view_p);
        for(int i = 0; i < v->polygons.size(); i++)
            render_shape(v->polygons[i]);
        painters_algorithm(t, v->r, view_p);
    }else if(view_sgn < -EPS){
        painters_algorithm(t, v->r, view_p);
        for(int i = 0; i < v->polygons.size(); i++)
            render_shape(v->polygons[i]);
        painters_algorithm(t, v->l, view_p);
    }else{
        painters_algorithm(t, v->l, view_p);
        for(int i = 0; i < v->polygons.size(); i++)
            render_shape(v->polygons[i]);
        painters_algorithm(t, v->r, view_p);
    }
}



void move(vector<bsp_tree::polygon> &pos, double dx, double dy, double dz){
    for(int i = selected_prism? 0: 5; i < (selected_prism? 5: pos.size()); i++){
        for(int j = 0; j < pos[i].points.size(); j++){
            if(pos[i].points[j].y + dy <= ground / 2 && dy < 0.0)
                return;
        }
    }

    if(selected_prism){
        centroid_prism[0] += dx; centroid_prism[1] += dy; centroid_prism[2] += dz;
    }
    else{
        centroid_cube[0] += dx; centroid_cube[1] += dy; centroid_cube[2] += dz;
    }

    for(int i = selected_prism? 0: 5; i < (selected_prism? 5: pos.size()); i++){
        for(int j = 0; j < pos[i].points.size(); j++){
            pos[i].points[j].x += dx;
            pos[i].points[j].y += dy;
            pos[i].points[j].z += dz;
        }
    }
}

void rotate(vector<bsp_tree::polygon> &pos, double angle_x, double angle_z, double angle_y){
    double *c = selected_prism ? centroid_prism: centroid_cube;

    for(int i = selected_prism? 0: 5; i < (selected_prism? 5: pos.size()); i++){
        for(int j = 0; j < pos[i].points.size(); j++){
            double c1 = cos(angle_x), c2 = cos(angle_y), c3 = cos(angle_z);
            double s1 = sin(angle_x), s2 = sin(angle_y), s3 = sin(angle_z);
            double x = pos[i].points[j].x - c[0], y = pos[i].points[j].y - c[1], z = pos[i].points[j].z - c[2];
            pos[i].points[j].x = y*(c3*s1*s2 - c1*s3) + z*(c1*c3*s2 + s1*s3) + c2*c3*x + c[0];
            pos[i].points[j].y = y*(c1*c3 + s1*s2*s3) + z*(c1*s2*s3 - c3*s1) + c2*s3*x + c[1];
            pos[i].points[j].z = c1*c2*z + c2*s1*y - s2*x + c[2];
        }
    }
}

bool process_keyboard(vector<bsp_tree::polygon> &pos){
    char ch = (char)getch();

    if(ch == 'w'){
        move(pos, 0.0, 0.1, 0.0); return true;
    }else if(ch == 'a'){
        move(pos, -0.1, 0.0, 0.0); return true;
    }else if(ch == 's'){
        move(pos, 0.0, -0.1, 0.0); return true;
    }else if(ch == 'd'){
        move(pos, 0.1, 0.0, 0.0); return true;
    }else if(ch == 'x'){
        selected_prism = selected_prism == 0? 1: 0; return false;
    }else if(ch == 'h'){
        rotate(pos, 0.0, 0.256, 0.0); return true;
    }else if(ch == 'j'){
        rotate(pos, 0.0, -0.256, 0.0); return true;
    }else if(ch == 'k'){
        rotate(pos, 0.0, 0.0, 0.256); return true;
    }else if(ch == 'l'){
        rotate(pos, 0.0, 0.0, -0.256); return true;
    }else if(ch == 'n'){
        rotate(pos, -0.256, 0.0, 0.0); return true;
    }else if(ch == 'm'){
        rotate(pos, 0.256, 0.0, 0.0); return true;
    }else if((int)ch == 75){
        light_dir[0] -= 0.1; return false;
    }else if((int)ch == 77){
        light_dir[0] += 0.1; return false;
    }
    return false;
}

void render_shadow(vector<bsp_tree::polygon> &polygons){

    for(int i = 0; i < polygons.size(); i++){
        bsp_tree::polygon poly;
        for(int j = 0; j < polygons[i].points.size(); j++){
            bsp_tree::point3d pt = polygons[i].points[j];

            double t = (ground - pt.y) / light_dir[1];

            pt.x += t*light_dir[0]; pt.y = ground; pt.z += t*light_dir[2];   
            poly.points.push_back(pt);
        }
        poly.color = RGB(0, 0, 0);
        render_shape(poly);
    }
}

int main(){
    srand(time(NULL));
    vector<bsp_tree::polygon> polygons =
    { // prism
     { {{-0.5, -0.5, -2.0}, {0.0, -0.5, -2.8}, {0.5, -0.5, -2.0}},  RGB(204, 36, 29) },
     { {{-0.5,  0.5, -2.0}, {0.5,  0.5, -2.0}, {0.0,  0.5, -2.8}},  RGB(204, 36, 29) },
     { {{-0.5, -0.5, -2.0}, {0.5, -0.5, -2.0}, {0.5,  0.5, -2.0}, {-0.5,  0.5, -2.0}}, RGB(152, 151, 26)},
     { {{0.0, -0.5, -2.8}, {-0.5, -0.5, -2.0}, {-0.5, 0.5, -2.0}, {0.0,  0.5, -2.8}},  RGB(215, 153, 33)},
     { {{0.5, -0.5, -2.0}, {0.0, -0.5, -2.8}, {0.0,  0.5, -2.8}, {0.5,  0.5, -2.0}},   RGB(69, 133, 136)},
        // cube (now it's actually a pyramid (it should be))
     { {{-2.0, 1.0, -2.0}, {-2.0, 0.0, -2.0}, {-1.0, 0.0, -2.0}, {-1.0, 1.0, -2.0}},   RGB(177, 98, 134)},
     { {{-1.0, 0.0, -2.0}, {-1.0, 0.0, -3.0}, {-1.0, 1.0, -3.0}, {-1.0, 1.0, -2.0}},   RGB(104, 157, 122)},
     { {{-2.0, 0.0, -2.0}, {-2.0, 1.0, -2.0}, {-2.0, 1.0, -3.0}, {-2.0, 0.0, -3.0}},   RGB(214, 93, 14)  },
     { {{-2.0, 0.0, -2.0}, {-2.0, 0.0, -3.0}, {-1.0, 0.0, -3.0}, {-1.0, 0.0, -2.0}},   RGB(146, 131, 116)},
     { {{-2.0, 1.0, -2.0}, {-1.0, 1.0, -2.0}, {-1.0, 1.0, -3.0}, {-2.0, 1.0, -3.0}},   RGB(251, 73, 53)  },
     { {{-2.0, 0.0, -3.0}, {-2.0, 1.0, -3.0}, {-1.0, 1.0, -3.0}, {-1.0, 0.0, -3.0}},   RGB(184, 187, 38) }
    };

    bsp_tree *t = new bsp_tree;

    bsp_tree::point3d view = {0.0, 0.0, 0.0};
    t->construct_tree(polygons, view);

    int win = initwindow(width, height, "3d");
    setbkcolor(RGB(24, 24, 24));

    while(1){
        painters_algorithm(t, t->root, view);
        if(kbhit()){
            bool needs_4update = process_keyboard(polygons);
            if(needs_4update){
                t->dealloc_nodes(t->root);
                t->construct_tree(polygons, view);
            }
        }
        render_shadow(polygons);
        swapbuffers();
        clearviewport();
        Sleep(1);
    }
    return 0;
}