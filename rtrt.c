/*
 *       _        _         
 *  _ __| |_ _ __| |_   ___ 
 * | '__| __| '__| __| / __|
 * | |  | |_| |  | |_ | (__ 
 * |_|   \__|_|   \__(_)___|
 *                          
 * rtrt.c -- a realtime raytracer written in C by Mark VandeWettering
 *
 * Copyright 2000, Mark T. VandeWettering
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <signal.h>

#include <GL/glut.h>


#define DEGTORAD(x)	((x)*M_PI/180.0)
#define XSIZE	(1280)
#define YSIZE	(720)

typedef float Vec[3] ;
typedef struct t_pixel {
    unsigned char r, g, b ;
} Color ;

float fov = DEGTORAD(45.0) ;

typedef struct t_sphere {
    Vec C ;
    float r2 ;
    Vec Q ;
} Sphere ;

typedef struct t_light {
    Vec L ;
} Light ;

Light light = {
    {-5.0f, 10.0f, 4.0f},
} ;

typedef struct t_ray {
    Vec P, D ;
} Ray ;

Sphere sphere = {
	{0.0, 0.0, 10.0},
	1.0,
	{-1.0, 0.0, -10.0}
} ;

float sint, cost ;

#define EPS	(.0001f)

static void
VecNormalize(Vec v)
{
    float l = v[0]*v[0]+v[1]*v[1]+v[2]*v[2] ;

    l = 1.0/sqrtf(l+0.0001) ;
    v[0] *= l ; v[1] *= l ; v[2] *= l ;
}

static void
VecComb(float a, Vec A, float b, Vec B, Vec R)
{
    R[0] = a*A[0] + b*B[0] ;
    R[1] = a*A[1] + b*B[1] ;
    R[2] = a*A[2] + b*B[2] ;
}

#define VecDot(a,b)	((a)[0]*(b)[0]+(a)[1]*(b)[1]+(a)[2]*(b)[2])
	
static void
SpecularDirection(Vec I, Vec N, Vec R)
{
	VecComb(1.0/fabs(VecDot(I,N)), I, 2.0, N, R);
	VecNormalize(R);
}

int SphereIntersect(Ray *r, Color *p, float *dist, int level);
int PlaneIntersect(Ray *r, Color *c, int level);

int
SphereIntersect(Ray *r, Color *p, float *dist, int level)
{
    float discrim ;
    float t, a, b, c ;
    Vec D ;

    if (level > 4) return 0 ;

    D[0] = r->D[0] ; D[1] = r->D[1] ; D[2] = r->D[2] ;
    sphere.Q[0] = r->P[0] - sphere.C[0] ;
    sphere.Q[1] = r->P[1] - sphere.C[1] ;
    sphere.Q[2] = r->P[2] - sphere.C[2] ;

    a = D[0]*D[0] + D[1]*D[1] + D[2]*D[2] ;

    b = 2.0 * (sphere.Q[0]*D[0] + 
	       sphere.Q[1]*D[1] +
	       sphere.Q[2]*D[2]) ;

    c = sphere.Q[0]*sphere.Q[0] + sphere.Q[1]*sphere.Q[1] +
	sphere.Q[2]*sphere.Q[2] - sphere.r2 ;


    discrim = b*b - 4.0f * a * c ;

    if (discrim < 0)
	return 0 ;

    discrim = sqrtf(discrim) ;
    a = 1.0f / (2.0f * a) ;

    if ((t = (-b - discrim)*a)  < EPS)
	if ((t = (-b + discrim)*a) < EPS)
	     return 0 ;

    
    *dist = t ;

    {
	Ray reflray ;
	Color tc ;
	Vec P, N, L ;
        float dp ;

	P[0] = r->P[0] + t * r->D[0] ;
	P[1] = r->P[1] + t * r->D[1] ;
	P[2] = r->P[2] + t * r->D[2] ;

        N[0] = P[0] - sphere.C[0] ;
        N[1] = P[1] - sphere.C[1] ;
        N[2] = P[2] - sphere.C[2] ;
	VecNormalize(N) ;

        L[0] = light.L[0] - P[0] ;
        L[1] = light.L[1] - P[1] ;
        L[2] = light.L[2] - P[2] ;
	VecNormalize(L) ;



	if ((dp = N[0]*L[0]+N[1]*L[1]+N[2]*L[2]) > 0) {
	    p->r = 16 + dp * 32 ;
	    p->g = 16 + dp * 16 ;
	    p->b = 16 + dp * 8 ;

	} else {
	    p->r = 16 ;
	    p->g = 8 ;
	    p->b = 4 ;
	}
#if 1
	{
 	    Vec H ;
	    H[0] = (L[0]-r->D[0]) ;
	    H[1] = (L[1]-r->D[1]) ;
	    H[2] = (L[2]-r->D[2]) ;
	    VecNormalize(H) ;

	    if ((dp = H[0]*N[0]+H[1]*N[1]+H[2]*N[2]) > 0) {
		dp *= dp ;
		dp *= dp ;
		dp *= dp ;
		dp *= dp ;
		dp *= dp ;
		dp *= dp ;
		dp *= dp ;
		dp *= 128 ;
		p->r += 32 + dp ;
		p->g += 32 + dp / 2 ;
		p->b += 32 + dp / 4 ;
	    }
	}
#endif

#if 1
	reflray.P[0] = P[0] ;
	reflray.P[1] = P[1] ;
	reflray.P[2] = P[2] ;

	SpecularDirection(r->D, N, reflray.D) ;
	
	if (PlaneIntersect(&reflray, &tc, level+1)) {
	    p->r += (tc.r >> 2) ;
	    p->g += (tc.g >> 3) ;
	    p->b += (tc.b >> 4) ;
	}
#endif
    }


    return 1 ;
}


int
PlaneIntersect(Ray *r, Color *c, int level)
{
    float t ;
    float yd = r->D[1] ;
    float x, z ;
    float nx, nz ;
    int ix, iz ;
    float dist ;
    Ray shadray, reflray ;
    Color tc ;

    if (yd > -EPS)
	return 0 ;

    t = (-1.0 - r->P[1]) / yd ;

    x = r->P[0] + t * r->D[0] ;
    z = r->P[2] + t * r->D[2] ;

    shadray.P[0] = reflray.P[0] = x ;
    shadray.P[1] = reflray.P[1] = -1 ;
    shadray.P[2] = reflray.P[2] = z ;

    reflray.D[0] = r->D[0] ;
    reflray.D[1] = -r->D[1] ;
    reflray.D[2] = r->D[2] ;

    shadray.D[0] = light.L[0] - x ;
    shadray.D[1] = light.L[1] + 1 ;
    shadray.D[2] = light.L[2] - z ;

    nx = x*cost+z*sint ;
    nz = z*cost-x*sint ;

    if (nx < 0)
	ix = (int) nx - 1 ;
    else
	ix = (int) nx ;

    if (nz < 0)
	iz = (int) nz - 1 ;
    else
	iz = (int) nz ;

    if (SphereIntersect(&shadray, &tc, &dist, level+1)) {
	if ((ix+iz)&1) {
	    c->r = 16 ;
	    c->g = 0 ;
	    c->b = 0 ;
	} else {
	    c->r = 0 ;
	    c->g = 16 ;
	    c->b = 0 ;
	}
    } else {
	if ((ix+iz)&1) {
	    c->r = 128 ;
	    c->g = 0 ;
	    c->b = 0 ;
	} else {
	    c->r = 0 ;
	    c->g = 128 ;
	    c->b = 0 ;
	}
    }

#if 0
    if (SphereIntersect(&reflray, &tc, &dist, level+1)) {
	c->r += tc.r / 2 ;
	c->g += tc.g / 2 ;
	c->b += tc.b / 2 ;
    }
#endif

    return 1 ;
}

unsigned char image[XSIZE*YSIZE*3] ;

int Frames = 0 ;
int Seconds = 0 ;

static void
updateseconds(void)
{
    Seconds ++ ;
    signal(SIGALRM, updateseconds) ;
    alarm(1) ;
    fprintf(stderr, "%f frames/second\n", (double) Frames/Seconds) ;
}

static void
redraw(void)
{
    glutPostRedisplay() ;
}

static 
void animate(void)
{
    Ray r ;
    Color C ;
    float hfov = tan(fov/2.0) ;
    float vfov = tan(fov/2.0) * YSIZE / XSIZE ;
    float hstep = 2.0 * hfov / XSIZE ;
    float vstep = 2.0 * vfov / YSIZE ; 
    int i, x, y ;
    unsigned char *bp ;
    static float s = 0.0, t = 0.0, u=0.0  ;
    static float dist ;
    Ray sray ;

    s += 0.02 ; if (s > 2.0 * M_PI) s = 0.0 ;
    light.L[0] = 50. * cosf(s) ;
    light.L[2] = 50. * sinf(s) ;

    t += 0.1 ; if (t > 2.0 * M_PI) t = 0.0 ;
    sphere.C[1] = fabsf(sinf(t)) ;

    u += 0.017 ; if (u > 2.0 * M_PI) u -= 2.0f*M_PI ;
    sint = sinf(u) ;
    cost = cosf(u) ;
    
    bp = image ;

    r.P[0] = 0.0f ; r.P[1] = 0.0f ; r.P[2] = 0.0f ;
    r.D[1] = -vfov ; r.D[2] = 1.0 ;

    for (y=0; y<YSIZE; y++, r.D[1] += vstep) {
	r.D[0] = -hfov ;
	for (x=0; x<XSIZE; x++, r.D[0] += hstep) {

	    if (SphereIntersect(&r, &C, &dist, 0)) {
		*bp++ = C.r ;
		*bp++ = C.g ;
		*bp++ = C.b ;
	    } else if (PlaneIntersect(&r, &C, 0)) {
		*bp++ = C.r ;
		*bp++ = C.g ;
		*bp++ = C.b ;
	    } else {
		*bp++ = 0 ;
		*bp++ = 0 ;
		*bp++ = 0 ;
	    }
	}
    }

    glDrawPixels(XSIZE, YSIZE, GL_RGB, GL_UNSIGNED_BYTE, (GLvoid *) image) ;
    Frames ++ ;

}

main(int argc, char *argv[])
{
    GLenum type ;

    glutInit(&argc, argv);

    type = GLUT_RGB ;

    glutInitDisplayMode(type) ;
    glutInitWindowSize(XSIZE, YSIZE) ;
    glutCreateWindow("Realtime Raytracing Demo") ;

    glutIdleFunc(redraw) ;
    glutDisplayFunc(animate) ;
    signal(SIGALRM, updateseconds) ;
    alarm(1) ;
    glutMainLoop() ;

    exit(0) ;
}
