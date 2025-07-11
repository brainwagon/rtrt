import taichi as ti
import math

## Converted by using github copilot functions
## Took about 30 minutes, and close to a dozen revisions.
## In the end, it is not actually faster, which is a bit
## surprising.

ti.init(arch=ti.gpu)

XSIZE = 1280
YSIZE = 720

Vec3 = ti.types.vector(3, float)
image = ti.Vector.field(3, float, shape=(XSIZE, YSIZE))

fov = math.radians(45.0)
hfov = math.tan(fov / 2.0)
vfov = -math.tan(fov / 2.0) * YSIZE / XSIZE

# Animation variables
s = ti.field(float, shape=())
t = ti.field(float, shape=())
u = ti.field(float, shape=())
s[None] = 0.0
t[None] = 0.0
u[None] = 0.0

# Sphere and light
sphere_C = ti.Vector.field(3, float, shape=())
sphere_C[None] = [0.0, 0.0, 10.0]
sphere_r2 = 1.0
light_L = ti.Vector.field(3, float, shape=())
light_L[None] = [-5.0, 10.0, 4.0]

sint = ti.field(float, shape=())
cost = ti.field(float, shape=())

EPS = 1e-4
MAX_DEPTH = 4

@ti.func
def VecNormalize(v):
    l = ti.sqrt(v.dot(v) + 1e-4)
    return v / l

@ti.func
def VecComb(a, A, b, B):
    return a * A + b * B

@ti.func
def SpecularDirection(I, N):
    R = VecComb(1.0 / ti.abs(I.dot(N)), I, 2.0, N)
    return VecNormalize(R)

@ti.func
def SphereIntersect(P, D):
    Q = P - sphere_C[None]
    a = D.dot(D)
    b = 2.0 * Q.dot(D)
    c = Q.dot(Q) - sphere_r2
    discrim = b * b - 4.0 * a * c

    hit = False
    t = 0.0
    hit_pos = Vec3([0.0, 0.0, 0.0])
    N = Vec3([0.0, 0.0, 0.0])

    if discrim >= 0:
        sqrt_discrim = ti.sqrt(discrim)
        a2 = 1.0 / (2.0 * a)
        t1 = (-b - sqrt_discrim) * a2
        t2 = (-b + sqrt_discrim) * a2

        t_candidate = 1e30
        if t1 >= EPS and t1 < t_candidate:
            t_candidate = t1
        if t2 >= EPS and t2 < t_candidate:
            t_candidate = t2

        if t_candidate < 1e30:
            hit = True
            t = t_candidate
            hit_pos = P + t * D
            N = VecNormalize(hit_pos - sphere_C[None])

    return hit, t, hit_pos, N

@ti.func
def PlaneIntersect(P, D):
    yd = D[1]
    hit = False
    t = 0.0
    hit_pos = Vec3([0.0, 0.0, 0.0])
    if yd < -EPS:
        t_val = (-1.0 - P[1]) / yd
        hit = True
        t = t_val
        hit_pos = P + t * D
    return hit, t, hit_pos

@ti.func
def Raytrace(P_init, D_init):
    # Iterative rather than recursive!
    color = Vec3([0.0, 0.0, 0.0])
    attenuation = Vec3([1.0, 1.0, 1.0])

    P = P_init
    D = D_init

    for depth in range(MAX_DEPTH + 1):
        # Try sphere
        hit_sphere, t_sphere, phit, N = SphereIntersect(P, D)
        if hit_sphere:
            # Diffuse shading
            L = VecNormalize(light_L[None] - phit)
            dp = N.dot(L)
            r = 16.0
            g = 8.0
            b = 4.0
            if dp > 0:
                r = 16.0 + dp * 32.0
                g = 16.0 + dp * 16.0
                b = 16.0 + dp * 8.0

            # Specular highlight
            H = VecNormalize(L - D)
            spec_dp = H.dot(N)
            if spec_dp > 0:
                temp = spec_dp
                for _ in range(7):
                    temp *= temp
                spec = temp * 128.0
                r += 32.0 + spec
                g += 32.0 + spec / 2.0
                b += 32.0 + spec / 4.0

            # Add local color
            color += attenuation * Vec3([r, g, b]) / 255.0

            # Reflection (sphere sees only the floor, not itself)
            refl_origin = phit
            refl_dir = SpecularDirection(D, N)

            # Check if reflection hits the plane
            hit_plane, t_plane, phit_plane = PlaneIntersect(refl_origin, refl_dir)
            # Predefine col for use in all branches
            col = Vec3([0.0, 0.0, 0.0])
            if hit_plane:
                # Update attenuation for next bounce (faint reflection)
                attenuation = attenuation * Vec3([0.25, 0.125, 0.0625])
                P = phit_plane
                D = refl_dir
                continue  # Next bounce: only reflect off the plane
            break  # Ray leaves scene or doesn't hit plane
        else:
            # Try plane
            hit_plane, t_plane, phit_plane = PlaneIntersect(P, D)
            col = Vec3([0.0, 0.0, 0.0])  # Predefine col before use
            if hit_plane:
                nx = phit_plane[0] * cost[None] + phit_plane[2] * sint[None]
                nz = phit_plane[2] * cost[None] - phit_plane[0] * sint[None]
                ix = ti.floor(nx)
                iz = ti.floor(nz)

                L = light_L[None] - phit_plane
                L = VecNormalize(L)
                shadow_origin = phit_plane + 1e-3 * L
                in_shadow, _, _, _ = SphereIntersect(shadow_origin, L)

                if in_shadow:
                    if (int(ix + iz) & 1):
                        col = Vec3([16, 0, 0])
                    else:
                        col = Vec3([0, 16, 0])
                else:
                    if (int(ix + iz) & 1):
                        col = Vec3([128, 0, 0])
                    else:
                        col = Vec3([0, 128, 0])
                color += attenuation * (col / 255.0)
            # If hit nothing or hit plane, we are done
            break

    # If nothing hit, color stays black
    return ti.min(color, 1.0)

@ti.kernel
def animate():
    for x, y in image:
        # Camera ray
        ucoord = -hfov + 2.0 * hfov * x / XSIZE
        vcoord = -vfov + 2.0 * vfov * y / YSIZE
        P = Vec3([0.0, 0.0, 0.0])
        D = Vec3([ucoord, -vcoord, 1.0])
        D = VecNormalize(D)
        image[x, y] = Raytrace(P, D)

def update_scene():
    # Animate light and sphere position, rotation
    s[None] += 0.02
    if s[None] > 2.0 * math.pi:
        s[None] = 0.0
    light_L[None][0] = 50.0 * math.cos(s[None])
    light_L[None][2] = 50.0 * math.sin(s[None])

    t[None] += 0.1
    if t[None] > 2.0 * math.pi:
        t[None] = 0.0
    sphere_C[None][1] = abs(math.sin(t[None]))

    u[None] += 0.017
    if u[None] > 2.0 * math.pi:
        u[None] -= 2.0 * math.pi
    sint[None] = math.sin(u[None])
    cost[None] = math.cos(u[None])

def main():
    gui = ti.GUI('Realtime Raytracing Demo', res=(XSIZE, YSIZE))
    while gui.running:
        update_scene()
        animate()
        gui.set_image(image.to_numpy())
        gui.show()

if __name__ == '__main__':
    main()
