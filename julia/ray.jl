# Translated from the Ray Tracer Language Comparison
# http://www.ffconsultancy.com/languages/ray_tracer/benchmark.html

import Printf
import StaticArrays
import LinearAlgebra: dot

const Vec = StaticArrays.SVector{3,Float64}
unitise(v::Vec) = v/sqrt(dot(v,v))

struct Hit
    l :: Float64
    d :: Vec
end

struct Ray
    orig :: Vec
    dir  :: Vec
end

struct Sphere
    center :: Vec
    radius :: Float64
end

struct Group
    bound :: Sphere
    child :: Vector{Union{Sphere,Group}}
end

function ray_sphere(r::Ray,s::Sphere)::Float64
    v = s.center - r.orig
    b = dot(v,r.dir)
    disc = b*b - dot(v,v) + s.radius*s.radius
    if disc < 0. return Inf end
    d = sqrt(disc)
    t2 = b + d
    if t2 < 0. return Inf end
    t1 = b - d
    return t1 > 0. ? t1 : t2
end

function inter(s::Sphere,hit::Hit,ray::Ray)::Hit
    λ = ray_sphere(ray,s)
    λ < hit.l ? Hit(λ,unitise(ray.orig + λ*ray.dir - s.center)) : hit
end

function inter(g::Group, hit::Hit,ray::Ray)::Hit
    l = ray_sphere(ray,g.bound)
    if (l < hit.l)
        for scene in g.child
            hit = inter(scene,hit,ray)
        end
    end
    hit
end

const zero3 = Vec(0.,0.,0.)
const hit0 = Hit(Inf,zero3)
const δ = sqrt(eps())


function raytrace(light :: Vec, ray :: Ray, scene::Union{Sphere,Group}) :: Float64
    hit = inter(scene,hit0,ray)
    if hit.l == Inf return 0. end
    g = dot(hit.d,light)
    if g >= 0.  return 0. end
    p = ray.orig + hit.l*ray.dir + δ*hit.d
    hit2 = hit0
    hit2 = inter(scene,hit2,Ray(p,-light))
    return hit2.l < Inf ? 0. : -g
end

function create(level,c::Vec,r::Float64)
    s = Sphere(c,r)
    if level == 1 return s else
        a = 3.0*r/sqrt(12.)
        aux(x, z) = create(level-1, c + Vec(x, a, z), r/2.)
        Group(Sphere(c, 3.0*r), [s, aux(-a, -a), aux(a, -a), aux(-a, a), aux(a, a)])
    end
end

const ss = 4

aux(x,n,d) = x - n/2 + d/ss
function antialiasedraytrace(x,y,n,light,scene,source)
    g = 0.
    for d in 0:ss^2-1
        dir = unitise(Vec(aux(x,n,d%ss), aux(y,n,d/ss), convert(Float64, n)))
        g += raytrace(light,Ray(source,dir),scene)
    end
    g = round(255.0*g/(ss*ss))
    isnan(g) ? 0x00 : convert(UInt8, g)
end

image(n,light,scene,source) = UInt8[antialiasedraytrace(x,y,n,light,scene,source) for x in 0:n-1 for y in n-1:-1:0]

function main(;n=512,level=9,light=unitise(Vec(-1., -3., 2.)),source=zero3,io=stdout)
    scene = create(level, Vec(0., -1., 4.), 1.)
    Printf.@printf(io,"P5\n%d %d\n255", n, n)
    write(io, image(n,light,scene, source))
end
main()


# For repl testing:
# import ImageView
# imagematrix(n,light,scene,source) = reshape(image(n,light,scene,source),(n,n))
# function viewRaytracedImage(;n=512,level=9,light=unitise(Vec(-1., -3., 2.)),source=zero3)
#     scene = create(level, Vec(0., -1., 4.), 1.)
#     ImageView.imshow(imagematrix(n,light,scene,source))
# end
# viewRaytracedImage()