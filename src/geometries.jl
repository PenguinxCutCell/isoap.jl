# -----------------------------------------------------------------------
# Predefined cell geometries
# -----------------------------------------------------------------------

"""
    cube() -> Polyhedron

Unit cube with vertices in [0,1]³.
Vertex numbering follows the isoap convention:
```
       7/----------/3
       /|         /|
      / |        / |
    8/__|______4/  |
     |  |       |  |
     |  /6------|--/2
     | /        | /
     |/_________|/
     5           1
```
"""
function cube()
    d0, d1 = 0.0, 1.0
    vertp = [
        (d1, d0, d1),  # 1
        (d1, d0, d0),  # 2
        (d1, d1, d0),  # 3
        (d1, d1, d1),  # 4
        (d0, d0, d1),  # 5
        (d0, d0, d0),  # 6
        (d0, d1, d0),  # 7
        (d0, d1, d1),  # 8
    ]
    ipv = [
        [1, 2, 3, 4],
        [2, 1, 5, 6],
        [3, 2, 6, 7],
        [4, 3, 7, 8],
        [1, 4, 8, 5],
        [6, 5, 8, 7],
    ]
    return Polyhedron(ipv, vertp)
end

"""
    tetrahedron() -> Polyhedron

A tetrahedron.
"""
function tetrahedron()
    vertp = [
        (0.0, 0.0, 0.0),
        (0.91, 0.24, 1.0),
        (0.72, 0.16, 0.07),
        (1.0, 1.0, 1.0),
    ]
    ipv = [
        [1, 3, 2],
        [3, 1, 4],
        [2, 3, 4],
        [1, 2, 4],
    ]
    return Polyhedron(ipv, vertp)
end

"""
    dodecahedron() -> Polyhedron

A regular dodecahedron.
"""
function dodecahedron()
    a = 1.0 / sqrt(3.0)
    b = sqrt((3.0 - sqrt(5.0)) / 6.0)
    c = sqrt((3.0 + sqrt(5.0)) / 6.0)

    vertp = [
        ( a,  a,  a),  # 1
        ( a,  a, -a),  # 2
        ( a, -a,  a),  # 3
        ( a, -a, -a),  # 4
        (-a,  a,  a),  # 5
        (-a,  a, -a),  # 6
        (-a, -a,  a),  # 7
        (-a, -a, -a),  # 8
        ( b,  c, 0.0), # 9
        (-b,  c, 0.0), # 10
        ( b, -c, 0.0), # 11
        (-b, -c, 0.0), # 12
        ( c, 0.0,  b), # 13
        ( c, 0.0, -b), # 14
        (-c, 0.0,  b), # 15
        (-c, 0.0, -b), # 16
        (0.0,  b,  c), # 17
        (0.0, -b,  c), # 18
        (0.0,  b, -c), # 19
        (0.0, -b, -c), # 20
    ]

    ipv = [
        [1, 9, 10, 5, 17],
        [1, 17, 18, 3, 13],
        [13, 3, 11, 4, 14],
        [10, 6, 16, 15, 5],
        [4, 20, 19, 2, 14],
        [8, 12, 7, 15, 16],
        [1, 13, 14, 2, 9],
        [9, 2, 19, 6, 10],
        [17, 5, 15, 7, 18],
        [7, 12, 11, 3, 18],
        [8, 16, 6, 19, 20],
        [8, 20, 4, 11, 12],
    ]
    return Polyhedron(ipv, vertp)
end

"""
    icosahedron() -> Polyhedron

A regular icosahedron.
"""
function icosahedron()
    t = (1.0 + sqrt(5.0)) / 2.0
    a = sqrt(1.0 + t^2)

    vertp = [
        ( t/a,  1.0/a,  0.0/a),   # 1
        (-t/a,  1.0/a,  0.0/a),   # 2
        ( t/a, -1.0/a,  0.0/a),   # 3
        (-t/a, -1.0/a,  0.0/a),   # 4
        ( 1.0/a,  0.0/a,  t/a),   # 5
        ( 1.0/a,  0.0/a, -t/a),   # 6
        (-1.0/a,  0.0/a,  t/a),   # 7
        (-1.0/a,  0.0/a, -t/a),   # 8
        ( 0.0/a,  t/a,  1.0/a),   # 9
        ( 0.0/a, -t/a,  1.0/a),   # 10
        ( 0.0/a,  t/a, -1.0/a),   # 11
        ( 0.0/a, -t/a, -1.0/a),   # 12
    ]

    ipv = [
        [1, 9, 5],
        [1, 6, 11],
        [3, 5, 10],
        [3, 12, 6],
        [2, 7, 9],
        [2, 11, 8],
        [4, 10, 7],
        [4, 8, 12],
        [1, 11, 9],
        [2, 9, 11],
        [3, 10, 12],
        [4, 12, 10],
        [5, 3, 1],
        [6, 1, 3],
        [7, 2, 4],
        [8, 4, 2],
        [9, 7, 5],
        [10, 5, 7],
        [11, 6, 8],
        [12, 8, 6],
    ]
    return Polyhedron(ipv, vertp)
end

"""
    complexcell() -> Polyhedron

A complex cell with 18 faces and 32 vertices.
"""
function complexcell()
    vertp = [
        (0.542827704131611521, 0.103161810115955432, 0.036812175979487702),
        (0.600855902479279003, 0.086471258131132003, 0.056436950730949877),
        (0.580277964506970667, 0.163039123851804080, 0.063539333285191152),
        (0.595353203042033874, 0.109843129311963592, 0.046191353733568030),
        (0.543760908872757964, 0.094099097938159251, 0.040511888762036471),
        (0.576819361018633514, 0.166401017225910497, 0.070695345095658793),
        (0.563208132688422847, 0.155406657743345611, 0.108888651199260125),
        (0.539375361775393247, 0.123540417024526422, 0.108841211139604044),
        (0.567496472062291590, 0.128649916760601390, 0.113035361390195654),
        (0.600028673063590978, 0.085143233872209040, 0.059111439263878407),
        (0.526213727798502173, 0.075636612108421930, 0.077638126137438007),
        (0.585556321506447985, 0.087406427365641859, 0.092250691458415329),
        (0.559954651117444135, 0.133408124832590513, 0.115142136143722096),
        (0.524564746594361031, 0.076122995267018601, 0.083047795609934638),
        (0.564747090248686634, 0.145410923059303421, 0.110587477743786994),
        (0.514395820566891815, 0.112406169868730282, 0.074847147170632553),
        (0.561075920034260656, 0.154023426807431363, 0.110265478973935432),
        (0.600086961628965465, 0.099753327716150253, 0.051150970749093465),
        (0.523529840650993616, 0.075758256340824406, 0.081652754048494577),
        (0.580306952541070453, 0.163300101039615730, 0.064077552617180802),
        (0.560784564570740995, 0.154556861208633767, 0.110001263915814051),
        (0.562766593882476296, 0.155209023571326599, 0.109311173701962666),
        (0.535792240519188834, 0.158071087967216972, 0.055846886029540778),
        (0.562461680790533380, 0.155100229495437003, 0.109460804789442229),
        (0.528789432867762255, 0.128164972986322317, 0.044046149723425451),
        (0.535717736632420283, 0.158416280041318303, 0.056565341810804984),
        (0.527966627100669772, 0.148427477131785279, 0.097455495587370988),
        (0.531197945823301376, 0.154922874892573614, 0.054094280101616717),
        (0.519703120392784435, 0.149093981287893751, 0.082950196610910618),
        (0.530093766604715744, 0.156612444958861563, 0.058283874375187901),
        (0.518229286434464531, 0.144572125378080146, 0.084271163208895244),
        (0.520436917533184662, 0.148367611160416940, 0.087663831740394063),
    ]

    ipv = [
        [5, 1, 4, 18, 2],
        [18, 20, 6, 7, 15, 9, 12, 10, 2],
        [10, 11, 5, 2],
        [23, 26, 6, 20, 3],
        [20, 18, 4, 3],
        [4, 1, 25, 28, 23, 3],
        [11, 19, 16, 25, 1, 5],
        [26, 30, 29, 32, 27, 21, 24, 22, 7, 6],
        [22, 15, 7],
        [27, 32, 31, 16, 19, 14, 8],
        [14, 12, 9, 13, 8],
        [13, 17, 21, 27, 8],
        [15, 22, 24, 17, 13, 9],
        [12, 14, 19, 11, 10],
        [31, 29, 30, 28, 25, 16],
        [24, 21, 17],
        [28, 30, 26, 23],
        [31, 32, 29],
    ]
    return Polyhedron(ipv, vertp)
end

"""
    distortedcube() -> Polyhedron

A distorted cube obtained by moving vertices 4 and 7 along the Y axis.
"""
function distortedcube()
    d0, d1 = 0.0, 1.0
    vertp = [
        (d1, d0, d1),       # 1
        (d1, d0, d0),       # 2
        (d1, d1, d0),       # 3
        (d1, 0.25, d1),     # 4
        (d0, d0, d1),       # 5
        (d0, d0, d0),       # 6
        (d0, 0.25, d0),     # 7
        (d0, d1, d1),       # 8
    ]
    ipv = [
        [1, 2, 3, 4],
        [2, 1, 5, 6],
        [3, 2, 6, 7],
        [4, 3, 7],
        [1, 4, 8, 5],
        [6, 5, 8, 7],
        [4, 7, 8],
    ]
    return Polyhedron(ipv, vertp)
end

"""
    pentapyramid() -> Polyhedron

A non-convex pentagonal pyramid.
"""
function pentapyramid()
    vertp = [
        (0.04, 0.77, 0.0),
        (0.0,  0.0,  0.0),
        (0.49, 0.22, 0.0),
        (1.0,  0.13, 0.0),
        (0.16, 1.0,  0.0),
        (0.1,  0.5,  1.0),
    ]
    ipv = [
        [1, 5, 4, 3, 2],
        [1, 2, 6],
        [2, 3, 6],
        [3, 4, 6],
        [4, 5, 6],
        [5, 1, 6],
    ]
    return Polyhedron(ipv, vertp)
end

"""
    cutcube() -> Polyhedron

A non-convex cell obtained by subtracting a pyramid from the cube.
"""
function cutcube()
    d0, d1 = 0.0, 1.0
    vertp = [
        (d1, d0, d1),    # 1
        (d1, d0, d0),    # 2
        (d1, d1, d0),    # 3
        (d1, d1, d1),    # 4
        (d0, d0, d1),    # 5
        (d0, d0, d0),    # 6
        (d0, d1, d0),    # 7
        (d0, d1, d1),    # 8
        (0.5, 0.5, 0.5), # 9
    ]
    ipv = [
        [1, 2, 3, 4],
        [2, 1, 5, 6],
        [3, 2, 6, 7],
        [4, 9, 8],
        [1, 4, 8, 5],
        [6, 5, 8, 7],
        [7, 8, 9],
        [3, 7, 9],
        [4, 3, 9],
    ]
    return Polyhedron(ipv, vertp)
end

"""
    stellatedcube() -> Polyhedron

A stellated cube.
"""
function stellatedcube()
    base = cube()
    d0, d1 = 0.0, 1.0
    a = 1.0

    # Outward face normals for the cube
    xns = [d1, d0, d0, d0, d0, -d1]
    yns = [d0, -d1, d0, d1, d0, d0]
    zns = [d0, d0, -d1, d0, d1, d0]

    new_vertp = collect(base.vertp)
    new_ipv = Vector{Vector{Int}}()
    nts0 = nfaces(base)
    ntp0 = nvertices(base)

    for is in 1:nts0
        xc, yc, zc = 0.0, 0.0, 0.0
        for ip in base.ipv[is]
            xc += base.vertp[ip][1]
            yc += base.vertp[ip][2]
            zc += base.vertp[ip][3]
        end
        n = length(base.ipv[is])
        ntp0 += 1
        push!(new_vertp, (xc / n + a * xns[is], yc / n + a * yns[is], zc / n + a * zns[is]))

        for iv in 1:n
            iv2 = iv == n ? 1 : iv + 1
            push!(new_ipv, [base.ipv[is][iv], base.ipv[is][iv2], ntp0])
        end
    end

    return Polyhedron(new_ipv, new_vertp)
end

"""
    nchexahedron() -> Polyhedron

A non-convex hexahedron.
"""
function nchexahedron()
    d0, d1 = 0.0, 1.0
    vertp = [
        (0.5, 0.75, d1),  # 1
        (0.5, 0.75, d0),  # 2
        (d1,  d1,   d0),  # 3
        (d1,  d1,   d1),  # 4
        (d0,  d0,   d1),  # 5
        (d0,  d0,   d0),  # 6
        (d0,  d1,   d0),  # 7
        (d0,  d1,   d1),  # 8
    ]
    ipv = [
        [1, 2, 3, 4],
        [2, 1, 5, 6],
        [3, 2, 6, 7],
        [4, 3, 7, 8],
        [1, 4, 8, 5],
        [6, 5, 8, 7],
    ]
    return Polyhedron(ipv, vertp)
end

"""
    _stellate(base::Polyhedron, a::Float64) -> Polyhedron

Helper: create a stellated version of any polyhedron by computing face normals
and pushing out new apex vertices.
"""
function _stellate(base::Polyhedron, a::Float64)
    new_vertp = collect(base.vertp)
    new_ipv = Vector{Vector{Int}}()
    nts = nfaces(base)
    ntp0 = nvertices(base)

    for is in 1:nts
        # Compute face normal from first three vertices
        ip1 = base.ipv[is][1]
        ip2 = base.ipv[is][2]
        ip3 = base.ipv[is][3]
        xv1 = base.vertp[ip2][1] - base.vertp[ip1][1]
        yv1 = base.vertp[ip2][2] - base.vertp[ip1][2]
        zv1 = base.vertp[ip2][3] - base.vertp[ip1][3]
        xv2 = base.vertp[ip3][1] - base.vertp[ip2][1]
        yv2 = base.vertp[ip3][2] - base.vertp[ip2][2]
        zv2 = base.vertp[ip3][3] - base.vertp[ip2][3]
        xn = yv1 * zv2 - zv1 * yv2
        yn = zv1 * xv2 - xv1 * zv2
        zn = xv1 * yv2 - yv1 * xv2
        dmod = sqrt(xn^2 + yn^2 + zn^2)
        xn /= dmod
        yn /= dmod
        zn /= dmod

        # Face centroid
        xc, yc, zc = 0.0, 0.0, 0.0
        n = length(base.ipv[is])
        for ip in base.ipv[is]
            xc += base.vertp[ip][1]
            yc += base.vertp[ip][2]
            zc += base.vertp[ip][3]
        end
        ntp0 += 1
        push!(new_vertp, (xc / n + a * xn, yc / n + a * yn, zc / n + a * zn))

        for iv in 1:n
            iv2 = iv == n ? 1 : iv + 1
            push!(new_ipv, [base.ipv[is][iv], base.ipv[is][iv2], ntp0])
        end
    end

    return Polyhedron(new_ipv, new_vertp)
end

"""
    stellateddodecahedron() -> Polyhedron

A stellated dodecahedron.
"""
stellateddodecahedron() = _stellate(dodecahedron(), 0.5)

"""
    stellatedicosahedron() -> Polyhedron

A stellated icosahedron.
"""
stellatedicosahedron() = _stellate(icosahedron(), 1.0)

"""
    hollowedcube() -> Polyhedron

A unit cube with a half-length cubic hollow in its center.
"""
function hollowedcube()
    d0, d1 = 0.0, 1.0
    d14 = 0.25

    # Outer cube vertices (1-8)
    outer = [
        (d1, d0, d1),       # 1
        (d1, d0, d0),       # 2
        (d1, d1, d0),       # 3
        (d1, d1, d1),       # 4
        (d0, d0, d1),       # 5
        (d0, d0, d0),       # 6
        (d0, d1, d0),       # 7
        (d0, d1, d1),       # 8
    ]
    # Inner cube vertices (9-16)
    inner = [
        (d1-d14, d0+d14, d1-d14), # 9
        (d1-d14, d0+d14, d0+d14), # 10
        (d1-d14, d1-d14, d0+d14), # 11
        (d1-d14, d1-d14, d1-d14), # 12
        (d0+d14, d0+d14, d1-d14), # 13
        (d0+d14, d0+d14, d0+d14), # 14
        (d0+d14, d1-d14, d0+d14), # 15
        (d0+d14, d1-d14, d1-d14), # 16
    ]
    vertp = vcat(outer, inner)

    # Outer faces (same winding as cube)
    outer_ipv = [
        [1, 2, 3, 4],
        [2, 1, 5, 6],
        [3, 2, 6, 7],
        [4, 3, 7, 8],
        [1, 4, 8, 5],
        [6, 5, 8, 7],
    ]
    # Inner faces (reversed winding)
    inner_ipv = [
        [12, 11, 10, 9],
        [14, 13, 9, 10],
        [15, 14, 10, 11],
        [16, 15, 11, 12],
        [13, 16, 12, 9],
        [15, 16, 13, 14],
    ]

    return Polyhedron(vcat(outer_ipv, inner_ipv), vertp)
end

"""
    drilledcube() -> Polyhedron

A non-simply connected polyhedron: a cube with a rectangular hole drilled through it.
"""
function drilledcube()
    d0, d1 = 0.0, 1.0
    d02, d12 = 0.25, 0.75

    vertp = [
        # Outer cube vertices (1-8)
        (d1, d0, d1),
        (d1, d0, d0),
        (d1, d1, d0),
        (d1, d1, d1),
        (d0, d0, d1),
        (d0, d0, d0),
        (d0, d1, d0),
        (d0, d1, d1),
        # Inner hole vertices (9-16)
        (d12, d0, d12),
        (d12, d0, d02),
        (d12, d1, d02),
        (d12, d1, d12),
        (d02, d0, d12),
        (d02, d0, d02),
        (d02, d1, d02),
        (d02, d1, d12),
    ]

    ipv = [
        # Outer cube faces
        [1, 2, 3, 4],
        [2, 1, 5, 6],
        [3, 2, 6, 7],
        [4, 3, 7, 8],
        [1, 4, 8, 5],
        [6, 5, 8, 7],
        # Inner hole faces (reversed winding)
        [12, 11, 10, 9],
        [14, 13, 9, 10],
        [15, 14, 10, 11],
        [16, 15, 11, 12],
        [13, 16, 12, 9],
        [15, 16, 13, 14],
    ]
    return Polyhedron(ipv, vertp)
end

"""
    zigzagcell(nzigs::Int=5) -> Polyhedron

A zig-zag prism.
"""
function zigzagcell(nzigs::Int=5)
    doff = 0.1
    d0, d1 = 0.0, 1.0
    ntp = 4 * nzigs
    nts = 2 * nzigs + 2

    vertp = Vector{NTuple{3,Float64}}(undef, ntp)
    for iv in 1:nzigs
        vertp[iv] = (d1 * (iv - 1), doff + d1 * mod(iv - 1, 2), d0)
        vertp[iv + nzigs] = (d1 * (nzigs - (iv - 1) - 1), d1 * mod(nzigs - (iv - 1) - 1, 2), d0)
        vertp[iv + 2*nzigs] = (d1 * (iv - 1), doff + d1 * mod(iv - 1, 2), d1)
        vertp[iv + 3*nzigs] = (d1 * (nzigs - (iv - 1) - 1), d1 * mod(nzigs - (iv - 1) - 1, 2), d1)
    end

    ipv = Vector{Vector{Int}}(undef, nts)

    # Face 1: bottom z=0 polygon
    ipv[1] = collect(1:2*nzigs)
    # Face 2: top z=1 polygon (reversed)
    ipv[2] = collect(4*nzigs:-1:2*nzigs+1)

    # Face 3: side face at start
    iv = 1
    ipv[3] = [iv, 2*nzigs - (iv - 1), 4*nzigs - (iv - 1), iv + 2*nzigs]
    # Face 4: side face at end
    iv = nzigs
    ipv[4] = [iv + 2*nzigs, 4*nzigs - (iv - 1), 2*nzigs - (iv - 1), iv]

    # Side faces
    for iv in 1:nzigs-1
        ipv[iv + 4] = [iv, iv + 2*nzigs, iv + 2*nzigs + 1, iv + 1]
        ipv[iv + 4 + (nzigs - 1)] = [
            2*nzigs - (iv - 1),
            2*nzigs - (iv - 1) - 1,
            2*nzigs - (iv - 1) - 1 + 2*nzigs,
            2*nzigs - (iv - 1) + 2*nzigs,
        ]
    end

    return Polyhedron(ipv, vertp)
end

# -----------------------------------------------------------------------
# Implicit functions for defining interface shapes
# -----------------------------------------------------------------------

"""
    func_sphere(x, y, z) -> Float64

Sphere with radius 0.325 centered at (0.5, 0.5, 0.5).
Exact volume = (4/3)π(0.325)³.
"""
func_sphere(x, y, z) = -((x - 0.5)^2 + (y - 0.5)^2 + (z - 0.5)^2 - 0.325^2)

"""
    func_torus(x, y, z) -> Float64

Torus with major radius 0.2, minor radius 0.1, centered at (0.5, 0.5, 0.5).
"""
func_torus(x, y, z) = 0.1^2 - (0.2 - sqrt((x - 0.5)^2 + (z - 0.5)^2))^2 - (y - 0.5)^2

"""
    func_orthocircle(x, y, z) -> Float64

Orthocircle surface centered at (1.25, 1.25, 1.25).
"""
function func_orthocircle(x, y, z)
    x = x - 1.25
    y = y - 1.25
    z = z - 1.25
    return ((x^2 + y^2 - 1)^2 + z^2) * ((y^2 + z^2 - 1)^2 + x^2) *
           ((z^2 + x^2 - 1)^2 + y^2) - 0.075^2 * (1 + 3 * (x^2 + y^2 + z^2))
end
