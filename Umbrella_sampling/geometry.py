from math import degrees, sqrt, acos


# Afstandsfunktion #
# Udregner afstanden mellem to atomer. Skal bruge to <x,y,z>.
def dist(v1, v2):
    dx = v2[0] - v1[0]
    dy = v2[1] - v1[1]
    dz = v2[2] - v1[2]
    result = dx ** 2 + dy ** 2 + dz ** 2
    return sqrt(result)

# Normaliser vektor


def normalize_3dvector(v1, v2, v3):
    v = [v1, v2, v3]
    v1 **= 2
    v2 **= 2
    v3 **= 2

    vector_sum = v1 + v2 + v3
    vector_length = sqrt(vector_sum)

    v[0] = v[0]/vector_length
    v[1] = v[1]/vector_length
    v[2] = v[2]/vector_length
    v_norm = v
    
    return v_norm



# Vinkelmåler #
def angle(v1, v2, v3):  # v2 er den relative orego
    # cos(theta) = dot(u*v)/(dist(v1, v2)*dist(v2*v3))
    # Tæller: dot-product
    ux = v1[0] - v2[0]
    uy = v1[1] - v2[1]
    uz = v1[2] - v2[2]
    vx = v3[0] - v2[0]
    vy = v3[1] - v2[1]
    vz = v3[2] - v2[2]
    dot = ux * vx + uy * vy + uz * vz
    # Nævner: de to længder ganget med hinanden
    dist_u = dist(v1, v2)
    dist_v = dist(v3, v2)
    theta = dot / (dist_u * dist_v)
    result = acos(theta)
    return degrees(result)


# Dihedral / Torsion #
def dihedral(p1, p2, p3, p4):
    # vektorer
    ux = p1[0] - p2[0]
    uy = p1[1] - p2[1]
    uz = p1[2] - p2[2]
    vx = p3[0] - p2[0]
    vy = p3[1] - p2[1]
    vz = p3[2] - p2[2]
    wx = p4[0] - p3[0]
    wy = p4[1] - p3[1]
    wz = p4[2] - p3[2]
    # Normalvektor
    n1 = [uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx]
    n2 = [(-1) * vy * wz - (-1) * vz * wy, (-1) * vz * wx - (-1) * vx * wz, (-1) * vx * wy - (-1) * vy * wx]
    dot = n1[0] * n2[0] + n1[1] * n2[1] + n1[2] * n2[2]  # Prikprodukt
    dist_n1 = sqrt(n1[0] ** 2 + n1[1] ** 2 + n1[2] ** 2)  # Afstand mellem normal vektor 1 og 2
    dist_n2 = sqrt(n2[0] ** 2 + n2[1] ** 2 + n2[2] ** 2)
    phi = dot / (dist_n1 * dist_n2)
    result = acos(phi)
    return degrees(result)


# Hydrogenbinding #
# For denne funktion vil p1 være O (hydrogenacceptor), p2 være det donerede H
# og p3 være det atom p2 er kovalent bundet til (altså H-donor atomet)
def is_hydrogenbond(p1, p2, p3):
    hydrogen_Dist = dist(p1, p2)
    hydrogen_Angle = angle(p1, p2, p3)
    if 1.5 < hydrogen_Dist < 3.5 and 110 < hydrogen_Angle < 180.0:
        value = True
    else:
        value = False
    return value