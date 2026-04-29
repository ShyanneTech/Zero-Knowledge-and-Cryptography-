# Confirming commutatitivity of points in elliptic circle under addition 
# var('y_t', 'y_u', 'x_t', 'x_u')
# lambda_p = (y_u - y_t)/(x_u - x_t)
# x_p = lambda_p^2 - x_t - x_u
# y_p = (lambda_p*(x_t - x_p) - y_t)

# lambda_q = (y_t - y_u)/(x_t - x_u)
# x_q = lambda_q^2 - x_u - x_t
# y_q = (lambda_q*(x_u - x_q) - y_u)

# y2=x3−x+1
def elliptic_circle_commutativity(x_t,y_t,x_u,y_u):
    # P = T + U ,  Q = U+T
    lambda_p = (y_u-y_t)/(x_u-x_t)
    x_p = lambda_p**2 - x_t - x_u
    y_p = (lambda_p*(x_t - x_p) - y_t)

    lambda_q = (y_t - y_u)/(x_t - x_u)
    x_q = lambda_q**2 - x_u - x_t
    y_q = (lambda_q*(x_u - x_q) - y_u)

    return (x_p == x_q, y_p== y_q)

print(elliptic_circle_commutativity(0,1,1,1))


# BN128 FORMULA USED BY ETHERUM

from py_ecc.bn128 import curve_order
# 21888242871839275222246405745257275088548364400416034343698204186575808495617
print(curve_order)


# Modular square root 
from libnum import has_sqrtmod_prime_power, sqrtmod_prime_power

# the functions take arguments# has_sqrtmod_prime_power(n, field_mod, k), where n**k,
# but we aren't interested in powers in modular fields, so we set k = 1
# check if sqrt(8) mod 11 exists
print(has_sqrtmod_prime_power(8, 11, 1))
print( has_sqrtmod_prime_power(7,11,1))
# False

# check if sqrt(5) mod 11 exists
print(has_sqrtmod_prime_power(5, 11, 1))
print(has_sqrtmod_prime_power(3, 11, 1))
# True

# compute sqrt(5) mod 11
print(list(sqrtmod_prime_power(5, 11, 1)))
print(list(sqrtmod_prime_power(3, 11, 1)))
# [4, 7]

assert (4 ** 2) % 11 == 5
assert (7 ** 2) % 11 == 5
assert (5 ** 2) % 11 == 3
assert (6 ** 2) % 11 == 3

# we expect 4 and 7 to be additive inverses of each other, because in "regular" math, the two solutions to a square root are sqrt and -sqrt
assert (4 + 7) % 11 == 0
assert (5 + 6) % 11 == 0


# Getting elliptic points satisfy y^2 = x^3 + 3 (mod 11)

# Getting myself familiar with matplotlib.pyplot package , it is used to draw graphs and charts etc in Python 

# import matplotlib.pyplot as plt 

# x= [1,2,3,4]
# y = [2,4,6,8]

# plt.plot(x,y)
# plt.title("SHARON THE ZK DEVELOPER ")
# plt.show()

# Now to solve the above equation in whatever mod 
import libnum
import matplotlib.pyplot as plt

def generate_points(mod):
    xs = []
    ys = []
    def y_squared(x):
        return (x**3 + 3) % mod

    for x in range(0, mod):
        if libnum.has_sqrtmod_prime_power(y_squared(x), mod, 1):
            square_roots = libnum.sqrtmod_prime_power(y_squared(x), mod, 1)

            # we might have two solutions
            for sr in square_roots:
                ys.append(sr)
                xs.append(x)
    return xs, ys


xs, ys = generate_points(11)
fig, (ax1) = plt.subplots(1, 1);
fig.suptitle('y^2 = x^3 + 3 (mod p)');
fig.set_size_inches(6, 6);
ax1.set_xticks(range(0,23));
ax1.set_yticks(range(0,23));
plt.grid()
plt.scatter(xs, ys)
plt.show();

# Doubling just means adding same point while additon is general addition between points 


# This doubling function would be referenced in the general function because we need to be sure it satisfies whatever elliptical curve equation is given y^2 = x^3 + 3 
def doubling_ponits(x,y,a,p):
    lambd = (((3*x**2 + a)% p) * pow(2*y,-1,p)) %p
    newx= (lambd**2 - 2 * x) % p
    newy= (lambd* x - lambd*newx - y) % p
    return (newx,newy)

def add_points(xp,yp,xq,yq,p, a=0):
    # None here represent the point of infinity which is the identity element
    if (xp == yp == None):
        return xq,yq
    if (xq == yq == None):
        return xp,yp
    
    # Whatever point given must satisfy the equation 
    assert (xq**3 + 3) % p == (yq ** 2) % p, "q not on curve"
    assert (xp**3 + 3) % p == (yp ** 2) % p, "p not on curve"

    if xq == xp and yq == yp:
        return doubling_ponits(xq, yq, a, p)
    elif xq == xp:
        return None, None

    lambd = ((yq - yp) * pow((xq - xp), -1, p) ) % p
    xr = (lambd**2 - xp - xq) % p
    yr = (lambd*(xp - xr) - yp) % p
    return xr, yr

print(add_points(1,2,1,2,11))    

# Every elliptic curve point in a cyclic group has a “number”

# for our purposes, (4, 10) is the generator point G
next_x, next_y = 4, 10
print(0, 4, 10)
points = [(next_x, next_y)]
for i in range(1, 13):
    # repeatedly add G to the next point to generate all the elements
    next_x, next_y = add_points(next_x, next_y, 4, 10, 11)
    print(i, next_x, next_y)
    points.append((next_x, next_y))

xs11, ys11 = generate_points(11)

fig, (ax1) = plt.subplots(1, 1);
fig.suptitle('y^2 = x^3 + 3 (mod 11)');
fig.set_size_inches(13, 6);

ax1.set_title("modulo 11")
ax1.scatter(xs11, ys11, marker='o');
ax1.set_xticks(range(0,11));
ax1.set_yticks(range(0,11));
ax1.grid()

for i in range(0, 11):
    plt.annotate(str(i+1), (points[i][0] + 0.1, points[i][1]), color="red");
plt.show()