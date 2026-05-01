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
# plt.show();

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
# plt.show()

# Python bn128 library

from py_ecc.bn128 import G1, multiply, add, eq, neg

print("THE TESTING BEGINS")

print(G1)
print(add(G1,G1))
print(multiply(G1,2))

# # 5G + 6G = 11G
assert eq(add(multiply(G1, 5), multiply(G1, 6)), multiply(G1, 11))

# numbers are large to make it difficult for attackers to brute force it 

import matplotlib.pyplot as plt
from py_ecc.bn128 import G1, multiply, neg

import numpy as np
xs = []
ys = []
for i in range(1,1000):
    xs.append(i)
    ys.append(int(multiply(G1, i)[1]))
    xs.append(i)
    ys.append(int(neg(multiply(G1, i))[1]))
plt.scatter(xs, ys, marker='.')
# plt.show()

# confirming if adding curve order to some point you get the point back and if it is the field modulus , it would give you a different point 
from py_ecc.bn128 import field_modulus


a= 7
print(eq(multiply(G1,a),multiply(G1,a+ curve_order)))
print(eq(multiply(G1, a), multiply(G1, a + field_modulus)))

# Encoding rational numbers 
five_over_two = (5 * pow(2, -1, curve_order)) % curve_order
one_half = pow(2, -1, curve_order)

# Essentially 5/2 = 2.5# 2.5 + 0.5 = 3
# but we are doing this in a finite field
assert eq(add(multiply(G1, five_over_two), multiply(G1, one_half)), multiply(G1, 3))


x = 6
y = 11
z = 17

# (x+y)+z = x+(y+z)
lhs = add(add(multiply(G1, x), multiply(G1, y)), multiply(G1, z))

rhs = add(multiply(G1, x), add(multiply(G1, y), multiply(G1, z)))

print(eq(lhs, rhs)) 


from py_ecc.bn128 import  neg, is_inf, Z1

# pick a field element
x = 4444444
# generate the point
p = multiply(G1, x)

# invert
p_inv = neg(p)
print("The inverse is :", p_inv)

# every element added to its inverse produces the identity element 
assert is_inf(add(p, p_inv))

# Z1 is just None, which is the point at infinity
assert Z1 is None

# special case: the inverse of the identity is itself
assert eq(neg(Z1), Z1)

print("testing the fact that elliptic curve points over real number , (x,y) the inverse is usually (x,INV(y) so INV(y)")

field_modulus = 21888242871839275222246405745257275088696311157297823662689037894645226208583
for i in range(1, 4):
    point = multiply(G1, i)
    print(point)
    print(neg(point))
    print('----')
    
    # x values are the same
    assert int(point[0]) == int(neg(point)[0])
    
    # y values are inverses of each other, we are adding y values
    # not ec points
    assert int(point[1]) + int(neg(point)[1]) == field_modulus


# Prover
secret_x = 5
secret_y = 10

x = multiply(G1, 5)
y = multiply(G1, 10)

proof = (x, y, 15)

# verifier
if multiply(G1, proof[2]) == add(proof[0], proof[1]):
    print("statement is true")
else:
    print("statement is false")
