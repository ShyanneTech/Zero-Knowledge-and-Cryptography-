import galois
import numpy as np
import random
from functools import reduce
p = 17
GFP = galois.GF(17)
x_values = GFP(np.array([1,2,3]))

def function_1(x):
    return x*x
def function_2(x):
    return x*x*x

v1= GFP(function_1(x_values))
print(v1)
v2= GFP(function_2(x_values))
print(v2)

def L(v):
    return galois.lagrange_poly(x_values,v)

lamda = random.randint(0,p)

print(L(v1))
print(L(v2))
print(L(v1+v2))
print(L(v1) + L(v2) == (L(v1+ v2)))
print(lamda, "The scalar field element")
print(lamda * L(v1) )
print(L(lamda * v1))

print("Another way for confirming scalar multiplication")

p = 17
GF = galois.GF(p)

xs = GF(np.array([1,2,3]))

# arbitrary vector
v =  GF(np.array([4,8,2]))

# arbitrary constant
lambda_ =  GF(15)

def L(v):
    return galois.lagrange_poly(xs, v)

print(L(v))
assert L(lambda_ * v) == lambda_ * L(v)

# Checking if Av1 = Bv2 ... A and B are matrices which are broken down to vectors 

p = 17 
GF = galois.GF(p)

x= GF(np.array([1,2]))
yp1 = GF(np.array([6,4]))
yp2 = GF(np.array([3,7]))
yq1 = GF(np.array([3,12]))
yq2 = GF(np.array([9,6]))

def L(v) :
    return galois.lagrange_poly(x,v)

p1= L(yp1)
p2= L(yp2)
q1= L(yq1)
q2= L(yq2)

# Successfuly converted matrix A and B to polynomials 
print('Successfuly converted matrix A and B to polynomials ')
print(p1)
print(p2)
print(q1)
print(q2)

u = random.randint(0, p)
tau = GF(u) # a random point

left_hand_side = p1(tau) * GF(2) + p2(tau) * GF(4)
right_hand_side = q1(tau) * GF(2) + q2(tau) * GF(2)

print(left_hand_side == right_hand_side) 


print('R1CS TO QAP over a Finite Field in Python')

# vector a = [1,out,x,y,v1,v2,v3]
L = np.array([
    [0, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0],
    [0, 0, 0, -5, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 1],
])

R = np.array([
    [0, 0, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0],
    [0, 0, 0, 1, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0],
])

O = np.array([
    [0, 0, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 0, 1, 0],
    [0, 0, 0, 0, 0, 0, 1],
    [0, 1, 0, 0, 0, -1, 0],
])

# to confirm R1CS was built well 
x = 4
y = -2
v1 = x * x
v2 = v1 * v1          
v3 = -5*y * y
z = v3*v1 + v2    
# witness
a = np.array([1, z, x, y, v1, v2, v3])

assert all(np.equal(np.matmul(L, a) * np.matmul(R, a), np.matmul(O, a))), "not equal"

# Procceds to convert it to field array 

GF = galois.GF(79)

L = (L + 79) % 79 
R = (R + 79) % 79
O = (O + 79) % 79

# for a negative number to field x= -2 , x= GF(-2+ 79)


L_galois = GF(L)
R_galois = GF(R)
O_galois = GF(O)

x = GF(4)
y = GF(-2 + 79) # we are using 79 as the field size, so 79 - 2 is -2
v1 = x * x  # not wrapped in GF because x is already a field element and multiplying two field element gives a field element  
v2 = v1 * v1         # x^4
v3 = GF(-5 + 79)*y * y
out = v3*v1 + v2  

witness = GF(np.array([1, out, x, y, v1, v2, v3]))

assert all(np.equal(np.matmul(L_galois, witness) * np.matmul(R_galois, witness), np.matmul(O_galois, witness))), "not equal"

def interpolate_column(col):
    xs = GF(np.array([1,2,3,4]))
    return galois.lagrange_poly(xs, col)

# axis 0 is the columns.
# apply_along_axis is the same as doing a for loop over the columns and collecting the results in an array
U_polys = np.apply_along_axis(interpolate_column, 0, L_galois)
V_polys = np.apply_along_axis(interpolate_column, 0, R_galois)
W_polys = np.apply_along_axis(interpolate_column, 0, O_galois)

print(U_polys)
print(V_polys)
print(W_polys)

# Now executing the QAP formula 

def inner_product_polynomials_with_witness(polys, witness):
    mul_ = lambda x, y: x * y
    sum_ = lambda x, y: x + y
    return reduce(sum_, map(mul_, polys, witness))

term_1 = inner_product_polynomials_with_witness(U_polys, witness)

term_2 = inner_product_polynomials_with_witness(V_polys, witness)

term_3 = inner_product_polynomials_with_witness(W_polys, witness)

print(term_1, "term 1")
print(term_2, "term 2")
print(term_3, "term 3")

# t = (x - 1)(x - 2)(x - 3)(x - 4)
t = galois.Poly([1, 78], field = GF) * galois.Poly([1, 77], field = GF) * galois.Poly([1, 76], field = GF) * galois.Poly([1, 75], field = GF)

h = (term_1 * term_2 - term_3) // t
print("h(x):",hex)

assert term_1 * term_2 == term_3 + h * t, "division has a remainder"


