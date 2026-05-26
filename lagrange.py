from scipy.interpolate import lagrange
x_values = [1, 2, 3, 4]
y_values = [4, 8, 2, 1]

print(lagrange(x_values, y_values))

x_values = [1, 4, 7, 4]
y_values = [4, 8, 2, 1]

print(lagrange(x_values, y_values))

import galois
import numpy as np
GF17 = galois.GF(17)

xs = GF17(np.array([1,2,3,4]))
ys = GF17(np.array([4,8,2,1]))

p = galois.lagrange_poly(xs, ys)


assert p(1) == GF17(4)
assert p(2) == GF17(8)
assert p(3) == GF17(2)
assert p(4) == GF17(1)

print(p)
import galois
import numpy as np

# Checking if two vectors are equal using SZL 
p = 103
GFP= galois.GF(p)

x= GFP(np.array([1,2,3,4,5]))

def poly(v):
    return galois.lagrange_poly(x,v)

v1=GFP(np.array([3,6,9,12,15]))
v2=GFP(np.array([3,6,9,12,15]))

pl = poly(v1)
p2= poly(v2)
import random

u = random.randint(0,p)

lhs = pl(u)
rhs= p2(u)

print(lhs==rhs)