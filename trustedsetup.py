from py_ecc.bn128 import G1, multiply, add
from functools import reduce

def inner_product(points, coeffs):
    return reduce(add, map(multiply, points, coeffs))

## Trusted Setup
tau = 88
degree = 3

# tau^3, tau^2, tau, 1
srs = [multiply(G1, tau**i) for i in range(degree,-1,-1)]

print(srs,"srs")

## Evaluate
# p(x) = 4x^2 + 7x + 8
coeffs = [0, 4, 7, 8]

poly_at_tau = inner_product(srs, coeffs)
print(poly_at_tau, "poly_at_tau")
