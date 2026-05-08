# 41 X103 = 4223 to R1CS

import numpy as np

O= np.matrix([[0,1,0,0]])
L= np.matrix([[0,0,1,0]])
R= np.matrix([[0,0,0,1]])

# Witness vector
a = np.array([1,4223,41,103])
 


# element wise for each constraint this works because it is just a constraint 
result = np.matmul(O,a) == np.matmul(L,a) * np.matmul(R,a) 
print(result)
print(result.all())

# assert shows the message if the condition is false
assert result.all(), "result contains an inequality"


# 20 X200 = 4000
# a=[1,4000,20,200]
# The matrices
O = np.matrix([[0,1,0,0]])
L = np.matrix([[0,0,1,0]])
R = np.matrix([[0,0,0,1]])
a = np.array([1,4000,20,200])
result = np.matmul(O,a) == np.matmul(L,a) * np.matmul(R,a)
print(result)

# EXAMPLE 2
#  r = x*y*z*u v1= xy v2=zu r = v1vu 
# order of equation  v1= xy v2=zu r = v1vu using random values for the variables to form my witness vector 

L = np.matrix([[0,0,1,0,0,0,0,0],[0,0,0,0,1,0,0,0],[0,0,0,0,0,0,1,0]])
R = np.matrix([[0,0,0,1,0,0,0,0],[0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,1]])
O = np.matrix([[0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,1],[0,1,0,0,0,0,0,0]])

# random values for x, y, z, and u
import random
x = random.randint(1,1000)
y = random.randint(1,1000)
z = random.randint(1,1000)
u = random.randint(1,1000)

# The algebraic circuits 
r = x * y * z * u
v1 = x*y
v2 = z*u

# The witness vector 
a = np.array([1,r,x,y,z,u,v1,v2])

result = np.matmul(O,a) == np.multiply(np.matmul(L,a) , np.matmul(R,a))

print('Major correction done here ')

if(result.all()):
    print(r,x,y,z,u,v1,v2)

print(result)
assert result.all(), "system contains an inequality"


# Example 3 z= x* y + 2  final equation -2 + z = x* y 

# Define the matrices
L = np.matrix([[0,0,1,0]])
R = np.matrix([[0,0,0,1]])
O = np.matrix([[-2,1,0,0]])

# pick random values to test the equation
x = random.randint(1,1000)
y = random.randint(1,1000)
z = x * y + 2 # witness vector
a = np.array([1, z, x, y])

result = np.matmul(O,a) == np.multiply(np.matmul(L,a) , np.matmul(R,a))
result2 = O.dot(a) == np.multiply(np.matmul(L, a), R.dot(a))

print('CHECKING IF IT IS DIFFERENT')
print(result)
print(result2)
assert result.all(), "result contains an inequality"


# Define the matrices
L = np.array([[0,0,3,0,0,0],
               [0,0,0,0,1,0],
               [0,0,5,0,0,0]])

R = np.array([[0,0,1,0,0,0],
               [0,0,0,1,0,0],
               [0,0,0,1,0,0]])

O = np.array([[0,0,0,0,1,0],
               [0,0,0,0,0,1],
               [-3,1,1,2,0,-1]])

# pick random values for x and y
x = random.randint(1,1000)
y = random.randint(1,1000)

# this is our orignal formula
out = 3 * x * x * y + 5 * x * y - x - 2 * y + 3 # the witness vector with the intermediate variables inside
v1 = 3*x*x
v2 = v1 * y
w = np.array([1, out, x, y, v1, v2])
print('Complex equation')
result = O.dot(w) == np.multiply(L.dot(w),R.dot(w))
result2 = np.matmul(O,w) == np.multiply(np.matmul(L,w) , np.matmul(R,w))
print(result)  
print(result2) 

p = 21888242871839275222246405745257275088548364400416034343698204186575808495617

# 1 - 2 = -1
print((1 - 2) % p)

# 21888242871839275222246405745257275088548364400416034343698204186575808495616
