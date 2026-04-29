# import galois
from libnum import has_sqrtmod_prime_power, sqrtmod_prime_power
p =11
print(pow(3,-1,p))

# assert (3*2) % p == 1

# Creaated by myself
def mul_inverse(a,p):
    inverse= pow(a,-1,p)
    return print(inverse)

mul_inverse(8,17)

# Turning fraction to field element

def fraction_to_field_element (num,den,p):
    inverse=  pow(den,-1,p)
    return print((num * inverse) % p )

fraction_to_field_element(1,2,7)
fraction_to_field_element(1,3,7)
fraction_to_field_element(5,6,7)

# Turning fraction to field element using galios
# def fraction_to_field_element_2 (num,den,p):
#     GFP = galois.GF(p)
#     inverse=  1/GFP(den)
#     return print((num * inverse) % p )

# fraction_to_field_element_2(1,2,7)
# fraction_to_field_element_2(1,3,7)
# fraction_to_field_element_2(5,6,7)

# Numbers that have roots in the field
numbers_with_roots = set()
prime= 11
for i in range(0,prime):
    numbers_with_roots.add((i*i) % prime)

print(numbers_with_roots)    

# Using libnum to verify sqr mod and list 
print(list(sqrtmod_prime_power(5,prime,1)))



def mod_sqrt(x,p):
     assert(p-3)%4 ==0
     exponent = (p+1)//4
     first_root = pow(x,exponent,p)
     second_root = pow(first_root,-1,p)
     return print(first_root,second_root)
     
mod_sqrt(7,19)


# Exercise: Convert the two equations to their finite field representation and see they are the same.
# x + 2 * y === y = 1/2 -x/2
# 4 * x + 8 * y === 1  y=1/8-x/2 it is obvius that the gradient is equal 

#  Primitive roots are elements that can generate all elements in a group( Cyclic group) under multiplication in this instance
# from galois import GF
# GF7 = GF(7)
# print(GF7.primitive_elements)
# # [3, 5]


