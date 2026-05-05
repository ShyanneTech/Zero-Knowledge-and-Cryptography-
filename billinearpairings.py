# from py_ecc.bn128 import G1, G2, pairing, add, multiply, eq,  curve_order

# print(G1)
# print(G2)

# x = 10 # chosen randomly
# assert eq(multiply(G2, x + curve_order), multiply(G2, x))
# assert eq(multiply(G1, x + curve_order), multiply(G1, x))

# # G2 similar behaviour with cycic group i.e G1
# print(eq(add(G1, G1), multiply(G1, 2)))

# print(eq(add(G2, G2), multiply(G2, 2)))

# # add(G1, G2)  Different groups .. A NO NO

# # BILLINEAR PAIRING IN PYTHON 
# P = multiply(G1,3)
# Q= multiply(G2,8)
# R = multiply(G1,24)

# # assert eq(pairing(Q, P), pairing(G2, R))
# # print(eq(pairing(P,Q), pairing(R,G2)))


# P_1 = multiply(G1, 3)
# P_2 = multiply(G2, 8)

# Q_1 = multiply(G1, 6)
# Q_2 = multiply(G2, 4)

# # lhs = pairing(P_2, P_1) 
# # rhs = pairing(Q_2, Q_1)

# # assert lhs == rhs
# # assert eq(pairing( P_2,P_1), pairing(Q_2, Q_1))

# # How to address elements in GT
# # 2 * 3 = 6
# P_1 = multiply(G1, 2)
# P_2 = multiply(G2, 3)

# # 4 * 5 = 20
# Q_1 = multiply(G1, 4)
# Q_2 = multiply(G2, 5)

# # 13 * 2 = 26
# R_1 = multiply(G1, 13)
# R_2 = multiply(G2, 2)


# lhs = pairing(P_2, P_1) * pairing(Q_2, Q_1)
# rhs = pairing(R_2, R_1)

# assert lhs == rhs
# print("pairing equation passed")

# # b ^ {2 * 3} * b ^ {4 * 5} = b ^ {13 * 2}
# # b ^ 6 * b ^ 20 = b ^ 26

# # assert eq(pairing(P_2, P_1) * pairing(Q_2, Q_1), pairing(R_2, R_1))





# import sys

# sys.setrecursionlimit(29793968203157093290)

# from py_ecc.bn128.bn128_curve import G1, G2, add, multiply, eq,  curve_order
# from py_ecc.bn128.bn128_pairing import  pairing
# print(G1)
# print(G2)

# x = 10 # chosen randomly
# assert eq(multiply(G2, x + curve_order), multiply(G2, x))
# assert eq(multiply(G1, x + curve_order), multiply(G1, x))

# # G2 similar behaviour with cycic group i.e G1
# print(eq(add(G1, G1), multiply(G1, 2)))

# print(eq(add(G2, G2), multiply(G2, 2)))

# # add(G1, G2)  Different groups .. A NO NO

# # BILLINEAR PAIRING IN PYTHON 
# P = multiply(G1,3)
# Q= multiply(G2,8)
# R = multiply(G1,24)

# assert pairing(Q, P)== pairing(G2, R)

# # print(eq(pairing(Q, P), pairing(G2, R)))


# P_1 = multiply(G1, 3)
# P_2 = multiply(G2, 8)

# Q_1 = multiply(G1, 6)
# Q_2 = multiply(G2, 4)

# lhs = pairing(P_2, P_1) 
# rhs = pairing(Q_2, Q_1)

# assert lhs == rhs
# # assert eq(pairing( P_2,P_1), pairing(Q_2, Q_1))

# # How to address elements in GT
# # 2 * 3 = 6
# P_1 = multiply(G1, 2)
# P_2 = multiply(G2, 3)

# # 4 * 5 = 20
# Q_1 = multiply(G1, 4)
# Q_2 = multiply(G2, 5)

# # 13 * 2 = 26
# R_1 = multiply(G1, 13)
# R_2 = multiply(G2, 2)


# lhs = pairing(P_2, P_1) * pairing(Q_2, Q_1)
# rhs = pairing(R_2, R_1)

# assert lhs == rhs
# print("pairing equation passed")


# # b ^ {2 * 3} * b ^ {4 * 5} = b ^ {13 * 2}
# # b ^ 6 * b ^ 20 = b ^ 26

# # assert eq(pairing(P_2, P_1) * pairing(Q_2, Q_1), pairing(R_2, R_1))


import sys
import threading

sys.setrecursionlimit(100000)
threading.stack_size(134217728)  # 128 MB

def main():
    from py_ecc.bn128 import G1, G2, pairing, add, multiply, eq, curve_order


    print(G1)
    print(G2)

    x = 10 # chosen randomly
    assert eq(multiply(G2, x + curve_order), multiply(G2, x))
    assert eq(multiply(G1, x + curve_order), multiply(G1, x))

    # G2 similar behaviour with cycic group i.e G1
    print(eq(add(G1, G1), multiply(G1, 2)))

    print(eq(add(G2, G2), multiply(G2, 2)))

    # add(G1, G2)  Different groups .. A NO NO

    # BILLINEAR PAIRING IN PYTHON 
    P = multiply(G1,3)
    Q= multiply(G2,8)
    R = multiply(G1,24)

    assert pairing(Q, P) == pairing(G2, R)
    print(pairing(Q, P) == pairing(G2, R))


    P_1 = multiply(G1, 3)
    P_2 = multiply(G2, 8)

    Q_1 = multiply(G1, 6)
    Q_2 = multiply(G2, 4)

    lhs = pairing(P_2, P_1) 
    rhs = pairing(Q_2, Q_1)

    assert lhs == rhs
    assert pairing(P_2, P_1) == pairing(Q_2, Q_1)

    # How to address elements in GT
    # 2 * 3 = 6
    P_1 = multiply(G1, 2)
    P_2 = multiply(G2, 3)

    # 4 * 5 = 20
    Q_1 = multiply(G1, 4)
    Q_2 = multiply(G2, 5)

    # 13 * 2 = 26
    R_1 = multiply(G1, 13)
    R_2 = multiply(G2, 2)


    lhs = pairing(P_2, P_1) * pairing(Q_2, Q_1)
    rhs = pairing(R_2, R_1)

    assert lhs == rhs
    print("pairing equation passed")

    # b ^ {2 * 3} * b ^ {4 * 5} = b ^ {13 * 2}
    # b ^ 6 * b ^ 20 = b ^ 26

    assert pairing(P_2, P_1) * pairing(Q_2, Q_1) == pairing(R_2, R_1)

    # result = pairing(G2, G1)
    # print(result)

t = threading.Thread(target=main)
t.start()
t.join()
