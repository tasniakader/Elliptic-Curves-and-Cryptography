import random
import matplotlib.pyplot as plt
import numpy as np

# modular inverse of x in F_p
def inverse(x, p):
    # Fermat's little theorem:
    # prime p -> inverse(x) = x^(p - 2) mod p
    return pow(x, p-2, p)

def inv(x, n):
    try:
        return pow(x, -1, n)
    except ValueError:
        return None

def elliptic_add(P, Q, a, p):
    x1 = P[0]
    y1 = P[1]
    x2 = Q[0]
    y2 = Q[1]
    
    # (None, None) represents the identity element
    if (P == (None, None)):
        return Q
    if (Q == (None, None)):
        return P

    
    if(P == Q):
        m = ((3 * x1 ** 2 + a) * inverse(2 * y1, p)) % p
    else:
        # if P and Q are inverses, return identity
        if (x1 == x2):
            return (None, None)
        
        m = ((y2 - y1) * inverse(x2 - x1, p)) % p
        
    x3 = (m ** 2 - x1 - x2) % p
    y3 = (m * (x1 - x3) - y1) % p
    
    return (x3, y3)

def elliptic_scalar_mul(k, P, a, p):
    if(k == 0):
        return (None, None)
    
    if(k == 1):
        return P
    
    R = elliptic_add(P, P, a, p)
    
    for i in range(0, k - 2):
        R = elliptic_add(P, R, a, p)
    
    return R

def elliptic_lin_combo(k, P, l, Q, a, p):
    R = elliptic_scalar_mul(k, P, a, p)
    Q = elliptic_scalar_mul(l, Q, a, p)
    
    return elliptic_add(R, Q, a, p)


def plot_points(path):
    plt.figure(figsize = (8, 8))
    
    x = []
    y = []
    
    for point in path:
        x.append(point[0])
        y.append(point[1])
        
    plt.scatter(x, y, color="blue", zorder=5)

    for i in range(len(path) - 1):
        plt.plot([x[i], x[i + 1]], [y[i], y[i + 1]], color="red")

    plt.title("Walk Used by Pollard's Algorithm")
    plt.grid(True)
    plt.axhline(0, color='black',linewidth=1)
    plt.axvline(0, color='black',linewidth=1)
    
    plt.show()

def pollard_rho(P, Q, a, p, n):
    # pick two random integers k, l
    r = random.randint(0, n - 1)
    s = random.randint(0, n - 1)
    print(f"r = {r}, s = {s}")
    
    R = elliptic_lin_combo(r, P, s, Q, a, p)
    print(f"R = rP ⊕ sQ = {r}{P} ⊕ {s}{Q} = {R}\n")
    
    # T: set of visited points
    T = [[], [], []]
    
    while(R not in T[0]):
        if(R == (None, None)):
            return -1
        
        T[0].append(R)
        T[1].append(r)
        T[2].append(s)
        
        print(f"Add R{R} to T:")
        print(f"T = {T}")
        
        if(R[0] % 3 == 0):
            print(f"{R[0]} is 0 mod 3 -> add Q{Q} to {R}")
            R = elliptic_add(Q, R, a, p)
            s = (s + 1) % n
        elif(R[0] % 3 == 1):
            print(f"{R[0]} is 1 mod 3 -> double {R}")
            R = elliptic_scalar_mul(2, R, a, p)
            r = (2 * r) % n
            s = (2 * s) % n
        else:
            print(f"{R[0]} is 2 mod 3 -> add P{P} to {R}")
            R = elliptic_add(P, R, a, p)
            r = (r + 1) % n
        
        print(f"R = {R}, r = {r}, s = {s}\n")
        
        if(R in T[0]):
            print("Collision!\n")

    index = T[0].index(R)
    r1 = T[1][index]
    s1 = T[2][index]

    r2 = r
    s2 = s
    
    print(f"Q = kP, k = (({r1} - {r2}) / ({s2} - {s1})) mod {n}")
    
    num = (r1 - r2) % n
    den = (s2 - s1) % n
    
    if(inv(den, n) is None):
        return -1
    
    T[0].append(R)
    plot_points(T[0])
    
    return (num * inv(den, n)) % n

while(True):
    print("(1) P ⊕ Q")
    print("(2) kP")
    print("(3) kP ⊕ lQ")
    print("(4) Pollard's rho-algorithm")
    print("(5) Quit")
    option = int(input("Choose an option: "))
    
    if(option > 5 or option < 1):
        print("invalid input")
        break;
    if(option == 5):
        break;
        
    a = int(input("a: "))
    p = int(input("p: "))
    P_str = input("P: ")
    P = tuple(map(int, P_str.split(',')))
    
    if(option == 1):
        Q_str = input("Q: ")
        Q = tuple(map(int, Q_str.split(',')))
        print(f"P ⊕ Q = {P} ⊕ {Q} = {elliptic_add(P, Q, a, p)}")
    elif(option == 2):
        k = int(input("k: "))
        print(f"kP = {k}{P} = {elliptic_scalar_mul(k, P, a, p)}")
    elif(option == 3):
        k = int(input("k: "))
        Q_str = input("Q: ")
        Q = tuple(map(int, Q_str.split(',')))
        l = int(input("l: "))
        print(f"kP ⊕ lQ = {k}{P} ⊕ {l}{Q} = {elliptic_lin_combo(k, P, l, Q, a, p)}")
    elif(option == 4):
        Q_str = input("Q: ")
        Q = tuple(map(int, Q_str.split(',')))
        n = int(input("n: "))
        print(f"k = {pollard_rho(P, Q, a, p, n)}")
    
    print()
