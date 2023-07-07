import numpy as np
import sys


def read_file(filename):
    f = open(sys.path[0] + filename)
    # Ignore line 1 and the brackets
    res = []
    i = 0
    for x in f:
        line_str = x.replace("[", " ").replace("]", " ").split()
        line_float = []
        # Converting string list to float list
        for x in line_str:
            line_float.append(float(x))

        if i == 0:
            n = int(float(x))
        if i > 0:
            res.append(line_float)
        i += 1

    f.close()
    return res, n


if __name__ == "__main__":
    p = 2**26 - 5
    A, n = read_file("/../data/Matrix_A.txt")
    B, n = read_file("/../data/Matrix_B.txt")

    A_np = np.array(A)
    B_np = np.array(B)
    print()
    print(A_np)
    print()
    print(B_np)

    C = A_np * B_np
    print()
    print(C)
    C = np.mod(C, p)
    print()
    print(p)
    print(C)
