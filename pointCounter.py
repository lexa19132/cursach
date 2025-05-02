import gmpy2 as gmp;
from typing import Callable;

def cycle_point_number(curve: Callable[[gmp.mpz, gmp.mpz], tuple[gmp.mpz, gmp.mpz]], p: gmp.mpz) -> int:
    counter = 0

    p_int = int(p)
    
    for i in range(p_int):
        for j in range(p_int):
            x = gmp.mpz(i)
            y = gmp.mpz(j)
            results = curve(x, y)
            if results[0] == results[1]:
                counter += 1
    return counter + 1

def legendre_point_number(curve: Callable[[gmp.mpz, gmp.mpz], tuple[gmp.mpz, gmp.mpz]], p: gmp.mpz) -> int:

    counter = 0
    p_int = int(p)
    
    for x in (gmp.mpz(i) for i in range(p_int)):
        y_squared = curve(x, gmp.mpz(0))[1]  # (x, y²)
        y_squared_mod = y_squared % p
        
        if y_squared_mod == 0:
            counter += 1
        else:
            ls = gmp.legendre(y_squared_mod, p)
            counter += 1 + ls
    
    return counter + 1

def extension_point_number(points_count: gmp.mpz, p: gmp.mpz, n: gmp.mpz) -> int:

    a = p + 1 - points_count
    a = gmp.mpz(a)
    p = gmp.mpz(p)
    n = gmp.mpz(n)

    def matrix_pow(n_power: gmp.mpz) -> gmp.mpz:
        matrix = [[a, -p], [gmp.mpz(1), gmp.mpz(0)]]
        result = [[gmp.mpz(1), gmp.mpz(0)], [gmp.mpz(0), gmp.mpz(1)]]

        while n_power > 0:
            if n_power % 2 == 1:
                new_row0 = [
                    result[0][0] * matrix[0][0] + result[0][1] * matrix[1][0],
                    result[0][0] * matrix[0][1] + result[0][1] * matrix[1][1]
                ]
                new_row1 = [
                    result[1][0] * matrix[0][0] + result[1][1] * matrix[1][0],
                    result[1][0] * matrix[0][1] + result[1][1] * matrix[1][1]
                ]
                result = [new_row0, new_row1]
            new_matrix00 = matrix[0][0] * matrix[0][0] + matrix[0][1] * matrix[1][0]
            new_matrix01 = matrix[0][0] * matrix[0][1] + matrix[0][1] * matrix[1][1]
            new_matrix10 = matrix[1][0] * matrix[0][0] + matrix[1][1] * matrix[1][0]
            new_matrix11 = matrix[1][0] * matrix[0][1] + matrix[1][1] * matrix[1][1]
            matrix = [[new_matrix00, new_matrix01], [new_matrix10, new_matrix11]]
            n_power //= 2
        return result[0][0] * a + result[0][1] * 2

    a_n = matrix_pow(n - 1)
    p_pow = gmp.powmod(p, n, gmp.mpz(2)**(67108864))
    Nn = p_pow + 1 - a_n
    return int(Nn)