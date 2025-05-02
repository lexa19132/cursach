import gmpy2 as gmp;
import pointCounter as pc

# Будем считать, что со старта у нас приведены все коэффициенты по модулю p.

a_1 = gmp.mpz(1)
a_2 = gmp.mpz(0)
a_3 = gmp.mpz(0)
a_4 = gmp.mpz(0)
a_6 = gmp.mpz(1)
p = gmp.mpz(2)
n = gmp.mpz(11)

def curve(x : gmp.mpz, y : gmp.mpz) -> tuple[gmp.mpz, gmp.mpz]:
   #y^2 = - a_1 * x * y + a_2 * x^2 - a_3 * y + a_4 * x + x^3 + a_6
   left_side = gmp.powmod(y, 2, p) # уже по модулю
   right_side = (- a_1 * x * y + a_2 * gmp.powmod(x, 2, p) - a_3 * y + a_4 * x + gmp.powmod(x, 3, p) + a_6) % p
   return [left_side, right_side]

number_of_points_cycle : int = pc.cycle_point_number(curve, p)
print(f"Number of points using cycle: {number_of_points_cycle}")

# number_of_points_legendre : int = pc.legendre_point_number(curve, p)
# print(f"Number of pounts using legenre symbol: {number_of_points_legendre}")

print(pc.extension_point_number(pc.cycle_point_number(curve, p), p, n))


