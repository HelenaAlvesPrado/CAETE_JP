import math

f0 = -0.5               # -90 to 90
f1 = f0  + 0.5


R_earth = 6371.0027 # km

t1 = abs(math.sin(f1) - math.sin(f0))

t2 = math.radians(0.5)

s = math.pow(R_earth,2) * t2 * t1

print(s)
