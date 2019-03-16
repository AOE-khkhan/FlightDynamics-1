import symengine

vars = symengine.symbols('x y')  # Define x and y variables
f = symengine.sympify(['y*x**2', '5*x + sin(y)']) # Define function
J = symengine.zeros(len(f), len(vars)) # Initiate Jacobian matrix

# Fill Jacobian Matrix with entries
for i, fi in enumerate(f):
    for j, s in enumerate(vars):
        J[i,j] = symengine.diff(fi, s)

print(J)

