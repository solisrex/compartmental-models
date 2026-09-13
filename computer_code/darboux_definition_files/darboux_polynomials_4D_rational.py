initialized_darbouxs = []

def log_init(init_function):
   def wrapper(self,*args,**kwargs):
       initialized_darbouxs.append(self)
       new_object = init_function(self,*args,**kwargs)
       return new_object
   return wrapper

from sympy import symbols, Function, Eq, sqrt, solve, diff, Derivative, simplify

As = [symbols('a'+str(i)) for i in range(6)]
Bs = [symbols('b'+str(i)) for i in range(12)]

t = symbols('t')

X1 = Function('X1')(t)
X2 = Function('X2')(t)   
X3 = Function('X3')(t)
X4 = Function('X4')(t)

X1_dot =  ( -As[0]*X1*X2 - As[1]*X1*X3 - As[2]*X1*X4)/(X1+X2+X3) - (Bs[3]+Bs[6]+Bs[9])*X1 + Bs[0]*X2 + Bs[1]*X3 + Bs[2]*X4;
X2_dot =  (  As[0]*X1*X2 - As[3]*X2*X3 - As[4]*X2*X4)/(X1+X2+X3) + Bs[3]*X1 - (Bs[0]+Bs[7]+Bs[10])*X2 + Bs[4]*X3 + Bs[5]*X4;
X3_dot =  (  As[1]*X1*X3 + As[3]*X2*X3 - As[5]*X3*X4)/(X1+X2+X3) + Bs[6]*X1 + Bs[7]*X2 - (Bs[1]+Bs[4]+Bs[11])*X3 + Bs[8]*X4;
X4_dot =  (  As[2]*X1*X4 + As[4]*X2*X4 + As[5]*X3*X4)/(X1+X2+X3) + Bs[9]*X1 + Bs[10]*X2 + Bs[11]*X3 - (Bs[2]+Bs[5]+Bs[8])*X4;

derivative_substitutions = {
    Derivative(X1, t):X1_dot,
    Derivative(X2, t):X2_dot,
    Derivative(X3, t):X3_dot,
    Derivative(X4, t):X4_dot
}

total_darboux = []

def track_items(cls):

    return cls


class DarbouxPolynomial:

    @log_init
    def __init__(self, name, constraints, darboux_list, cofactor,validate=True): 
        self.name = name
        self.constraints = constraints
        self.darboux_list = darboux_list
        self.cofactor = cofactor
        self.validate = validate

    def __repr__(self):
        return self.name

    def __str__(self):
        return self.name
        
def check_darboux(darboux_polynomial):
    if not darboux_polynomial.validate:
            return False
    
    if len(darboux_polynomial.darboux_list) < 1:
        return False
    

    combined_constraints = solve(darboux_polynomial.constraints,As+Bs,dict=True)
    if len(combined_constraints) == 0:
        return False
    
    cofactor = darboux_polynomial.cofactor
    for darboux in darboux_polynomial.darboux_list:
        for constraint_element in combined_constraints:
            darboux_dot = diff(darboux,t).subs(derivative_substitutions).subs(constraint_element)
            product = (darboux*cofactor).subs(constraint_element)
            difference = simplify(darboux_dot-product)
            if not difference.equals(0):
                return False
            
    return True


from sympy import Eq

# Generated from 4D_rational_2_solved(1).txt
# lambda has been fixed to 1; nonzero conditions are intentionally omitted.
# Source entries of the form `nan = 0` are omitted as invalid constraints.


darboux_2 = DarbouxPolynomial(
    "darboux_2",
    [Eq((As[2] - As[4])*(As[2] - As[5] + Bs[11]), 0), Eq((As[2] - As[5] + Bs[11])*(As[2] - Bs[2] - Bs[5] - Bs[8]), 0), Eq(Bs[9], 0), Eq(Bs[10], 0)],
    [(-Bs[11]/(As[2] - As[5]))*X1 + (-Bs[11]/(As[2] - As[5]))*X2 + (-Bs[11]/(As[2] - As[5]))*X3 + (1)*X4],
    (-As[2] + As[5] - Bs[11])*X3/(X1 + X2 + X3),
)



darboux_3 = DarbouxPolynomial(
    "darboux_3",
    [Eq(-As[2] + Bs[2] + Bs[5] + Bs[8], 0), Eq(-As[4] + Bs[2] + Bs[5] + Bs[8], 0), Eq(Bs[9], 0), Eq(Bs[10], 0)],
    [(Bs[11]/(As[5] - Bs[2] - Bs[5] - Bs[8]))*X1 + (Bs[11]/(As[5] - Bs[2] - Bs[5] - Bs[8]))*X2 + (Bs[11]/(As[5] - Bs[2] - Bs[5] - Bs[8]))*X3 + (1)*X4],
    (As[5] - Bs[11] - Bs[2] - Bs[5] - Bs[8])*X3/(X1 + X2 + X3),
)


darboux_4 = DarbouxPolynomial(
    "darboux_4",
    [Eq((As[2] - As[4])*(As[2] - As[5] + Bs[11]), 0), Eq((As[2] - As[5] + Bs[11])*(As[2] - Bs[2] - Bs[5] - Bs[8]), 0), Eq(Bs[9], 0), Eq(Bs[10], 0)],
    [(-Bs[11]/(As[2] - As[5]))*X1 + (-Bs[11]/(As[2] - As[5]))*X2 + (-Bs[11]/(As[2] - As[5]))*X3 + (1)*X4],
    (-As[2] + As[5] - Bs[11])*X3/(X1 + X2 + X3),
)


darboux_5 = DarbouxPolynomial(
    "darboux_5",
    [Eq(As[2] - As[3] - As[5] + Bs[11], 0), Eq(As[2] - As[4], 0), Eq(-As[2] + Bs[2] + Bs[5] + Bs[8], 0), Eq(Bs[9], 0), Eq(Bs[10], 0)],
    [(1)*X1 + (1)*X2 + (1)*X3 + (-(As[2] - As[5])/Bs[11])*X4],
    (-As[2] + As[5] - Bs[11])*X3/(X1 + X2 + X3),
)


darboux_6 = DarbouxPolynomial(
    "darboux_6",
    [Eq(As[1], 0), Eq(As[3], 0)],
    [(1)*X1 + (1)*X2 + (1)*X3 + (1)*X4],
    (0)*X3/(X1 + X2 + X3),
)

darboux_7 = DarbouxPolynomial(
    "darboux_7",
    [Eq(-As[1] - As[5] + Bs[1] + Bs[11] + Bs[4], 0), Eq(-As[3] - As[5] + Bs[1] + Bs[11] + Bs[4], 0), Eq(Bs[6], 0), Eq(Bs[7], 0), Eq(Bs[8], 0)],
    [(As[5]/(Bs[1] + Bs[11] + Bs[4]))*X1 + (As[5]/(Bs[1] + Bs[11] + Bs[4]))*X2 + (1)*X3 + (As[5]/(Bs[1] + Bs[11] + Bs[4]))*X4],
    (As[5] - Bs[1] - Bs[11] - Bs[4])*X3/(X1 + X2 + X3),
)

darboux_8 = DarbouxPolynomial(
    "darboux_8",
    [Eq(As[1] - As[3], 0), Eq(-As[1] - As[5] + Bs[1] + Bs[11] + Bs[4], 0), Eq(Bs[6], 0), Eq(Bs[7], 0), Eq(Bs[8], 0)],
    [(1)*X1 + (1)*X2 + ((As[1] + As[5])/As[5])*X3 + (1)*X4],
    (-As[1])*X3/(X1 + X2 + X3),
)

darboux_9 = DarbouxPolynomial(
    "darboux_9",
    [Eq(As[0], 0), Eq(As[3], 0)],
    [(1)*X1 + (1)*X2 + (1)*X3 + (1)*X4],
    (0)*X3/(X1 + X2 + X3),
)

darboux_10 = DarbouxPolynomial(
    "darboux_10",
    [Eq(As[0], 0), Eq(As[3]*(As[2]*Bs[3] - Bs[2]*Bs[3] - Bs[3]*Bs[5] - Bs[3]*Bs[8] - Bs[5]*Bs[9]), 0), Eq(As[3]*(As[4]*Bs[3] + As[4]*Bs[9] - Bs[2]*Bs[3] - Bs[3]*Bs[5] - Bs[3]*Bs[8] - Bs[5]*Bs[9]), 0), Eq(As[3]*(Bs[11]*Bs[3] - Bs[4]*Bs[9])*(As[3]*Bs[3] + As[5]*Bs[3] - Bs[11]*Bs[3] - Bs[2]*Bs[3] - Bs[3]*Bs[5] - Bs[3]*Bs[8] + Bs[4]*Bs[9] - Bs[5]*Bs[9]), 0), Eq(As[3]*(Bs[0]*Bs[9] + Bs[10]*Bs[3] + Bs[10]*Bs[9] + Bs[7]*Bs[9]), 0)],
    [(1)*X1 + ((As[3]*Bs[9] + Bs[11]*Bs[3] - Bs[4]*Bs[9])/(Bs[11]*Bs[3] - Bs[4]*Bs[9]))*X2 + (1)*X3 + (-(As[3]*Bs[3] - Bs[11]*Bs[3] + Bs[4]*Bs[9])/(Bs[11]*Bs[3] - Bs[4]*Bs[9]))*X4],
    (-As[3])*X3/(X1 + X2 + X3),
    validate=False
)

darboux_11 = DarbouxPolynomial(
    "darboux_11",
    [Eq(As[0], 0), Eq(As[1]*(As[2]*Bs[0]*Bs[6] + As[2]*Bs[10]*Bs[6] + As[2]*Bs[3]*Bs[7] + As[2]*Bs[6]*Bs[7] - Bs[0]*Bs[2]*Bs[6] - Bs[0]*Bs[5]*Bs[6] - Bs[0]*Bs[6]*Bs[8] - Bs[0]*Bs[8]*Bs[9] - Bs[10]*Bs[2]*Bs[6] - Bs[10]*Bs[3]*Bs[8] - Bs[10]*Bs[6]*Bs[8] - Bs[10]*Bs[8]*Bs[9] - Bs[2]*Bs[3]*Bs[7] - Bs[2]*Bs[6]*Bs[7] - Bs[3]*Bs[5]*Bs[7] - Bs[3]*Bs[7]*Bs[8] - Bs[5]*Bs[6]*Bs[7] - Bs[5]*Bs[7]*Bs[9] - Bs[6]*Bs[7]*Bs[8] - Bs[7]*Bs[8]*Bs[9]), 0), Eq((As[1] - As[3])*(Bs[0]*Bs[1]*Bs[9] + Bs[0]*Bs[11]*Bs[6] + Bs[0]*Bs[11]*Bs[9] + Bs[0]*Bs[4]*Bs[9] + Bs[1]*Bs[10]*Bs[3] + Bs[1]*Bs[10]*Bs[9] + Bs[1]*Bs[7]*Bs[9] + Bs[10]*Bs[11]*Bs[3] + Bs[10]*Bs[11]*Bs[6] + Bs[10]*Bs[11]*Bs[9] + Bs[10]*Bs[3]*Bs[4] + Bs[10]*Bs[4]*Bs[6] + Bs[10]*Bs[4]*Bs[9] + Bs[11]*Bs[3]*Bs[7] + Bs[11]*Bs[6]*Bs[7] + Bs[11]*Bs[7]*Bs[9]), 0), Eq(As[1]*(As[4]*Bs[0]*Bs[6] + As[4]*Bs[3]*Bs[7] + As[4]*Bs[6]*Bs[7] + As[4]*Bs[7]*Bs[9] - Bs[0]*Bs[2]*Bs[6] - Bs[0]*Bs[5]*Bs[6] - Bs[0]*Bs[6]*Bs[8] - Bs[0]*Bs[8]*Bs[9] - Bs[10]*Bs[2]*Bs[6] - Bs[10]*Bs[3]*Bs[8] - Bs[10]*Bs[6]*Bs[8] - Bs[10]*Bs[8]*Bs[9] - Bs[2]*Bs[3]*Bs[7] - Bs[2]*Bs[6]*Bs[7] - Bs[3]*Bs[5]*Bs[7] - Bs[3]*Bs[7]*Bs[8] - Bs[5]*Bs[6]*Bs[7] - Bs[5]*Bs[7]*Bs[9] - Bs[6]*Bs[7]*Bs[8] - Bs[7]*Bs[8]*Bs[9]), 0), Eq(As[1]*(Bs[0]*Bs[1]*Bs[9] + Bs[0]*Bs[11]*Bs[6] + Bs[0]*Bs[11]*Bs[9] + Bs[0]*Bs[4]*Bs[9] + Bs[1]*Bs[10]*Bs[3] + Bs[1]*Bs[10]*Bs[9] + Bs[1]*Bs[7]*Bs[9] + Bs[10]*Bs[11]*Bs[3] + Bs[10]*Bs[11]*Bs[6] + Bs[10]*Bs[11]*Bs[9] + Bs[10]*Bs[3]*Bs[4] + Bs[10]*Bs[4]*Bs[6] + Bs[10]*Bs[4]*Bs[9] + Bs[11]*Bs[3]*Bs[7] + Bs[11]*Bs[6]*Bs[7] + Bs[11]*Bs[7]*Bs[9])*(As[1]*Bs[0]*Bs[6] + As[1]*Bs[0]*Bs[9] + As[1]*Bs[10]*Bs[3] + As[1]*Bs[10]*Bs[6] + As[1]*Bs[10]*Bs[9] + As[1]*Bs[3]*Bs[7] + As[1]*Bs[6]*Bs[7] + As[1]*Bs[7]*Bs[9] + As[5]*Bs[0]*Bs[6] + As[5]*Bs[0]*Bs[9] + As[5]*Bs[10]*Bs[3] + As[5]*Bs[10]*Bs[6] + As[5]*Bs[10]*Bs[9] + As[5]*Bs[3]*Bs[7] + As[5]*Bs[6]*Bs[7] + As[5]*Bs[7]*Bs[9] - Bs[0]*Bs[1]*Bs[9] - Bs[0]*Bs[11]*Bs[6] - Bs[0]*Bs[11]*Bs[9] - Bs[0]*Bs[2]*Bs[6] - Bs[0]*Bs[4]*Bs[9] - Bs[0]*Bs[5]*Bs[6] - Bs[0]*Bs[6]*Bs[8] - Bs[0]*Bs[8]*Bs[9] - Bs[1]*Bs[10]*Bs[3] - Bs[1]*Bs[10]*Bs[9] - Bs[1]*Bs[7]*Bs[9] - Bs[10]*Bs[11]*Bs[3] - Bs[10]*Bs[11]*Bs[6] - Bs[10]*Bs[11]*Bs[9] - Bs[10]*Bs[2]*Bs[6] - Bs[10]*Bs[3]*Bs[4] - Bs[10]*Bs[3]*Bs[8] - Bs[10]*Bs[4]*Bs[6] - Bs[10]*Bs[4]*Bs[9] - Bs[10]*Bs[6]*Bs[8] - Bs[10]*Bs[8]*Bs[9] - Bs[11]*Bs[3]*Bs[7] - Bs[11]*Bs[6]*Bs[7] - Bs[11]*Bs[7]*Bs[9] - Bs[2]*Bs[3]*Bs[7] - Bs[2]*Bs[6]*Bs[7] - Bs[3]*Bs[5]*Bs[7] - Bs[3]*Bs[7]*Bs[8] - Bs[5]*Bs[6]*Bs[7] - Bs[5]*Bs[7]*Bs[9] - Bs[6]*Bs[7]*Bs[8] - Bs[7]*Bs[8]*Bs[9]), 0)],
    [(-(As[1]*Bs[0]*Bs[9] + As[1]*Bs[10]*Bs[3] + As[1]*Bs[10]*Bs[9] + As[1]*Bs[7]*Bs[9] - Bs[0]*Bs[1]*Bs[9] - Bs[0]*Bs[11]*Bs[6] - Bs[0]*Bs[11]*Bs[9] - Bs[0]*Bs[4]*Bs[9] - Bs[1]*Bs[10]*Bs[3] - Bs[1]*Bs[10]*Bs[9] - Bs[1]*Bs[7]*Bs[9] - Bs[10]*Bs[11]*Bs[3] - Bs[10]*Bs[11]*Bs[6] - Bs[10]*Bs[11]*Bs[9] - Bs[10]*Bs[3]*Bs[4] - Bs[10]*Bs[4]*Bs[6] - Bs[10]*Bs[4]*Bs[9] - Bs[11]*Bs[3]*Bs[7] - Bs[11]*Bs[6]*Bs[7] - Bs[11]*Bs[7]*Bs[9])/(Bs[0]*Bs[1]*Bs[9] + Bs[0]*Bs[11]*Bs[6] + Bs[0]*Bs[11]*Bs[9] + Bs[0]*Bs[4]*Bs[9] + Bs[1]*Bs[10]*Bs[3] + Bs[1]*Bs[10]*Bs[9] + Bs[1]*Bs[7]*Bs[9] + Bs[10]*Bs[11]*Bs[3] + Bs[10]*Bs[11]*Bs[6] + Bs[10]*Bs[11]*Bs[9] + Bs[10]*Bs[3]*Bs[4] + Bs[10]*Bs[4]*Bs[6] + Bs[10]*Bs[4]*Bs[9] + Bs[11]*Bs[3]*Bs[7] + Bs[11]*Bs[6]*Bs[7] + Bs[11]*Bs[7]*Bs[9]))*X1 + (-(As[1]*Bs[0]*Bs[9] + As[1]*Bs[10]*Bs[3] + As[1]*Bs[10]*Bs[6] + As[1]*Bs[10]*Bs[9] - Bs[0]*Bs[1]*Bs[9] - Bs[0]*Bs[11]*Bs[6] - Bs[0]*Bs[11]*Bs[9] - Bs[0]*Bs[4]*Bs[9] - Bs[1]*Bs[10]*Bs[3] - Bs[1]*Bs[10]*Bs[9] - Bs[1]*Bs[7]*Bs[9] - Bs[10]*Bs[11]*Bs[3] - Bs[10]*Bs[11]*Bs[6] - Bs[10]*Bs[11]*Bs[9] - Bs[10]*Bs[3]*Bs[4] - Bs[10]*Bs[4]*Bs[6] - Bs[10]*Bs[4]*Bs[9] - Bs[11]*Bs[3]*Bs[7] - Bs[11]*Bs[6]*Bs[7] - Bs[11]*Bs[7]*Bs[9])/(Bs[0]*Bs[1]*Bs[9] + Bs[0]*Bs[11]*Bs[6] + Bs[0]*Bs[11]*Bs[9] + Bs[0]*Bs[4]*Bs[9] + Bs[1]*Bs[10]*Bs[3] + Bs[1]*Bs[10]*Bs[9] + Bs[1]*Bs[7]*Bs[9] + Bs[10]*Bs[11]*Bs[3] + Bs[10]*Bs[11]*Bs[6] + Bs[10]*Bs[11]*Bs[9] + Bs[10]*Bs[3]*Bs[4] + Bs[10]*Bs[4]*Bs[6] + Bs[10]*Bs[4]*Bs[9] + Bs[11]*Bs[3]*Bs[7] + Bs[11]*Bs[6]*Bs[7] + Bs[11]*Bs[7]*Bs[9]))*X2 + (1)*X3 + (-(As[1]*Bs[0]*Bs[6] + As[1]*Bs[0]*Bs[9] + As[1]*Bs[10]*Bs[3] + As[1]*Bs[10]*Bs[6] + As[1]*Bs[10]*Bs[9] + As[1]*Bs[3]*Bs[7] + As[1]*Bs[6]*Bs[7] + As[1]*Bs[7]*Bs[9] - Bs[0]*Bs[1]*Bs[9] - Bs[0]*Bs[11]*Bs[6] - Bs[0]*Bs[11]*Bs[9] - Bs[0]*Bs[4]*Bs[9] - Bs[1]*Bs[10]*Bs[3] - Bs[1]*Bs[10]*Bs[9] - Bs[1]*Bs[7]*Bs[9] - Bs[10]*Bs[11]*Bs[3] - Bs[10]*Bs[11]*Bs[6] - Bs[10]*Bs[11]*Bs[9] - Bs[10]*Bs[3]*Bs[4] - Bs[10]*Bs[4]*Bs[6] - Bs[10]*Bs[4]*Bs[9] - Bs[11]*Bs[3]*Bs[7] - Bs[11]*Bs[6]*Bs[7] - Bs[11]*Bs[7]*Bs[9])/(Bs[0]*Bs[1]*Bs[9] + Bs[0]*Bs[11]*Bs[6] + Bs[0]*Bs[11]*Bs[9] + Bs[0]*Bs[4]*Bs[9] + Bs[1]*Bs[10]*Bs[3] + Bs[1]*Bs[10]*Bs[9] + Bs[1]*Bs[7]*Bs[9] + Bs[10]*Bs[11]*Bs[3] + Bs[10]*Bs[11]*Bs[6] + Bs[10]*Bs[11]*Bs[9] + Bs[10]*Bs[3]*Bs[4] + Bs[10]*Bs[4]*Bs[6] + Bs[10]*Bs[4]*Bs[9] + Bs[11]*Bs[3]*Bs[7] + Bs[11]*Bs[6]*Bs[7] + Bs[11]*Bs[7]*Bs[9]))*X4],
    (-As[1])*X3/(X1 + X2 + X3),
    validate=False
)

darboux_12 = DarbouxPolynomial(
    "darboux_12",
    [Eq(As[0], 0), Eq(As[1], 0)],
    [(1)*X1 + (1)*X2 + (1)*X3 + (1)*X4],
    (0)*X3/(X1 + X2 + X3),
)

darboux_13 = DarbouxPolynomial(
    "darboux_13",
    [Eq(As[1], 0), Eq(As[3], 0)],
    [(1)*X1 + (1)*X2 + (1)*X3 + (1)*X4],
    (0)*X3/(X1 + X2 + X3),
)

darboux_14 = DarbouxPolynomial(
    "darboux_14",
    [Eq(As[1] - As[3], 0), Eq(-As[1] - As[5] + Bs[1] + Bs[11] + Bs[4], 0), Eq(Bs[6], 0), Eq(Bs[7], 0), Eq(Bs[8], 0)],
    [(1)*X1 + (1)*X2 + ((As[1] + As[5])/As[5])*X3 + (1)*X4],
    (-As[1])*X3/(X1 + X2 + X3),
)

darboux_15 = DarbouxPolynomial(
    "darboux_15",
    [Eq(As[1] - As[3], 0), Eq(-As[2]*Bs[6] + Bs[2]*Bs[6] + Bs[5]*Bs[6] + Bs[6]*Bs[8] + Bs[8]*Bs[9], 0), Eq(-As[4]*Bs[6] + Bs[2]*Bs[6] + Bs[5]*Bs[6] + Bs[6]*Bs[8] + Bs[8]*Bs[9], 0), Eq(-As[5]*Bs[6] - As[5]*Bs[9] + Bs[2]*Bs[6] + Bs[5]*Bs[6] + Bs[6]*Bs[8] + Bs[8]*Bs[9], 0), Eq(-As[3]*Bs[6] - As[3]*Bs[9] + Bs[1]*Bs[9] + Bs[11]*Bs[6] + Bs[11]*Bs[9] + Bs[4]*Bs[9], 0), Eq(Bs[10]*Bs[6] - Bs[7]*Bs[9], 0)],
    [(1)*X1 + (1)*X2 + ((Bs[6] + Bs[9])/Bs[6])*X3 + (0)*X4],
    (-As[3])*X3/(X1 + X2 + X3),
    validate=False
)

darboux_16 = DarbouxPolynomial(
    "darboux_16",
    [Eq(-As[2] + Bs[1] + Bs[4], 0), Eq(As[1] - As[3], 0), Eq(-As[4] + Bs[1] + Bs[4], 0), Eq(Bs[1] - Bs[2] + Bs[4] - Bs[5], 0), Eq(Bs[10] + Bs[7], 0), Eq(Bs[6] + Bs[9], 0)],
    [(-(As[1] - Bs[1] - Bs[4])/(Bs[1] + Bs[4]))*X1 + (-(As[1] - Bs[1] - Bs[4])/(Bs[1] + Bs[4]))*X2 + (1)*X3 + (1)*X4],
    (-As[1])*X3/(X1 + X2 + X3),
)

darboux_17 = DarbouxPolynomial(
    "darboux_17",
    [Eq(As[1] - As[3], 0), Eq(As[5], 0), Eq(-As[3] + Bs[1] + Bs[11] + Bs[4], 0), Eq(Bs[6], 0), Eq(Bs[7], 0), Eq(Bs[8], 0)],
    [(0)*X1 + (0)*X2 + (1)*X3 + (0)*X4],
    (-As[3])*X3/(X1 + X2 + X3),
)

darboux_18 = DarbouxPolynomial(
    "darboux_18",
    [Eq(As[0], 0), Eq(-As[2]*Bs[3] + Bs[2]*Bs[3] + Bs[3]*Bs[5] + Bs[3]*Bs[8] + Bs[5]*Bs[9], 0), Eq(-As[4]*Bs[3] - As[4]*Bs[9] + Bs[2]*Bs[3] + Bs[3]*Bs[5] + Bs[3]*Bs[8] + Bs[5]*Bs[9], 0), Eq(-As[5]*Bs[3] + Bs[2]*Bs[3] + Bs[3]*Bs[5] + Bs[3]*Bs[8] + Bs[5]*Bs[9], 0), Eq(Bs[0]*Bs[9] + Bs[10]*Bs[3] + Bs[10]*Bs[9] + Bs[7]*Bs[9], 0), Eq(As[3]*Bs[3] - Bs[11]*Bs[3] + Bs[4]*Bs[9], 0)],
    [(1)*X1 + ((Bs[3] + Bs[9])/Bs[3])*X2 + (1)*X3 + (0)*X4],
    (-As[3])*X3/(X1 + X2 + X3),
    validate=False
)

darboux_19 = DarbouxPolynomial(
    "darboux_19",
    [Eq(As[0], 0), Eq(As[2]*(As[1] - As[3]), 0), Eq(As[1]*(As[2]*As[4] - As[2]*Bs[5] - As[4]*Bs[2]), 0), Eq(As[1]*(As[2]*Bs[0] + As[2]*Bs[10] + As[2]*Bs[7] - Bs[0]*Bs[2] - Bs[0]*Bs[5] - Bs[10]*Bs[2] - Bs[2]*Bs[7]), 0), Eq(As[1]*(As[2]*Bs[4] - As[2]*Bs[5] + Bs[1]*Bs[5] - Bs[2]*Bs[4]), 0), Eq(As[1]*(As[2]*Bs[3] - Bs[2]*Bs[3] - Bs[3]*Bs[5] - Bs[5]*Bs[6] - Bs[5]*Bs[9]), 0)],
    [(-(As[1] - As[2])/As[2])*X1 + (-(As[1]*As[2] - As[1]*Bs[2] - As[2]*Bs[5])/(As[2]*Bs[5]))*X2 + (1)*X3 + (1)*X4],
    (-As[1])*X3/(X1 + X2 + X3),
    validate=False
)

darboux_20 = DarbouxPolynomial(
    "darboux_20",
    [Eq(As[0], 0), Eq(-As[2]*Bs[0] - As[2]*Bs[10] + Bs[0]*Bs[2] + Bs[0]*Bs[5] + Bs[0]*Bs[8] + Bs[10]*Bs[2], 0), Eq(-As[4]*Bs[0] + Bs[0]*Bs[2] + Bs[0]*Bs[5] + Bs[0]*Bs[8] + Bs[10]*Bs[2], 0), Eq(-As[5]*Bs[0] + Bs[0]*Bs[2] + Bs[0]*Bs[5] + Bs[0]*Bs[8] + Bs[10]*Bs[2], 0), Eq(As[1]*Bs[0] - Bs[0]*Bs[11] + Bs[1]*Bs[10], 0), Eq(Bs[0]*Bs[9] + Bs[10]*Bs[3] + Bs[10]*Bs[6] + Bs[10]*Bs[9], 0)],
    [((Bs[0] + Bs[10])/Bs[0])*X1 + (1)*X2 + (1)*X3 + (0)*X4],
    (-As[1])*X3/(X1 + X2 + X3),
    validate=False
)

darboux_21 = DarbouxPolynomial(
    "darboux_21",
    [Eq(As[2] - As[4], 0), Eq(As[3] - As[4] + As[5], 0), Eq(-As[4] + Bs[2] + Bs[5] + Bs[8], 0), Eq(Bs[9], 0), Eq(Bs[10], 0), Eq(Bs[11], 0)],
    [(0)*X1 + (0)*X2 + (0)*X3 + (1)*X4],
    (-As[3])*X3/(X1 + X2 + X3),
)

darboux_22 = DarbouxPolynomial(
    "darboux_22",
    [Eq(-As[1]*Bs[1] - As[1]*Bs[4] - As[5]*Bs[1] - As[5]*Bs[4] + Bs[1]*Bs[2] + Bs[1]*Bs[5] + Bs[1]*Bs[8] + Bs[11]*Bs[2] + Bs[11]*Bs[5] + Bs[2]*Bs[4] + Bs[4]*Bs[5] + Bs[4]*Bs[8], 0), Eq(-As[2]*Bs[1] - As[2]*Bs[11] - As[2]*Bs[4] + Bs[1]*Bs[2] + Bs[1]*Bs[5] + Bs[1]*Bs[8] + Bs[11]*Bs[2] + Bs[11]*Bs[5] + Bs[2]*Bs[4] + Bs[4]*Bs[5] + Bs[4]*Bs[8], 0), Eq(-As[3]*Bs[1] - As[3]*Bs[4] - As[5]*Bs[1] - As[5]*Bs[4] + Bs[1]*Bs[2] + Bs[1]*Bs[5] + Bs[1]*Bs[8] + Bs[11]*Bs[2] + Bs[11]*Bs[5] + Bs[2]*Bs[4] + Bs[4]*Bs[5] + Bs[4]*Bs[8], 0), Eq(-As[4]*Bs[1] - As[4]*Bs[11] - As[4]*Bs[4] + Bs[1]*Bs[2] + Bs[1]*Bs[5] + Bs[1]*Bs[8] + Bs[11]*Bs[2] + Bs[11]*Bs[5] + Bs[2]*Bs[4] + Bs[4]*Bs[5] + Bs[4]*Bs[8], 0), Eq(Bs[1]*Bs[9] + Bs[11]*Bs[6] + Bs[11]*Bs[9] + Bs[4]*Bs[9], 0), Eq(Bs[1]*Bs[10] + Bs[10]*Bs[11] + Bs[10]*Bs[4] + Bs[11]*Bs[7], 0)],
    [(-Bs[11]/(Bs[1] + Bs[4]))*X1 + (-Bs[11]/(Bs[1] + Bs[4]))*X2 + (0)*X3 + (1)*X4],
    ((As[5]*Bs[1] + As[5]*Bs[4] - Bs[1]*Bs[2] - Bs[1]*Bs[5] - Bs[1]*Bs[8] - Bs[11]*Bs[2] - Bs[11]*Bs[5] - Bs[2]*Bs[4] - Bs[4]*Bs[5] - Bs[4]*Bs[8])/(Bs[1] + Bs[4]))*X3/(X1 + X2 + X3),
    validate=False
)

darboux_23 = DarbouxPolynomial(
    "darboux_23",
    [Eq(As[1] - As[2] + As[5] - Bs[11], 0), Eq(As[1] - As[3], 0), Eq(As[1] - As[4] + As[5] - Bs[11], 0), Eq(-As[1] - As[5] + Bs[11] + Bs[2] + Bs[5] + Bs[8], 0), Eq(Bs[9], 0), Eq(Bs[10], 0)],
    [(1)*X1 + (1)*X2 + (1)*X3 + (-(As[1] - Bs[11])/Bs[11])*X4],
    (-As[1])*X3/(X1 + X2 + X3),
)

darboux_24 = DarbouxPolynomial(
    "darboux_24",
    [Eq(As[1] - As[3], 0), Eq(As[3] + As[5], 0), Eq(Bs[1] + Bs[11] + Bs[4], 0), Eq(Bs[6], 0), Eq(Bs[7], 0), Eq(Bs[8], 0)],
    [(1)*X1 + (1)*X2 + (0)*X3 + (1)*X4],
    (-As[3])*X3/(X1 + X2 + X3),
)

darboux_25 = DarbouxPolynomial(
    "darboux_25",
    [Eq(As[0], 0), Eq(As[2], 0), Eq(Bs[0] + Bs[7], 0), Eq(As[3] + As[5] - Bs[11] - Bs[4], 0), Eq(Bs[2] + Bs[8], 0), Eq(Bs[3] + Bs[9], 0)],
    [((As[3] + As[5])/As[5])*X1 + (1)*X2 + ((As[3] + As[5])/As[5])*X3 + (1)*X4],
    (-As[3])*X3/(X1 + X2 + X3),
)

darboux_26 = DarbouxPolynomial(
    "darboux_26",
    [Eq(As[0], 0), Eq(As[2], 0), Eq(As[3] + As[5] - Bs[11] - Bs[4], 0), Eq(Bs[0] + Bs[7], 0), Eq(Bs[2] + Bs[8], 0), Eq(Bs[3] + Bs[9], 0)],
    [(1)*X1 + (As[5]/(Bs[11] + Bs[4]))*X2 + (1)*X3 + (As[5]/(Bs[11] + Bs[4]))*X4],
    (As[5] - Bs[11] - Bs[4])*X3/(X1 + X2 + X3),
)

darboux_27 = DarbouxPolynomial(
    "darboux_27",
    [Eq(As[0], 0), Eq(As[1], 0), Eq(As[2], 0), Eq(As[3], 0)],
    [(1)*X1 + (1)*X2 + (1)*X3 + (1)*X4],
    (0)*X3/(X1 + X2 + X3),
)

darboux_28 = DarbouxPolynomial(
    "darboux_28",
    [Eq(As[0], 0), Eq(As[2], 0), Eq(Bs[0], 0), Eq(Bs[1], 0), Eq(Bs[2], 0), Eq(Bs[3] + Bs[6] + Bs[9], 0)],
    [(1)*X1 + (0)*X2 + (0)*X3 + (0)*X4],
    (-As[1])*X3/(X1 + X2 + X3),
)

darboux_29 = DarbouxPolynomial(
    "darboux_29",
    [Eq(As[0], 0), Eq(As[1], 0), Eq(As[3], 0), Eq(As[4], 0)],
    [(1)*X1 + (1)*X2 + (1)*X3 + (1)*X4],
    (0)*X3/(X1 + X2 + X3),
)

darboux_30 = DarbouxPolynomial(
    "darboux_30",
    [Eq(As[0], 0), Eq(As[4], 0), Eq(As[1] + As[5] - Bs[1] - Bs[11], 0), Eq(Bs[3] + Bs[6], 0), Eq(Bs[5] + Bs[8], 0), Eq(Bs[0] + Bs[10], 0)],
    [(1)*X1 + ((As[1] + As[5])/As[5])*X2 + ((As[1] + As[5])/As[5])*X3 + (1)*X4],
    (-As[1])*X3/(X1 + X2 + X3),
)

darboux_31 = DarbouxPolynomial(
    "darboux_31",
    [Eq(As[0], 0), Eq(-As[2]*Bs[3] + Bs[2]*Bs[3] + Bs[3]*Bs[5] + Bs[3]*Bs[8] + Bs[5]*Bs[9], 0), Eq(-As[3]*Bs[3] - As[5]*Bs[3] + Bs[2]*Bs[3] + Bs[3]*Bs[5] + Bs[3]*Bs[8] + Bs[5]*Bs[9], 0), Eq(-As[4]*Bs[3] - As[4]*Bs[9] + Bs[2]*Bs[3] + Bs[3]*Bs[5] + Bs[3]*Bs[8] + Bs[5]*Bs[9], 0), Eq(Bs[0]*Bs[9] + Bs[10]*Bs[3] + Bs[10]*Bs[9] + Bs[7]*Bs[9], 0), Eq(Bs[11]*Bs[3] - Bs[4]*Bs[9], 0)],
    [(0)*X1 + (-Bs[9]/Bs[3])*X2 + (0)*X3 + (1)*X4],
    ((As[5]*Bs[3] - Bs[2]*Bs[3] - Bs[3]*Bs[5] - Bs[3]*Bs[8] - Bs[5]*Bs[9])/Bs[3])*X3/(X1 + X2 + X3),
    validate=False
)

darboux_32 = DarbouxPolynomial(
    "darboux_32",
    [Eq(As[0], 0), Eq(-As[1]*Bs[0]*Bs[1] - As[1]*Bs[0]*Bs[4] - As[1]*Bs[1]*Bs[10] - As[1]*Bs[1]*Bs[7] - As[5]*Bs[0]*Bs[1] - As[5]*Bs[0]*Bs[4] - As[5]*Bs[1]*Bs[10] - As[5]*Bs[1]*Bs[7] + Bs[0]*Bs[1]*Bs[2] + Bs[0]*Bs[1]*Bs[5] + Bs[0]*Bs[1]*Bs[8] + Bs[0]*Bs[11]*Bs[2] + Bs[0]*Bs[11]*Bs[5] + Bs[0]*Bs[2]*Bs[4] + Bs[0]*Bs[4]*Bs[5] + Bs[0]*Bs[4]*Bs[8] + Bs[1]*Bs[10]*Bs[2] + Bs[1]*Bs[10]*Bs[8] + Bs[1]*Bs[2]*Bs[7] + Bs[1]*Bs[5]*Bs[7] + Bs[1]*Bs[7]*Bs[8] + Bs[10]*Bs[11]*Bs[2] + Bs[10]*Bs[2]*Bs[4] + Bs[11]*Bs[2]*Bs[7], 0), Eq(-As[2]*Bs[0]*Bs[1] - As[2]*Bs[0]*Bs[11] - As[2]*Bs[0]*Bs[4] - As[2]*Bs[1]*Bs[10] - As[2]*Bs[1]*Bs[7] - As[2]*Bs[10]*Bs[11] - As[2]*Bs[10]*Bs[4] - As[2]*Bs[11]*Bs[7] + Bs[0]*Bs[1]*Bs[2] + Bs[0]*Bs[1]*Bs[5] + Bs[0]*Bs[1]*Bs[8] + Bs[0]*Bs[11]*Bs[2] + Bs[0]*Bs[11]*Bs[5] + Bs[0]*Bs[2]*Bs[4] + Bs[0]*Bs[4]*Bs[5] + Bs[0]*Bs[4]*Bs[8] + Bs[1]*Bs[10]*Bs[2] + Bs[1]*Bs[10]*Bs[8] + Bs[1]*Bs[2]*Bs[7] + Bs[1]*Bs[5]*Bs[7] + Bs[1]*Bs[7]*Bs[8] + Bs[10]*Bs[11]*Bs[2] + Bs[10]*Bs[2]*Bs[4] + Bs[11]*Bs[2]*Bs[7], 0), Eq(-As[3]*Bs[0]*Bs[1] - As[3]*Bs[0]*Bs[4] - As[3]*Bs[1]*Bs[10] - As[3]*Bs[1]*Bs[7] - As[5]*Bs[0]*Bs[1] - As[5]*Bs[0]*Bs[4] - As[5]*Bs[1]*Bs[10] - As[5]*Bs[1]*Bs[7] + Bs[0]*Bs[1]*Bs[2] + Bs[0]*Bs[1]*Bs[5] + Bs[0]*Bs[1]*Bs[8] + Bs[0]*Bs[11]*Bs[2] + Bs[0]*Bs[11]*Bs[5] + Bs[0]*Bs[2]*Bs[4] + Bs[0]*Bs[4]*Bs[5] + Bs[0]*Bs[4]*Bs[8] + Bs[1]*Bs[10]*Bs[2] + Bs[1]*Bs[10]*Bs[8] + Bs[1]*Bs[2]*Bs[7] + Bs[1]*Bs[5]*Bs[7] + Bs[1]*Bs[7]*Bs[8] + Bs[10]*Bs[11]*Bs[2] + Bs[10]*Bs[2]*Bs[4] + Bs[11]*Bs[2]*Bs[7], 0), Eq(-As[4]*Bs[0]*Bs[1] - As[4]*Bs[0]*Bs[11] - As[4]*Bs[0]*Bs[4] - As[4]*Bs[1]*Bs[7] + Bs[0]*Bs[1]*Bs[2] + Bs[0]*Bs[1]*Bs[5] + Bs[0]*Bs[1]*Bs[8] + Bs[0]*Bs[11]*Bs[2] + Bs[0]*Bs[11]*Bs[5] + Bs[0]*Bs[2]*Bs[4] + Bs[0]*Bs[4]*Bs[5] + Bs[0]*Bs[4]*Bs[8] + Bs[1]*Bs[10]*Bs[2] + Bs[1]*Bs[10]*Bs[8] + Bs[1]*Bs[2]*Bs[7] + Bs[1]*Bs[5]*Bs[7] + Bs[1]*Bs[7]*Bs[8] + Bs[10]*Bs[11]*Bs[2] + Bs[10]*Bs[2]*Bs[4] + Bs[11]*Bs[2]*Bs[7], 0), Eq(Bs[0]*Bs[1]*Bs[9] + Bs[0]*Bs[11]*Bs[6] + Bs[0]*Bs[11]*Bs[9] + Bs[0]*Bs[4]*Bs[9] + Bs[1]*Bs[10]*Bs[3] + Bs[1]*Bs[10]*Bs[9] + Bs[1]*Bs[7]*Bs[9] + Bs[10]*Bs[11]*Bs[3] + Bs[10]*Bs[11]*Bs[6] + Bs[10]*Bs[11]*Bs[9] + Bs[10]*Bs[3]*Bs[4] + Bs[10]*Bs[4]*Bs[6] + Bs[10]*Bs[4]*Bs[9] + Bs[11]*Bs[3]*Bs[7] + Bs[11]*Bs[6]*Bs[7] + Bs[11]*Bs[7]*Bs[9], 0)],
    [(-(Bs[0]*Bs[11] + Bs[10]*Bs[11] + Bs[10]*Bs[4] + Bs[11]*Bs[7])/(Bs[0]*Bs[1] + Bs[0]*Bs[4] + Bs[1]*Bs[10] + Bs[1]*Bs[7]))*X1 + (-(Bs[0]*Bs[11] - Bs[1]*Bs[10])/(Bs[0]*Bs[1] + Bs[0]*Bs[4] + Bs[1]*Bs[10] + Bs[1]*Bs[7]))*X2 + (0)*X3 + (1)*X4],
    ((As[5]*Bs[0]*Bs[1] + As[5]*Bs[0]*Bs[4] + As[5]*Bs[1]*Bs[10] + As[5]*Bs[1]*Bs[7] - Bs[0]*Bs[1]*Bs[2] - Bs[0]*Bs[1]*Bs[5] - Bs[0]*Bs[1]*Bs[8] - Bs[0]*Bs[11]*Bs[2] - Bs[0]*Bs[11]*Bs[5] - Bs[0]*Bs[2]*Bs[4] - Bs[0]*Bs[4]*Bs[5] - Bs[0]*Bs[4]*Bs[8] - Bs[1]*Bs[10]*Bs[2] - Bs[1]*Bs[10]*Bs[8] - Bs[1]*Bs[2]*Bs[7] - Bs[1]*Bs[5]*Bs[7] - Bs[1]*Bs[7]*Bs[8] - Bs[10]*Bs[11]*Bs[2] - Bs[10]*Bs[2]*Bs[4] - Bs[11]*Bs[2]*Bs[7])/(Bs[0]*Bs[1] + Bs[0]*Bs[4] + Bs[1]*Bs[10] + Bs[1]*Bs[7]))*X3/(X1 + X2 + X3),
    validate=False
)

darboux_33 = DarbouxPolynomial(
    "darboux_33",
    [Eq(As[0], 0), Eq(As[4], 0), Eq(Bs[3], 0), Eq(Bs[4], 0), Eq(Bs[5], 0), Eq(Bs[0] + Bs[10] + Bs[7], 0)],
    [(0)*X1 + (1)*X2 + (0)*X3 + (0)*X4],
    (-As[3])*X3/(X1 + X2 + X3),
)

darboux_34 = DarbouxPolynomial(
    "darboux_34",
    [Eq(As[1] - As[3], 0), Eq(As[2], 0), Eq(As[4], 0), Eq(Bs[1] + Bs[4], 0), Eq(Bs[2] + Bs[5], 0), Eq(Bs[6] + Bs[9], 0), Eq(Bs[10] + Bs[7], 0)],
    [(1)*X1 + (1)*X2 + (0)*X3 + (0)*X4],
    (-As[3])*X3/(X1 + X2 + X3),
)

darboux_35 = DarbouxPolynomial(
    "darboux_35",
    [Eq(As[1] - As[3], 0), Eq(As[2], 0), Eq(As[4], 0), Eq(Bs[1] + Bs[4], 0), Eq(Bs[2] + Bs[5], 0), Eq(Bs[10] + Bs[7], 0), Eq(Bs[6] + Bs[9], 0)],
    [(1)*X1 + (1)*X2 + (0)*X3 + (0)*X4],
    (-As[3])*X3/(X1 + X2 + X3),
)

darboux_36 = DarbouxPolynomial(
    "darboux_36",
    [Eq(As[0], 0), Eq(As[2], 0), Eq(As[5], 0), Eq(Bs[0] + Bs[7], 0), Eq(-As[3] + Bs[11] + Bs[4], 0), Eq(Bs[2] + Bs[8], 0), Eq(Bs[3] + Bs[9], 0)],
    [(1)*X1 + (0)*X2 + (1)*X3 + (0)*X4],
    (-As[3])*X3/(X1 + X2 + X3),
)

darboux_37 = DarbouxPolynomial(
    "darboux_37",
    [Eq(As[0], 0), Eq(As[4], 0), Eq(As[5], 0), Eq(-As[1] + Bs[1] + Bs[11], 0), Eq(Bs[3] + Bs[6], 0), Eq(Bs[5] + Bs[8], 0), Eq(Bs[0] + Bs[10], 0)],
    [(0)*X1 + (1)*X2 + (1)*X3 + (0)*X4],
    (-As[1])*X3/(X1 + X2 + X3),
)

darboux_38 = DarbouxPolynomial(
    "darboux_38",
    [Eq(As[0], 0), Eq(As[1] - As[3], 0), Eq(As[2], 0), Eq(As[4], 0), Eq(Bs[0]*Bs[1] + Bs[0]*Bs[4] + Bs[1]*Bs[10] + Bs[1]*Bs[7], 0), Eq(Bs[1]*Bs[5] - Bs[2]*Bs[4], 0), Eq(Bs[1]*Bs[3] + Bs[3]*Bs[4] + Bs[4]*Bs[6] + Bs[4]*Bs[9], 0)],
    [(1)*X1 + (-Bs[1]/Bs[4])*X2 + (0)*X3 + (0)*X4],
    (-As[3])*X3/(X1 + X2 + X3),
)

darboux_39 = DarbouxPolynomial(
    "darboux_39",
    [Eq(As[0], 0), Eq(As[2], 0), Eq(As[3] + As[5], 0), Eq(Bs[0] + Bs[7], 0), Eq(Bs[11] + Bs[4], 0), Eq(Bs[2] + Bs[8], 0), Eq(Bs[3] + Bs[9], 0)],
    [(0)*X1 + (1)*X2 + (0)*X3 + (1)*X4],
    (-As[3])*X3/(X1 + X2 + X3),
)

darboux_40 = DarbouxPolynomial(
    "darboux_40",
    [Eq(As[0], 0), Eq(As[1] - As[3], 0), Eq(As[2], 0), Eq(As[3] + As[5], 0), Eq(Bs[0]*Bs[1] + Bs[0]*Bs[11] + Bs[0]*Bs[4] + Bs[1]*Bs[7], 0), Eq(Bs[0]*Bs[8] - Bs[2]*Bs[7], 0), Eq(Bs[0]*Bs[6] + Bs[3]*Bs[7] + Bs[6]*Bs[7] + Bs[7]*Bs[9], 0)],
    [(1)*X1 + (Bs[0]/(Bs[0] + Bs[7]))*X2 + (0)*X3 + (Bs[0]/(Bs[0] + Bs[7]))*X4],
    (-As[3])*X3/(X1 + X2 + X3),
)

darboux_41 = DarbouxPolynomial(
    "darboux_41",
    [Eq(As[0], 0), Eq(As[1] - As[3], 0), Eq(As[4], 0), Eq(As[3] + As[5], 0), Eq(Bs[0]*Bs[1] + Bs[0]*Bs[11] + Bs[0]*Bs[4] + Bs[1]*Bs[10] + Bs[1]*Bs[7] + Bs[10]*Bs[11] + Bs[10]*Bs[4] + Bs[11]*Bs[7], 0), Eq(Bs[1]*Bs[3] + Bs[11]*Bs[3] + Bs[3]*Bs[4] + Bs[4]*Bs[6], 0), Eq(Bs[1]*Bs[5] + Bs[11]*Bs[5] + Bs[4]*Bs[5] + Bs[4]*Bs[8], 0)],
    [(1)*X1 + (-(Bs[1] + Bs[11])/Bs[4])*X2 + (0)*X3 + (1)*X4],
    (-As[3])*X3/(X1 + X2 + X3),
)

for darboux in initialized_darbouxs:
    print(darboux.name,check_darboux(darboux))