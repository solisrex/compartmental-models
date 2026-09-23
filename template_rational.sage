def form_monomial(variablelist,powers):
    if len(variablelist) != len(powers):
        raise Exception()
    return prod(a**b for a,b in zip(variablelist,powers))

t = var('t');

X1 = function('X1')(t);
X2 = function('X2')(t);
X3 = function('X3')(t);

dependent_vars = [X1,X2,X3]

As = [var('A_'+'%d'%i) for i in range(3)]
Bs = [var('B_'+'%d'%i) for i in range(6)]

X1_quad = (-As[0]*X1*X2 - As[1]*X1*X3)/(X1+X2)
X1_lin  = -(Bs[2]+Bs[4])*X1 + Bs[0]*X2 + Bs[1]*X3

X2_quad = (As[0]*X1*X2 - As[2]*X2*X3)/(X1+X2)
X2_lin  = Bs[2]*X1 - (Bs[0]+Bs[5])*X2 + Bs[3]*X3

X3_quad = (As[1]*X1*X3 + As[2]*X2*X3)/(X1+X2)
X3_lin  = Bs[4]*X1 + Bs[5]*X2 - (Bs[1]+Bs[3])*X3

X1_dot = X1_quad + X1_lin
X2_dot = X2_quad + X2_lin
X3_dot = X3_quad + X3_lin

derivatives = [
        diff(X1,t) == X1_dot,
        diff(X2,t) == X2_dot,
        diff(X3,t) == X3_dot]

Cs = [var('C_'+'%d'%i) for i in range(4)]
Ds = [var('D_'+'%d'%i) for i in range(4)]

darboux  = Ds[0] + Ds[1]*X1 + Ds[2]*X2 + Ds[3]*X3;
cofactor = (Cs[0] + Cs[1]*X1 + Cs[2]*X2 + Cs[3]*X3)/(X1+X2);

product = (darboux*cofactor);

darboux_dot = (diff(darboux,t).subs(derivatives))

difference = darboux_dot - product;

difference_numerator = difference.simplify_rational().numerator()

powers = []

for k in [1,2]:
    powers += WeightedIntegerVectors(k, [1, 1, 1]).list()

elim_vars = [X1 == 0, X2 == 0, X3 == 0]

eqns = [difference_numerator.subs(elim_vars) == 0]

for power in powers:
    monomial = form_monomial(dependent_vars,power)

    eqns += [expand(difference_numerator).coefficient(monomial).subs(elim_vars) == 0]
