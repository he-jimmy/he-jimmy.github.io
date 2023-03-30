#Symbolic check of Lemma 2.15
#variables corresponding to parameters in deformed boson model, with n1 and n2 being the number of arrows leaving
#from the top and he0, he1, ho0, ho1 being the Rogers-Szego polynomials (even and odd), viewed as formal variables.
var('q t v a n1 n2 he0 he1 ho0 ho1')

#boundary vertex weights
def p_diag(a,t,v):
    return (1/a^2-1)/(1/a-v*t)/(1/a+1/v)


def bd_weight(i,j,a,t,v):
    if i==0 and j==0:
        return 1-t*p_diag(a,t,v)
    elif i==0 and j==1:
        return t*p_diag(a,t,v)
    elif i==1 and j==0:
        return p_diag(a,t,v)
    elif i==1 and j==1:
        return 1-p_diag(a,t,v)
    else:
        return 0

#vertex weights in the first deformed boson model
def vtx_weight_1(i,j,m,n,a,q):
    if i+m==j+n:
        if i==0 and j==0:
            return 1
        elif i==0 and j==1:
            return a
        elif i==1 and j==0:
            return 1-q^(m+1)
        elif i==1 and j==1:
            return a
        else:
            return 0
    else:
        return 0

#vertex weights in the second deformed boson model
def vtx_weight_2(i,j,m,n,a,q):
    if i+m==j+n:
        if i==0 and j==0:
            return a
        elif i==0 and j==1:
            return 1
        elif i==1 and j==0:
            return a*(1-q^(m+1))
        elif i==1 and j==1:
            return 1
        else:
            return 0
    else:
        return 0

#odd Rogers-Szego factor, with recurrence applied if m1=n1+1
def H_odd(m,a,q,t,v):
    if m==n1-1:
        return ho0
    elif m==n1:
        return (-t*v)*ho1
    elif m==n1+1:
        return (-t*v)^2*((1-1/(t*v^2))*ho1-1/(t*v^2)*(q^n1-1)*ho0)

#even Rogers-Szego factor, with recurrence applied if m2=n2+1
def H_even(m,a,q,t,v):
    if m==n2-1:
        return he0
    elif m==n2:
        return he1
    elif m==n2+1:
        return ((1-t)*he1-t*(q^n2-1)*he0)

#Left hand side of equality
def part_fn_1(i,j):
    return sum([bd_weight(i,k,a,t,v)*vtx_weight_1(k,l,m2,n2,a,q)*vtx_weight_1(l,j,m1,n1,a,q)*H_odd(m1,a,q,t,v)*H_even(m2,a,q,t,v) for k in range(0,2) for l in range(0,2) for m1 in [n1-1,n1,n1+1] for m2 in [n2-1,n2,n2+1]])


#Right hand side of equality
def part_fn_2(i,j):
    return sum([vtx_weight_2(i,k,m2,n2,a,q)*vtx_weight_2(k,l,m1,n1,a,q)*bd_weight(l,j,a,t,v)*H_odd(m1,a,q,t,v)*H_even(m2,a,q,t,v) for k in range(0,2) for l in range(0,2) for m1 in [n1-1,n1,n1+1] for m2 in [n2-1,n2,n2+1]])

#Check that the four cases, corresponding to i,j=0,1, all evaluate to 0
(part_fn_1(0,0)-part_fn_2(0,0)).full_simplify()
(part_fn_1(0,1)-part_fn_2(0,1)).full_simplify()
(part_fn_1(1,0)-part_fn_2(1,0)).full_simplify()
(part_fn_1(1,1)-part_fn_2(1,1)).full_simplify()

