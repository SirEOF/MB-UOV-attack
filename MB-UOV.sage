sage_server.MAX_STDOUT_SIZE=sage_server.MAX_STDOUT_SIZE*10

#---------------------basic parameters
q=11
v=20
d=2
op=5
o=d*op
n=v+o
K.<a>=GF(q)
P=PolynomialRing(K,'x',n)
x=P.gens()
x_vec=vector(x)



#---------------------generating public key
A00=[0 for i in range(o)]
A01=[0 for i in range(o)]
B0=[0 for i in range(o)]
B=[[] for i in range(o)]
A=[0 for i in range(o)]
D=[0 for i in range(o)]
D[0]=random_matrix(K,v,v)

rotate_matrix=[]
I=matrix.identity(K,v,v)
rotate_matrix.append(I.row(v-1).list())
for i in range(v-1):
    rotate_matrix.append(I.row(i).list())
rotate_matrix=matrix(K,v,v,rotate_matrix)



for i in range(o):#-------------------B[i] is vv part
    D[i]=rotate_matrix^i*D[0]

Ar=[0 for i in range(op)]
for r in range(op):
    A01[r]=[]
    for i in range(v*op):
        A01[r].append(K.random_element())
    Ar[r]=matrix(K,op,v,A01[r]).transpose()
    for i in range(v*(o-op)):
        A01[r].append(K(0))
    A01[r]=matrix(K,o,v,A01[r]).transpose()

rotate_matrix=[]
for i in range(op):
    for j in range(o):
        if(j==(o-op+i)):
            rotate_matrix.append(K(1))
        else:
            rotate_matrix.append(K(0))
for i in range(o-op):
    for j in range(o):
        if(j==i):
            rotate_matrix.append(K(1))
        else:
            rotate_matrix.append(K(0))
rotate_matrix=matrix(K,o,o,rotate_matrix).transpose()
for i in range(o):
    A01[i]=A01[i%op]*rotate_matrix^(i//op)

for i in range(o):
    A00[i]=D[i]
for i in range(o):
    A[i]=matrix.block([[A00[i],A01[i]],[matrix(K,o,v),matrix(K,o,o)]])

#---------------omit linear and constant parts for simplicity

F=[0 for i in range(o)]
for i in range(o):
    F[i]=A[i]#-------------------central quadratic matrix
while true:
    T=random_matrix(K,n)
    if T.is_invertible():
        break
Tt=T.transpose()
P=[0 for i in range(o)]
for i in range(o):
    P[i]=Tt*F[i]*T#------------public key of MB_UOV

for i in range(o):
    P[i]=P[i]+P[i].transpose()

#--------------------------------------Constructing an equivalent key
VV=VectorSpace(K,n)
TO=matrix(K,o,n)
TO.set_block(0,0,P[0].kernel().matrix())
TO.set_block(op-1,0,P[op+1].kernel().matrix())
TO=VV.span(TO).matrix()
BB=VV.span(TO).complement().matrix()
TT=block_matrix(2,1,[BB,TO])


print "An equivalent key:"
for i in range(o):
    print (TT*P[i]*TT.transpose()).str()
    print " "
