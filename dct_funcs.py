import numpy as np

def prep_coef_cos(type_cc, Nd=8):
    coef_cos = np.zeros([Nd,Nd])
    coef_cos.astype(type_cc)
    for i in range(Nd):
        for j in range(Nd):
            coef_cos[i,j] = np.cos(np.pi*i*(2*j+1)/(2*Nd))

    coef_cos.astype(type_cc)

    return coef_cos

def dct_W(coef, u, v, type_k, type_o, N=8):
    coef = coef.astype(type_k)

    a_u = 1./np.sqrt(N*4)
    if u>0:
        a_u=1./np.sqrt(N*2)
    
    a_v = 1./np.sqrt(N*4)
    if v>0:
        a_v=1./np.sqrt(N*2)

    
    aau=np.array([a_u])
    aau=aau.astype(type_k)
    a_u = aau[0]

    aav=np.array([a_v])
    aav=aav.astype(type_k)
    a_v = aav[0]
    
    
    G=np.zeros([N,N,1,1])
    G=G.astype(type_o)

    for i in range(N):
        for j in range(N):            
            G[i,j,0,0] = 4.*a_u*a_v*coef[u,i]*coef[v,j]


    G=G.astype(type_o)

    return G


def idct_W(coef, u, v, type_k, type_o, N=8):
    coef = coef.astype(type_k)

    G=np.zeros([N,N,1,1])
    G=G.astype(type_o)

    for i in range(N):
        a_i=1./np.sqrt(N)
        if i>0:
            a_i=np.sqrt(2./N)
        for j in range(N):
            a_j=1./np.sqrt(N)
            if j>0:
                a_j=np.sqrt(2./N)
            #G = G + a_i*a_j*image_in[i,j]*coef[i,u]*coef[j,v]
            G[i,j,0,0] = a_i*a_j*coef[i,u]*coef[j,v]


    G=G.astype(type_o)

    return G
