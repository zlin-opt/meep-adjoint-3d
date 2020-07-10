import numpy as np

def stepfunc(zz,eta,beta):

    b = beta if beta>1e-2 else 1e-2

    r1 = np.tanh(b*eta) + np.tanh(b*(zz-eta))
    r2 = np.tanh(b*eta) + np.tanh(b*(1.-eta))

    return 1. - np.divide(r1,r2)

def stepgrad(zz,eta,beta):

    b = beta if beta>1e-2 else 1e-2

    sech = np.divide(1.,np.cosh(b*(zz-eta)))
        
    return b * (np.sinh(b*zz)*np.sinh(b-b*zz)/np.sinh(b)) * np.multiply(sech,sech)

def varh_expand(hdof, nx,ny, nlay, mz, beta):

    edof=[]
    
    for i in range(nlay):
        tmp=np.zeros((nx,ny,mz[i]))
        for iz in range(mz[i]):
            zz = float(iz)/float(mz[i])
            tmp[:,:,iz]=stepfunc(zz,hdof[i][:,:],beta)
        edof.append(tmp)

    return edof

def varh_contract(egrad, hdof, nx,ny, nlay, mz, beta):

    hgrad=[]

    for i in range(nlay):
        tmp=np.zeros((nx,ny))
        for iz in range(mz[i]):
            zz = float(iz)/float(mz[i])
            tmp[:,:] += stepgrad(zz,hdof[i][:,:],beta) * egrad[i][:,:,iz]
        hgrad.append(tmp)

    return hgrad

def density_filter(nx,ny, rad,sig, normalized):

    n=nx*ny
    Q=np.identity(n)

    if rad>0.9:

        for iy in range(ny):
            for ix in range(nx):
                i=ix + nx *iy
                jx0=int(ix-np.ceil(rad)+1)
                jx1=int(ix+np.ceil(rad))
                jy0=int(iy-np.ceil(rad)+1)
                jy1=int(iy+np.ceil(rad))
                norm=0.
                for jy in range(jy0,jy1):
                    for jx in range(jx0,jx1):

                        jjx=(jx+nx)%nx
                        jjy=(jy+ny)%ny
                        dist2=float((ix-jx)**2 + (iy-jy)**2)
                        sig2=sig**2
                        j=jjx+nx*jjy
                        wt=np.exp(-dist2/sig2)
                        norm+=wt
                        if dist2<rad**2:
                            Q[i,j]=wt
                        
                if normalized==1:
                    Q[i,:]=Q[i,:]/norm

    return Q

# nx=200
# ny=200
# Q = density_filter(nx,ny, 10,10, 1.)
# v = np.zeros((nx,ny))
# v[nx/2,ny/2]=1.0
# v=np.reshape( np.dot(Q,v.flatten(order='F')), (nx,ny), order='F' )

# import h5py as hp
# fid=hp.File('v.h5','w')
# fid.create_dataset('data',data=v)
# fid.close()


                
        
