import numpy as np
import meep as mp
import solver3d
import filters

fcen=1.0
df=0.1
courant=0.5
lbda=1./fcen
du=0.05
npix=int(lbda/du)
npml=int(npix/2)
nx=6*npix+2*npml
ny=6*npix+2*npml
nz=5*npix+2*npml

#make sure nx, ny and nz are even
if (nx%2) != 0:
    print 'WARNING: nx is not even.\n'
if (ny%2) != 0:
    print 'WARNING: ny is not even.\n'
if (nz%2) != 0:
    print 'WARNING: nz is not even.\n'
             
eps=np.ones((nx,ny,nz))
nlay=2
mz=[20,20]
gap=[3,3]
pml2struct=18
beta=60.
epsbkg=1.0
epsdiff=1.37

z0=[npml+pml2struct]
for i in range(1,nlay):
    z0.append(z0[i-1]+mz[i-1]+gap[i-1])
print 'z0 : {}'.format(z0)
#eps[:,:,z0[0]-3:z0[nlay-1]]=epsbkg+epsdiff

np.random.seed(1234)
hdof=[]
for i in range(nlay):
    hdof.append(np.random.rand(nx,ny))
    hdof[i][0:npml+3,:]=0.01
    hdof[i][:,0:npml+3]=0.01
    hdof[i][nx-npml-3:nx,:]=0.01
    hdof[i][:,ny-npml-3:ny]=0.01

edof=filters.varh_expand(hdof, nx,ny, nlay, mz, beta)
for il in range(nlay):
    eps[:,:,z0[il]:z0[il]+mz[il]] = epsbkg + epsdiff * edof[il][:,:,:]

print 'npix,npml,nx,ny,nz: {},{},{},{},{}'.format(npix,npml,nx,ny,nz)

jx=np.zeros((nx,ny,nz),dtype=complex)
jx[nx/2,ny/2,nz-npml-10]=1./(du*du*du)

jy=np.zeros((nx,ny,nz),dtype=complex)

jz=np.zeros((nx,ny,nz),dtype=complex)

px=np.zeros((nx,ny,nz),dtype=complex)
px[npml+3:nx-npml-3,npml+3:ny-npml-3,npml+10]=1.

py=np.zeros((nx,ny,nz),dtype=complex)
py[npml+3:nx-npml-3,npml+3:ny-npml-3,npml+10]=1.

pz=np.zeros((nx,ny,nz),dtype=complex)
pz[npml+3:nx-npml-3,npml+3:ny-npml-3,npml+10]=1.

niter=200
dh=0.001
ix=nx/2+5
iy=ny/2+4
il=nlay/2
hdof[il][ix,iy]=0.4
for i in range(niter):
    
    hdof[il][ix,iy]+=dh
    edof=filters.varh_expand(hdof, nx,ny, nlay, mz, beta)
    for jl in range(nlay):
        eps[:,:,z0[jl]:z0[jl]+mz[jl]] = epsbkg + epsdiff * edof[jl][:,:,:]

    omega=2.*np.pi*fcen
    ex, ey, ez, dt = solver3d.fdtd(eps,jx,jy,jz,
                                   fcen,df,courant,
                                   nx,ny,nz,du,npml)
    print 'shape of ex, ey, ez: {}, {}, {}'.format(np.shape(ex), np.shape(ey), np.shape(ez))

    
    obj=np.real(np.sum(   np.multiply(px,np.multiply(np.conj(ex),ex)) +
                          np.multiply(py,np.multiply(np.conj(ey),ey)) +
                          np.multiply(pz,np.multiply(np.conj(ez),ez)) ))

    iw = (1.0 - np.exp(-1j*omega*dt)) * (1.0/dt)
    adjx=np.multiply(px,np.conj(ex)) * (df/iw)
    adjy=np.multiply(py,np.conj(ey)) * (df/iw)
    adjz=np.multiply(pz,np.conj(ez)) * (df/iw)
    ux, uy, uz, dt = solver3d.fdtd(eps,adjx,adjy,adjz,
                                   fcen,df,courant,
                                   nx,ny,nz,du,npml)

    w2 = omega*omega
    grad = 2. * np.real( w2 * ( np.multiply(ux,ex)+np.multiply(uy,ey)+np.multiply(uz,ez) ))

    egrad=[]
    for jl in range(nlay):
        egrad.append( epsdiff*grad[:,:,z0[jl]:z0[jl]+mz[jl]] )
    hgrad=filters.varh_contract(egrad, hdof, nx,ny, nlay, mz, beta)
        
    print 'obj: {}, {}, {}'.format(hdof[il][ix,iy],obj,hgrad[il][ix,iy])
                             
