import numpy as np
import meep as mp

def dat2pos(r,data,nx,ny,nz,du,val):

    decplace=5

    x0=np.around(nx*du/2.,decplace)
    y0=np.around(ny*du/2.,decplace)
    z0=np.around(nz*du/2.,decplace)
    rx=np.around(r.x,decplace)
    ry=np.around(r.y,decplace)
    rz=np.around(r.z,decplace)
    jx=np.around( (rx+x0)/du, decplace)
    jy=np.around( (ry+y0)/du, decplace)
    jz=np.around( (rz+z0)/du, decplace)
    ix = int(np.floor( jx ))
    iy = int(np.floor( jy ))
    iz = int(np.floor( jz ))
    
    if 0<=ix<nx and 0<=iy<ny and 0<=iz<nz:
        return data[ix,iy,iz]
    else:
        return val

def fdtd(eps,jx,jy,jz,
         fcen,df,courant,
         nx,ny,nz,du,npml):

    dpml=npml*du
    Lx=nx*du
    Ly=ny*du
    Lz=nz*du

    resolution=1./du

    mtr_dt=10./fcen
    mtr_c=mp.Ex
    mtr_r=mp.Vector3(0,0,0)
    mtr_tol=1e-8
    
    def epsfun(r):
        return dat2pos(r,eps,nx,ny,nz,du,1.)

    def jxfun(r):
        return dat2pos(r,jx,nx,ny,nz,du,0.)

    def jyfun(r):
        return dat2pos(r,jy,nx,ny,nz,du,0.)

    def jzfun(r):
        return dat2pos(r,jz,nx,ny,nz,du,0.)

    cell = mp.Vector3(Lx,Ly,Lz)
    pml_layers = [mp.PML(dpml)]
    sources = [ mp.Source(mp.GaussianSource(fcen,fwidth=df),
                          component=mp.Ex,
                          center=mp.Vector3(0,0,0),
                          size=cell,
                          amp_func=jxfun),
                mp.Source(mp.GaussianSource(fcen,fwidth=df),
                          component=mp.Ey,
                          center=mp.Vector3(0,0,0),
                          size=cell,
                          amp_func=jyfun),
                mp.Source(mp.GaussianSource(fcen,fwidth=df),
                          component=mp.Ez,
                          center=mp.Vector3(0,0,0),
                          size=cell,
                          amp_func=jzfun) ]
    sim = mp.Simulation(cell_size=cell,
                        boundary_layers=pml_layers,
                        epsilon_func=epsfun,
                        eps_averaging=False,
                        sources=sources,
                        resolution=resolution,
                        force_complex_fields=True,
                        Courant=courant)

    sim.init_sim()

    yeeL=sim.gv.little_corner()
    yeeU=sim.gv.big_corner()
    print 'Yee lower left coords: {}'.format([yeeL.x()*du/2.,yeeL.y()*du/2.,yeeL.z()*du/2.])
    print 'Yee upper right coords: {}'.format([yeeU.x()*du/2.,yeeU.y()*du/2.,yeeU.z()*du/2.])
             
    dft_vol = mp.Volume(center=mp.Vector3(0,0,0), size=cell)
    dft_obj = sim.add_dft_fields([mp.Ex,mp.Ey,mp.Ez], fcen, fcen,1, where=dft_vol, yee_grid=True)

    sim.run( until_after_sources=mp.stop_when_fields_decayed(mtr_dt, mtr_c, mtr_r, mtr_tol) )

    ex = sim.get_dft_array(dft_obj, mp.Ex, 0)
    ey = sim.get_dft_array(dft_obj, mp.Ey, 0)
    ez = sim.get_dft_array(dft_obj, mp.Ez, 0)

    dt = sim.fields.dt
    return ex, ey, ez, dt

