"""Illustrative magnetopause shielding currents and R1/R2 boundary plasma flows.

The magnetopause geometry is the sigma=1.08 surface in geopack.t96.t96.
K direction uses (B_inside - B_outside) cross n with B_outside=0 and an
inside centered dipole plus T96 dipole shielding. No current magnitude,
conductance, reconnection or substorm dynamics is inferred.
"""
import numpy as np
from scipy.integrate import solve_ivp
from geopack import t96


def boundary_systems(pressure, tilt):
    ps=np.deg2rad(tilt)
    rot=np.array([[np.cos(ps),0,np.sin(ps)],[0,1,0],[-np.sin(ps),0,np.cos(ps)]])
    axis=rot[:,2]
    scale=(pressure/2.)**.14
    center=np.array([(5.48-70)/scale,0.,0.])
    axes=np.array([70*1.08,70*np.sqrt(1.08**2-1),70*np.sqrt(1.08**2-1)])/scale
    theta_max=np.arccos((-25-center[0])/axes[0])
    theta=np.linspace(0,theta_max,40)
    phi=np.linspace(0,2*np.pi,65)
    surface=np.array([center+axes*np.array([np.cos(t),np.sin(t)*np.cos(f),np.sin(t)*np.sin(f)]) for t in theta for f in phi])

    def direction(p):
        n=(p-center)/axes**2;n/=np.linalg.norm(n)
        r=np.linalg.norm(p)
        dipole=30500/r**3*(axis-3*np.dot(axis,p)*p/r**2)
        shielding=np.array(t96.dipshld(ps,*(p*scale)))*scale**3
        current=np.cross(dipole+shielding,n)
        return current/max(np.linalg.norm(current),1e-12)

    paths=[]
    for z in [-10,-7,-4,0,4,7,10]:
        seed=np.array([center[0]+axes[0]*np.sqrt(1-(z/axes[2])**2),0.,z])
        def stop(s,p):return p[0]+25
        stop.terminal=True;stop.direction=-1
        halves=[]
        for sign in [-1,1]:
            departure=np.sign(sign*direction(seed)[1])
            def meridian(s,p):return departure if s<.01 else p[1]
            meridian.terminal=True;meridian.direction=-departure
            sol=solve_ivp(lambda s,p:sign*direction(p),(0,95),seed,events=[stop,meridian],
                          rtol=2e-7,atol=1e-9,max_step=.18)
            halves.append(sol.y.T)
        points=np.concatenate([halves[0][:0:-1],halves[1]])
        paths.append(dict(kind='tail',points=np.round(points,6).tolist(),
            label='Chapman-Ferraro / tail return current on magnetopause',
            quantity='conventional-current',termination='truncated streamline'))

    # Narrow illustrative flow channels, not a solution for E or conductivity.
    # Opposite R1/R2 closure E fields imply opposite zonal ExB drift senses.
    for hemisphere in [1,-1]:
        for start,end,side in [(23,18,'premidnight westward'),(1,6,'postmidnight eastward')]:
            for latitude in [66.3,66.5,66.7]:
                lat=np.deg2rad(hemisphere*latitude)
                lon=(np.linspace(start,end,90)-12)*np.pi/12
                points=1.02*np.column_stack([np.cos(lat)*np.cos(lon),np.cos(lat)*np.sin(lon),np.full(len(lon),np.sin(lat))])@rot.T
                paths.append(dict(kind='bcbf',points=np.round(points,6).tolist(),hemisphere=hemisphere,
                    quantity='plasma-flow',label=f'BCBF: {side} ExB plasma drift, schematic R1/R2 boundary channel'))
    return dict(points=np.round(surface,6).tolist(),rows=len(theta),columns=len(phi),
                model='T96 sigma=1.08 magnetopause; clipped at GSM X=-25 RE',
                center=center.tolist(),axes=axes.tolist()),paths
