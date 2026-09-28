"""Generate GSM T96 + centered-dipole traces; schematic closures are separate.
Run: conda run -n base python build_traces.py
Scientific arrays: traces.h5. Browser transport: traces.json.
"""
from pathlib import Path
import json
import importlib.metadata
import numpy as np
import h5py
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator
from geopack import t96
from boundary_systems import boundary_systems

ROOT = Path(__file__).resolve().parent
PRESETS = [('quiet', 2., -10., 0., 0., 0.), ('active', 4., -60., 3., -5., 10.)]
RE_METERS = 6_371_000.0
MU0 = 4e-7*np.pi

def current_density_xz(field, magnetopause):
    """Return Jy=curl(B)y/mu0 on GSM Y=0, including the boundary sheet.

    The model field is set to zero outside the T96 sigma=1.08 surface.  Its
    jump at that surface is therefore part of the numerical curl rather than
    being replaced by a disconnected illustrative current.
    """
    x=np.linspace(-24.5,12.0,147);z=np.linspace(-15.0,15.0,121)
    bx=np.zeros((len(z),len(x)));bz=np.zeros_like(bx)
    center=np.asarray(magnetopause['center']);axes=np.asarray(magnetopause['axes'])
    inside=np.zeros_like(bx,dtype=bool)
    for iz,zz in enumerate(z):
        for ix,xx in enumerate(x):
            p=np.array([xx,0.,zz]);r=np.linalg.norm(p)
            within=np.sum(((p-center)/axes)**2)<=1
            if within and r>=1.05:
                b=field(p);bx[iz,ix],bz[iz,ix]=b[0],b[2];inside[iz,ix]=True
    # B is in nT and distance in RE.  nT/RE -> nA/m^2 contributes 1/(mu0*RE).
    jy=(np.gradient(bx,z,axis=0)-np.gradient(bz,x,axis=1))/(MU0*RE_METERS)
    earth=np.hypot(z[:,None],x[None,:])<2.0
    vmax=float(np.percentile(np.abs(jy[~earth]),99.0))
    jy[earth]=0
    return dict(x=np.round(x,5).tolist(),z=np.round(z,5).tolist(),
                jy=np.round(jy,5).tolist(),units='nA/m^2',vmax=round(vmax,5),
                definition='Jy = (dBx/dz - dBz/dx) / mu0; B=0 outside T96 magnetopause')

def generate(name, pressure, dst, by, bz, tilt):
    ps = np.deg2rad(tilt)
    rot = np.array([[np.cos(ps), 0, np.sin(ps)], [0, 1, 0], [-np.sin(ps), 0, np.cos(ps)]])
    axis = rot[:, 2]
    par = [pressure, dst, by, bz, 0, 0, 0, 0, 0, 0]
    def field(p):
        r = np.linalg.norm(p)
        internal = 30500 / r**3 * (axis - 3 * np.dot(axis, p) * p/r**2)
        return internal + np.array(t96.t96(par, ps, *p))
    def seed(lat, mlt, radius=1.02):
        lat, lon = np.deg2rad(lat), (mlt-12)*np.pi/12
        return rot @ (radius*np.array([np.cos(lat)*np.cos(lon), np.cos(lat)*np.sin(lon), np.sin(lat)]))
    def trace(lat, mlt, half=False):
        direction = -np.sign(lat)
        def rhs(s, p):
            b = field(p)
            return direction*b/np.linalg.norm(b)
        def surface(s,p): return np.linalg.norm(p)-1.019
        def limit(s,p): return 25-np.linalg.norm(p)
        def equator(s,p): return np.dot(axis,p)
        for ev in (surface,limit,equator): ev.terminal=True
        surface.direction=limit.direction=-1
        events=[surface,limit,equator] if half else [surface,limit]
        sol=solve_ivp(rhs,[0,85],seed(lat,mlt),events=events,rtol=2e-6,atol=1e-8,max_step=.16)
        p=sol.y.T
        assert np.isfinite(p).all() and len(p)>5
        return p, ('equator' if half and len(sol.t_events[2]) else 'earth' if len(sol.t_events[0]) else 'domain-limit')
    paths=[]
    def add(kind, points, label, **extra):
        paths.append(dict(kind=kind, points=np.round(points,6).tolist(),label=label,**extra))
    for lat in [55,65,72,79]:
        for mlt in range(0,24,3):
            p,stop=trace(lat,mlt)
            add('field',p,'T96 magnetic field line',termination=stop)
    magnetopause,extra=boundary_systems(pressure,tilt)
    bc=np.asarray(magnetopause['center']);ba=np.asarray(magnetopause['axes'])
    def smooth_route(way, samples=241):
        """Round a schematic circuit without moving its physical waypoints."""
        way=np.asarray(way,float)
        s=np.r_[0,np.cumsum(np.linalg.norm(np.diff(way,axis=0),axis=1))]
        ss=np.linspace(0,s[-1],samples)
        route=np.column_stack([PchipInterpolator(s,way[:,i])(ss) for i in range(3)])
        # Interpolation may overshoot very slightly even though every waypoint
        # is inside the T96 boundary.  Keep the drawing within that boundary.
        for ir in range(1,len(route)-1):
            q=(route[ir]-bc)/ba;u=np.linalg.norm(q)
            if u>1:route[ir]=bc+q/u*.995*ba
        route[0]=way[0];route[-1]=way[-1]
        return route
    def mp_point(x,phi,hem,scale=.985):
        """Point just inside the T96 magnetopause at fixed GSM X."""
        section=np.sqrt(max(0.,1-((x-bc[0])/ba[0])**2))
        return np.array([x,scale*ba[1]*section*np.cos(phi),
                         hem*scale*ba[2]*section*np.sin(phi)])
    # R2 joins the partial ring current.  R1 continues through a schematic
    # boundary/plasma-sheet generator route instead of stopping at z=0.
    for hem in [1,-1]:
        for dusk in [18]:
            dawn=24-dusk
            legs={}
            for region,lat in [('r1',70),('r2',63)]:
                for side,mlt in [('dusk',dusk),('dawn',dawn)]:
                    p,stop=trace(hem*lat,mlt,True)
                    assert stop=='equator', (name,hem,region,mlt,stop)
                    legs[region,side]=p
                    upward=(region=='r1')==(side=='dusk')
                    # One representative pair is intentional: both Figure 4
                    # alternatives attach to the same R1 endpoints, and every
                    # displayed FAC remains part of a complete circuit.
                    add(region,p if upward else p[::-1],f'{region.upper()} · {side} · '+('upward' if upward else 'downward')+' conventional current',termination=stop)
            # Smooth radial/polar interpolation in SM equatorial plane.
            for kind,start,end,via in [('partial',legs['r2','dawn'][-1],legs['r2','dusk'][-1],'night')]:
                a,b=rot.T@start,rot.T@end
                th0,th1=np.arctan2(a[1],a[0]),np.arctan2(b[1],b[0])
                if via=='night' and th1>th0: th1-=2*np.pi
                t=np.linspace(0,1,120); th=th0+(th1-th0)*t
                rr=(1-t)*np.linalg.norm(a)+t*np.linalg.norm(b)
                pts=np.column_stack([rr*np.cos(th),rr*np.sin(th),np.zeros(len(t))])@rot.T
                add(kind,pts,'Schematic partial ring-current closure')
            # Figure 4 is drawn as two *alternative complete magnetospheric
            # branches from the same representative R1 sheet*.
            if dusk==18:
                # Dusk R1 is upward (Earth -> magnetosphere), dawn R1 downward.
                start,end=legs['r1','dusk'][-1],legs['r1','dawn'][-1]

                # Short/open branch: each R1 leg reaches the dayside boundary;
                # the connection between them lies directly on the high-latitude
                # magnetopause, as in the small loop at the left of Figure 4.
                direct=[start,mp_point(5.,.72,hem),mp_point(6.,np.pi/2,hem),
                        mp_point(5.,np.pi-.72,hem),end]
                add('r1dayside',smooth_route(direct),
                    'R1 Figure 4 route: direct closure through the dayside magnetopause / solar-wind generator',
                    hemisphere=hem,segment='direct-to-magnetopause')

                # Long/closed branch: current first runs antisunward in the
                # near-equatorial plasma sheet, turns around the far-tail edge,
                # and only then returns earthward on the high-latitude boundary.
                # This is deliberately not a Y-directed bridge at one X.
                plasma_y=.72*start[1]
                far=[start,[-10.,plasma_y,hem*.25],[-18.,plasma_y*.82,hem*.20],
                     mp_point(-24.,0.,hem)]
                far += [mp_point(-24.,p,hem) for p in (.35,.75,1.15,1.55,1.95,2.30)]
                far += [mp_point(x,p,hem) for x,p in [(-20.,2.30),(-16.,2.27),
                                                       (-12.,2.23),(-8.,2.18)]]
                far += [[-6.,.90*end[1],hem*7.],end]
                add('r1tail',smooth_route(far,321),
                    'R1 Figure 4 route: tailward in the plasma sheet, around the far tail, earthward on the magnetopause',
                    hemisphere=hem,segment='far-tail-plasma-sheet-and-boundary-return')

                # Shared ionospheric part of either alternative circuit.
                dawn_foot=legs['r1','dawn'][0];dusk_foot=legs['r1','dusk'][0]
                pole=hem*axis*1.02;cap=[]
                for a,b in [(dawn_foot,pole),(pole,dusk_foot)]:
                    for tt in np.linspace(0,1,45,endpoint=False):
                        q=(1-tt)*a+tt*b;cap.append(1.02*q/np.linalg.norm(q))
                cap.append(dusk_foot)
                add('r1pedersen',cap,'Shared ionospheric closure of the two Figure 4 R1 alternatives: dawn to dusk',
                    hemisphere=hem,segment='shared-ionospheric-closure')
            for mlt,l0,l1 in [(dawn,70,63),(dusk,63,70)]:
                add('pedersen',[seed(hem*l,mlt) for l in np.linspace(l0,l1,35)],'Ionospheric Pedersen closure (schematic)')
    # Equivalent substorm wedge: upward west/premidnight, downward east/
    # postmidnight, westward ionospheric and eastward equatorial closure.
    # T96 only supplies the FAC geometry, not substorm dynamics or amplitudes.
    for hem in [1,-1]:
        west,sw=trace(hem*66,22,True)
        east,se=trace(hem*66,2,True)
        assert sw==se=='equator', (name,hem,sw,se)
        meta=dict(hemisphere=hem)
        add('wedge',west,'Substorm wedge: upward FAC, western / premidnight edge',segment='upward',**meta)
        add('wedge',east[::-1],'Substorm wedge: downward FAC, eastern / postmidnight edge',segment='downward',**meta)
        add('electrojet',[seed(hem*66,h) for h in np.linspace(26,22,100)],
            'Substorm westward auroral electrojet (equivalent current closure)',segment='electrojet',**meta)
        a,b=rot.T@west[-1],rot.T@east[-1]
        t=np.linspace(0,1,120)
        th0=np.arctan2(a[1],a[0]);th1=np.arctan2(b[1],b[0])
        if th1<th0:th1+=2*np.pi
        th=th0+(th1-th0)*t;rr=(1-t)*np.linalg.norm(a)+t*np.linalg.norm(b)
        add('wedge',np.column_stack([rr*np.cos(th),rr*np.sin(th),np.zeros(len(t))])@rot.T,
            'Substorm wedge: eastward equatorial closure (schematic)',segment='tail-closure',**meta)
    for radius in [3.5,4.3]:
        th=np.linspace(0,-2*np.pi,180)
        add('ring',np.column_stack([radius*np.cos(th),radius*np.sin(th),np.zeros(len(th))])@rot.T,'Symmetric westward ring current (schematic)')
    # Cross-tail current: dawn to dusk; its return lies on the same T96
    # magnetopause used by the Chapman-Ferraro shielding streamlines.
    center=np.asarray(magnetopause['center']);axes=np.asarray(magnetopause['axes'])
    for x in [-12,-17,-22]:
        section=np.sqrt(max(0.,1-((x-center[0])/axes[0])**2))
        yr,zr=axes[1]*section,axes[2]*section
        for sign in [-1,1]:
            th=np.linspace(0,np.pi,100)
            cross=np.column_stack([np.full(70,x),np.linspace(-yr,yr,70),np.zeros(70)])
            back=np.column_stack([np.full(100,x),yr*np.cos(th),sign*zr*np.sin(th)])
            add('tail',np.concatenate([cross,back]),
                'Cross-tail current with Chapman-Ferraro magnetopause return')
    paths.extend(extra)
    return dict(name=name,pressure=pressure,dst=dst,by=by,bz=bz,tilt=tilt,
                paths=paths,magnetopause=magnetopause,current_density_xz=current_density_xz(field,magnetopause))

if __name__=='__main__':
    result={'model':'T96 external + 30500 nT centered dipole', 'coordinates':'GSM', 'units':'Earth radii', 'geopack':importlib.metadata.version('geopack'), 'presets':[]}
    for preset in PRESETS:
        print('Tracing',preset[0],flush=True)
        result['presets'].append(generate(*preset))
    with h5py.File(ROOT/'traces.h5','w') as h:
        for key in ['model','coordinates','units','geopack']: h.attrs[key]=result[key]
        h.attrs['generator']='build_traces.py'
        for preset in result['presets']:
            g=h.create_group(preset['name'])
            for key in ['pressure','dst','by','bz','tilt']:g.attrs[key]=preset[key]
            surface=g.create_dataset('magnetopause',data=preset['magnetopause']['points'],compression='gzip')
            for key in ['rows','columns','model','center','axes']:surface.attrs[key]=preset['magnetopause'][key]
            for i,path in enumerate(preset['paths']):
                ds=g.create_dataset(str(i),data=path['points'],compression='gzip')
                for key in ['kind','label']:ds.attrs[key]=path[key]
                for key in ['hemisphere','quantity','segment','termination']:
                    if key in path:ds.attrs[key]=path[key]
            current=preset['current_density_xz']
            ds=g.create_dataset('current_density_xz',data=np.asarray(current['jy']),compression='gzip')
            ds.attrs['x']=current['x'];ds.attrs['z']=current['z'];ds.attrs['units']=current['units']
            ds.attrs['definition']=current['definition'];ds.attrs['vmax']=current['vmax']
    (ROOT/'traces.json').write_text(json.dumps(result,separators=(',',':')))
    print('Saved traces.h5 and traces.json',flush=True)
