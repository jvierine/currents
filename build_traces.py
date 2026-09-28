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
    # Cover the full displayed T96 boundary/tail scene at 0.25 RE spacing.
    x=np.linspace(-27.0,13.0,161);z=np.linspace(-30.0,30.0,241)
    bx=np.zeros((len(z),len(x)));bz=np.zeros_like(bx)
    center=np.asarray(magnetopause['center']);axes=np.asarray(magnetopause['axes'])
    inside=np.zeros_like(bx,dtype=bool)
    for iz,zz in enumerate(z):
        for ix,xx in enumerate(x):
            p=np.array([xx,0.,zz]);r=np.linalg.norm(p)
            within=np.sum(((p-center)/axes)**2)<=1
            if within and r>=5.0:
                b=field(p);bx[iz,ix],bz[iz,ix]=b[0],b[2];inside[iz,ix]=True
    # B is in nT and distance in RE.  nT/RE -> nA/m^2 contributes 1/(mu0*RE).
    jy=(np.gradient(bx,z,axis=0)-np.gradient(bz,x,axis=1))/(MU0*RE_METERS)
    # The internal centered-dipole/IGRF-like contribution is not a
    # magnetospheric current system.  Exclude the near-Earth volume entirely.
    earth=np.hypot(z[:,None],x[None,:])<5.0
    signal=(np.abs(jy)>0)&~earth
    vmax=float(np.percentile(np.abs(jy[signal]),99.0))
    jy[earth]=0
    return dict(x=np.round(x,5).tolist(),z=np.round(z,5).tolist(),
                jy=np.round(jy,5).tolist(),units='nA/m^2',vmax=round(vmax,5),
                definition='Jy = (dBx/dz - dBz/dx) / mu0; B=0 outside T96 magnetopause; r<5 RE excluded')

def magnetopause_current_zero_crossings(current, magnetopause):
    """Locate dayside Chapman-Ferraro Jy reversals on the X-Z cut."""
    x=np.asarray(current['x']);z=np.asarray(current['z']);jy=np.asarray(current['jy'])
    center=np.asarray(magnetopause['center']);axes=np.asarray(magnetopause['axes'])
    found={}
    for hem in (-1,1):
        keep=(hem*z>=4.)&(hem*z<=min(15.,.9*axes[2]))
        zz=z[keep];values=[]
        for zi in zz:
            xb=center[0]+axes[0]*np.sqrt(max(0.,1-(zi/axes[2])**2))
            # Sample just inside the boundary sheet; interpolation avoids
            # locking the zero to a 0.25 RE pixel center.
            values.append(np.interp(xb-.20,x,jy[np.argmin(np.abs(z-zi))]))
        values=np.asarray(values)
        crossings=np.where(values[:-1]*values[1:]<=0)[0]
        assert len(crossings), (hem,'no magnetopause Jy reversal')
        i=crossings[np.argmin(np.abs(zz[crossings]-hem*9.))]
        f=-values[i]/(values[i+1]-values[i])
        z0=zz[i]+f*(zz[i+1]-zz[i])
        x0=center[0]+axes[0]*np.sqrt(max(0.,1-(z0/axes[2])**2))
        found[hem]=np.array([x0,0.,z0])
    return found

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
    # Build geometry first, then use the same X-Z Jy diagnostic shown in the
    # browser to anchor the Chapman-Ferraro cusp reversal points.
    magnetopause,_=boundary_systems(pressure,tilt,geometry_only=True)
    current_xz=current_density_xz(field,magnetopause)
    cusp_points=magnetopause_current_zero_crossings(current_xz,magnetopause)
    magnetopause,extra=boundary_systems(pressure,tilt,cusp_points=cusp_points)
    magnetopause['current_zero_crossings']={str(h):np.round(p,6).tolist() for h,p in cusp_points.items()}
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
                    # R1 is drawn below as one continuous far-tail circuit that
                    # reaches the ionosphere itself.  R2 retains its T96 FACs.
                    if region=='r2':
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
            # Retain the Figure 4 far-tail route only.  It is a single complete
            # R1 path between ionospheric footpoints, with no separate FAC seam
            # at the magnetic equator and no imposed dawn/dusk displacement in
            # the magnetosphere.  Small mirrored Y offsets near Earth make the
            # two ionospheric connections readable; both sweep tailward equally.
            if dusk==18:
                dusk_foot=legs['r1','dusk'][0]
                dawn_foot=legs['r1','dawn'][0]
                def top(x,scale=.985): return mp_point(x,np.pi/2,hem,scale)
                far=[dusk_foot,[-1.8,.45*dusk_foot[1],hem*2.2],
                     [-4.5,0.,hem*.85],[-10.,0.,hem*.22],[-18.,0.,hem*.12],
                     [-24.,0.,0.],[-24.,0.,hem*6.],top(-24.)]
                far += [top(x) for x in (-20.,-16.,-12.,-8.,-4.)]
                far += [[-2.8,0.,hem*4.2],[-1.8,.45*dawn_foot[1],hem*2.2],dawn_foot]
                route=smooth_route(far,321)
                route[route[:,0]<-3.,1]=0. # exact noon-midnight symmetry outside the entry bends
                add('r1tail',route,
                    'R1 Figure 4 far-tail route: ionosphere to plasma sheet, around the far tail, boundary return to ionosphere',
                    hemisphere=hem,segment='ionosphere-plasma-sheet-far-tail-boundary-ionosphere')

                # Ionospheric closure of the single far-tail R1 circuit.
                pole=hem*axis*1.02;cap=[]
                for a,b in [(dawn_foot,pole),(pole,dusk_foot)]:
                    for tt in np.linspace(0,1,45,endpoint=False):
                        q=(1-tt)*a+tt*b;cap.append(1.02*q/np.linalg.norm(q))
                cap.append(dusk_foot)
                add('r1pedersen',cap,'R1 ionospheric closure: dawn to dusk',
                    hemisphere=hem,segment='ionospheric-closure')
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
                paths=paths,magnetopause=magnetopause,current_density_xz=current_xz)

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
