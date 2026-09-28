"""Data QA: finite traces, Earth exclusion, FAC alignment, circuit continuity."""
import json
from pathlib import Path
import numpy as np
from geopack import t96
from build_traces import current_density_xz, MU0, RE_METERS

# Analytic sign check: Bx=2z and Bz=-3x gives curl(B)y=+5.
analytic=current_density_xz(lambda p:np.array([2*p[2],0.,-3*p[0]]),
    {'center':[0,0,0],'axes':[100,100,100]})
ix=np.argmin(np.abs(np.asarray(analytic['x'])+10));iz=np.argmin(np.abs(analytic['z']))
assert abs(analytic['jy'][iz][ix]-5/(MU0*RE_METERS))<1e-4

d=json.loads(Path(__file__).with_name('traces.json').read_text())
for preset in d['presets']:
    ps=np.deg2rad(preset['tilt']);m=np.array([np.sin(ps),0,np.cos(ps)])
    par=[preset[k] for k in ['pressure','dst','by','bz']]+[0]*6
    def b(p):
        r=np.linalg.norm(p)
        return 30500/r**3*(m-3*np.dot(m,p)*p/r**2)+np.array(t96.t96(par,ps,*p))
    errors=[]
    for path in preset['paths']:
        p=np.array(path['points']);assert np.isfinite(p).all()
        assert np.linalg.norm(p,axis=1).min()>.999
        if (path['kind'] in ['r1','r2'] and path.get('segment')!='magnetospheric-closure') or path.get('segment') in ['upward','downward']:
            for i in range(0,len(p)-1,8):
                tangent=p[i+1]-p[i];field=b((p[i]+p[i+1])/2);cos=abs(np.dot(tangent,field))/(np.linalg.norm(tangent)*np.linalg.norm(field));errors.append(np.rad2deg(np.arccos(np.clip(cos,-1,1))))
    assert max(errors)<2, max(errors)
    circuit=[p for p in preset['paths'] if p['kind'] not in ['field','ring','tail','chapman','bcbf']]
    # Every schematic circuit, including R1, is continuous.
    starts=np.array([p['points'][0] for p in circuit]);ends=np.array([p['points'][-1] for p in circuit])
    assert not any(p['kind']=='outer' for p in preset['paths'])
    gap=max(np.linalg.norm(starts-end,axis=1).min() for end in ends)
    assert gap<3e-6,gap
    assert max(np.linalg.norm(ends-start,axis=1).min() for start in starts)<3e-6
    for hem in [1,-1]:
        r1paths=[p for p in preset['paths'] if p.get('hemisphere')==hem]
        assert not any(p['kind']=='r1pedersen' for p in r1paths)
        r1={p['side']:np.array(p['points']) for p in r1paths if p['kind']=='r1'}
        r2={p['side']:np.array(p['points']) for p in r1paths if p['kind']=='r2'}
        ped={p['side']:np.array(p['points']) for p in r1paths if p['kind']=='pedersen'}
        boundary=[np.array(p['points']) for p in r1paths if p['kind']=='r1boundary']
        partial=[np.array(p['points']) for p in r1paths if p['kind']=='partial']
        assert set(r1)==set(r2)==set(ped)=={'dawn','dusk'}
        assert len(boundary)==len(partial)==1
        boundary,partial=boundary[0],partial[0]
        # One complete R1-R2 loop, with no separate polar-cap R1 closure.
        sequence=[r1['dusk'],boundary,r1['dawn'],ped['dawn'],r2['dawn'],partial,r2['dusk'],ped['dusk']]
        for first,second in zip(sequence,sequence[1:]+sequence[:1]):
            assert np.linalg.norm(first[-1]-second[0])<3e-6
        assert boundary[0,1]>0 and boundary[-1,1]<0
        assert np.max(boundary[:,0])>5 # dayside high-latitude Figure 7 branch
        assert np.max(hem*boundary[:,2])>8
        wedge={p['segment']:np.array(p['points']) for p in preset['paths'] if p.get('hemisphere')==hem and p.get('segment') in ['upward','downward','electrojet','tail-closure']}
        assert set(wedge)=={'upward','downward','electrojet','tail-closure'}
        up,down,jet,tail=[wedge[k] for k in ['upward','downward','electrojet','tail-closure']]
        for end,start in [(up[-1],tail[0]),(tail[-1],down[0]),(down[-1],jet[0]),(jet[-1],up[0])]:
            assert np.linalg.norm(end-start)<3e-6
        assert np.linalg.norm(up[0])<np.linalg.norm(up[-1])
        assert np.linalg.norm(down[0])>np.linalg.norm(down[-1])
        assert jet[0,1]<0 and jet[-1,1]>0 # postmidnight to premidnight
        assert np.max(np.abs(np.linalg.norm(jet,axis=1)-1.02))<2e-6
        assert np.max(tail[:,0])<0 # nightside, never a dayside shortcut
    boundary=preset['magnetopause'];center=np.array(boundary['center']);axes=np.array(boundary['axes'])
    for path in preset['paths']:
        p=np.array(path['points'])
        if path.get('quantity')=='conventional-current' and 'magnetopause' in path['label']:
            assert path['quantity']=='conventional-current'
            assert np.max(np.abs(np.sum(((p-center)/axes)**2,axis=1)-1))<1e-5
            assert p[:,0].min()>=-25.00001
            nose=np.argmax(p[:,0])
            if 0<nose<len(p)-1:assert p[nose+1,1]>p[nose-1,1], 'dayside current must be dawn to dusk'
        if path['kind']=='bcbf':
            assert path['quantity']=='plasma-flow'
            assert np.max(np.abs(np.linalg.norm(p,axis=1)-1.02))<2e-6
            assert np.sign(p[len(p)//2,2])==path['hemisphere']
        if path['kind']=='r1boundary':
            u=np.sum(((p-center)/axes)**2,axis=1)
            assert np.max(u)<=1.001
            assert np.max(p[:,0])>5 and np.max(path['hemisphere']*p[:,2])>8
    assert sum(p['label'].startswith('Chapman-Ferraro /') for p in preset['paths'])==6
    current=preset['current_density_xz'];x=np.asarray(current['x']);z=np.asarray(current['z']);jy=np.asarray(current['jy'])
    assert np.isfinite(jy).all() and current['units']=='nA/m^2'
    assert (x[0],x[-1],z[0],z[-1])==(-27.,13.,-30.,30.)
    rr=np.hypot(z[:,None],x[None,:]);assert np.max(np.abs(jy[rr<5]))==0
    zeros=preset['magnetopause']['current_zero_crossings']
    assert set(zeros)=={'-1','1'}
    for hem in (-1,1):
        q=np.asarray(zeros[str(hem)]);assert np.sign(q[2])==hem and 5<q[0]<9
    # Cross-tail conventional current is dawn-to-dusk (+Y) near midnight.
    tail=jy[np.argmin(abs(z)),np.argmin(abs(x+15))]
    assert tail>0,(preset['name'],tail)
    print(preset['name'],len(preset['paths']),'paths; max FAC tangent error',round(max(errors),3),'deg; max closure gap',gap,'RE')
