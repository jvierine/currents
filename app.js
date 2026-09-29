import * as THREE from 'three';
import {OrbitControls} from 'three/addons/controls/OrbitControls.js';

const $=s=>document.querySelector(s), viewport=$('#viewport');
const defs={field:['Magnetic field lines','#426477'],r1:['Region 1 FAC','#ff5b63'],r1boundary:['R1 · high-latitude boundary closure','#ff5b63'],r2:['Region 2 FAC','#64dfce'],pedersen:['R1 ↔ R2 Pedersen closure','#f6e8a5'],partial:['Region 2 + partial ring closure','#64dfce'],ring:['Symmetric ring current','#ef779d'],tail:['Tail + high-latitude Chapman-Ferraro return','#73a9ff']};
defs.wedge=['Substorm wedge FAC + tail','#ff835c'];
defs.chapman=['Chapman–Ferraro · dayside shielding','#39cd7d'];
defs.electrojet=['Auroral westward electrojet','#a5f575'];
defs.magnetopause=['Magnetopause surface','#76969f'];
defs.bcbf=['BCBF · plasma flow, not current','#f4f8ff'];
const scene=new THREE.Scene();scene.background=new THREE.Color('#061017');
const camera=new THREE.PerspectiveCamera(40,1,.01,250);camera.up.set(0,0,1);
const renderer=new THREE.WebGLRenderer({antialias:true});renderer.setPixelRatio(Math.min(devicePixelRatio,2));viewport.append(renderer.domElement);
const controls=new OrbitControls(camera,renderer.domElement);controls.enableDamping=true;controls.minDistance=1.4;controls.maxDistance=240;
scene.add(new THREE.AmbientLight(0x92b4c7,1.2));const sun=new THREE.DirectionalLight(0xe8f4ff,2);sun.position.set(20,1,4);scene.add(sun);
const earth=new THREE.Mesh(new THREE.SphereGeometry(1,64,40),new THREE.MeshPhongMaterial({color:0x153c50,shininess:18}));scene.add(earth);
function line(points,color,opacity=.5){const g=new THREE.BufferGeometry().setFromPoints(points.map(p=>new THREE.Vector3(...p)));return new THREE.Line(g,new THREE.LineBasicMaterial({color,transparent:true,opacity}));}
// Coordinate graticule: GSM, not geographic longitude or a geographic map.
for(let lat=-75;lat<=75;lat+=15){const a=THREE.MathUtils.degToRad(lat),pts=[];for(let j=0;j<=100;j++){const t=j/100*2*Math.PI;pts.push([1.004*Math.cos(a)*Math.cos(t),1.004*Math.cos(a)*Math.sin(t),1.004*Math.sin(a)]);}scene.add(line(pts,0x70a5ac,.28));}
for(let lon=0;lon<360;lon+=30){const t=THREE.MathUtils.degToRad(lon),pts=[];for(let j=0;j<=80;j++){const a=-Math.PI/2+j/80*Math.PI;pts.push([1.004*Math.cos(a)*Math.cos(t),1.004*Math.cos(a)*Math.sin(t),1.004*Math.sin(a)]);}scene.add(line(pts,0x70a5ac,.28));}
const labels=[];
function label(text,p,small=false){const el=document.createElement('span');el.className='axis-label'+(small?' small':'');el.textContent=text;viewport.append(el);labels.push({el,p:new THREE.Vector3(...p)});}
for(const [text,p] of [['SUN · +X',[13,0,0]],['DUSK · +Y',[0,13,0]],['DAWN',[0,-13,0]],['NORTH · +Z',[0,0,11]],['TAIL',[-25,0,0]]]){scene.add(line([[0,0,0],p],0x426071,.23));label(text,p);}
for(const r of [5,10,15,20]){const pts=[];for(let j=0;j<=160;j++)pts.push([r*Math.cos(j/160*2*Math.PI),r*Math.sin(j/160*2*Math.PI),0]);scene.add(line(pts,0x294551,.22));label(`${r} Rᴇ`,[-r,0,0],true);}
const groups={},visible={};let animated=[],pickables=[],data,currentPlane,playing=!matchMedia('(prefers-reduced-motion: reduce)').matches,phase=0;
for(const [kind,[name,color]] of Object.entries(defs)){const row=document.createElement('label');row.className='row';row.innerHTML=`<input type="checkbox" checked data-layer="${kind}"><i class="swatch" style="background:${color}"></i><span>${name}</span>`;$('#layers').append(row);visible[kind]=true;row.querySelector('input').addEventListener('change',e=>{visible[kind]=e.target.checked;updateVisibility();});}
function updateVisibility(){for(const [kind,g] of Object.entries(groups)){g.visible=visible[kind];for(const o of g.children)if(o.userData.south)o.visible=$('#south').checked;}}
$('#south').addEventListener('change',updateVisibility);
function seismic(t){
  t=Math.max(-1,Math.min(1,t));
  if(t<0){const u=t+1;return [Math.round(255*u),Math.round(255*u),255];}
  return [255,Math.round(255*(1-t)),Math.round(255*(1-t))];
}
function buildCurrentPlane(c,magnetopause){
  if(currentPlane){scene.remove(currentPlane);currentPlane.geometry.dispose();currentPlane.material.map.dispose();currentPlane.material.dispose();}
  const nx=c.x.length,nz=c.z.length,pixels=new Uint8Array(nx*nz*4),vmax=c.vmax||1;
  const center=magnetopause.center,axes=magnetopause.axes;
  for(let iz=0;iz<nz;iz++)for(let ix=0;ix<nx;ix++){
    const xx=c.x[ix],zz=c.z[iz],k=4*(iz*nx+ix),[r,g,b]=seismic(c.jy[iz][ix]/vmax),radius=Math.hypot(xx,zz);
    const inside=((xx-center[0])/axes[0])**2+(zz/axes[2])**2<=1;
    pixels[k]=r;pixels[k+1]=g;pixels[k+2]=b;pixels[k+3]=inside&&radius>=5?195:0;
  }
  const texture=new THREE.DataTexture(pixels,nx,nz,THREE.RGBAFormat);texture.magFilter=THREE.NearestFilter;texture.minFilter=THREE.NearestFilter;texture.needsUpdate=true;
  const geometry=new THREE.PlaneGeometry(c.x.at(-1)-c.x[0],c.z.at(-1)-c.z[0]);
  currentPlane=new THREE.Mesh(geometry,new THREE.MeshBasicMaterial({map:texture,transparent:true,side:THREE.DoubleSide,depthWrite:false}));
  currentPlane.rotation.x=Math.PI/2;currentPlane.position.set((c.x[0]+c.x.at(-1))/2,0,(c.z[0]+c.z.at(-1))/2);currentPlane.renderOrder=-1;
  currentPlane.visible=$('#current-density').checked;scene.add(currentPlane);
  const f=n=>(n<0?'−':'')+Math.abs(n).toPrecision(2);$('#current-min').textContent=f(-vmax);$('#current-max').textContent=f(vmax);
}
$('#current-density').addEventListener('change',e=>{if(currentPlane)currentPlane.visible=e.target.checked;$('#current-colorbar').hidden=!e.target.checked;if(e.target.checked)view('side');});
function build(index){
  for(const g of Object.values(groups)){scene.remove(g);g.traverse(o=>{o.geometry?.dispose();if(o.material)o.material.dispose();});}animated=[];pickables=[];
  for(const k of Object.keys(defs)){groups[k]=new THREE.Group();scene.add(groups[k]);}
  const p=data.presets[index];$('#params').textContent=`Pdyn ${p.pressure} nPa · Dst ${p.dst} nT\nIMF By ${p.by} / Bz ${p.bz} nT · tilt ${p.tilt}°`;
  if(p.magnetopause){
    const {points,rows,columns}=p.magnetopause,indices=[];
    for(let i=0;i<rows-1;i++)for(let j=0;j<columns-1;j++){
      const a=i*columns+j,b=a+columns;indices.push(a,b,a+1,b,b+1,a+1);
    }
    const geometry=new THREE.BufferGeometry();geometry.setAttribute('position',new THREE.Float32BufferAttribute(points.flat(),3));geometry.setIndex(indices);geometry.computeVertexNormals();
    const shell=new THREE.Mesh(geometry,new THREE.MeshBasicMaterial({color:defs.magnetopause[1],transparent:true,opacity:.055,side:THREE.DoubleSide,depthWrite:false}));
    groups.magnetopause.add(shell);
    for(let j=0;j<columns-1;j+=8)groups.magnetopause.add(line(Array.from({length:rows},(_,i)=>points[i*columns+j]),defs.magnetopause[1],.22));
  }
  for(const path of p.paths){if(!defs[path.kind])continue;const color=defs[path.kind][1],pts=path.points.map(a=>new THREE.Vector3(...a));
    const southern=path.hemisphere===-1||(['r1','r1boundary','r2','pedersen'].includes(path.kind)&&pts[Math.floor(pts.length/2)].z<0);
    const object=line(path.points,color,path.kind==='field'?.30:.9);object.userData={label:path.label,south:southern};groups[path.kind].add(object);
    if(path.kind==='bcbf'){object.material.dispose();object.material=new THREE.LineDashedMaterial({color,dashSize:.022,gapSize:.014});object.computeLineDistances();}
    if(path.kind==='field')continue;pickables.push(object);
    const lengths=[0];for(let i=1;i<pts.length;i++)lengths.push(lengths[i-1]+pts[i].distanceTo(pts[i-1]));
    const ionospheric=['pedersen','electrojet','bcbf'].includes(path.kind);
    const r1=['r1','r1boundary'].includes(path.kind);
    const count=ionospheric?(path.kind==='pedersen'?1:3):Math.max(2,Math.ceil(lengths.at(-1)/(r1?8:5)));
    for(let i=0;i<count;i++){const arrow=new THREE.Mesh(new THREE.ConeGeometry(ionospheric?.018:r1?.24:.11,ionospheric?.07:r1?.8:.36,10),new THREE.MeshBasicMaterial({color}));arrow.userData.south=southern;groups[path.kind].add(arrow);animated.push({arrow,pts,lengths,offset:i/count});}
  }
  // A current sheet is visible as a ribbon between neighboring current
  // paths, as in Ganushkina Figure 7. Width denotes schematic extent only.
  for(const hem of [1,-1])for(const kind of ['r1','r1boundary','chapman'])for(const side of kind==='r1'?['dawn','dusk']:[undefined]){
    const paths=p.paths.filter(q=>q.kind===kind&&(kind==='chapman'?Math.sign(q.points.reduce((s,v)=>s+v[2],0))===hem:q.hemisphere===hem)&&q.side===side).sort((a,b)=>(a.sheet_index??Math.max(...a.points.map(v=>hem*v[2])))-(b.sheet_index??Math.max(...b.points.map(v=>hem*v[2]))));
    const rows=paths.map(q=>{
      const pts=q.points.map(a=>new THREE.Vector3(...a)),s=[0];for(let i=1;i<pts.length;i++)s.push(s.at(-1)+pts[i].distanceTo(pts[i-1]));
      return Array.from({length:201},(_,i)=>{const v=s.at(-1)*i/200;let k=1;while(k<s.length-1&&s[k]<v)k++;return pts[k-1].clone().lerp(pts[k],(v-s[k-1])/(s[k]-s[k-1]));});
    });
    const verts=rows.flat().flatMap(v=>v.toArray()),idx=[];
    for(let i=0;i<rows.length-1;i++)for(let j=0;j<200;j++){const a=i*201+j,b=a+201;idx.push(a,b,a+1,b,b+1,a+1);}
    const g=new THREE.BufferGeometry();g.setAttribute('position',new THREE.Float32BufferAttribute(verts,3));g.setIndex(idx);
    const mesh=new THREE.Mesh(g,new THREE.MeshBasicMaterial({color:defs[kind][1],transparent:true,opacity:.38,side:THREE.DoubleSide,depthWrite:false}));mesh.userData.south=kind!=='chapman'&&hem<0;groups[kind].add(mesh);
  }
  buildCurrentPlane(p.current_density_xz,p.magnetopause);
  updateVisibility();$('#selected').textContent='T96 geometry + illustrative current circuits and plasma flows · click a path';
}
const up=new THREE.Vector3(0,1,0),direction=new THREE.Vector3();
function positionArrow(a){const s=((phase+a.offset)%1)*a.lengths.at(-1);let lo=0,hi=a.lengths.length-1;while(hi-lo>1){const mid=(lo+hi)>>1;if(a.lengths[mid]<s)lo=mid;else hi=mid;}const t=(s-a.lengths[lo])/(a.lengths[hi]-a.lengths[lo]||1);a.arrow.position.lerpVectors(a.pts[lo],a.pts[hi],t);direction.subVectors(a.pts[hi],a.pts[lo]).normalize();a.arrow.quaternion.setFromUnitVectors(up,direction);const scale=Math.min(1,Math.max(.12,a.arrow.position.length()/3));a.arrow.scale.setScalar(scale);}
function view(name){const mobile=innerWidth<760;const settings={global:[[-5,0,0],[24,32,23]],north:[[0,0,.3],[.7,-1.4,4.5]],side:[[-5,0,0],[0,-42,3]],tail:[[-12,0,0],[-37,-26,17]]};const [target,pos]=settings[name];controls.target.set(...target);camera.position.set(...pos);if(name==='global'){const aspect=mobile?innerWidth/innerHeight:Math.max(300,innerWidth-320)/innerHeight;const half=Math.atan(Math.tan(THREE.MathUtils.degToRad(camera.fov/2))*Math.min(1,aspect));camera.position.sub(controls.target).normalize().multiplyScalar(33/Math.sin(half)).add(controls.target);}controls.update();for(const b of document.querySelectorAll('[data-view]'))b.setAttribute('aria-pressed',String(b.dataset.view===name));}
document.querySelectorAll('[data-view]').forEach(b=>b.addEventListener('click',()=>view(b.dataset.view)));$('#reset').onclick=()=>view('global');$('#preset').onchange=e=>build(+e.target.value);
function playLabel(){$('#play').textContent=playing?'Pause arrows':'Play arrows';}playLabel();$('#play').onclick=()=>{playing=!playing;playLabel();};
const ray=new THREE.Raycaster();ray.params.Line.threshold=.12;let down;
renderer.domElement.addEventListener('pointerdown',e=>down=[e.clientX,e.clientY]);renderer.domElement.addEventListener('pointerup',e=>{if(!down||Math.hypot(e.clientX-down[0],e.clientY-down[1])>5)return;const r=renderer.domElement.getBoundingClientRect();ray.setFromCamera(new THREE.Vector2((e.clientX-r.left)/r.width*2-1,-(e.clientY-r.top)/r.height*2+1),camera);const hit=ray.intersectObjects(pickables.filter(o=>o.visible&&o.parent.visible))[0];if(hit)$('#selected').textContent=hit.object.userData.label;});
function resize(){const w=viewport.clientWidth,h=viewport.clientHeight;renderer.setSize(w,h);camera.aspect=w/h;camera.updateProjectionMatrix();}addEventListener('resize',resize);resize();view('global');
let last=performance.now();function frame(now){const dt=Math.min((now-last)/1000,.05);last=now;if(playing)phase=(phase+dt*.045)%1;controls.update();for(const a of animated)positionArrow(a);for(const l of labels){const p=l.p.clone().project(camera);l.el.hidden=p.z>1||p.z< -1||Math.abs(p.x)>1||Math.abs(p.y)>1;l.el.style.left=`${(p.x+1)*viewport.clientWidth/2}px`;l.el.style.top=`${(1-p.y)*viewport.clientHeight/2}px`;}renderer.render(scene,camera);requestAnimationFrame(frame);}requestAnimationFrame(frame);
$('#figure7').onclick=()=>{for(const k of Object.keys(visible)){visible[k]=['r1','r1boundary','chapman','magnetopause'].includes(k);document.querySelector(`[data-layer="${k}"]`).checked=visible[k];}$('#south').checked=false;updateVisibility();view('global');controls.target.set(-5,0,0);camera.position.set(20,72,30);controls.update();$('#selected').textContent='Figure 7 comparison · red Region 1 sheet · green Chapman–Ferraro · northern hemisphere';};
fetch('./traces.json?v=r1-sheet-2').then(r=>{if(!r.ok)throw Error(`Trace data HTTP ${r.status}`);return r.json();}).then(d=>{data=d;build(0);if(new URLSearchParams(location.search).get('view')==='figure7')$('#figure7').click();}).catch(e=>{$('#selected').textContent=`Unable to load traces: ${e.message}`;$('#selected').classList.add('error');console.error(e);});
