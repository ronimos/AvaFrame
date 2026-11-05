import json, numpy as np, rasterio
from scipy.ndimage import map_coordinates
from shapely.geometry import LineString
from shapely.ops import linemerge
import pathlib; S=str(pathlib.Path(__file__).resolve().parent)
W,H=1300,975
dem=rasterio.open(f'{S}/usgs10.tif'); Z=dem.read(1); DT=dem.transform
def zAt(x,y):
    c,r=~DT*(np.asarray(x),np.asarray(y)); return map_coordinates(Z,[r-.5,c-.5],order=1,mode='nearest')
# OSM CO 133 as one line, split at the hairpin (southernmost point) into the west (upper) and east (lower) legs
segs=json.load(open(f'{S}/co133_utm.json'))['co133']
road=linemerge([LineString(s) for s in segs])
if road.geom_type!='LineString': road=max(road.geoms,key=lambda g:g.length)
xy=np.array(road.coords); k=np.argmin(xy[:,1])
legA,legB=xy[:k+1],xy[k:]
upper,lower=(legA,legB) if legA[:,0].mean()<legB[:,0].mean() else (legB,legA)
def dens(a,step=5):
    L=LineString(a); d=np.arange(0,L.length,step); p=np.array([L.interpolate(s).coords[0] for s in d]); return p
def clip(p): return p[(p[:,1]<4333200)&(p[:,1]>4330300)]
UP=clip(dens(upper)); LO=clip(dens(lower))
UP3=np.c_[UP,zAt(*UP.T)+1]; LO3=np.c_[LO,zAt(*LO.T)+1]
# traced pixels (x right, y down) - upper leg incl. hairpin top, and lower leg incl. hairpin bottom
tUP=np.array([(450,700),(500,690),(550,660),(600,640),(650,622),(700,600),(750,585),(800,570),(850,555),(900,540),(950,525),(1000,515),(1050,502),(1100,490),(1150,470),(1200,450),(1250,435),(1300,422),
              (440,705),(425,720),(410,735),(390,750),(370,760),(340,770),(310,775),(290,778)],float)
tLO=np.array([(275,790),(275,810),(280,830),(300,843),(330,848),(360,848),(400,843),(430,835),(460,822),(490,808),(510,795),(530,780),(550,765),(570,748),(600,738),(650,738),(700,742),(750,730),
              (800,712),(850,710),(900,713),(940,722),(1000,730),(1050,730),(1100,728),(1150,728),(1200,725),(1250,725)],float)
def R(yaw,pitch,roll):
    # camera looks along azimuth yaw (deg from north, clockwise), tilted down by pitch, rolled by roll
    a,p,r=np.radians([yaw,pitch,roll])
    fwd=np.array([np.sin(a)*np.cos(p),np.cos(a)*np.cos(p),-np.sin(p)])
    right=np.array([np.cos(a),-np.sin(a),0.0])
    up=np.cross(right,fwd)
    right,up=np.cos(r)*right+np.sin(r)*up, -np.sin(r)*right+np.cos(r)*up
    return np.vstack([right,up,fwd])
def project(P,prm):
    X,Y,Zc,yaw,pitch,roll,f=prm
    q=(P-np.array([X,Y,Zc]))@R(yaw,pitch,roll).T
    u=W/2+f*q[:,0]/q[:,2]; v=H/2-f*q[:,1]/q[:,2]
    return np.c_[u,v],q[:,2]
def segdist(pts,poly):
    a=poly[:-1][None]; ab=(poly[1:]-poly[:-1])[None]; p=pts[:,None]
    t=np.clip(((p-a)*ab).sum(2)/np.maximum((ab*ab).sum(2),1e-9),0,1)
    d=np.hypot(*np.moveaxis(a+ab*t[...,None]-p,2,0)); return d.min(1)
def resid(prm):
    pu,du=project(UP3,prm); pl,dl=project(LO3,prm)
    if (du<=0).any() or (dl<=0).any(): return np.full(len(tUP)+len(tLO),1e4)
    return np.r_[segdist(tUP,pu),segdist(tLO,pl)]
