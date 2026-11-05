import sys; import pathlib; S=str(pathlib.Path(__file__).resolve().parent); sys.path.insert(0,S); sys.path.insert(0,str(pathlib.Path(S).parents[3]/'siteTools'))
from cam import *
import avasite as b, geopandas as gpd
from PIL import Image, ImageDraw
from shapely.geometry import shape
from shapely.ops import unary_union
from rasterio.features import shapes, rasterize
from rasterio.transform import from_origin
from scipy import ndimage
POLY=[(495,240),(503,270),(517,300),(537,327),(547,343),(550,370),(557,400),(567,420),(577,440),(587,465),(600,507),(610,540),(630,573),(660,600),(677,613),(693,623),(703,640),(708,667),(707,683),(670,685),(640,667),(620,643),(613,627),(603,590),(597,550),(588,505),(575,470),(560,430),(540,400),(517,373),(500,347),(483,333),(480,323),(463,320),(462,300),(463,287),(470,270)]
def cast(prm):
    X0,Y0,Z0,yaw,pitch,roll,f=prm
    m=Image.new('L',(W,H),0); ImageDraw.Draw(m).polygon(POLY,fill=1); m=np.asarray(m).astype(bool)
    v,u=np.nonzero(m); u=u+.5; v=v+.5
    D=np.c_[(u-W/2)/f,-(v-H/2)/f,np.ones_like(u)]@R(yaw,pitch,roll); D/=np.linalg.norm(D,axis=1)[:,None]
    C=np.array([X0,Y0,Z0]); hit=np.full((len(u),3),np.nan); t=np.arange(800,8000,2.0)
    for i in range(0,len(u),2000):
        d=D[i:i+2000]; P=C+t[None,:,None]*d[:,None,:]
        zt=zAt(P[...,0].ravel(),P[...,1].ravel()).reshape(P.shape[:2]); below=P[...,2]<zt; k=np.argmax(below,1); ok=below.any(1)
        hit[i:i+2000][ok]=P[np.arange(len(d))[ok],k[ok]]
    ok=~np.isnan(hit[:,0]); return np.c_[hit[ok],v[ok]]
def topoly(P,c=5):
    x0,y1=P[:,0].min()-20,P[:,1].max()+20; nx=int((P[:,0].max()+20-x0)/c)+1; ny=int((y1-P[:,1].min()+20)/c)+1
    g=np.zeros((ny,nx),'uint8'); g[((y1-P[:,1])/c).astype(int),((P[:,0]-x0)/c).astype(int)]=1
    g=ndimage.binary_fill_holes(ndimage.binary_closing(g,iterations=2)).astype('uint8')
    gg=unary_union([shape(s) for s,v in shapes(g,transform=from_origin(x0,y1,c,c)) if v==1])
    return (max(gg.geoms,key=lambda q:q.area) if gg.geom_type=='MultiPolygon' else gg).simplify(3)
def stats(g, z, T, slope, asp):
    m=rasterize([g],out_shape=z.shape,transform=T).astype(bool); a=np.radians(asp[m])
    return dict(area=g.area, zmin=np.nanmin(z[m]), zmax=np.nanmax(z[m]), slope=np.median(slope[m]),
                aspect=(np.degrees(np.arctan2(np.sin(a).mean(),np.cos(a).mean()))+360)%360)
