#!/usr/bin/env python
"""
 shapefile workaround and network analysis in SciPython
"""
__author__ = """Antoine Lucas (alucas@caltech.edu), Caltech 2011"""
#-----------------------------------------------------------
from os import *
import os
#import Sci/Num
import numpy as np
from numpy import *
import scipy as Sci
import scipy.linalg
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, show
from matplotlib.patches import Ellipse
from pylab import *

#import Geo/GDAL/NetWorkX
import shapefile
from  pylab import *
from numpy import linalg as LA
import networkx as nx
from networkx.algorithms import bipartite
import nx_spatial as ns
import pygraphviz as pgv
import osgeo.gdal 
import osgeo.ogr
import osgeo.osr
import osgeo.gdalnumeric
import osgeo.gdalconst
import sys
from osgeo import ogr
import random


#from shapely.geometry import Point
from shapely.geometry import MultiPoint
#from descartes.patch import PolygonPatch


#-----------------------------------------------------------
#USER PARAMETERS
usrhome=os.environ['HOME']

shpfile=usrhome + '/Documents/NetWorkX/db/fluvial/hbetts_densified.shp'
shpfile=usrhome + '/Documents/NetWorkX/db/fluvial/allchannels_densified_v1.shp'
RTitan = 2575.0;    # Titan's radius [km]
pixres = 0.351;     # Pixel Resolution [km] (c.f. SAR in the ArcGIS Project)
radius = 4*pixres;  # Threshold for the search radius [km]
dist_min=radius/RTitan; # search radius [radians]


print "Working on", shpfile
sf = shapefile.Reader(shpfile)
shapes = sf.shapes()
N_channels = len(shapes)

N_channels = 100
channels = ((xrange(N_channels)))


#-----------------------------------------------------------
# FUNCTIONS

# Computes great arc distance on a sphere
def gaa(lon1,lat1,lon2,lat2):   
     lat1=(pi/2)-lat1; lat2=(pi/2)-lat2
     EE=2*pi
     aa1 = multiply(sin(EE-lat1),cos(lon1)) 
     aa2 = multiply(sin(EE-lat2),cos(lon2))
     aa  = multiply(aa1,aa2)
     bb1 = multiply(sin(EE-lat1),sin(lon1))
     bb2 = multiply(sin(EE-lat2),sin(lon2))
     bb  = multiply(bb1,bb2)
     cc  = multiply(cos(EE-lat1),cos((EE-lat2)))
     dradians = arccos( aa + bb + cc)
     dradians = real(dradians)
     nani = (argwhere(isnan(dradians)))
     
     if (len(nani)>0):
         nani = nani[0,0,0]
         if (lon1[nani] == lon2): dradians[nani] = 0.0
         if (lon1[nani] != lon2): dradians[nani] = 1e+308
     return dradians





def lbtoangle(l1,b1,l2,b2):
     p1 = cos((l1-l2))
     p2 = cos((b1-b2))
     p3 = cos((b1+b2))
     return (acos(((p1*(p2+p3))+(p2-p3))/2))

def interpgeo(lat,lon,maxdist):
    
    dlat = abs(diff(lat))
    dlon = abs(diff(lon))
    dist = zeros(len(lon) -1)  + maxdist

    for i in arange(len(lon)-1):
     dist[i] = gaa(lon[i],lat[i],lon[i+1],lat[i+1])  

    indx  = argwhere( dist > maxdist);
  
    print len(indx), len(dist)

    if (len(indx)>0):
       steps = ceil(dist[indx]/maxdist);      #%  No points added each location
       totalpts = sum(steps)-len(steps);   #%  Total points to be added
       lastpt  = len(lat);                 #%  Current last point in data set
       lato = zeros(len(lat)+totalpts)
       lono = zeros(len(lon)+totalpts)
       print steps, totalpts       


    return lon,lat



def gcd(long_1,lat_1,long_2,lat_2):
        import math  
#a = sin((lat2-lat1)/2).^2 + cos(lat1) .* cos(lat2) .* sin((lon2-lon1)/2).^2;
#rng = 2 * atan2(sqrt(a),sqrt(1 - a));
        dlong = array(long_2 - long_1)
        dlat = array(lat_2 - lat_1)
        a = array(multiply(sin(dlat/2),sin(dlat/2)))
        b = array(multiply(cos(lat_1),cos(lat_2)))
        c = array(multiply(sin(dlong/2),sin(dlong/2)))    
        c = array(multiply(b,c))  
        dist = zeros(len(c))
        for i in xrange(0,len(c)):
           dist[i] = 2 * math.atan2(sqrt(c[i]),sqrt(1-c[i]))  
        return vstack(dist)


def min_dist(n1,n2,loc):
    d=gaa(loc[n1,0],loc[n1,1],loc[n2,0],loc[n2,1])
    return  d

def getfieldinfo(lyr, feature, flds):
    f = feature
    return [f.GetField(f.GetFieldIndex(x)) for x in flds]

def restart_line():
    sys.stdout.write('\r')
    sys.stdout.flush()

def wraptopi(lon):
   bb = pi*((abs(lon)/pi) - 2*ceil(((abs(lon)/pi)-1)/2))
   lon = multiply(bb,sign(lon))
   return lon

def NorthAzimuth(x1,y1,x2,y2):
     degBearing = math.degrees(math.atan2((x2 - x1),(y2 - y1)))
      if (degBearing < 0):
          degBearing += 360.0
     return degBearing



def gcaz(lat1,lon1,lat2,lon2):
    import math
    aa = array(multiply(cos(lat2),sin(lon2-lon1)))
    bb = array(multiply(cos(lat1),sin(lat2)))
    cc = array(multiply( sin(lat1) , cos(lat2) , cos(lon2-lon1)))
    bb = bb - cc
    az = zeros(len(aa))
    
    for i in xrange(0,len(aa)):
        az[i] = math.atan2(aa[i],bb[i])

    ll = argwhere(lat1 <= -pi/2)
    if (ll.size>0): az[ll] = 0
    ll = argwhere(lat1 >= pi/2)
    if (ll.size>0): az[ll] = 0
    ll = argwhere(lat2 <= -pi/2)
    if (ll.size>0): az[ll] = pi
    ll = argwhere(lat2 >= pi/2)
    if (ll.size>0): az[ll] = pi   

    return vstack(az)



def poly2spharea(lat,lon,radius):
  lat = vstack([lat[:],lat[0]])
  lon = vstack([lon[:],lon[0]])
  lat0 = 0
  lon0 = 0
  colat = gaa(lat0,lon0,lat,lon)
  az = gcaz(lat0,lon0,lat,lon)
  daz=diff(az, n=1, axis=0)
  daz=wraptopi(daz)
  deltas=diff(colat, n=1, axis=0)
  deltas=deltas/2
  t = colat[0:-1]
  colat=t+deltas
  integrands=multiply((1-cos(colat)),daz)
  area = abs(sum(integrands))/(4*pi)
  area = min(area,1-area);
  area = area*4*pi*radius*radius;
  return area


def spherical_angle(a,b,c):
    import math
    aa = (cos(c) - cos(a)*cos(b))
    #aa = sin(c/2)**2 - sin((a-b)/2)**2
    #hav_C = aa /  (sin(a)*sin(b))
    #d = max(a,b)
    #sin_C = sin(c)*sin(D)/sin(d)
    aa = aa / (sin(a)*sin(b))
    CC = math.acos(aa)
    #CC = math.acos(1 - hav_C * 2) 

    return CC


def mean_endmember_angle(loc_em):
    ex=(loc_em[:,0])
    ey=(loc_em[:,1])
    mpx = mean(ex)
    mpy = mean(ey)
    a=gaa(mpx,mpy,ex[0],ey[0])
    cc= zeros(len(ex))
    for i in xrange(1,len(ex)):
        b=gaa(mpx,mpy,ex[i],ey[i])
        c=gaa(ex[0],ey[0],ex[i],ey[i])
        cc[i] = spherical_angle(a,b,c)
        cc[i]=gcaz(mpy,mpx,ey[i],ex[i]) -gcaz(mpy,mpx,ey[0],ex[0]) 

    return cc


#-----------------------------------------------------------
# MAIN PROCEDURE

#A = zeros((N_channels,N_channels));
#A = mat(A);



close('all')
G=nx.DiGraph()
###G = nx.path_graph(N_channels)
 

maxlen = zeros((N_channels))

print "Populating Nodes/channels"

xxx=0
loc=[];

shp = ogr.Open(shpfile)
lyrcount = shp.GetLayerCount() #multiple layers indicate a directory

lyr = shp.GetLayerByIndex(0)
flds = [x.GetName() for x in lyr.schema]

w_con  =  1000
w_chan =  -1000

for findex in channels:     #xrange(lyr.GetFeatureCount()):
           f = lyr.GetFeature(findex)
           flddata = getfieldinfo(lyr, f, flds)
           g = f.geometry()
           last = g.GetPointCount()  - 1
           ptsi=mat(shapes[findex].points)
           xx = ptsi[:,0]
           yy = ptsi[:,1]

           nanx=argwhere(isnan(xx))          
           nany=argwhere(isnan(yy))
           
           if (nanx == nany and len(nanx) != 0):
              xx[nanx] = []
              yy[nany] = []
           elif (len(nanx) != 0 or len(nany) != 0):
              xx = []
              yy = []

           if (len(xx) > 0):
              [G.add_edge(x,x+1,weight=w_chan) for x in range(xxx,xxx+last)]    

           xxx=xxx+last+1


loc=[((shapes[x].points)) for x in range(0,len(shapes))];
loc = vstack(loc);
 

loc = loc*pi/180;


loc_bad = argwhere(abs(loc)>pi);
loc_bad = loc_bad[:,0]
loc[loc_bad,:] = [0,0];

attributes=dict([(x,loc[x,:]) for x in range(0,len(loc))]);
pos = attributes

print "Populating Edges between channels"
populate=1
if (populate>0):
     xxx=0
     for i in channels:
          ptsi=mat(shapes[i].points)
          xx = ptsi[:,0]
          yy = ptsi[:,1]

          nanx=argwhere(isnan(xx))          
          nany=argwhere(isnan(yy))
          
          if (nanx == nany and len(nanx) != 0):
              xx[nanx] = []
              yy[nany] = []
          elif (len(nanx) != 0 or len(nany) != 0):
              xx = []
              yy = []

          if (len(xx) > 0):          
           yyy = 0 
           dd = str(i*100/N_channels) + '%\r\r'
           os.write(1,dd)
           dist = zeros(len(xx)-1)
           #for dd in arange(len(xx)-1):
           #         dist[dd] = gaa(xx[dd],yy[dd],xx[dd+1],yy[dd+1])  

           #dist = max(dist_min,mean(dist))
           dist = dist_min
           #print "dist", dist
           #print str(i*100/N_channels) + '%\r\r'
           for j in channels:
               pts=mat(shapes[j].points)
               xx2 = pts[:,0]
               yy2 = pts[:,1]
               
               nanx=argwhere(isnan(xx2))          
               nany=argwhere(isnan(yy2))
          
               if (nanx == nany and len(nanx) != 0):
                xx2[nanx] = []
                yy2[nany] = []
               elif (len(nanx) != 0 or len(nany) != 0):
                xx2 = []
                yy2 = []
               if (len(xx2) > 0):  
                #densify geopath 
                if (i != j): 
                   #print "i,j",i,j
                   ds = gaa(xx,yy,xx2[0],yy2[0])
                   de = gaa(xx,yy,xx2[-1],yy2[-1])
                   #print max(ds), max(de)
                   nande = ~isnan(de)
                   nands = ~isnan(ds)
                   fe = np.argwhere(de<=dist)
                   fs = np.argwhere(ds<=dist)
                   nozeroe = array(where(de != 0))
                   nozeroe = nozeroe[0]
                   nozeros = array(where(ds != 0))
                   nozeros = nozeros[0] 
                   if(fs.size >0):
                        iis=array(np.argwhere(array(ds==min(ds))))

                        #iis = np.argwhere(ds == min(ds[nands]))
                        iis = iis.flatten()
                        iis = iis[0]
                        p1 = xxx + iis 
                        p2 = yyy
                        #print "i,j,p1,p2",i,j, p1,p2
                        #print "xx,yy,xx2,yy2,",xx[iis],yy[iis], xx2[0],yy2[0]
                        G.add_edge(p2,p1,weight=w_con)

                   elif (fe.size > 0):
                        iie=array(np.argwhere(array(de==min(de))))
                        iie = iie.flatten()
                        iie = iie[0]
                        p1 = xxx + iie
                        p2 = yyy  + len(pts) -1
                        #print "i,j,p1,p2",i,j, p1,p2
                        #print "xx,yy,xx2,yy2,",xx[iie],yy[iie], xx2[-1],yy2[-1]
                        G.add_edge(p2,p1,weight=w_con)

               yyy = yyy + len(pts)

          xxx = xxx + len(ptsi)



G=G.to_undirected() 
G=nx.minimum_spanning_tree(G)


fig=1


ggg=0
C=nx.connected_component_subgraphs(G)
print "Number of Networks found", len(C)

#file = open("fluvial_network.txt", "w")
newLine="#     <Lon>          <Lat>          A          Length_max          Density"
print newLine
#file.write(newLine)
for g in C[1:len(C)]:
#if (ggg==0):
        figure(ggg)
        
        endmember_indx=argwhere(array(g.degree().values()) == 1)
        #c=[random.random()]*nx.number_of_nodes(da) # random color...
        ce=[random.random()]*(nx.number_of_nodes(g)-1) # random color...
        #nx.draw(g,pos,node_color=c,vmin=0.0,vmax=1.0,with_labels=False,node_size=10,alpha=0.45)
        #nx.draw(g,pos,edge_color=ce,vmin=0.0,vmax=1.0,with_labels=False,alpha=0.45)
        
        loc_em=zeros((len(endmember_indx),2))
        degree_sequence=sorted(nx.degree(g).values(),reverse=True) # degree sequence
        dmax=max(degree_sequence)
        plt.loglog(degree_sequence,'b-',marker='o')
            #plt.title(["Degree rank plot",ggg])
        plt.title('Fluvial Network #'+str(ggg))   #+' Area:'+str(area)+' density:'+str(ddensity))
        plt.ylabel("degree")
        plt.xlabel("rank")
        plt.axes([0.45,0.45,0.45,0.45])
        plt.axis('off')
            #nx.draw_networkx_nodes(g,pos,node_size=20)
        nx.draw_networkx_edges(g,pos,edge_color=ce) 



        if len(endmember_indx)>1:
            
            for i in arange(len(endmember_indx)):
                endmember=g.nodes()[endmember_indx[i]]
                loc_em[i] = loc[endmember]

            loc_em=mat(loc_em)
            ex=(loc_em[:,0])
            ey=(loc_em[:,1])
            #drainage mean point location (mean of endmembers) 
            mpx = mean(ex)
            mpy = mean(ey)
   
       
            endmember_angle = mean_endmember_angle(loc_em)
            #print "angle", endmember_angle*180/pi    
            em_order = endmember_angle.argsort()
            ex = ex[em_order]
            ey = ey[em_order]
            area = poly2spharea(ey,ex,RTitan)
  
            gnodes=g.nodes()
            g_ex=loc[gnodes,0]
            g_ey=loc[gnodes,1]

            #print "g: min/max:", min(g_ex),max(g_ex),min(g_ey),max(g_ey)

            exd = [max(g_ex), max(g_ex), min(g_ex),min(g_ex)]
            eyd = [max(g_ey), min(g_ey), min(g_ey),max(g_ey)]
            #print ex,exd,ey,eyd
            area = poly2spharea(vstack(eyd),vstack(exd),RTitan)



            #points1 = MultiPoint([(0, 0), (1, 1), (0, 2), (2, 2), (3, 1), (1, 0),(2,3),(3,6),(5,1)])
            points1 = MultiPoint(loc[g.nodes()])

            hull1 = points1.convex_hull
            #patch1 = PolygonPatch(hull1, facecolor='#6699cc', edgecolor='#6699cc', alpha=0.5, zorder=10)

            envelop = vstack(list(hull1.exterior.coords))
            area = poly2spharea(vstack(envelop[:,1]),vstack(envelop[:,0]),RTitan)


            d = zeros(len(g.edges()))
            for i in xrange(0,len(g.edges())):
                tt = array(g.edges()[i])
                d[i] = gaa(loc[tt[0],0],loc[tt[0],1],loc[tt[1],0],loc[tt[1],1])
          
            wnan = argwhere(isnan(d))
            d[wnan] = 0;
            length_max = max(cumsum(d))*RTitan
            ddensity = length_max/area
     
            newline=ggg,mean(loc[g.nodes(),0]),mean(loc[g.nodes(),1]), area, length_max, ddensity
            print newline
            #file.write(newLine)
            #endmember_angle=matrix(vstack(endmember_angle))
            #bmat('endmember_angle,ex,ey')
            endmember_indx = endmember_indx[em_order]

            da=nx.Graph()

            #for i in arange(len(endmember_indx)-1):
            #    da.add_edges_from([(g.nodes()[endmember_indx[i]],g.nodes()[endmember_indx[i+1]])])

            #da.add_edges_from([(g.nodes()[endmember_indx[-1]],g.nodes()[endmember_indx[0]])])
          
            for i in arange(len(envelop)-1):
                da.add_edges_from([(i,i+1)])

            da.add_edges_from([(len(envelop)-1,0)])

            #da.add_edges_from([(0,1),(1,2),(2,3),(3,0)])
            da.to_undirected()

            #dapos=dict([(x,[exd[x],eyd[x]])  for x in range(0,len(exd))]);
            dapos=dict([(x,[envelop[x,0],envelop[x,1]])  for x in range(0,len(envelop))]);
          
            
            nx.draw_networkx_edges(da,dapos,alpha=0.45)
            plot(ex,ey,'ro')
            #plot(mpx,mpy,'g*')
            #title('Fluvial Network #'+str(ggg))
        plt.savefig("networkx_"+str(ggg)+".eps")
        ##ns.write_shp(da,"networkx_"+str(ggg)+".shp")
        close('all')
        ggg=ggg+1


test8=0
if (test8>0):
    g=C[3]
    g2=nx.Graph()
    g2.add_nodes_from(g.edges())
    for e1 in g2.nodes():
       for e2 in g2.nodes():
          if e1!=e2 and (np.in1d(array(e1),array(e2))).any():
             g2.add_edge(e1,e2)

    fig=fig+1
    figure(fig)
    subplot(121)
    nx.draw_networkx_edges(g) 
    title("graphe des noeuds")
    subplot(122)
    nx.draw_networkx_edges(g2) 
    title("graphe des liens")
    #savefig("directedgraph")


deg_seq=0
if (deg_seq==1):
    g=C[3]
    g2=nx.Graph()
    g2.add_nodes_from(g.edges())
    for e1 in g2.nodes():
        for e2 in g2.nodes():
            if e1!=e2 and (np.in1d(array(e1),array(e2))).any():
                g2.add_edge(e1,e2)

    degree_sequence=sorted(nx.degree(g2).values(),reverse=True)


fig2=0
if (fig2>0):
    ggg=1
    for g in C[0:len(C)]:       
        gg=ggg 
        fig=fig+1
        figure(fig)
        degree_sequence=sorted(nx.degree(g).values(),reverse=True) # degree sequence
        dmax=max(degree_sequence)
        plt.loglog(degree_sequence,'b-',marker='o')
        plt.title(["Degree rank plot",ggg])
        plt.ylabel("degree")
        plt.xlabel("rank")
        plt.axes([0.45,0.45,0.45,0.45])
        plt.axis('off')
        nx.draw_networkx_nodes(g,pos,node_size=20)
        nx.draw_networkx_edges(g,pos,alpha=0.4)
        ggg =ggg + 1
        plt.savefig("rank_network"+str(gg)+".eps")
       
show()

#file.close()


