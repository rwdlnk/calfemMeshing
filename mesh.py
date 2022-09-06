
def psi(xi, eta):  # Element basis functions on 4-node quad element
  return [(1.-xi)*(1.-eta)/4., (1.+xi)*(1.-eta)/4.,(1+xi)*(1+eta)/4,(1-xi)*(1+eta)/4]
# END OF psi

def dpsidxi(xi,eta):  # xi and eta derivatives of element basis functions
  ddxi = [[-(1-eta)/4, (1-eta)/4, (1+eta)/4, -(1+eta)/4],
          [-(1-xi)/4, -(1+xi)/4,  (1+xi)/4,  (1-xi)/4]]
  return ddxi
# END OF dpsidxi

def CGmetric(rx, ry, dpsi_xi):  # Element metric tensor 
  # Two dimensional coarse-grained metric
  gab = [[0,0],[0,0]]  # Covariant metric tensor

  for a in range(2): # Components in xi_i = (xi, eta) 
    for b in range(2):
      for m in range(4):  # 4 nodes per quad
        for n in range(4):  # 4 nodes per quad
          gab[a][b] += dpsi_xi[a][m]*dpsi_xi[b][n]*(rx[m]*rx[n] + ry[m]*ry[n])

  g = gab[0][0]*gab[1][1] - gab[1][0]*gab[0][1]  # determinate of gab

  gAB = [[gab[1][1], -gab[0][1]],[-gab[1][0], gab[0][0]]]/g  # Contravariant metric tensor

  return g, gAB

# END OF CGmetric

def targetMetric(elmx, elmy):  # Use target metric per Hansen, Douglass, Zardecki 9.2.3
  # Two dimensional metric
  gab = [[0,0],[0,0]] # Covariant metric tensor approx.

  dpsi_dxii = dpsidxi(0.,0.); # target metric is evaluated at the element center (xi, eta) = (0,0)

  # Using Eq 9.19 (2-D form) to compute gab for each element     
  for a in range(2):
    for b in range(2):
      for i in range(4):  # dot product sum
        for j in range(4):  # dot product sum
          gab[a][b] += (elmx[i]*elmx[j] + elmy[i]*elmy[j])*dpsi_dxii[a][i]*dpsi_dxii[b][j]

  g = gab[0][0]*gab[1][1] - gab[1][0]*gab[0][1]  # determinate of gab

  gAB = [[gab[1][1]/g, -gab[0][1]/g],[-gab[1][0]/g, gab[0][0]/g]]  # Contravariant metric tensor

  return g, gAB
# END OF targetMetric

def metric(elmx, elmy, dpsi_xi):  # Element metric tensor 
  # Two dimensional element metric tensor

  gab = [[0,0],[0,0]]  # Covariant metric tensor
  g = 0. # determinant of gab
  
  for a in range(2): # Components in xi_i = (xi, eta) 
    for b in range(2):
      for m in range(4):  # 4 nodes per quad
        for n in range(4):  # 4 nodes per quad
          gab[a][b] += dpsi_xi[a][m]*dpsi_xi[b][n]*(elmx[m]*elmx[n] + elmy[m]*elmy[n])

  g = gab[0][0]*gab[1][1] - gab[1][0]*gab[0][1]  # determinate of gab

  gAB = [[gab[1][1]/g, -gab[0][1]]/g,[-gab[1][0]/g, gab[0][0]]/g]  # Contravariant metric tensor

  return g, gAB

# END OF metric

def NSF(ex, ey, E):  # H,D,Z Sect. 9.3

  from math import sqrt

#  N = 2 # second order integration
#  sq3 = sqrt(3) 
#  xi =  [-sq3, sq3]
#  eta = [-sq3, sq3]
#  W =   [1., 1.]

#  N = 3 # second order integration
#  sq35 = sqrt(3/5) 
#  xi =  [-sq35, 0., sq35]
#  eta = [-sq35, 0., sq35]
#  W =   [5./9, 8./9., 5./9.]

  N = 4 # second order integration
  sq65 = sqrt(6/5)
  c37 = 3/7
  c27 = 2/7
  cc = c27*sq65
  xi =   [-sqrt(c37 - cc), -sqrt(c37 + cc), sqrt(c37-cc), sqrt(c37+cc)]
  eta =  [-sqrt(c37 - cc), -sqrt(c37 + cc), sqrt(c37-cc), sqrt(c37+cc)]
  c12 = 1/2
  cs = sqrt(30)/36
  W =   [c12+cs, c12-cs, c12+cs, c12-cs]

  I = 0.  # Total smoothness funtional

  for e in range(E):
    Ie = 0. # Smoothness function integral for element e

    for i in range(N):   # Integration outer loop
      for j in range(N): # Integration inner loop
        dPsi_xi = dpsidxi(xi[i],eta[j])  # compute derivative of psi at integration points
        g, gAB = metric(ex[e], ey[e], dPsi_xi)  # compute metric tensor at integration points
        Ie += W[i]*W[j]*(gAB[0][0] + gAB[1][1]) * sqrt(g) # compute element integral

    I += Ie/2.

  return I

#  END OF NSF

#def Kmn(ex, ey, coords, nodeCon):  # Element stiffness matrix
#def Kmn(ex, ey, rx, ry):  # Element stiffness matrix
def Kmn(rx, ry):  # Element stiffness matrix
  from math import sqrt

#  N = 9 # ninth order integration
#  sq35 = sqrt(3/5) 
#  xi =  [-sq35, 0, sq35, -sq35, 0, sq35, -sq35, 0, sq35]
#  eta = [-sq35, -sq35,-sq35, 0, 0, 0, sq35, sq35, sq35]
#  W =   [25/81, 40/81, 25/81, 40/81, 64/81, 40/81, 25/81, 40/81, 25/81]

#  N = 2 # second order integration
#  sq3 = sqrt(3) 
#  xi =  [-sq3, sq3]
#  eta = [-sq3, sq3]
#  W =   [1, 1]

#  N = 3 # second order integration
#  sq35 = sqrt(3/5) 
#  xi =  [-sq35, 0., sq35]
#  eta = [-sq35, 0., sq35]
#  W =   [5./9, 8./9., 5./9.]

  N = 4 # second order integration
  sq65 = sqrt(6/5)
  c37 = 3/7
  c27 = 2/7
  cc = c27*sq65
  xi =   [-sqrt(c37 - cc), -sqrt(c37 + cc), sqrt(c37-cc), sqrt(c37+cc)]
  eta =  [-sqrt(c37 - cc), -sqrt(c37 + cc), sqrt(c37-cc), sqrt(c37+cc)]
  c12 = 1/2
  cs = sqrt(30)/36
  W =   [c12+cs, c12-cs, c12+cs, c12-cs]

  # Integrate the stiffness matrix
  
  K = [[0,0,0,0],
       [0,0,0,0],
       [0,0,0,0],
       [0,0,0,0]]

  #g, gAB = targetMetric(ex,ey)  # compute metric tensor at element center

  for i in range(N):   # Integration outer loop
    for j in range(N): # Integration inner loop
      dPsi_xi = dpsidxi(xi[i],eta[j])  # compute derivative of psi at integration points
      #g, gAB = metric(ex, ey, dPsi_xi)  # compute metric tensor at integration points
      g, gAB = CGmetric(rx, ry, dPsi_xi)  # compute metric tensor at integration points

      for a in range(2):  # alpha loop
        for b in range(2): # beta loop
          for m in range(4): # node outer loop
            for n in range(4): # node inner loop
              K[m][n] += W[i]*W[j]*sqrt(g)*dPsi_xi[a][m]*gAB[a][b]*dPsi_xi[b][n]
  return K

# END OF Kmn

#def Kmatrix(ex, ey, enodes, N):  # list of element x,y coodrs, list of element nodes, list of boundary nodes, num nodes
def Kmatrix(ex, ey, enodes, N, rbar):  # list of element x,y coodrs, list of element nodes, list of boundary nodes, num nodes
  # Setup the global stiffness matrix 
  import numpy as np
  E = len(ex)  # number of elements

  K = np.zeros([N,N])

  rx = []
  ry = []

  for e in range(E):
    x = ex[e]
    y = ey[e]
    for n in range(4):
      rx.append(rbar[enodes[e][n]-1][0])
      ry.append(rbar[enodes[e][n]-1][1])

    #Ktmp = Kmn(x, y, coords)
    Ktmp = Kmn(rx, ry)

# Assemble global matrix
    for m in range(4):
      I = enodes[e][m]-1 # enodes, bnodes are Fortran-like 1 based; subtract 1 to get 0 based
      for n in range(4):
        J = enodes[e][n]-1
        K[I][J] += Ktmp[m][n]
  return K  
# END OF Kmatrix

def elementConnectivity(edof, E):
# Compute element to element connectivity
# Assume mesh consists of only 4-node quads

  import numpy as np

  elmCon = np.zeros([E,4],dtype=int)
  for element in range(E): 
    nodes =  edof[element]
    sides = [[nodes[0],nodes[1]],
             [nodes[1],nodes[2]],
             [nodes[2],nodes[3]],
             [nodes[3],nodes[0]]]
    for s in range(4):
      for e in range(E):
        if( e == element):  # skip current element
          continue 
        Z = set(edof[e]) & set(sides[s])
        if(len(Z) < 2):
          elmCon[element][s] = -1  # Side s is on an external boundary
        else:
          elmCon[element][s] = e + 1
          break   # Found the neighboring element for side s of element

  return elmCon
        
# END OF elementConnectivity

def nodeConnectivity(elmCon, edof, E, N, bdyNodes, BN):
  # Find neighbor nodes connected to each element node
  # Search for neighbor sides/nodes containing node nn in element e

  nodeCon = [[] for i in range(N)]

# Do boundary nodes first
  for n in range(BN):     # For each element

    item = bdyNodes[n]

    if(n == 0):
      nds = [bdyNodes[n+1],bdyNodes[-1]]
    elif(n == BN-1):
      nds = [bdyNodes[n-1], bdyNodes[0]]
    else:
      nds = [bdyNodes[n+1], bdyNodes[n-1]]
    
    nodeCon[item-1].append(nds[0])  
    nodeCon[item-1].append(nds[1])  

# Now, insert interior nodes 
  for e in range(E):     # For each element
    #print("e: ", e)

    neighborEs = elmCon[e] # Get list of neighboring elements to e
    #print("neighborEs: ", neighborEs)

    for n in range(4):  # for each node in element e
    #  print("n: ", n)

      item = edof[e][n]  # node at center of connectivity
    #  print("item: ", item)

      nList = nodeCon[item-1]  
    #  print("nList: ", nList)
      for i in range(4):  # For each neighbor element
        ne = neighborEs[i]
    #    print("neighborE: ", ne, " i: ", i, "  edof[ne-1]: ", edof[ne-1])

        if(ne == -1):   # Neighbor "element" is a boundary, skip
          continue 
        try:
          index = list(edof[ne-1]).index(item)
    #      print("index: ", index)
          if(index == 0):
            nList.append(edof[ne-1][1])
            nList.append(edof[ne-1][3])
          elif(index == 3):
            nList.append(edof[ne-1][0])
            nList.append(edof[ne-1][2])
          else:
            nList.append(edof[ne-1][index+1])
            nList.append(edof[ne-1][index-1])
    #      print("nList before set: ", nList)
          tmp = set(nList)  # eliminates duplicates
    #      print("  tmp: ",   tmp)
          nodeCon[item-1] = list(tmp)
    #      print("\n")
        except ValueError:
          continue
  return nodeCon

# END OF nodeConnectivity

def elementCentroid(ex, ey, E):

  elCent = np.zeros([E,2],dtype=float)
  Atot = 0.

  for e in range(E):

    x = [ex[e][0],ex[e][1],ex[e][2],ex[e][3],ex[e][0]]
    y = [ey[e][0],ey[e][1],ey[e][2],ey[e][3],ey[e][0]]
    
    A = 0.

    for i in range(4):
      xy = x[i]*y[i+1] - x[i+1]*y[i]
      A += xy
      elCent[e][0] += (x[i]+x[i+1])*xy
      elCent[e][1] += (y[i]+y[i+1])*xy

    A /= 2. 

    elCent[e][0] /= (6.*A)
    elCent[e][1] /= (6.*A)
    Atot += A

  return elCent, Atot

# END OF elementCentroid

def extractNodeCoords(ex, ey, edof, E, N):

  X = [[] for i in range(N)]
  Y = [[] for i in range(N)]

  for e in range(E):
    enodes = edof[e]
    for i in range(4):
      nd = edof[e][i] -1
      if(len(X[nd]) > 0):
        continue
      xn = ex[e][i]
      yn = ey[e][i]
      X[nd].append(xn)
      Y[nd].append(yn)

  return X, Y

# END OF extractNodeCoords

def rbar(nodeCon, N, coords):

  rb = [[] for i in range(N)]

  for n in range(N):
    xsum = 0.
    ysum = 0.

    nn = len(nodeCon[n])
    #print("nn: ", nn, " nodeCon: ", nodeCon[n])
    for j in range(nn):
      node = nodeCon[n][j]
      #print("nodeCon[", n, "][", j, "]: ", node)
      #print("   len(xn): ", len(xn))
      xsum += coords[node-1][0]
      ysum += coords[node-1][1]
    rb[n].append(xsum/(float(nn))) 
    rb[n].append(ysum/(float(nn))) 

  return rb

# END OF rbar

def LBsmooth(coords, edof, nodes, bdyNodes):
  import numpy as np
  from math import sqrt
  from scipy.linalg import lu_factor, lu_solve

  E = len(edof) # number of elements
  N = np.size(nodes) # number of nodes
  #print("nNodes: ", N, "  nElements: ", E)

  elmCon = elementConnectivity( edof, E)
  #print("\n elmCon: ", elmCon[:15,:])

  BN = np.size(bdyNodes) # number of boundary nodes
  #print("nBdyNodes: ", BN)

  nodeCon = nodeConnectivity(elmCon, edof, E, N, bdyNodes, BN)
  #print("\n nodeCon: ", nodeCon[:15])

  bdyCoords = coords[:BN] # boundary node x,y pairs
  bdyX = bdyCoords[:,0] # boundary node x coordinates
  bdyY = bdyCoords[:,1] # boundary node y coordinates

  ex, ey = getElementXY(edof, coords, dofs) # Element x,y coordinates; for each node in the element

  I0 = NSF(ex, ey, E) # Average 
  print("Initial NSF I0: ", I0)

  epsilon = 1.e-6
  error = 1.
  iter = 1
  maxIter = 50

  while error >= epsilon and iter <= maxIter:

    #nx, ny = extractNodeCoords(ex, ey, edof, E, N)

    rBar = rbar(nodeCon, N, coords)

    #elmCent, Atot = elementCentroid(ex, ey, E)
    #print("\n elmCent: ", elmCent[:5,:])

    intCoords = coords[BN:] # interior node x,y pairs

    Xint = intCoords[:,0] # interior node x coordinates
    Yint = intCoords[:,1] # interior node y coordinates

    # Construct stiffness matrix
    #K = Kmatrix(ex, ey, edof, N)
    K = Kmatrix(ex, ey, edof, N, rBar)

    # Partition K to apply fixed node boundary conditions.
    # Solving
    #  Q X = - P bdyX
    #  Q Y = - P bdyY where

    P = K[BN:,0:BN]
    Q = K[BN:,BN:]
    # and x and y are partitioned as

    BCnodes = nodes[0:BN] # boundary node numbers
    INTnodes = nodes[BN:]

    RHSx = -np.matmul(P, bdyX)
    RHSy = -np.matmul(P, bdyY)

    lu, piv = lu_factor(Q)

    X = lu_solve((lu, piv), RHSx)
    if(not np.allclose(Q @ X - RHSx, np.zeros(N-BN,), rtol=1.e-4)):
      print("Solution for X failed!\n")

    dx = 0.
    XXint = zip(X, Xint)
    for x1,x2 in XXint:
      dx += abs(x1-x2) 

    Y = lu_solve((lu, piv), RHSy)
    if(not np.allclose(Q @ Y - RHSy, np.zeros(N-BN,), rtol=1.e-4)):
      print("Solution for Y failed!\n")

    dy = 0.
    YYint = zip(Y, Yint)
    for y1,y2 in YYint:
      dy += abs(y1-y2) 

    error = (dx*dx + dy*dy)/(N-BN)

    j = 0
    for i in INTnodes:
      coords[i-1][0] = X[j]
      coords[i-1][1] = Y[j]
      j += 1

    ex, ey = getElementXY(edof, coords, dofs) # Element x,y coordinates; for each node in the element

    I = NSF(ex, ey, E) 
    
    print("Iteration: ", iter, ", average nodal correction: ", error, ", NSF: ", I/I0)

    iter += 1

  writeData(N, BN, coords, ex, ey, edof, bdofs, bdyNodes, nodes, bdyX, bdyY, elmCon, nodeCon)

# END OF LBsmooth

def getElementXY(edof, coords, dofs):
  import calfem.core as cfc
  ex, ey = cfc.coordxtr(edof, coords, dofs)
  return ex, ey
# END OF getElementXY

def writeData(N, BN, coords, ex, ey, edof, bdofs, bdyNodes, nodes, bdyX, bdyY, elmCon, nodeCon):

  # Write out the element, node and boundary information
  exy = zip(ex,ey)
  with open("elementCoords.dat", 'w+') as fp:
    for X,Y in exy:
      print(X, " ", Y, file=fp)
  fp.close()

  with open("elementNodes.dat", 'w+') as fp:
     for nd in edof:
      print(nd, file=fp)
  fp.close()

  with open("elementConnectivity.dat", 'w+') as fp:
     for nd in elmCon:
      print(nd, file=fp)
  fp.close()

  with open("nodeConnectivity.dat", 'w+') as fp:
     for nd in nodeCon:
      print(nd, file=fp)
  fp.close()

  with open("bdryNodes.dat", 'w+') as fp:
    print(bdyNodes, file=fp)
  fp.close()

  with open("Nodes.dat", 'w+') as fp:
    for n in nodes: 
      print(n, file=fp)
  fp.close()

  #print("bdyX: ", bdyX)
  #print("bdyY: ", bdyY)

  #print("Bndy coords: ", coords[:BN,:][:][1])
  #print("Bndy coords: ", coords[:BN,:][0][1])
  #print("coords[0]: ", coords[0])
  #print("Bndy coords: ", coords[:BN])
  #print("Nbdynodes: ", BN)
  #print("Nbdynodes: ", len(bdofs[0]))
  #print("Node: ", nodes[0], "  coords: ", coords[nodes[0]-1])
  #print("Node: ", nodes[253], "  coords: ", coords[nodes[253]-1])
  #print("Node: ", nodes[N-1], "  coords: ", coords[nodes[N-1]-1])
  #print("length of coords: ", len(coords), "  length of edof: ", len(edof), " len of bdyNodes: ", len(bdyNodes))
  #print("nodes of element[0]: ", edof[0],"  coords of element-node[0]: ", ex[0][0], ", ", ey[0][0])
  #print("nodes of element[0]: ", edof[0],"  coords of element-node[0]: ", coords[edof[0][0]-1][0], ", ", coords[edof[0][0]-1][1])
  #print("nodes of element[0]: ", edof[0],"  coords of element-node[0]: ", coords[edof[0][0]-1])
  #print("nodes of element[0]: ", edof[0],"  coords of element-node[1]: ", coords[edof[0][1]-1])
  #print("nodes of element[0]: ", edof[0],"  coords of element-node[2]: ", coords[edof[0][2]-1])
  #print("nodes of element[0]: ", edof[0],"  coords of element-node[3]: ", coords[edof[0][3]-1])

  #print("BDY coords[bdofs[5][:]]:\n", coords[bdofs[5][:]])
  #for i in range(15):
  #  print("bdy node coords: Node ", bdyNodes[i]-1, coords[bdyNodes[i]-1])
  #  print("bdy node coords: Node ", bdyNodes[253-i]-1, coords[bdyNodes[253-i]-1])
# END OF writeData
 
# MAIN PROGRAM
import sys
sys.path.insert(0,"/home/rod/ResearchII/work/WFpy/src")

from WilsonFowler import *
#from math import atan, sin, sqrt, cos

import numpy as np
import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.vis as cfv
import calfem.utils as cfu

print("\nMesh Boundary\n")

x = [0.,0.,0.,0.,1.,2.,2.,2., 3.,4.,4.,4.,4.,4.,4.,5.,6.,7.,8.,8.,8.,8.,8.,8.,8.,7.,6.,5.,4.,3.,2.,1.,0.]
y = [0.,1.,2.,3.,3.,3.,2.,1.,1.,1.,2.,3.,4.,5.,6.,6.,6.,6.,6.,5.,4.,3.,2.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.]

N = 33

#sbs = True
#Tx0 = 1./sqrt(5.)
#Ty0 = 2./sqrt(5.)

#ses = True
#TxL = 1.
#TyL = 0.

#wf = WF(x,y,maxIter=10,setBeginSlope=sbs,beginTx = Tx0,beginTy = Ty0, \
#        setEndSlope = ses, endTx = TxL, endTy = TyL,useLineSearch=True)
wf = WF(x,y,maxIter=20,useLineSearch=True)

wf.writeXY("mesh.XY")

xn = []
yn = []

#newN = 66
newN = 128

xn,yn, = wf.Distrib(newN)
#xn,yn, = wf.equalArcLengthDistrib(newN)
wf.writeXY("mesh_new.XY", (xn,yn))

wfn = WF(xn, yn, maxIter=30,useLineSearch=True)
#print(wfn)

# Begin meshing

g = cfg.Geometry()

# Add points
N = newN-1
xy = zip(xn[:-1],yn[:-1])
for X,Y in xy:
  g.point([X,Y],marker=5)

# Add curves
for i in range(N-1):
  g.spline([i, i+1])
g.spline([N-1,0])

g.surface([i for i in range(N)])

cfv.drawGeometry(g)
cfv.showAndWait()

# Create mesh
mesh = cfm.GmshMesh(g)

mesh.el_type = 3
mesh.dofs_per_node = 1  # Degrees of freedom per node.
#mesh.el_size_factor = 0.5  # Factor that changes element sizes.

coords, edof, dofs, bdofs, element_markers = mesh.create()
#cfu.save_mesh(mesh, name='Non-box-quad-domain.cfm')

cfv.figure()

# Draw the mesh.

cfv.drawMesh( 
    coords=coords, 
    edof=edof, 
    dofs_per_node=mesh.dofsPerNode,
    el_type=mesh.elType,
    filled=True,
    title="Gaming Non-box Domain 01")
cfv.showAndWait()
print(dofs)
LBsmooth(coords, edof, dofs[:,0], bdofs[0])

# Draw the smoothed mesh.

cfv.drawMesh( 
    coords=coords, 
    edof=edof, 
    dofs_per_node=mesh.dofsPerNode,
    el_type=mesh.elType,
    filled=True,
    title="Gaming Non-box Domain 01")
cfv.showAndWait()

del wf
del wfn
# END OF main
