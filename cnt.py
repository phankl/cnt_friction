# CNT and Unit Cell Class

import constants as const
import numpy as np
from scipy import spatial

class UnitCell:
  
  axis = np.array([0.0, 0.0, 1.0])

  def __init__(self, n, m):
    
    self.n = n
    self.m = m
 
    # Find atoms in unit cell on graphene sheet

    nR = np.gcd(2*n+m, 2*m+n)
    nC = 4 * (n*n + n*m + m*m) // nR
    t1 = (2*m+n) // nR
    t2 = -(2*n+m) // nR

    chiralVector = n*const.A1 + m*const.A2
    axisVector = t1*const.A1 + t2*const.A2

    chiralVectorLength = np.linalg.norm(chiralVector)
    axisVectorLength = np.linalg.norm(axisVector)

    self.latticeVector = axisVectorLength * self.axis

    chiralVectorNormalised = chiralVector / chiralVectorLength
    axisVectorNormalised = axisVector / axisVectorLength

    radius = 0.5 * chiralVectorLength / np.pi

    self.radius = radius
    self.length = axisVectorLength
   
    nMin = min(0, n, t1, n+t1)
    nMax = max(0, n, t1, n+t1)
    mMin = min(0, m, t2, m+t2)
    mMax = max(0, m, t2, m+t2)

    nRange = range(nMin, nMax+1)
    mRange = range(mMin, mMax+1)

    siteCandidatesA = np.array([n_*const.A1 + m_*const.A2 for n_ in nRange for m_ in mRange])
    siteCandidatesB = np.array([n_*const.A1 + m_*const.A2 + const.D for n_ in nRange for m_ in mRange])

    siteCandidates = np.concatenate((siteCandidatesA, siteCandidatesB))

    # Check if candidates are in unit cell and compute 3D coordinates in CNT

    circ = siteCandidates @ chiralVectorNormalised
    z = siteCandidates @ axisVectorNormalised
    circMask = (circ > -const.EPS) & (circ < chiralVectorLength - const.EPS) 
    zMask = (z > -const.EPS) & (z < axisVectorLength - const.EPS)
    mask = circMask & zMask
    
    self.sites = np.zeros((np.sum(mask), 3))
    circ = circ[mask] / radius
    self.sites[:, 0] = np.sin(circ) * radius
    self.sites[:, 1] = np.cos(circ) * radius
    self.sites[:, 2] = z[mask]

class CNT:
  
  def __init__(self, n, m, length, cellMode=True, rot=0.0, origin=(0.0, 0.0, 0.0), axis=(0.0, 0.0, 1.0)):

    self.n = n
    self.m = m
    
    unitCell = UnitCell(n, m)
    self.unitCell = unitCell

    if cellMode:
      self.length = length * unitCell.length
    else:
      self.length = length
    
    self.origin = np.array(origin)
    self.axis = np.array(axis / np.linalg.norm(axis))

    self.radius = unitCell.radius

    cellLength = unitCell.length
    cellSites = unitCell.sites
    cellAxis = unitCell.axis
   
    # Rotate unit cell along origin axis

    rotationAxis = cellAxis
    rotationAngleCos = np.cos(rot)
    rotationAngleSin = np.sin(rot)

    cellSites = [rotationAngleCos*site + rotationAngleSin*np.cross(rotationAxis, site) + (1.0-rotationAngleCos)*np.dot(rotationAxis, site)*rotationAxis for site in cellSites]

    # Rotate unit cell into axis

    rotationAxis = np.cross(cellAxis, self.axis)
    rotationAngleCos = np.dot(cellAxis, self.axis)
    if rotationAngleCos + 1.0 < const.EPS:
      cellSites = -cellSites
    else:
      cellSites = [rotationAngleCos*site + np.cross(rotationAxis, site) + np.dot(rotationAxis, site)/(1.0+rotationAngleCos)*rotationAxis for site in cellSites]

    latticeVector = cellLength * np.array(self.axis)
    
    cellNumber = length
    if not cellMode:
      cellNumber = np.ceil(length/cellLength).astype('int')

    cells = [cellSites + i*latticeVector + origin for i in range(cellNumber-1)]
    lastCell = cellSites + (cellNumber-1)*latticeVector + origin
    if not cellMode:
      lastCell = lastCell[np.dot(lastCell - origin, self.axis) < length + const.EPS]
    
    if cellNumber > 1:
      cells.append(lastCell)
      self.sites = np.concatenate(cells)
    else:
      self.sites = lastCell
