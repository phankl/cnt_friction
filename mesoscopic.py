import numpy as np

import cnt

def dataFile(n, m, lTube, nSegment, fileName):
  tube = cnt.CNT(n, m, lTube, cellMode=False)
  radius = tube.radius
  mSegment = 12 * len(tube.sites) / nSegment
  lSegment = lTube / nSegment

  with open(fileName, 'w') as f:
    f.write("LAMMPS CNT Data File\n\n")
    
    f.write(f"  {2*nSegment} atoms\n")
    f.write(f"  {2*nSegment} bonds\n")
    f.write(f"  {2*nSegment} angles\n")
    f.write("  0 dihedrals\n")
    f.write("  0 impropers\n\n")

    f.write("  3 atom types\n")
    f.write("  1 bond types\n")
    f.write("  1 angle types\n")
    f.write(f"  0 {lTube} xlo xhi\n")
    f.write(f"  -100 100 ylo yhi\n")
    f.write(f"  -100 100 zlo zhi\n\n")

    f.write("Masses\n\n")

    f.write(f"  1 {mSegment}\n")
    f.write(f"  2 {mSegment}\n")
    f.write(f"  3 {mSegment}\n\n")

    f.write("Atoms\n\n")

    # atoms
    
    atomID = 1
    for i in range(nSegment):
      if i == 0:
        f.write(f"{atomID} 1 2 0 {i*lSegment} 0 0\n")
      else:
        f.write(f"{atomID} 1 1 0 {i*lSegment} 0 0\n")
      atomID += 1

    for i in range(nSegment):
      if i == 0:
        f.write(f"{atomID} 2 3 0 {i*lSegment} {2*radius+3} 0\n")
      else:
        f.write(f"{atomID} 2 1 0 {i*lSegment} {2*radius+3} 0\n")
      atomID += 1
        
    # bonds

    f.write("\nBonds\n\n")

    bondID = 1
    for i in range(nSegment-1):
      f.write(f"{bondID} 1 {i+1} {i+2}\n")
      bondID += 1
    f.write(f"{bondID} 1 {nSegment} 1\n")
    bondID += 1

    for i in range(nSegment-1):
      f.write(f"{bondID} 1 {nSegment+i+1} {nSegment+i+2}\n")
      bondID += 1
    f.write(f"{bondID} 1 {2*nSegment} {nSegment+1}\n")

    # angles

    f.write("\nAngles\n\n")
    
    angleID = 1
    for i in range(nSegment-2):
      f.write(f"{angleID} 1 {i+1} {i+2} {i+3}\n")
      angleID += 1
    f.write(f"{angleID} 1 {nSegment-1} {nSegment} 1\n")
    angleID += 1
    f.write(f"{angleID} 1 {nSegment} 1 2\n")
    angleID += 1
    
    for i in range(nSegment-2):
      f.write(f"{angleID} 1 {nSegment+i+1} {nSegment+i+2} {nSegment+i+3}\n")
      angleID += 1
    f.write(f"{angleID} 1 {2*nSegment-1} {2*nSegment} {nSegment+1}\n")
    angleID += 1
    f.write(f"{angleID} 1 {2*nSegment} {nSegment+1} {nSegment+2}\n")

dataFile(10, 10, 250, 25, "meso.data")
