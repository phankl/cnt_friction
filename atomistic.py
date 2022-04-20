import numpy as np

import cnt

def dataFile(cnt1, cnt2, cellMin, cellMax, fileName):
  sites1 = cnt1.sites
  sites2 = cnt2.sites

  sites = np.concatenate((sites1, sites2))

  with open(fileName, 'w') as f:
    f.write("LAMMPS CNT Data File\n\n")
    
    f.write(f"  {len(sites1)+len(sites2)} atoms\n")
    f.write("  0 bonds\n")
    f.write("  0 angles\n")
    f.write("  0 dihedrals\n")
    f.write("  0 impropers\n\n")

    f.write("  1 atom types\n")
    f.write(f"  {cellMin[0]} {cellMax[0]} xlo xhi\n")
    f.write(f"  {cellMin[1]} {cellMax[1]} ylo yhi\n")
    f.write(f"  {cellMin[2]} {cellMax[2]} zlo zhi\n\n")

    f.write("Masses\n\n")

    f.write("  1 12.0\n\n")

    f.write("Atom\n\n")

    atomID = 1
    for site in sites1:
      f.write(f"{atomID} 1 1 0 {site[0]} {site[1]} {site[2]}\n")
      atomID += 1
    for site in sites2:
      f.write(f"{atomID} 2 1 0 {site[0]} {site[1]} {site[2]}\n")
      atomID += 1

n = 10
m = 10

angle = 0.0
angle *= np.pi / 180.0

cells = 102

dist = 2.0 * 6.785 + 5

cnt1 = cnt.CNT(n, m, cells, cellMode=True, axis=(1.0, 0.0, 0.0))
cnt2 = cnt.CNT(n, n, cells, cellMode=True, origin=(0.0, dist, 0.0), axis=(1.0, 0.0, 0.0))

dataFile(cnt1, cnt2, (0.0, -100, -100), (cnt1.length, 100, 100), "cnt.data")
