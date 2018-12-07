import MDAnalysis
import MDAnalysis.lib.mdamath

# open DL_POLY trajectory and associated .psf file
u = MDAnalysis.Universe('data.psf','HISTORY',format='HISTORY',dt=0.1)
# MDAnalysis.core.flags['use_pbc'] = True

# prepare output for the original trajectory and one centred on the solute
t = MDAnalysis.Writer("traj.xtc", u.trajectory.n_atoms)
w = MDAnalysis.Writer("init.gro", multiframe = False)

# write first frame, that is joined, centred and wrapped
for residue in u.residues:
    MDAnalysis.lib.mdamath.make_whole(residue)
w.write(u.atoms)
u.trajectory.rewind()
# run through trajectory and:
# 1) 'join' (make whole) all residues, to remove split bonds across box edges
for ts in u.trajectory:
    for residue in u.residues:
        MDAnalysis.lib.mdamath.make_whole(residue)
    # write original trajectory in xtc format
    u.atoms.wrap(compound='residues',center='com')
    t.write(u.atoms)
t.close()
