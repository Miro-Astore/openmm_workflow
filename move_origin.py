import MDAnalysis as mda
import numpy as np

# Load your solvated system
u = mda.Universe('step3_input.psf','step3_input.pdb')

# Assuming the PBC metadata is stored in the trajectory's dimensions
# We can extract the box dimensions
dimensions = u.trajectory.ts.dimensions
if dimensions == None:
    print('hi')
    x_span = np.max(u.atoms.positions[:,0]) - np.min(u.atoms.positions[:,0])
    y_span = np.max(u.atoms.positions[:,1]) - np.min(u.atoms.positions[:,1])
    z_span = np.max(u.atoms.positions[:,2]) - np.min(u.atoms.positions[:,2])
    dimensions = np.array([x_span, y_span, z_span, 90, 90, 90] )

# We can define the corner of the box as the minimum x, y, and z
# which are the first three elements in dimensions
print(dimensions)
box_dimensions = dimensions[0:3]

current_origin = u.atoms.center_of_mass()

# Then we can shift all atom positions by subtracting the corner coordinates
u.atoms.positions -= (current_origin + box_dimensions/2)

# Now the system's corner corresponds to the origin.
# You can now save your new system with the corner at the origin
u.atoms.write('ionized.pdb')
