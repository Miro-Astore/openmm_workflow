import MDAnalysis as mda
import numpy as np

# Load your solvated system
u = mda.Universe('step3_input.psf','step3_input.pdb')
out_file = 'out.rst7'


#Assuming ihe PBC metadata is stored in the trajectory's dimensions
# We can extract the box dimensions
dimensions = u.trajectory.ts.dimensions
if dimensions is None:
    x_span = np.max(u.atoms.positions[:,0]) - np.min(u.atoms.positions[:,0])
    y_span = np.max(u.atoms.positions[:,1]) - np.min(u.atoms.positions[:,1])
    z_span = np.max(u.atoms.positions[:,2]) - np.min(u.atoms.positions[:,2])
    dimensions = np.array([x_span, y_span, z_span, 90, 90, 90] )
    u.trajectory.ts.dimensions = dimensions


def write_rst7(filename, title, natom, coordinates):
    with open(filename, 'w') as file:
        file.write(title + '\n')  # write the title
        file.write(f'{natom:6d}\n')  # write the number of atoms
        # Write the coordinates
        for i in range(0, natom):

            for j in range(3):
                file.write(f'{coordinates[i, j]:12.7f}')
            if (i % 2) == 1: 
                file.write('\n')
        if (i % 2) == 0: 
            file.write('\n')
        for j in range(6):
            file.write(f'{dimensions [j]:12.7f}')


        file.close()

# We can define the corner of the box as the minimum x, y, and z
# which are the first three elements in dimensions
box_dimensions = dimensions[0:3]

current_origin = u.atoms.center_of_mass()

# Then we can shift all atom positions by subtracting the corner coordinates
u.atoms.positions -= (current_origin + box_dimensions/2)

running_number = 0 

write_rst7('out.rst7','written',u.atoms.n_atoms,u.atoms.positions)
# Now the system's corner corresponds to the origin.
# You can now save your new system with the corner at the origin
#u.atoms.write('ionized.pdb')
