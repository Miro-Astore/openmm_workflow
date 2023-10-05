from __future__ import print_function
import argparse
import threading
import sys
import os
import numpy as np

import pdb
import rlcompleter

pdb.Pdb.complete=rlcompleter.Completer(locals()).complete

from openmm_workflow.utilities.omm_readinputs import *
from openmm_workflow.utilities.omm_readparams import *
from openmm_workflow.utilities.omm_vfswitch import *
from openmm_workflow.utilities.omm_barostat import *
from openmm_workflow.utilities.omm_restraints import *
from openmm_workflow.utilities.omm_rewrap import *

from simtk.unit import *
from openmm import *
from openmm.app import *


parser = argparse.ArgumentParser()
parser.add_argument('--platform', nargs=1, help='OpenMM platform (default: CUDA or OpenCL)')
parser.add_argument('-i', dest='inpfile', help='Input parameter file', required=True)
parser.add_argument('-p', dest='topfile', help='Input topology file', required=True)
parser.add_argument('-c', dest='crdfile', help='Input coordinate file')
parser.add_argument('-ipdb', dest='input_pdbfile', help='Input coordinate pdb file')
parser.add_argument('-t', dest='toppar', help='Input CHARMM-GUI toppar stream file (optional)')
parser.add_argument('-b', dest='sysinfo', help='Input CHARMM-GUI sysinfo stream file (optional)')
parser.add_argument('-ff', dest='fftype', help='Input force field type (default: CHARMM)', default='CHARMM')
parser.add_argument('-icrst', metavar='RSTFILE', dest='icrst', help='Input CHARMM RST file (optional)')
parser.add_argument('-irst', metavar='RSTFILE', dest='irst', help='Input restart file (optional)')
parser.add_argument('-ichk', metavar='CHKFILE', dest='ichk', help='Input checkpoint file (optional)')
parser.add_argument('-opdb', metavar='PDBFILE', dest='opdb', help='Output PDB file (optional)')
parser.add_argument('-orst', metavar='RSTFILE', dest='orst', help='Output restart file (optional)')
parser.add_argument('-ochk', metavar='CHKFILE', dest='ochk', help='Output checkpoint file (optional)')
parser.add_argument('-odcd', metavar='DCDFILE', dest='odcd', help='Output trajectory file (optional)')
parser.add_argument('-rewrap', dest='rewrap', help='Re-wrap the coordinates in a molecular basis (optional)', action='store_true', default=False)
parser.add_argument('--restart-timer', dest='restart_timer', help='Choose whether to restart the timer from zero or not.', action='store_true', default=False)
args = parser.parse_args()

# Load parameters
print("Loading parameters")
inputs = read_inputs(args.inpfile)


print('making top')
top = read_top(args.topfile, args.fftype.upper())
print('done')

crd = read_crd(args.crdfile, args.fftype.upper())

if args.input_pdbfile is not None : 
    crd = PDBFile(args.input_pdbfile, args.fftype.upper())
    params = read_params(args.toppar)
    top = read_box(top, args.sysinfo) if args.sysinfo else gen_box(top, crd)

# Build system
print("building system")
nboptions = dict(nonbondedMethod=inputs.coulomb,
                 nonbondedCutoff=inputs.r_off*nanometers,
                 constraints=inputs.cons,
                 ewaldErrorTolerance=inputs.ewald_Tol,
                 hydrogenMass = 4 * amu)
if inputs.vdw == 'Switch': nboptions['switchDistance'] = inputs.r_on*nanometers
if inputs.vdw == 'LJPME':  nboptions['nonbondedMethod'] = LJPME
if inputs.implicitSolvent:
    nboptions['implicitSolvent'] = inputs.implicitSolvent
    nboptions['implicitSolventSaltConc'] = inputs.implicit_salt*(moles/liter)
    nboptions['temperature'] = inputs.temp*kelvin
    nboptions['soluteDielectric'] = inputs.solut_diele
    nboptions['solventDielectric'] = inputs.solve_diele
    nboptions['gbsaModel'] = inputs.gbsamodel

if   args.fftype.upper() == 'CHARMM': system = top.createSystem(params, **nboptions)
elif args.fftype.upper() == 'AMBER':  
    system = top.createSystem(**nboptions)

# Define your two groups of atoms
print (inputs.pulling_group1)
print (inputs.pulling_group2)
group1_indices = inputs.pulling_group1
group2_indices = inputs.pulling_group2
restraint_force = openmm.CustomCentroidBondForce(2, "0.5*pull_strength*(distance(g1,g2)-equil)^2")
restraint_force.addGroup(group1_indices)
restraint_force.addGroup(group2_indices)

restraint_force.addGlobalParameter('equil', inputs.start_pull)

print(inputs.pull_strength)
restraint_force.addGlobalParameter('pull_strength', inputs.pull_strength)
restraint_force.addBond([0, 1], [])
system.addForce(restraint_force)

if inputs.vdw == 'Force-switch': system = vfswitch(system, top, inputs)
if inputs.lj_lrc == 'yes':
    for force in system.getForces():
        if isinstance(force, NonbondedForce): force.setUseDispersionCorrection(True)
        if isinstance(force, CustomNonbondedForce) and force.getNumTabulatedFunctions() != 1:
            force.setUseLongRangeCorrection(True)
if inputs.e14scale != 1.0:
    for force in system.getForces():
        if isinstance(force, NonbondedForce): nonbonded = force; break
    for i in range(nonbonded.getNumExceptions()):
        atom1, atom2, chg, sig, eps = nonbonded.getExceptionParameters(i)
        nonbonded.setExceptionParameters(i, atom1, atom2, chg*inputs.e14scale, sig, eps)

if inputs.pcouple == 'yes':      system = barostat(system, inputs)
if inputs.rest == 'yes':         system = restraints(system, crd, inputs)

print('done')

# Set platform
DEFAULT_PLATFORMS = 'CUDA', 'OpenCL', 'CPU'
enabled_platforms = [Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())]
if args.platform:
    if not args.platform[0] in enabled_platforms:
        print("Unable to find OpenMM platform '{}'; exiting".format(args.platform[0]), file=sys.stderr)
        sys.exit(1)

    platform = Platform.getPlatformByName(args.platform[0])
else:
    for platform in DEFAULT_PLATFORMS:
        if platform in enabled_platforms:
            platform = Platform.getPlatformByName(platform)
            break
    if isinstance(platform, str):
        print("Unable to find any OpenMM platform; exiting".format(args.platform[0]), file=sys.stderr)
        sys.exit(1)

print("Using platform:", platform.getName())
prop = dict(CudaPrecision='single') if platform.getName() == 'CUDA' else dict()

# Build simulation context
print('building system')

integrator = LangevinMiddleIntegrator(inputs.temp*kelvin, inputs.fric_coeff/picosecond, inputs.dt*picoseconds)
simulation = Simulation(top.topology, system, integrator, platform, prop)

if args.icrst:
    print('restarting from restart file ' + str (args.icrst))
    charmm_rst = read_charmm_rst(args.icrst)
    simulation.context.setPositions(charmm_rst.positions)
    simulation.context.setVelocities(charmm_rst.velocities)
    simulation.context.setPeriodicBoxVectors(charmm_rst.box[0], charmm_rst.box[1], charmm_rst.box[2])
else:
    simulation.context.setPositions(crd.positions)

if args.irst:
    print('restarting from restart file ' + str (args.irst))
    with open(args.irst, 'r') as f:
        simulation.context.setState(XmlSerializer.deserialize(f.read()))
if args.ichk:
    print('restarting from restart file' + str (args.ichk))
    with open(args.ichk, 'rb') as f:
        simulation.context.loadCheckpoint(f.read())

if args.restart_timer == True:
    print('resetting timer')
    simulation.context.setTime (0 * unit.picoseconds)

# Re-wrap
if args.rewrap:
    simulation = rewrap(simulation)

# Calculate initial system energy
print("\nInitial system energy")
print(simulation.context.getState(getEnergy=True).getPotentialEnergy())

# Energy minimization
if inputs.mini_nstep > 0:
    print("\nEnergy minimization: %s steps" % inputs.mini_nstep)
    simulation.minimizeEnergy(tolerance=inputs.mini_Tol*kilojoule/mole, maxIterations=inputs.mini_nstep)
    print(simulation.context.getState(getEnergy=True).getPotentialEnergy())

# Generate initial velocities
if inputs.gen_vel == 'yes':
    print("\nGenerate initial velocities")
    if inputs.gen_seed:
        simulation.context.setVelocitiesToTemperature(inputs.gen_temp, inputs.gen_seed)
    else:
        simulation.context.setVelocitiesToTemperature(inputs.gen_temp)


def checkpoint_dcd  (simulation):
    if inputs.nstdcd > 0:
        if not args.odcd: args.odcd = 'output.dcd'
        simulation.reporters.append(DCDReporter(args.odcd, inputs.nstdcd))
    simulation.reporters.append(
        StateDataReporter(sys.stdout, inputs.nstout, step=True, time=True, potentialEnergy=True, temperature=True, progress=True,
                          remainingTime=True, speed=True, totalSteps=inputs.nstep, separator='\t')
    )

pulling_file_name = inputs.pulling_out_file
print('Writing pulling coordinates to ' + str(pulling_file_name))

if args.restart_timer == True:
    with open (str(pulling_file_name),'w') as f : 
        f.write('')
    f.close()
        
def calculate_centroid_distance(indices1, indices2, state):

    positions = state.getPositions()
    centroid1 = np.mean([positions[i].value_in_unit(nanometers) for i in indices1], axis=0)
    centroid2 = np.mean([positions[i].value_in_unit(nanometers) for i in indices2], axis=0)
    distance = np.linalg.norm(centroid1 - centroid2)
    return distance

def write_pulling_coords (filename, *args): 
    with open (str(filename),'a') as f : 
        write_str = [repr (float(x)) for x in args] 
        for x in args:
            float_x = float(x)
            f.write(f"{float_x:12.8f}   ")
        f.write('\n')
        #f.write(f"{num1:10.3f} {num2:10.3f}\n"
        #write_str = ' '.join(write_str) + '\n'


        #f.write(write_str)
    f.close()
    
def write_wrapper (indices1, indices2, simulation, filename): 
    state = simulation.context.getState(getPositions=True)
    current_time = state.getTime().value_in_unit(nanoseconds)
    distance = calculate_centroid_distance (indices1,indices2, state) 
    write_pulling_coords (filename, current_time, distance)
        

def write_restart(simulation):
    # Write restart file
    print ('writing restart')
    if not (args.orst or args.ochk): args.orst = 'output.rst'
    if args.orst:
        state = simulation.context.getState( getPositions=True, getVelocities=True )
        with open(args.orst, 'w') as f:
            f.write(XmlSerializer.serialize(state))
    if args.ochk:
        with open(args.ochk, 'wb') as f:
            f.write(simulation.context.createCheckpoint())
    if args.opdb:
        crd = simulation.context.getState(getPositions=True).getPositions()
        PDBFile.writeFile(top.topology, crd, open(args.opdb, 'w'))


start = inputs.start_pull
end = inputs.end_pull
pulling_steps = inputs.pull_step
reporting_interval = inputs.reporting_interval
print(reporting_interval)
restraint_steps =  np.linspace(start,end,int(pulling_steps/reporting_interval))
simulation.context.setParameter('equil', start)

checkpoint_dcd(simulation)

print("\nMD run: %s steps" % inputs.nstep)
i = 0

start_steps = simulation.context.getStepCount() 

if inputs.nstdcd % inputs.reporting_interval != 0 :
    raise ValueError('dcd reporting must be a multiple of umbrella sampling cv reporting')

run = False

while simulation.context.getStepCount() - start_steps < pulling_steps : 
    run = True
    simulation.step(inputs.reporting_interval)
    simulation.context.setParameter('equil', restraint_steps[i])
    #t = threading.Thread(target=write_wrapper, args=(group1_indices, group2_indices, simulation, pulling_file_name))
    #t.start()
    write_wrapper(group1_indices, group2_indices, simulation, pulling_file_name)

    if (simulation.context.getStepCount() - start_steps) % inputs.nstdcd == 0  :
        #t = threading.Thread(target=write_restart,args=(simulation,))
        #t.start ()
        write_restart(simulation)
    i = i + 1

if run == True:
    print('pulling finished')
else:
    print('no pulling needed')
simulation.context.setParameter('equil', end)

while (simulation.context.getStepCount() - start_steps) < inputs.nstep : 
    simulation.step(inputs.reporting_interval)
    #t = threading.Thread(target=write_wrapper, args=(group1_indices, group2_indices, simulation, pulling_file_name))
    #t.start()
    write_wrapper(group1_indices, group2_indices, simulation, pulling_file_name)
    if (simulation.context.getStepCount() - start_steps) % inputs.nstdcd == 0  :
        #t = threading.Thread(target=write_restart,args=(simulation,))
        #t.start ()
        write_restart(simulation)
 
