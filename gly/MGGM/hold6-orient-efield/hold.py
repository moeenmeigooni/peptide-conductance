from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout, exit, stderr, argv
import mdtraj as md
import numpy as np

###################

input_prefix  = 'MGGM'
output_prefix = 'hold'
_, rep = argv

###################
# make sure n_run_steps is divisible by chk_freq
n_run_steps = 50000000 # 200 ns
chk_freq    = 100000 
dcd_freq    = 1000
state_freq  = 1000

###################

# load CHARMM files
psf_filename = '../build/%s_water_%s.psf'%(input_prefix, rep)
pdb_filename = '../build/%s_water_%s.pdb'%(input_prefix, rep)
psf = CharmmPsfFile(psf_filename)
pdb = PDBFile(pdb_filename)

toppar_dir = '/Scr/meigoon2/polymer/toppar'
param_files = ['par_all36m_prot.prm',
               'par_all36_na.prm',
               'par_all36_carb.prm',
               'par_all36_lipid.prm',
               'par_all36_cgenff.prm',
               'toppar_water_ions.str']
param_list = ['%s/%s'%(toppar_dir, param_file) for param_file in param_files]
params = CharmmParameterSet(*param_list)

# get relevant atom indices for custom forces
t = md.load(pdb_filename) # box info is preserved if pdb is written with mdtraj
sequence = [residue.name for residue in t.top.residues if residue.chain.index == 0]
n_residues = len(sequence)
n_atoms = t.n_atoms
protein_ind      = [int(ind) for ind in t.top.select('protein')] # indices must be int, not np.int64
s_ind1, s_ind2   = [int(ind) for ind in t.top.select('name SD and resid 0 %i'%(n_residues - 1))]
cg_ind1, cg_ind2 = [int(ind) for ind in t.top.select('name CG and resid 0 %i'%(n_residues - 1))]
ce_ind1, ce_ind2 = [int(ind) for ind in t.top.select('name CE and resid 0 %i'%(n_residues - 1))]

# get box information
boxX, boxY, boxZ = t.unitcell_lengths[0] 
psf.setBox(boxX*nanometer, boxY*nanometer, boxZ*nanometer)

# create system
system = psf.createSystem(params, nonbondedMethod=PME, 
                          nonbondedCutoff=1.2*nanometer, switchDistance=1.*nanometer,
                          constraints=HBonds, rigidWater=True)

# add S-S force via CustomCompoundBondForce
sspull_constant    = 418.4 # units: kJ/(mol*nm^2)
sspull_target_dist = 0.6   # units: nm
sspull = CustomCompoundBondForce(2, '0.5*k*((z2 - z1) - targetdist)^2')
sspull.addGlobalParameter('k', sspull_constant)
sspull.addGlobalParameter('targetdist', sspull_target_dist)
sspull.addBond([s_ind1, s_ind2], [])
system.addForce(sspull)


# add force orienting sulfur lone pair to z axis 
ssorient_constant = 4.184 * 10. # 10 kcal/mol when unaligned, -10 kcal/mol when aligned
ssorient = CustomCompoundBondForce(3, 'dir*k*(z1 - (z2+z3)/2)/pointdistance(x1, y1, z1, (x2+x3)/2, (y2+y3)/2, (z2+z3)/2)')
ssorient.addPerBondParameter('dir')
ssorient.addGlobalParameter('k', ssorient_constant)
ssorient.addBond([s_ind1, cg_ind1, ce_ind1], [1.])
ssorient.addBond([s_ind2, cg_ind2, ce_ind2], [-1.])
system.addForce(ssorient)

# add electric field
applied_potential   = 0.1       # Volts
separation_distance = sspull_target_dist  # nm
au_s_bond_length    = 0.25      # nm
avogadro_number     = 6.022e23
electron_volt       = 1.602e-19 # C/e-
J_to_kJ             = 0.001
efield_strength     = applied_potential * electron_volt * avogadro_number * J_to_kJ / (separation_distance + 2. * au_s_bond_length)
efield = CustomExternalForce('q*Ez*z')
efield.addPerParticleParameter('q')
efield.addGlobalParameter('Ez', efield_strength) # kJ/(mol*nm*e-)
nonbonded = [f for f in system.getForces() if isinstance(f, NonbondedForce)][0]
for i in range(n_atoms):
    charge, sigma, epsilon = nonbonded.getParticleParameters(i)
    efield.addParticle(i, [charge])
system.addForce(efield)

# NPT stuff
system.setDefaultPeriodicBoxVectors(*psf.boxVectors)
system.addForce(MonteCarloBarostat(1*atmosphere, 300*kelvin))

# set up for GPU
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 4*femtoseconds)
platform = Platform.getPlatformByName('CUDA')
properties = {'DeviceIndex': rep, 'Precision': 'mixed'}
simulation = Simulation(psf.topology, system, integrator, platform, properties)
simulation.context.setPositions(pdb.positions)

# do quick minimization
simulation.minimizeEnergy()#maxIterations=10000)

# set up reporters
simulation.reporters.append(DCDReporter('%s_%s.dcd'%(output_prefix, rep), dcd_freq, enforcePeriodicBox=True))
simulation.reporters.append(StateDataReporter(stdout, state_freq, step=True, potentialEnergy=True, 
                                              temperature=True, speed=True, elapsedTime=True))
simulation.reporters.append(CheckpointReporter('%s_%s.chk'%(output_prefix, rep), chk_freq))

# run
n_step_per_cycle = n_run_steps // chk_freq
n_cycle = n_run_steps // n_step_per_cycle
for i in range(n_cycle):
    cycle_finished = False 
    while(not cycle_finished):
        try:
            simulation.step(n_step_per_cycle)
            cycle_finished = True
        except (ValueError, OpenMMException):
            simulation.loadCheckpoint('%s_%s.chk'%(output_prefix, rep))


