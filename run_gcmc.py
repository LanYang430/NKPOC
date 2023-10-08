#!/usr/bin/env python
import sys
import re
from sys import stdout
import numpy as np
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
# from heat_flux import *
from run_md import *
import copy
import scipy
from scipy import special
import pickle

# this function places a molecule into a box in random position and random velocity
def get_rand_pos(positions, box):
    r_center = np.average(positions, axis=0)
    r_relat = positions - r_center
    s = np.random.random(3)
    R = scipy.stats.special_ortho_group.rvs(3)
    return R.dot(r_rela) + box.dot(s)

# volume should carry unit of nanometers**3
def pressure_correction(vol, n, n_max, kT):
    n_ghost = int(n_max - n)
    P = n_ghost*kT/vol
    return np.ones(3) * P

# note that the ghost gas molecules still contribute to the pressure
# need to deduct the pressure
def update_P(simulation, P, n, n_max, kT):
    state = simulation.context.getState()
    vol = state.getPeriodicBoxVolume()
    P_corr = pressure_correction(vol, n, n_max, kT)
    barostat = simulation.system.getForces()[-1]
    barostat.setDefaultPressure(P + P_corr)
    return

def get_gas_params(system_gas):
    params = {'q':[], 'sigma':[], 'epsil':[]}
    nbforce = system_gas.getForce(int(np.where([('NonbondedForce' in f.__repr__()) for f in system_gas.getForces()])[0][0]))
    for i in range(system_gas.getNumParticles()):
        p = nbforce.getParticleParameters(i)
        params['q'].append(p[0])
        params['sigma'].append(p[1])
        params['epsil'].append(p[2])
    return params


def modify_particle_param(qqforce, ljforce, context, i, q, sigma, epsil, lambd_sc):
    qqforce.setParticleParameters(i, q, 1.0, 0.0)
    sigma = sigma._value
    epsil = epsil._value
    param = [np.sqrt(4*epsil*sigma**6), np.sqrt(4*epsil*sigma**12), lambd_sc]
    ljforce.setParticleParameters(i, param)
    qqforce.updateParametersInContext(context)
    ljforce.updateParametersInContext(context)
    return 

# a lambda-MC step
# n_atoms = (n_atoms_cages, n_atoms_per_gas
# it will update simulation, system, and return new value of n, and mol_indices
def monte_carlo_lambda(simulation, params_gas, mu, therm_wavelength, n, dn, n_max, kT, P_ext, n_atoms, mol_indices, dic_forces):
    s = np.random.randint(2)
    state = simulation.context.getState(getPositions=True, getEnergy=True)
    vol = state.getPeriodicBoxVolume()
    pot0 = state.getPotentialEnergy()
    qqforce = dic_forces['NonbondedForce']
    ljforce = dic_forces['CustomNonbondedForce']
    n_atoms_cages, n_atoms_per_gas = n_atoms
    fac = vol/therm_wavelength**3
    beta = 1 / (kT*AVOGADRO_CONSTANT_NA)
    # if we happen to have integer number of molecule, shuffle the molecule order
    if abs(n - int(n)) < 1e-8:
        np.random.shuffle(mol_indices[:int(n)])
        np.random.shuffle(mol_indices[int(n):])
    # decrease n
    if s == 0:
        # zero particle, reject
        if abs(n - 0.0) < 1e-8:
            return 0, n, mol_indices
        i_partial = mol_indices[int(np.ceil(n)) - 1]
        lambd = n - (np.ceil(n) - 1)
        n_new = n - dn
        lambd_new = lambd - dn
        ib = n_atoms_cages + i_partial*n_atoms_per_gas
        ie = ib + n_atoms_per_gas
        delta_n = -dn
    # increase n
    else:
        if abs(n - n_max) < 1e-8:
            print('Warning: exceed maximum number of particles')
            return 0, n, mol_indices
        i_partial = mol_indices[int(np.floor(n))]
        lambd = n - np.floor(n)
        n_new = n + dn
        lambd_new = lambd + dn
        ib = n_atoms_cages + i_partial*n_atoms_per_gas
        ie = ib + n_atoms_per_gas
        delta_n = dn
    # in principle the lambda should be non-negative, but due to the machine precision, 
    # it can also be a very tiny negative number, which causes issues
    if lambd_new < 0:
        lambd_new = 0.0
    for ii, i in enumerate(range(ib, ie)):
        lambd_qq, lambd_lj = get_lambds(lambd_new)
        q = params_gas['q'][ii] * lambd_qq
        sigma = params_gas['sigma'][ii]
        epsil = params_gas['epsil'][ii] * lambd_lj
        lambd_sc = lambd_lj
        modify_particle_param(qqforce, ljforce, simulation.context, i, q, sigma, epsil, lambd_sc)
    state = simulation.context.getState(getPositions=True, getEnergy=True)
    pot1 = state.getPotentialEnergy()
    dU = pot1 - pot0
    prob = np.exp(special.loggamma(n+1)-special.loggamma(n_new+1)) * (fac**(delta_n)) * np.exp(beta * (mu*delta_n - dU))
    # success
    if np.random.random() < prob:
        # reset pressure correction if necessary (if we completely deleted a molecule, or we start to add a new particle)
        if (s == 0 and abs(n_new - int(n_new)) < 1e-8) or (s == 1 and abs(n - int(n)) < 1e-8): 
            update_P(simulation, P_ext, n, n_max, kT)
        return 1, n_new, mol_indices
    else:
        # reset the system
        for ii, i in enumerate(range(ib, ie)):
            lambd_qq, lambd_lj = get_lambds(lambd)
            q = params_gas['q'][ii] * lambd_qq
            sigma = params_gas['sigma'][ii]
            epsil = params_gas['epsil'][ii] * lambd_lj
            lambd_sc = lambd_lj
            modify_particle_param(qqforce, ljforce, simulation.context, i, q, sigma, epsil, lambd_sc)
        return 0, n, mol_indices

def load_chk(ifn):
    with open(ifn, 'rb') as f:
        data = pickle.load(f)
    pos = data['pos']*nanometers
    box = data['box'] * nanometers
    n = data['n']
    mol_indices = data['mol_indices']
    n_max = len(mol_indices)
    return pos, box, n, mol_indices, n_max

if __name__ == '__main__':
    #-------- the system parameter --------
    n_max = 128
    restart = None  # a fresh calculation, otherwise the restart chk file name 
    n_atoms_per_gas = 5
    n_cages = 32
    n_atoms_per_cage = 122
    n_atoms_cage = n_cages * n_atoms_per_cage
    n = 0
    dn = 1.0           # time step for n
    P = 0.1 * bar      # pressure of CO2, used to determine \mu
    P_ext = 1.0 * bar  # total external pressure
    T = 298 * kelvin
    kT = BOLTZMANN_CONSTANT_kB * T
    freq_mc = 50
    n_cycles = 10000000
    f_cycle_nout = 10  # print out number of particles every 10 cycles
    f_cycle_nchk = 100 # checkpoint every 100 cycles
    platform_name = 'CUDA'
    n_eq = 10000  
    dt = 1.0 # fs
    nout_log = 10000
    nout_traj = 10000
    mol_indices = np.array(range(n_max))
    np.random.shuffle(mol_indices)
    #-------- the co2 molecule --------
    mol_gas = PDBFile('co2.pdb')
    ff_gas = ForceField('co2.xml')
    system_gas = ff_gas.createSystem(mol_gas.topology, nonbondedMethod=PME, nonbondedCutoff=1.0*nanometer, constraints=None, removeCMMotion=True)
    params_gas = get_gas_params(system_gas)
    system_gas.addConstraint(0, 1, 0.195948) # rigid epm2 model for co2
    m_gas = 0.0 * dalton
    for i in range(system_gas.getNumParticles()):
        m_gas += system_gas.getParticleMass(i)
    PLANCK_CONSTANT_h = 6.62607004e-34 * joules * seconds
    # lambda = h / sqrt{2*pi*m*kT}
    therm_wavelength = np.sqrt(PLANCK_CONSTANT_h**2 / 2 / np.pi / (m_gas/AVOGADRO_CONSTANT_NA) / kT)
    # chemical potential computed as idea gas: \mu - \mu_rot
    # -kT * ln(V/N/lambda^3)
    mu = -kT * np.log(kT/P/therm_wavelength**3)
    mu = mu*AVOGADRO_CONSTANT_NA
    #---------- create systems and simulations -----------
    # mol, top, system = create_system('supercell.pdb', 'supercell.top', nbparams=None)
    mol, top, system = create_system('supercell.pdb', 'supercell.top', nbparams='cage_nb_redefined.itp')
    # if it is a restart, then overwrite previous geometry and configuration settings with the settings in the restart file
    if restart is not None:
        pos, box, n, mol_indices, n_max = load_chk(restart)
        print('Restart file %s is read!'%restart)
        stdout.flush()
    # insert gas molecules, and read the init geom
    n_floor = int(np.floor(n))
    lambd = n - n_floor
    lambdas = list(np.ones(n_atoms_cage))
    for i in range(n_max):
        ip = np.where(mol_indices==i)[0][0]
        if ip < n_floor:
            system = insert_mol_to_sys(system, system_gas)
            lambdas.extend(list(np.ones(n_atoms_per_gas)*1.0))
        elif ip == n_floor:
            system = insert_mol_to_sys(system, system_gas, lambd=lambd)
            lambdas.extend(list(np.ones(n_atoms_per_gas)*lambd))
        else:
            system = insert_mol_to_sys(system, system_gas, lambd=0)
            lambdas.extend(list(np.ones(n_atoms_per_gas)*0.0))
    system = set_soft_lj(system, type=1, lambdas=lambdas)
    label = padding(n_max)
    mol = PDBFile('geometries/cage_%s.pdb'%label)
    # set up simulation context
    platform = Platform.getPlatformByName(platform_name)
    if platform_name == 'CPU':
        properties = {}
    else:
        properties = {'CudaPrecision': 'single'}
    integrator = LangevinIntegrator(T, 1.0/picosecond, dt*femtosecond)
    # integrator.setRandomNumberSeed(12345)
    # correct the pressure from the ghost molecules
    system.addForce(MonteCarloAnisotropicBarostat(np.array([P_ext, P_ext, P_ext])*bar, T, True, True, True, 5))
    simulation = Simulation(mol.topology, system, integrator, platform, properties)
    if restart is None:
        simulation.context.setPositions(mol.positions)
    else:
        simulation.context.setPositions(pos)
        simulation.context.setPeriodicBoxVectors(box[0], box[1], box[2])
    simulation.context.setVelocitiesToTemperature(T)
    update_P(simulation, P_ext, n, n_max, kT)
    # dictionary of forces
    dic_forces = set_force_groups(simulation.system)
    # ---------- outputs -----------
    # report n and mol_indices
    class GCMCStatusReporter(object):
        def __init__(self, file, reportInterval):
            self._out = open(file, 'w')
            self._reportInterval = reportInterval
    
        def __del__(self):
            self._out.close()
    
        def describeNextReport(self, simulation):
            steps = self._reportInterval - simulation.currentStep%self._reportInterval
            return (steps, False, False, True, False)
    
        def report(self, simulation, state):
            t = state.getTime()._value
            self._out.write('%f: %g\n'%(t, n))
            for i in mol_indices:
                self._out.write('%d '%i)
            self._out.write('\n')
            self._out.flush()

    simulation.reporters.append(StateDataReporter(stdout, nout_log, time=True, 
       potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
       temperature=True, density=True))
    simulation.reporters.append(DCDReporter('traj.dcd', nout_traj))
    simulation.reporters.append(GCMCStatusReporter('gcmc_status.dat', nout_traj))

    if restart is None:
        simulation.minimizeEnergy()
        simulation.step(n_eq)
    #---------- start running -----------
    n_atoms = (n_atoms_cage, n_atoms_per_gas)
    ofile = open('npart.xvg', 'w')
    for i_cycle in range(n_cycles):
        simulation.step(freq_mc)
        succ, n, mol_indices = monte_carlo_lambda(simulation, params_gas, mu, therm_wavelength, n, dn, n_max, kT, P_ext, n_atoms, mol_indices, dic_forces)
        if i_cycle%f_cycle_nout == 0:
            print(i_cycle, n, file=ofile)
            ofile.flush()
        if i_cycle%f_cycle_nchk == 0:
            chk = open('chk.pickle', 'wb')
            state = simulation.context.getState(getPositions=True, getEnergy=True)
            if state.getPotentialEnergy()._value >  1e8:
                sys.exit('Energy exploded ...')
            data = {}
            data['pos'] = np.array(state.getPositions()._value)
            data['box'] = np.array(state.getPeriodicBoxVectors()._value)
            data['n'] = n
            data['mol_indices'] = mol_indices
            pickle.dump(data, chk)
            chk.close()
    ofile.close()



