#!/usr/bin/env python
import sys
import re
from sys import stdout
import numpy as np
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
# from heat_flux import *

def get_connectivity(top):
    n_atoms = top.topology.getNumAtoms()
    connectivity = np.zeros((n_atoms, n_atoms))
    for bond in top.topology.bonds():
        i0 = bond.atom1.index
        i1 = bond.atom2.index
        connectivity[i0, i1] = 1
        connectivity[i1, i0] = 1
    return connectivity

# 1-2 and 1-3 exclusions
def get_exclusion_lists(top):
    connectivity = get_connectivity(top)
    n_atoms = top.topology.getNumAtoms()
    exclusions = []
    for i in range(n_atoms):
        first_neighbours = np.where(connectivity[i]==1)[0]
        for j in first_neighbours:
            if j != i:
                exclusions.append((i, j))
                second_neighbours = np.where(connectivity[j]==1)[0]
                for k in second_neighbours:
                    if k != i:
                        exclusions.append((i, k))
    return exclusions

def read_nb_params(ifn):
    params = {}
    with open(ifn, 'r') as f:
        for line in f:
            words = line.split()
            if len(words) == 0:
                continue
            if line.startswith(';') or words[0] == ';' or '[' in line:
                continue
            label = (words[0], words[1])
            sigma = float(words[3])
            epsil = float(words[4])
            C = 4 * epsil * sigma**6
            A = 4 * epsil * sigma**12
            params[label] = (C, A)
    return params

def get_lambds(lambd):
    if lambd >= 0.5:
        lambd_qq = (lambd - 0.5) * 2.0
        lambd_lj = 1.0
        if lambd_qq > 1.0:
            lambd_qq = 1.0
        elif lambd_qq < 0.0:
            lambd_qq = 0.0
    else:
        lambd_qq = 0.0
        lambd_lj = lambd * 2
        if lambd_lj < 0.0:
            lambd_lj = 0.0
        elif lambd_lj > 1.0:
            lambd_lj = 1.0
    return lambd_qq, lambd_lj

# set force groups and get a dictionary of forces
def set_force_groups(system, verbose=False):
    dic_forces = {}
    forces = system.getForces()
    for i, f in enumerate(forces):
        name = re.search('openmm.([A-Za-z]+);', f.__repr__()).group(1)
        if verbose:
            print(i, f)
        f.setForceGroup(i)
        dic_forces[name] = f
    return dic_forces

def create_system(pdb_file, top_file, nbparams='cage_nb_redefined.itp'):
    # set up system and topology
    mol = PDBFile(pdb_file)
    top = GromacsTopFile(top_file, periodicBoxVectors=mol.topology.getPeriodicBoxVectors())
    exclusions = get_exclusion_lists(top)
    system = top.createSystem(nonbondedMethod=PME, nonbondedCutoff=1.0*nanometer)
    dic_forces = set_force_groups(system)
    # print(dic_forces)
    # -------------------------------------------------------------------
    # ----------------- start my adjust to the forces -------------------
    # -------------------------------------------------------------------
    qqforce = dic_forces['NonbondedForce']
    ljforce = dic_forces['CustomNonbondedForce']
    lj14force = dic_forces['CustomBondForce']
    bondforce = dic_forces['HarmonicBondForce']
    angleforce = dic_forces['HarmonicAngleForce']
    dihedforce = dic_forces['RBTorsionForce']
    if nbparams is not None:
        nbtype_list = []
        nbtype_dic = {}
        index = 0
        # construct the nonbonding atomtype list
        for name_mol, n_mol in top._molecules:
            for i_mol in range(n_mol):
                for atom in top._moleculeTypes[name_mol].atoms:
                    nbtype = atom[1]
                    nbtype_list.append(nbtype)
                    if nbtype in nbtype_dic:
                        nbtype_dic[nbtype].append(index)
                    else:
                        nbtype_dic[nbtype] = [index]
                    index += 1
        nbtype_list = np.array(nbtype_list)
        nb_params = read_nb_params(nbparams) # nb parameters defined by user
        involved_atypes = []
        for label in nb_params:
            involved_atypes.append(label[0])
            involved_atypes.append(label[1])
        # intra-cage pairs
        nb_params_intra = {}
        for label in nb_params.keys():
            at0, at1 = label
            imol0 = int(at0.split('_')[-1])
            imol1 = int(at1.split('_')[-1])
            if imol0 == imol1:  # within the same cage, possibly involving 1-4 pairs
                nb_params_intra[label] = nb_params[label]
        # deal with modified 1-4 pair, also keep a record which pairs are already excluded
        pairs_excluded = []
        for ip in range(lj14force.getNumBonds()):
            i0, i1, p = lj14force.getBondParameters(ip)
            at0 = nbtype_list[i0]
            at1 = nbtype_list[i1]
            label = (at0, at1)
            if label in nb_params_intra:
                p = nb_params_intra[label]
                p_scaled = (p[0]/2, p[1]/2)
                lj14force.setBondParameters(ip, i0, i1, p_scaled)
                pairs_excluded.append(label)
                pairs_excluded.append(label[::-1])
            elif label[::-1] in nb_params_intra:
                p = nb_params_intra[label[::-1]]
                p_scaled = (p[0]/2, p[1]/2)
                lj14force.setBondParameters(ip, i0, i1, p_scaled)
                pairs_excluded.append(label)
                pairs_excluded.append(label[::-1])
        for i in range(ljforce.getNumExclusions()):
            i0, i1 = ljforce.getExclusionParticles(i)
            at0 = nbtype_list[i0]
            at1 = nbtype_list[i1]
            if at0 in involved_atypes and at1 in involved_atypes:
                pairs_excluded.append((i0, i1))
                pairs_excluded.append((i1, i0))
        # deal with the true nonbonded pairs, that has not been excluded
        ##############################################################################
        # WARNING: Our current approach to deal with pairwise LJ modifications       #
        # does not account for cutoff scheme. That is, the correction is calculated  #
        # for all LJ interactions, regardless of their distance.                     #
        # This leads to ~3 kJ/mol energy difference per cage with gromacs            #
        # I assume this problem is not serious for a somewhat qualitative calculation#
        # for NKPOC                                                                  #
        ##############################################################################
        for label in nb_params.keys():
            at0, at1 = label
            p = nb_params[label]
            for i0 in nbtype_dic[at0]:
                p0_std = ljforce.getParticleParameters(i0)
                for i1 in nbtype_dic[at1]:
                    p1_std = ljforce.getParticleParameters(i1)
                    if at0[-1] == at1[-1]:
                        # if it has been already taken cared by the 1-4 part
                        if (i0, i1) in pairs_excluded: 
                            continue
                        # same atomtype, then halve our task
                        elif at0 == at1 and i0 >= i1:
                            continue
                    # give a special treatment
                    p_corr = (p[0]-p0_std[0]*p1_std[0], p[1]-p0_std[1]*p1_std[1])
                    # p_corr = (-p0_std[0]*p1_std[0], -p0_std[1]*p1_std[1])
                    lj14force.addBond(i0, i1, p_corr)
    # -------------------------------------------------------------------
    # --------------------- end of my adjustment ------------------------
    # -------------------------------------------------------------------
    qqforce.setEwaldErrorTolerance(1e-4)
    ljforce.setUseLongRangeCorrection(False)
    ljforce.setUseSwitchingFunction(False)
    ljforce.setSwitchingDistance(0.9*nanometer)
    lj14force.setUsesPeriodicBoundaryConditions(True)
    bondforce.setUsesPeriodicBoundaryConditions(True)
    angleforce.setUsesPeriodicBoundaryConditions(True)
    dihedforce.setUsesPeriodicBoundaryConditions(True)

    return mol, top, system

# add a new gas molecule to the existing topology
def insert_mol_to_sys(system, system_mol, lambd=1.0):
    n0 = system.getNumParticles()
    dn = system_mol.getNumParticles()
    index_map_o2n = {}     # old index to new index
    index_map_n2o = {}     # new index to old index
    # adding particles
    for i in range(dn):
        m = system_mol.getParticleMass(i)
        i_new = system.addParticle(m)
        index_map_o2n[i] = i_new
        index_map_n2o[i_new] = i
    # print(index_map_o2n)
    # define the new virtual sites, if any
    for i in range(dn):
        if system_mol.isVirtualSite(i):
            i_new = index_map_o2n[i]
            virt = system_mol.getVirtualSite(i)
            if not (isinstance(virt, TwoParticleAverageSite) or isinstance(virt, ThreeParticleAverageSite)):
                sys.exit('ERROR: only supports 2/3 particle average sites')
            indices_new = []
            weights = []
            for i_virt in range(virt.getNumParticles()):
                indices_new.append(index_map_o2n[virt.getParticle(i_virt)])
                weights.append(virt.getWeight(i_virt))
            if isinstance(virt, TwoParticleAverageSite):
                virt_new = TwoParticleAverageSite(indices_new[0], indices_new[1], weights[0], weights[1])
            elif isinstance(virt, ThreeParticleAverageSite):
                virt_new = ThreeParticleAverageSite(indices_new[0], indices_new[1], indices_new[2], weights[0], weights[1], weights[2])
            # print(i_new, virt_new.getParticle(0), virt_new.getParticle(1), virt_new.getWeight(0), virt_new.getWeight(1))
            system.setVirtualSite(i_new, virt_new)
    # set atomic charges, and lj parameters
    dic_forces = set_force_groups(system)
    dic_forces_mol = set_force_groups(system_mol)
    qqforce = dic_forces['NonbondedForce']
    ljforce = dic_forces['CustomNonbondedForce']
    lj14force = dic_forces['CustomBondForce']
    bondforce = dic_forces['HarmonicBondForce']
    angleforce = dic_forces['HarmonicAngleForce']
    dihedforce = dic_forces['RBTorsionForce']
    nbforce_mol = dic_forces_mol['NonbondedForce']
    for i in range(dn):
        i_new = index_map_o2n[i]
        q, sigma, epsil = nbforce_mol.getParticleParameters(i)
        # insert partial particles
        lambd_qq, lambd_lj = get_lambds(lambd)
        q = q * lambd_qq
        sigma = sigma._value
        epsil = epsil._value * lambd_lj
        i_new_check = qqforce.addParticle(q, 1.0, 0.0)
        if i_new_check != i_new:
            sys.exit('ERROR: particle number issue in nb force')
        # note the energy function is: 'A1*A2/r^12-C1*C2/r^6'
        # and the parameters are (C, A)
        param_lj = [np.sqrt(4*epsil*sigma**6), np.sqrt(4*epsil*sigma**12)]
        i_new_check = ljforce.addParticle(param_lj)
        if i_new_check != i_new:
            sys.exit('ERROR: particle number issue in nb force')
    # set exclusions and pair exceptions
    for i_except in range(nbforce_mol.getNumExceptions()):
        i0, i1, qq, sigma, epsil = nbforce_mol.getExceptionParameters(i_except)
        sigma = sigma._value
        epsil = epsil._value
        i0_new = index_map_o2n[i0]
        i1_new = index_map_o2n[i1]
        # qqforce exception deals with only q-q interactions
        qqforce.addException(i0_new, i1_new, qq, sigma, 0.0)
        ljforce.addExclusion(i0_new, i1_new)
        # print(i0_new, i1_new, qq, sigma, epsil)
        if epsil > 1e-8:
            # parameter is C, A, '-C/r^6+A/r^12'
            param = (4*epsil*sigma**6, 4*epsil*sigma**12)
            lj14force.addBond(i0_new, i1_new, param)
    # set bonds and angles and dihedrals
    if 'HarmonicBondForce' in dic_forces_mol:
        bondforce_mol = dic_forces_mol['HarmonicBondForce']
        for ib in range(bondforce_mol.getNumBonds()):
            i0, i1, r0, k = bondforce_mol.getBondParameters(ib)
            i0_new = index_map_o2n[i0]
            i1_new = index_map_o2n[i1]
            bondforce.addBond(i0_new, i1_new, r0._value, k._value)
    # WARNING
    # All the angle, dihedral stuff is to be implemented ...
    # add constraints
    for i in range(system_mol.getNumConstraints()):
        i0, i1, r0 = system_mol.getConstraintParameters(i)
        i0_new = index_map_o2n[i0]
        i1_new = index_map_o2n[i1]
        system.addConstraint(i0_new, i1_new, r0)
    return system


def set_soft_lj(system, type=0, lambdas=None):
    dic_forces = set_force_groups(system)
    ljforce = dic_forces['CustomNonbondedForce']
    if type == -1:
        pass
    elif type == 0:
        ljforce.setEnergyFunction('A1*A2/r^3/r^3/r^3/r^3-C1*C2/r^3/r^3')
    elif type == 1: # the soft-core from gromacs, need to specify scaling factors for each atom
        ljforce.setEnergyFunction('A1*A2/(r^6+0.000001*(1-S1*S2))^2 - C1*C2/(r^6+0.000001*(1-S1*S2))')
        ljforce.addPerParticleParameter('S')
        n_atoms = ljforce.getNumParticles()
        for i_atom in range(n_atoms):
            lambd_qq, lambd_lj = get_lambds(lambdas[i_atom])
            lambd_sc = lambd_lj
            p_old = ljforce.getParticleParameters(i_atom)
            p_new = (p_old[0], p_old[1], lambd_lj)
            ljforce.setParticleParameters(i_atom, p_new)
    elif type == 2:
        ljforce.setEnergyFunction('A1*A2/(r+0.001)^12-C1*C2/(r+0.001)^6')
    return system


def padding(n):
    s = '%d'%n
    while len(s) < 3:
        s = '0' + s
    return s

if __name__ == '__main__':
    # create system
    # mol, top, system = create_system('supercell.pdb', 'supercell.top', nbparams=None)
    mol, top, system = create_system('supercell.pdb', 'supercell.top', nbparams='cage_nb_redefined.itp')
    mol_co2 = PDBFile('co2.pdb')
    ff_co2 = ForceField('co2.xml')
    system_co2 = ff_co2.createSystem(mol_co2.topology, nonbondedMethod=PME, nonbondedCutoff=1.0*nanometer, constraints=None, removeCMMotion=True)
    system_co2.addConstraint(0, 1, 0.195948) # rigid epm2 model for co2
    # insert gas molecules, and read the init geom
    n_co2 = 0
    label = padding(n_co2)
    mol_cage_co2 = PDBFile('geometries/cage_%s.pdb'%label)
    for i in range(n_co2):
        system = insert_mol_to_sys(system, system_co2)
    # run an NPT simulation, anisotropic coupling
    system.addForce(MonteCarloAnisotropicBarostat(np.array([1.0, 1.0, 1.0])*bar, 300.*kelvin, True, True, True, 10))
   
    # set up simulation
    platform = Platform.getPlatformByName('CUDA')
    properties = {'CudaPrecision': 'mixed'}
#     platform = Platform.getPlatformByName('CPU')
#     properties = {}
    integrator = LangevinIntegrator(300*kelvin, 1.0/picosecond, 1.0*femtosecond)
    simulation = Simulation(mol_cage_co2.topology, system, integrator, platform, properties)
    simulation.context.setPositions(mol_cage_co2.positions)

    # check energies
    energies = []
    for i in range(system.getNumForces()):
        state = simulation.context.getState(getEnergy=True, getPositions=True, groups=2**i)
        energies.append(state.getPotentialEnergy().value_in_unit(kilojoules_per_mole))
        print(energies[i])
    print('total energy:', np.sum(np.array(energies)))

    # outputs
    simulation.reporters.append(StateDataReporter(stdout, 100, time=True, 
       potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
       temperature=True, density=True))
    simulation.reporters.append(DCDReporter('traj.dcd', 100))

    # run
    simulation.step(100000)


