# SPDX-License-Identifier: BSD-3-Clause-LBNL
# Derived from: https://github.com/fbpic/fbpic/blob/dev/fbpic/picmi/__init__.py
# and from :    https://github.com/fbpic/fbpic/blob/dev/fbpic/picmi/simulation.py
# Copyright (c) 2016–2025 FBPIC contributors
# Additional U.S. DOE/LBNL notice may apply; see upstream LICENSE.txt/LEGAL.txt.
# Modifications (c) 2025 Stefan Costea, LeCAD-PEG

"""
This file is part of the Heterogenous EXAscale Particle-In-Cell code (HEXAPIC).
It defines the picmi standard interface.
"""

# Define general variables that each PICMI code should define
codename = 'hexapic'

# Check that the `picmistandard` package has been installed
try:
    from picmistandard.base import register_codename, register_constants
    register_codename(codename)
except ImportError:
    raise ImportError(
        "In order to use HEXAPIC with PICMI, you should install the \n"
        "`picmistandard` package, e.g. with: `pip install picmistandard`")

from numpy import min as np_min, sqrt, array

from scipy import constants
class constants:
    # Put the constants in their own namespace
    c = constants.c
    ep0 = constants.epsilon_0
    mu0 = constants.mu_0
    q_e = constants.e
    m_e = constants.m_e
    m_p = constants.m_p
    m_u = constants.atomic_mass
register_constants(constants)

# Create dictionaries with species_type defined in openPMD2
particle_charge = { 'electron': -constants.q_e,
                    'positron': constants.q_e,
                    'H+': constants.q_e,
                    'H': 0.0                   }

particle_mass = {   'electron': constants.m_e,
                    'positron': constants.m_e,
                    'H+': constants.m_p,
                    'H':1.00784*constants.m_u  }

# Import picmi objects
from picmistandard import PICMI_Simulation
from picmistandard import PICMI_ElectrostaticSolver
from picmistandard import PICMI_Cartesian2DGrid
from picmistandard import PICMI_ConstantAppliedField
from picmistandard import PICMI_PseudoRandomLayout
from picmistandard import PICMI_UniformDistribution
from picmistandard import PICMI_UniformFluxDistribution
from picmistandard import PICMI_Species, PICMI_MultiSpecies

cskw = {} # code-specific keywords
#kw = {}

class ConstantAppliedField(PICMI_ConstantAppliedField):
    # Redefine the `init` method, as required by the picmi `_ClassWithInit`
    def init(self, kw):
        self.V_xmin = cskw.pop('V_xmin', None)
        self.V_xmax = cskw.pop('V_xmax', None)
        self.V_ymin = cskw.pop('V_ymin', None)
        self.V_ymax = cskw.pop('V_ymax', None)

# Add injection distribution to Species class
class Species(PICMI_Species):
    def add_injection_distribution(self, distribution=None):
        self.injection_distribution = distribution

# Custom write input
class Simulation(PICMI_Simulation):
    def write_input_file(self, file_name='input_file.inp', description=''):

        assert len(self.applied_fields)<2 # only one append is allowed

        content = "HEXAPIC \n{\t "+description+" \n}\n"
        
        content += "Region \n{\n"
        
        massmin = array([self.species[i].mass for i in range(len(self.species))])
        massmin = np_min(massmin)

        for specie in self.species:
            content += "Species\n{\n"
            content += "\t name = "+ specie.name + "\n"
            content += "\t m = "+ str(specie.mass) + "\n"
            content += "\t q = "+ str(specie.charge) + "\n"
            for i in range(len(self.SRC)):
                if specie.name in self.SRC[i]['sp_source_f'].keys():
                    f = self.SRC[i]['sp_source_f'][specie.name]
                    content += "\t source" + str(i)+ " = "+ str(f) + "\n"
                if specie.name in self.SRC[i]['sp_heat_f'].keys():
                    f = self.SRC[i]['sp_heat_f'][specie.name]
                    content += "\t heat" + str(i)+ " = "+ str(f) + "\n"
            subcycle = max(1, int(sqrt(specie.mass/massmin)/4 + 0.5))
            content += "\t subcycle = "+ str(subcycle) + "\n}\n"
        
        grid = self.solver.grid
        content += "Grid\n{\n"
        content += "\t J = "+ str(grid.number_of_cells[0]) + "\n"
        content += "\t x1s = "+ str(grid.lower_bound[0]) + "\n"
        content += "\t x1f = "+ str(grid.upper_bound[0]) + "\n"
        content += "\t n1 = 1.0 \n"
        content += "\t K = "+ str(grid.number_of_cells[1]) + "\n"
        content += "\t x2s = "+ str(grid.lower_bound[1]) + "\n"
        content += "\t x2f = "+ str(grid.upper_bound[1]) + "\n"
        content += "\t n2 = 1.0 \n"
        content += "\t Geometry = 1\n"

        PeriodicFlagX1, PeriodicFlagX2 = False, False
        BoundarycX = [grid.lower_boundary_conditions[0],
                      grid.upper_boundary_conditions[0]]
        BoundarycY = [grid.lower_boundary_conditions[1],
                      grid.upper_boundary_conditions[1]]

        if 'periodic' in BoundarycX:
            assert not('dirichlet' in BoundarycX), \
            Exception('For PBC bc_xmin and bc_xmax should be \'periodic\'')
            grid.lower_boundary_conditions[0] = 'periodic'
            grid.upper_boundary_conditions[0] = 'periodic'
            if self.applied_fields:
                field = self.applied_fields[0]
                assert field.Ex == None, \
                Exception('For PBC in x direction Ex should be \'None\'')
            content += "\t PeriodicFlagX1 = 1 \n"
            PeriodicFlagX1 = True
        if 'periodic' in BoundarycY:
            assert not('dirichlet' in BoundarycY), \
            Exception('For PBC bc_ymin and bc_ymax should be \'periodic\'')
            grid.lower_boundary_conditions[1] = 'periodic'
            grid.upper_boundary_conditions[1] = 'periodic'
            if self.applied_fields:
                field = self.applied_fields[0]
                assert field.Ey == None, \
                Exception('For PBC in y direction Ey should be \'None\'')
            content += "\t PeriodicFlagX2 = 1 \n"
            PeriodicFlagX2 = True
        content += "}\n"

        content += "Control\n{\n"
        content += "\t dt = "+ str(self.time_step_size) + "\n"
        if not type(self.solver) == PICMI_ElectrostaticSolver:
            raise ValueError('HEXAPIC only supports electrostatic solver.')
        if self.applied_fields:
            field = self.applied_fields[0]
            if field.Bx !=None:
                content += "\t B01 = "+ str(field.Bx) + "\n"
            if field.By !=None:
                content += "\t B02 = "+ str(field.By) + "\n"
            if field.Bz !=None:
                content += "\t B03 = "+ str(field.Bz) + "\n"
        if ('periodic' in BoundarycX) or ('periodic' in BoundarycY):
            content += "\t ElectrostaticFlag = 5 \n}\n"
            # Periodic DADI field solver for XOOPIC
        else:
            content += "\t ElectrostaticFlag = 1 \n}\n"
        
        if self.applied_fields:
            field = self.applied_fields[0]
            # there is nothing on the periodic boundary
            if not(PeriodicFlagX1):
                if BoundarycX[0] == 'dirichlet':
                    content += "Equipotential\n{\n"
                    content += "\t C = "+ str(field.V_xmin) + "\n"
                    content += "\t reflection = "
                    content += str(grid.reflection[0]) + "\n"
                elif BoundarycX[0] == 'neumann':
                    content += "Dielectric\n{\n"
                    content += "\t QuseFlag = 0" + "\n"
                    content += "\t reflection = "
                    content += str(grid.reflection[0]) + "\n"
                else:
                    # open poundary as dielectric with no reflection
                    content += "Dielectric\n{\n"
                    content += "\t QuseFlag = 0" + "\n"
                    content += "\t reflection = 0" + "\n"
                content += "\t j1 = 0 \n"
                content += "\t j2 = 0 \n"
                content += "\t k1 = 0 \n"
                content += "\t k2 = "+ str(grid.number_of_cells[1]) + "\n"
                content += "\t normal = 1 \n}\n"
                if BoundarycX[1] == 'dirichlet':
                    content += "Equipotential\n{\n"
                    content += "\t C = "+ str(field.V_xmax) + "\n"
                    content += "\t reflection = "
                    content += str(grid.reflection[1]) + "\n"
                elif BoundarycX[1] == 'neumann':
                    content += "Dielectric\n{\n"
                    content += "\t QuseFlag = 0" + "\n"
                    content += "\t reflection = "
                    content += str(grid.reflection[1]) + "\n"
                else:
                    content += "Dielectric\n{\n"
                    content += "\t QuseFlag = 0" + "\n"
                    content += "\t reflection = 0" + "\n"
                content += "\t j1 = "+ str(grid.number_of_cells[0]) + "\n"
                content += "\t j2 = "+ str(grid.number_of_cells[0]) + "\n"
                content += "\t k1 = 0 \n"
                content += "\t k2 = "+ str(grid.number_of_cells[1]) + "\n"
                content += "\t normal = -1 \n}\n"
            if not(PeriodicFlagX2):
                if BoundarycY[0] == 'dirichlet':
                    content += "Equipotential\n{\n"
                    content += "\t C = "+ str(field.V_ymin) + "\n"
                    content += "\t reflection = "
                    content += str(grid.reflection[2]) + "\n"
                elif BoundarycY[0] == 'neumann':
                    content += "Dielectric\n{\n"
                    content += "\t QuseFlag = 0" + "\n"
                    content += "\t reflection = "
                    content += str(grid.reflection[2]) + "\n"
                else:
                    content += "Dielectric\n{\n"
                    content += "\t QuseFlag = 0" + "\n"
                    content += "\t reflection = 0" + "\n"
                content += "\t j1 = 0 \n"
                content += "\t j2 = "+ str(grid.number_of_cells[0]) + "\n"
                content += "\t k1 = 0 \n"
                content += "\t k2 = 0 \n"
                content += "\t normal = 1 \n}\n"
                if BoundarycY[1] == 'dirichlet':
                    content += "Equipotential\n{\n"
                    content += "\t C = "+ str(field.V_ymax) + "\n"
                    content += "\t reflection = "
                    content += str(grid.reflection[3]) + "\n"
                elif BoundarycY[1] == 'neumann':
                    content += "Dielectric\n{\n"
                    content += "\t QuseFlag = 0" + "\n"
                    content += "\t reflection = "
                    content += str(grid.reflection[3]) + "\n"
                else:
                    content += "Dielectric\n{\n"
                    content += "\t QuseFlag = 0" + "\n"
                    content += "\t reflection = 0" + "\n"
                content += "\t j1 = 0 \n"
                content += "\t j2 = "+ str(grid.number_of_cells[0]) + "\n"
                content += "\t k1 = "+ str(grid.number_of_cells[1]) + "\n"
                content += "\t k2 = "+ str(grid.number_of_cells[1]) + "\n"
                content += "\t normal = -1 \n}\n"

        for i,specie in enumerate(self.species):
            idist = specie.initial_distribution
            if idist is not None:
                # np2c = int(idist.density*Vol/self.layouts[i].n_macroparticles)
                np2c = np2c_global
                content += "Load\n{\n"
                content += "\t speciesName = "+ specie.name + "\n"
                content += "\t x1MinMKS = "+ str(idist.lower_bound[0])+"\n"
                content += "\t x1MaxMKS = "+ str(idist.upper_bound[0])+"\n"
                content += "\t x2MinMKS = "+ str(idist.lower_bound[1])+"\n"
                content += "\t x2MaxMKS = "+ str(idist.upper_bound[1])+"\n"
                content += "\t density = "+ str(idist.density)+"\n"
                content += "\t np2c = "+ "{:.1e}".format(np2c)+"\n"
                content += "\t v1thermal = "+ str(idist.rms_velocity[0])+"\n"
                content += "\t v2thermal = "+ str(idist.rms_velocity[1])+"\n"
                content += "\t v3thermal = "+ str(idist.rms_velocity[2])+"\n"
                content += "\t v1drift = "+ str(idist.directed_velocity[0])+"\n"
                content += "\t v2drift = "+ str(idist.directed_velocity[1])+"\n"
                content += "\t v3drift = "+ str(idist.directed_velocity[2])+"\n"
                content += "\t LoadMethodFlag = 0\n}\n"

            inj = specie.injection_distribution
            if inj is not None:
                np2c = np2c_global
                content += "EmitPort\n{\n"
                content += "\t speciesName = "+ specie.name + "\n"
                current = float(inj.flux)*self.time_step_size
                if inj.flux_normal_axis == 'x':
                    current *= inj.upper_bound[1]-inj.lower_bound[1]
                elif inj.flux_normal_axis == 'y':
                    current *= inj.upper_bound[0]-inj.lower_bound[0]
                else:
                    raise ValueError('Injection plane normal not supported.')
                inj_dx = (grid.upper_bound[0]-grid.lower_bound[0])/ \
                    grid.number_of_cells[0]
                inj_dy = (grid.upper_bound[1]-grid.lower_bound[1])/ \
                    grid.number_of_cells[1]
                j1 = int(inj.lower_bound[0]/inj_dx)
                j2 = int(inj.upper_bound[0]/inj_dx)
                k1 = int(inj.lower_bound[1]/inj_dy)
                k2 = int(inj.upper_bound[1]/inj_dy)
                content += "\t I = "+ str(current)+"\n"
                content += "\t j1 = "+ str(j1)+"\n"
                content += "\t j2 = "+ str(j2)+"\n"
                content += "\t k1 = "+ str(k1)+"\n"
                content += "\t k2 = "+ str(k2)+"\n"
                content += "\t normal = "+ str(inj.flux_direction)+"\n"
                content += "\t np2c = "+ "{:.1e}".format(np2c)+"\n"
                content += "\t v1thermal = "+ str(inj.rms_velocity[0])+"\n"
                content += "\t v2thermal = "+ str(inj.rms_velocity[1])+"\n"
                content += "\t v3thermal = "+ str(inj.rms_velocity[2])+"\n"
                content += "\t v1drift = "+str(inj.directed_velocity[0])+"\n"
                content += "\t v2drift = "+str(inj.directed_velocity[1])+"\n"
                content += "\t v3drift = "+str(inj.directed_velocity[2])+"\n}\n"


        SRC = self.SRC
        for i in range(len(SRC)):
            content += "Source\n{\n"
            content += "\t index = "+ str(i) + "\n"
            content += "\t n_inj = "+ str(SRC[i]['n_inj']) + "\n"
            content += "\t shape_x = "+ str(SRC[i]['shape_x']) + "\n"
            content += "\t shape_y = "+ str(SRC[i]['shape_y']) + "\n"
            content += "\t j1 = "+ str(SRC[i]['j1']) + "\n"
            content += "\t j2 = "+ str(SRC[i]['j2']) + "\n"
            content += "\t k1 = "+ str(SRC[i]['k1']) + "\n"
            content += "\t k2 = "+ str(SRC[i]['k2']) + "\n"
            content += "\t Tpar = "+ str(SRC[i]['Tpar']) + "\n"
            content += "\t Tper = "+ str(SRC[i]['Tper']) + "\n"
            content += "\t heat_freq = "+ "{:.3e}".format(SRC[i]['heat_freq']) + "\n"
            content += "\t t_dep = "+ str(SRC[i]['t_dep']) + "\n"
            if SRC[i]['t_dep'] == 1:
                content += "\t t0 = "+ str(SRC[i]['t0']) + "\n"
                content += "\t t1 = "+ str(SRC[i]['t1']) + "\n"
                content += "\t t2 = "+ str(SRC[i]['t2']) + "\n"
                content += "\t t3 = "+ str(SRC[i]['t3']) + "\n}\n"
        

        for i in range(len(self.MCC)):
            content += "MCC\n{\n"
            content += "\t collName = "+ self.MCC[i][0] + "\n"
            content += "\t specie1Name = "+ self.MCC[i][1] + "\n"
            content += "\t specie2Name = "+ self.MCC[i][2] + "\n"
            content += "\t reactionIdx = "+ self.MCC[i][3] + "\n"
            content += "\t collXSectionFile = "+ self.MCC[i][4] + "\n}\n"


        PSI = self.PSI
        for i in range(len(PSI)):
            psi = PSI[i]
            for j in range(len(psi["boundary"])):
                for k in range(len(psi["secondary1"])):
                    content += "PSI\n{\n"
                    content += "\t boundary = " + psi["boundary"][j] + "\n"
                    content += "\t specie = " + psi["specie"] + "\n"
                    content += "\t secondary1 = " + psi["secondary1"][k] + "\n"
                    if len(psi["secondary2"][k]) > 0:
                        content += "\t secondary2 = " + psi["secondary2"][k] + "\n"
                    if len(psi["secondary3"][k]) > 0:
                        content += "\t secondary3 = " + psi["secondary3"][k] + "\n"
                    content += "\t sec_type = " + str(psi["sec_types"][k]) + "\n"
                    content += "\t inj_type = " + str(psi["inj_types"][k]) + "\n"
                    content += "\t Ew = " + str(psi["Ew"][k]) + "\n"
                    content += "\t E0 = " + str(psi["E0"][k]) + "\n"
                    content += "\t gamma0 = " + str(psi["gamma0"][k]) + "\n"
                    content += "\t T0 = " + str(psi["T0"][k]) + "\n"
                    content += "\t Ep = " + str(psi["Ep"][k]) + "\n}\n"


        content += "}\n" # End of Region

        f = open(file_name, "w")
        f.write(content)
        f.close()

# ---- test case----
NX = 100 # number of cells along X
NY = 100 # number of cells along Y
dx = 1e-5 # cell size in m 
dy = 1e-5 # cell size in m
Lx = NX*dx # the spatial extent of the simulation along X
Ly = NY*dy # the spatial extent of the simulation along Y
Vol = Lx*Ly*1
bc_xmin = 'open'
bc_xmax = 'open'
bc_ymin = 'dirichlet'
bc_ymax = 'dirichlet'
reflection = [0.0, 0.0, 0.0, 0.0] # boundary reflections
#xmin, xmax, ymin, ymax
cskw['V_ymin'] =  -10 # boundary potential in V
cskw['V_ymax'] =  10 # boundary potential in V
cskw['V_xmax'] =  0 # boundary potential in V
cskw['V_xmin'] =  0
# Warning: if Ex or Ey = 0.0, it is assumed that the corresponding boundaries 
# are grounded. Use None for open boundaries.
Ex = None # external electric field in V/m
Ey = (cskw['V_ymax'] - cskw['V_ymin']) / Ly # external electric field in V/m
Ez = None # external electric field in V/m
Bx = 0.1 # external magnetic field in T 
By = 0.0 # external magnetic field in T 
Bz = None # external magnetic field in T 
dt = 1.0e-12 # time step in s
max_steps = 1000 # maximum number of simulation time steps
solver_method = 'Multigrid'
solver_max_iter = 1000
particle_names = ['electron', 'H+', 'H']
specie_density = [1.0e18, 1.0e18, 0.0] # particles/m^3
specie_temperature = [5.0, 5.0, 1.0] # eV
np2c_global = 1e6
np2c = []
specie_mp = []
for i,_ in enumerate(particle_names):
    np2c.append(np2c_global)
    specie_mp.append(int(specie_density[i]*Vol/np2c[i]))
    assert (specie_mp[i]>1 or specie_mp[i]==0)
injection_plane_axis = 'x'
injection_plane_position = 0.0
injection_plane_lower_bound = [0.0, Ly/2-Ly/4, None]
injection_plane_upper_bound = [0.0, Ly/2+Ly/4, None]
injection_plane_normal = 1
particle_currents = [5e5/dt, 5e5/dt, 0.0/dt] # particles/(m^2*s)


# MCC list of collision types
# [[name, specie1, specie2, reaction_index, cross-section_file], ...]
# MCC = [['Null', 'electron', 'positron', '1', 'e_p_test.txt']]
MCC = [['Elastic', 'electron', 'H', '1', 'eH_el.txt'],
       ['Exciation', 'electron', 'H', '2', 'eH_exc.txt'],
       ['Ionisation', 'electron', 'H', '3', 'eH_ion.txt'],
       ['Charge_excange', 'H+', 'H', '5', 'HH+_cx.txt'],
       ['Elastic', 'H+', 'H', '1', 'HH+_el.txt'],
       ['Elastic', 'H', 'H', '1', 'HH_el.txt']]


# PSI for each of the boundaries xmin, xmax, ymin, ymax
psi_e = {'boundary' : ['xmin', 'xmax', 'ymin', 'ymax'],
         'specie' : 'electron',
         'secondary1' : ['electron'], # possible secondaries1
         'secondary2' : [''], # possible secondaries2 if '' no additional sp
         'secondary3' : [''], # possible secondaries3 if '' no additional sp
         'sec_types' : [0], # types of secondary injections
         'inj_types' : [0], # types of velocity distributions for injections
         # parameters for injection
         # effect of different coefficients explained inside PSI routines
         # it is not necessary that E0 is allways energy for non 0 sec_type
         # something else of similar importance might be saved here
         # what is actualy saved is explained in PSI routines
         'Ew' : [4.0], # [eV]
         'E0' : [4.0], # [eV]
         'gamma0' : [1.0], 
         'T0' : [4.0], # [eV]
         'Ep' : [5.0]  # [eV]
        }

psi_i = {'boundary' : ['xmin', 'xmax', 'ymin', 'ymax'],
         'specie' : 'H+',
         'secondary1' : ['electron', 'H'], # possible secondaries1
         'secondary2' : ['', ''], # possible secondaries2 if '' no additional sp
         'secondary3' : ['', ''], # possible secondaries3 if '' no additional sp
         'sec_types' : [1, 3], # types of secondary injections
         'inj_types' : [0, 1], # types of velocity distributions for injections
         'Ew' : [4.0, 4.0], # [eV]
         'E0' : [1e4, 4.0], # [m/s] for sec type 1
         'gamma0' : [1.0, 1.0], 
         'T0' : [4.0, 4.0], # [eV]
         'Ep' : [5.0, 5.0]  # [eV]
        }

PSI = [psi_e, psi_i] # put in all posible species reacting with wall
# for each boundary you can do its own, then specify just one (['xmin'])


source0 = {'n_inj' : 0.0, # number of injected particles [s^-1 m^-3]
           'shape_x' : 0,  # uniform/cos/normal
           'shape_y' : 1,
           'j1' : 30,     # bottom-left corner
           'k1' : 30,
           'j2' : 50,     # top-right corner
           'k2' : 50,
           'Tpar' : 7.0,   # source paralel temperature [eV]
           'Tper' : 6.0,   # source perpendicular temperature [eV]
           'heat_freq' : 0.0, # fraquency of source heating per particle
           't_dep' : 1,
           't0' : 0.0,
           't1' : 1e-11,
           't2' : 3e-11,
           't3' : 4e-11,
           'sp_source_f' : {'electron' : 0.5, 'H+' : 0.5},
           'sp_heat_f' : {'electron' : 0.5, 'H+' : 0.5}
          }
source1 = {'n_inj' : 0.0, # number of injected particles [s^-1 m^-3]
           'shape_x' : 2,  # uniform/cos/normal
           'shape_y' : 2,
           'j1' : 20,     # bottom-left corner
           'k1' : 20,
           'j2' : 50,     # top-right corner
           'k2' : 60,
           'Tpar' : 5.0,   # source paralel temperature [eV]
           'Tper' : 5.0,   # source perpendicular temperature [eV]
           'heat_freq' : 0.0, # fraquency of source heating per particle
           't_dep' : 0,
           'sp_source_f' : {'H' : 1.0},
           'sp_heat_f' : {'H' : 1.0}
          }

SRC = [source0, source1]


case_description =  "Uniformly distributed electron-proton plasma "
case_description += "between two biased plates with wall injection. "
case_description += "MCC + PSI."
file_name='input_file.inp'
# ------------------  

velth = {} # thermal velocity in any direction
veld = {}  # drift velocity in any direction
for i,specie in enumerate(particle_names):
    velth[specie] = [0.0, 0.0, 0.0]
    for j in range(3):
        velth[specie][j] = \
            sqrt(specie_temperature[i]*constants.q_e/particle_mass[specie])
    veld[specie] = [0.0, 0.0, 0.0]


grid = PICMI_Cartesian2DGrid(number_of_cells=None, 
    lower_bound=None, upper_bound=None,
    lower_boundary_conditions=[bc_xmin,bc_ymin], 
    upper_boundary_conditions=[bc_xmax,bc_ymax], nx=NX, ny=NY, xmin=0.0,
    xmax=Lx, ymin=0.0, ymax=Ly, bc_xmin=None, bc_xmax=None, 
    bc_ymin=None, bc_ymax=None, moving_window_velocity=None, 
    refined_regions=[], lower_bound_particles=None, upper_bound_particles=None,
    xmin_particles=None, xmax_particles=None, 
    ymin_particles=None, ymax_particles=None, 
    lower_boundary_conditions_particles=None, 
    upper_boundary_conditions_particles=None, 
    bc_xmin_particles=None, bc_xmax_particles=None, 
    bc_ymin_particles=None, bc_ymax_particles=None, 
    guard_cells=None, pml_cells=None)

grid.reflection = reflection

solver = PICMI_ElectrostaticSolver(grid, 
    method=solver_method, 
    required_precision=None, 
    maximum_iterations=solver_max_iter)

field = ConstantAppliedField(Ex=Ex, Ey=Ey, Ez=Ez, Bx=Bx, By=By, Bz=Bz, 
    lower_bound=[grid.lower_bound[0], grid.lower_bound[1], None], 
    upper_bound=[grid.upper_bound[0], grid.upper_bound[1], None])

sim = Simulation(solver=solver, time_step_size=dt, max_steps=max_steps, 
    max_time=None, verbose=None, particle_shape='linear', gamma_boost=None, 
    load_balancing=None)

sim.add_applied_field(field)

for i,p in enumerate(particle_names):
    init_dist = PICMI_UniformDistribution(specie_density[i], 
        lower_bound=[grid.lower_bound[0], grid.lower_bound[1], None], 
        upper_bound=[grid.upper_bound[0], grid.upper_bound[1], None], 
        rms_velocity=velth[p], directed_velocity=veld[p], fill_in=None)
    particle_layout = PICMI_PseudoRandomLayout(n_macroparticles=specie_mp[i], 
        n_macroparticles_per_cell=None, seed=None, grid=None)
    specie = Species(particle_type=None, name=p, 
        charge_state=None, charge=particle_charge[p], mass=particle_mass[p], 
        initial_distribution=init_dist, 
        particle_shape=None, density_scale=None, method='Boris')
    inject_dist = PICMI_UniformFluxDistribution(flux=particle_currents[i], 
        flux_normal_axis=injection_plane_axis,
        surface_flux_position=injection_plane_position, 
        flux_direction=injection_plane_normal,
        lower_bound=injection_plane_lower_bound, 
        upper_bound=injection_plane_upper_bound,
        rms_velocity=velth[p], directed_velocity=veld[p], 
        flux_tmin = None, flux_tmax = None, 
        gaussian_flux_momentum_distribution = True)
    specie.add_injection_distribution(inject_dist)

    sim.add_species(species=specie, layout=particle_layout,
        initialize_self_field=None)

sim.MCC = MCC
sim.PSI = PSI
sim.SRC = SRC
    
sim.write_input_file(file_name, case_description)
