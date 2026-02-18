/**
 * SPDX-License-Identifier: EUPL-1.2
 * SPDX-FileCopyrightText: 2025 Stefan Costea, Miha Radez, LeCAD-PEG
 *
 * @file collisions.cpp
 * @brief Monte-Carlo Collisions (MCC): elastic, ionization, excitation, and
 * charge-exchange processes.
 */

/* HEXAPIC particle collisions section */

// Coulomb collisions should be separated from other collisions
// not yet implemented...

#include "hexapic.hpp"

void collisions_init();
void read_add_collision(int i);
void test_sp_reaction_compatibility(int sp1, int sp2, int reaction_index);
void add_new_empty_collision(int sp1, int sp2);
reactionFunction get_reaction(int reaction_index);
void add_sp_out(int coll_i, int reaction_index);
void read_add_cross_sections(int coll_i, std::string filename);
void set_maximal_sigma_v(int coll_i);
REAL read_cross_section(int coll_i, int j, REAL E);

void MCC();
void collision_in_cell(int coll_i, int cell_i);
REAL relative_velocity(int cell_i, int sp1, int sp2, int p1, int p2);
REAL reduced_mass(REAL m1, REAL m2);
REAL density_computed_particles(int cell_i, int sp);
void COM_velocity(REAL v_com[3], int cell_i, int sp_1, int sp_2, int p_1,
                  int p_2);
void isotropic_scatter(REAL v_scatt[3]);
void isotropic_scatter(REAL v_scatt[3], REAL v);
void null_collision(int coll_i, int react_i, int cell_i, int sp_1, int sp_2,
                    int p_1, int p_2, REAL E);
void elastic_scattering_collision(int coll_i, int react_i, int cell_i, int sp_1,
                                  int sp_2, int p_1, int p_2, REAL E);
void excitation_collision(int coll_i, int react_i, int cell_i, int sp_1,
                          int sp_2, int p_1, int p_2, REAL E);
void ionisation_collision(int coll_i, int react_i, int cell_i, int sp_1,
                          int sp_2, int p_1, int p_2, REAL E);
void charge_exchange_collision(int coll_i, int react_i, int cell_i, int sp_1,
                               int sp_2, int p_1, int p_2, REAL E);

/*********************************/
/*    initialization functions   */
/*********************************/

void collisions_init() {
  collisions = {};
  for (int i = 0; i < MCC_init_data.size(); i++)
    read_add_collision(i);
  for (int i = 0; i < collisions.size(); i++)
    set_maximal_sigma_v(i);

  // Clear the data and release allocated memory
  for (auto &inner : MCC_init_data) {
    inner.clear();
    inner.shrink_to_fit();
  }
  MCC_init_data.clear();
  MCC_init_data.shrink_to_fit();
}

void read_add_collision(int i) {
  // reads cross-section file and adds collision (or just a new reaction)
  std::vector<std::string> MCC_init = MCC_init_data[i];
  std::string name, specie1, specie2, filename;
  int reaction_index, sp1, sp2;

  name = MCC_init[0];
  specie1 = MCC_init[1];
  specie2 = MCC_init[2];
  reaction_index = std::stoi(MCC_init[3]);
  filename = MCC_init[4];

  bool new_collision = true;
  int coll_i;

  sp1 = specie_index(specie1);
  sp2 = specie_index(specie2);
  if (species[sp1].mass > species[sp2].mass)
    std::swap(sp1, sp2);
  test_sp_reaction_compatibility(sp1, sp2, reaction_index);

  for (coll_i = 0; coll_i < collisions.size(); coll_i++) {
    if ((collisions[coll_i].sp[0] == sp1 && collisions[coll_i].sp[1] == sp2) ||
        (collisions[coll_i].sp[1] == sp1 && collisions[coll_i].sp[0] == sp2)) {
      new_collision = false;
      break;
    }
  }

  if (new_collision)
    add_new_empty_collision(sp1, sp2);
  Collision &coll = collisions[coll_i];
  coll.reaction_names.push_back(name);
  coll.n_reactions += 1;
  coll.processCollision.push_back(get_reaction(reaction_index));
  add_sp_out(coll_i, reaction_index);
  read_add_cross_sections(coll_i, filename);
}

int specie_index(std::string specie_name) {
  for (int i = 0; i < species.size(); i++) {
    if (species[i].name == specie_name)
      return i;
  }
  std::cout << "Specie name was not found in used species vector!\n";
  return species.size();
}

void test_sp_reaction_compatibility(int sp1, int sp2, int reaction_index) {
  // TO DO: add errors, stop run
  bool error = false;
  // ...

  if (error)
    std::cout << "Species " << species[sp1].name << ", " << species[sp2].name
              << " not compatible with reaction " << reaction_index << "\n";
}

void add_new_empty_collision(int sp1, int sp2) {
  Collision new_collision;
  new_collision.reaction_names = {};
  new_collision.sp[0] = sp1;
  new_collision.sp[1] = sp2;
  new_collision.n_reactions = 0;
  new_collision.processCollision = {};
  new_collision.sp_out = {};
  new_collision.sigma = {};
  new_collision.E_min = {};
  collisions.push_back(new_collision);
}

reactionFunction get_reaction(int reaction_index) {
  // return pointer to proper reaction function between two particles
  // void f(CellID, sp1ID, sp2ID, part1ID, part2ID)
  if (reaction_index == 0)
    return &null_collision;
  else if (reaction_index == 1)
    return &elastic_scattering_collision;
  else if (reaction_index == 2)
    return &excitation_collision;
  else if (reaction_index == 3)
    return &ionisation_collision;
  else if (reaction_index == 5)
    return &charge_exchange_collision;
  else {
    std::cout << "Reaction not jet specified!\n";
    return &null_collision;
  }
}

void add_sp_out(int coll_i, int reaction_index) {
  // index of species that come out of reaction
  int sp1 = collisions[coll_i].sp[0], sp2 = collisions[coll_i].sp[1];

  if (reaction_index == 0 || reaction_index == 1 || reaction_index == 2)
    sp1 = sp1;

  else if (reaction_index == 3) { // ionisation
    // find specie with name(sp2) + "+"
    sp2 = specie_index(species[sp2].name + "+");
    if (sp2 == collisions[coll_i].sp[1] && sp2 == species.size()) {
      std::cout << "ERROR: Proper ion specie of " << species[sp2].name
                << " not specified!\n";
      sp2 = collisions[coll_i].sp[1];
    }
  }

  else if (reaction_index == 5) { // charge exchange
    std::string name1 = species[sp1].name, name2 = species[sp1].name;
    if (name1.back() != '+')
      sp1 = specie_index(name1 + "+");
    else {
      name1.pop_back();
      sp1 = specie_index(name1);
    }

    if (name2.back() != '+')
      sp2 = specie_index(name2 + "+");
    else {
      name2.pop_back();
      sp2 = specie_index(name2);
    }

  }

  else {
    std::cout << "Sp_out not specified jet!\n";
  }

  collisions[coll_i].sp_out.push_back({sp1, sp2});
}

void read_add_cross_sections(int coll_i, std::string filename) {
  // reads E and sigma from file

  int i, N = 0;
  REAL E_min;
  std::string line;
  std::set<std::pair<REAL, REAL>> E_sigma;
  std::ifstream file(filename);

  if (!file) {
    std::cerr << "Unable to open file: " << filename << std::endl;
    collisions[coll_i].sigma.push_back(E_sigma);
  }

  // Read the first line to extract N
  if (getline(file, line)) {
    std::istringstream iss(line);

    if (!(iss >> N >> i >> E_min)) { // Read N and E_min
      std::cerr << "Failed to read the number of data pairs.\n";
      collisions[coll_i].sigma.push_back(E_sigma);
    }
    collisions[coll_i].E_min.push_back(E_min);
  }

  for (i = 0; i < N && getline(file, line); i++) {
    std::istringstream iss(line);
    double E, sigma;
    if (iss >> E >> sigma) {
      E_sigma.insert({E, sigma});
    }
  }

  file.close();
  collisions[coll_i].sigma.push_back(E_sigma);
}

void set_maximal_sigma_v(int coll_i) {
  // check for maximum of total cross-section times v

  REAL sigma_t, sigma_v_max;
  Collision &coll = collisions[coll_i];
  std::set<REAL> E_all = {};
  REAL m = reduced_mass(species[coll.sp[0]].mass, species[coll.sp[1]].mass);

  // all tabulated energies
  for (int j = 0; j < coll.n_reactions; j++)
    for (auto E_sigma : coll.sigma[j])
      E_all.insert(E_sigma.first);

  sigma_v_max = 0.0;
  for (auto E : E_all) {
    // total cross-section
    sigma_t = 0.0;
    for (int j = 0; j < coll.n_reactions; j++)
      sigma_t += read_cross_section(coll_i, j, E);

    sigma_v_max = std::max(sigma_v_max, sigma_t * sqrt(2 * E * q_e / m));
  }
  coll.sigma_v_max = sigma_v_max;
}

REAL read_cross_section(int coll_i, int j, REAL E) {
  // collision_i, j -> which collision reaction
  // linear interpolation between tabulated values
  // outside tabulated domain -> sigma = 0
  REAL x;
  std::set<std::pair<REAL, REAL>> &sigma = collisions[coll_i].sigma[j];

  auto it_next = sigma.upper_bound({E, 0.0});
  if (it_next == sigma.end())
    return 0.0;
  else if (it_next == sigma.begin())
    return 0.0;
  else {
    auto it_prev = std::prev(it_next);
    x = (it_next->first - E) / (it_next->first - it_prev->first);
    return x * it_prev->second + (1 - x) * it_next->second;
  }
}

/*********************************/
/*      colliding functions      */
/*********************************/

void MCC() {
  // For now done for equal particle weights (np2c_global)
  for (int coll_i = 0; coll_i < collisions.size(); coll_i++) {
    for (int cell_i = 0; cell_i < grid.Ncx * grid.Ncy; cell_i++)
      collision_in_cell(coll_i, cell_i);
  }
}

void collision_in_cell(int coll_i, int cell_i) {
  REAL P_tot, nu_max, nu_cum, sp2_density, R, E, v, m;
  int N_coll, i, j, sp1, sp2, p1, p2;

  Collision &coll = collisions[coll_i];
  sp1 = coll.sp[0];
  sp2 = coll.sp[1];
  Cell &cell = cells[cell_i];

  // target density
  sp2_density = density_computed_particles(cell_i, sp2) * np2c_global;
  nu_max = sp2_density * coll.sigma_v_max;
  P_tot = 1 - exp(-nu_max * dt); ///< maximal probabillity for collision

  // N_coll_actual = Ptot * (#computed_particles * np2c_global)
  // N_coll_computed = N_coll_actual / np2c_global
  // adding frand() not to underestimate number of collisions each time
  N_coll = cell.part[sp1].size() * P_tot + frand();

  for (i = 0; i < N_coll; i++) {
    // Choose random particles
    p1 = cell.part[sp1].size() * frand();
    p2 = cell.part[sp2].size() * frand();

    v = relative_velocity(cell_i, sp1, sp2, p1, p2);
    m = reduced_mass(species[sp1].mass, species[sp2].mass);
    E = m * v * v / (2 * q_e); ///< Energy for cross-section [eV]

    nu_cum = 0.0; ///< cumulative likelihood between 0 and nu_max;
    R = nu_max * frand();
    for (j = 0; j < coll.n_reactions; j++) {
      // null collision event if nu_cum remains under R
      nu_cum += sp2_density * v * read_cross_section(coll_i, j, E);
      if (nu_cum > R) { // collision j actually takes place
        // coll_i and j required for additional data and sp_out
        coll.processCollision[j](coll_i, j, cell_i, sp1, sp2, p1, p2, E);
        break;
      }
    }
  }
}

REAL density_computed_particles(int cell_i, int sp) {
  // density of computed paticles = density / np2c
  REAL dV = grid.dx * grid.dy * 1;
  return cells[cell_i].part[sp].size() / dV;
}

REAL relative_velocity(int cell_i, int sp1, int sp2, int p1, int p2) {
  Particle &part1 = cells[cell_i].part[sp1][p1];
  Particle &part2 = cells[cell_i].part[sp2][p2];
  REAL v = (part1.vx - part2.vx) * (part1.vx - part2.vx);
  v += (part1.vy - part2.vy) * (part1.vy - part2.vy);
  v += (part1.vz - part2.vz) * (part1.vz - part2.vz);
  return sqrt(v);
}

REAL reduced_mass(REAL m1, REAL m2) { return m1 * m2 / (m1 + m2); }

void COM_velocity(REAL v_com[3], int cell_i, int sp_1, int sp_2, int p_1,
                  int p_2) {
  REAL m1 = species[sp_1].mass;
  REAL m2 = species[sp_2].mass;
  REAL M = m1 + m2;
  Particle &p1 = cells[cell_i].part[sp_1][p_1];
  Particle &p2 = cells[cell_i].part[sp_2][p_2];

  v_com[0] = (p1.vx * m1 + p2.vx * m2) / M;
  v_com[1] = (p1.vy * m1 + p2.vy * m2) / M;
  v_com[2] = (p1.vz * m1 + p2.vz * m2) / M;
}

void isotropic_scatter(REAL v_scatt[3]) {
  // function for isotriopic scattering in COM
  // takes |v_scatt| as input
  REAL cos_th, sin_th, phi, v;

  cos_th = 1 - 2 * frand();
  sin_th = sqrt(1 - cos_th * cos_th);
  phi = TWOPI * frand();
  v = sqrt(v_scatt[0] * v_scatt[0] + v_scatt[1] * v_scatt[1] +
           v_scatt[2] * v_scatt[2]);

  v_scatt[2] = v * cos_th;
  v_scatt[0] = v * sin_th * cos(phi);
  v_scatt[1] = v * sin_th * sin(phi);
}

void isotropic_scatter(REAL v_scatt[3], REAL v) {
  // function for isotriopic scattering in COM
  // takes v = |v_scatt| as input
  REAL cos_th, sin_th, phi;

  cos_th = 1 - 2 * frand();
  sin_th = sqrt(1 - cos_th * cos_th);
  phi = TWOPI * frand();

  v_scatt[2] = v * cos_th;
  v_scatt[0] = v * sin_th * cos(phi);
  v_scatt[1] = v * sin_th * sin(phi);
}

// Different reaction functions
/*
        V. Vahedi, M. Surendra,
        A Monte Carlo collision model for the particle-in-cell method:
        applications to argon and oxygen discharges
*/

void null_collision(int coll_i, int react_i, int cell_i, int sp_1, int sp_2,
                    int p_1, int p_2, REAL E) {
  // reaction_index = 0
  // nothing happens
}

void elastic_scattering_collision(int coll_i, int react_i, int cell_i, int sp_1,
                                  int sp_2, int p_1, int p_2, REAL E) {
  // reaction_index = 1
  // e + A → e + A or (B+) + A → (B+) + A
  // isotropic scattering

  REAL mass_r;
  REAL v_com[3], v_scatt[3]; ///< center of mass frame, scattered v
  Particle &p1 = cells[cell_i].part[sp_1][p_1];
  Particle &p2 = cells[cell_i].part[sp_2][p_2];

  COM_velocity(v_com, cell_i, sp_1, sp_2, p_1, p_2);
  v_scatt[0] = p1.vx - v_com[0];
  v_scatt[1] = p1.vy - v_com[1];
  v_scatt[2] = p1.vz - v_com[2];
  isotropic_scatter(v_scatt);
  mass_r = species[sp_1].mass / species[sp_2].mass;

  // new velocities
  p1.vx = v_com[0] + v_scatt[0];
  p1.vy = v_com[1] + v_scatt[1];
  p1.vz = v_com[2] + v_scatt[2];
  // m1v1 + mvv2 in COM remains equal to 0
  p2.vx = v_com[0] - v_scatt[0] * mass_r;
  p2.vy = v_com[1] - v_scatt[1] * mass_r;
  p2.vz = v_com[2] - v_scatt[2] * mass_r;
}

void excitation_collision(int coll_i, int react_i, int cell_i, int sp_1,
                          int sp_2, int p_1, int p_2, REAL E) {
  // reaction_index = 2
  // e + A → e + A*
  // isotropic scattering

  /*
  for now deexcitation happens momentarily A* → A + γ so it works just as
  an energy sink (as we do not track photons)
  there are no excitated species in simulation
  there can be different excitations for each atom (different γ energies)
  */

  REAL m1, m2, mass_r, v;
  REAL v_com[3], v_scatt[3]; ///< center of mass frame, scatter velocity
  Particle &p1 = cells[cell_i].part[sp_1][p_1];
  Particle &p2 = cells[cell_i].part[sp_2][p_2];

  m1 = species[sp_1].mass;
  m2 = species[sp_2].mass;
  mass_r = m1 / m2;

  COM_velocity(v_com, cell_i, sp_1, sp_2, p_1, p_2);
  // in COM reduce energy by excitational energy
  v = (E - collisions[coll_i].E_min[react_i]) * 2 * q_e; ///< [J]
  v /= reduced_mass(m1, m2);
  v = sqrt(v); ///< this is relative velocity
  v = m2 * v / (m2 + m1);
  isotropic_scatter(v_scatt, v);

  // new velocities
  p1.vx = v_com[0] + v_scatt[0];
  p1.vy = v_com[1] + v_scatt[1];
  p1.vz = v_com[2] + v_scatt[2];
  p2.vx = v_com[0] - v_scatt[0] * mass_r;
  p2.vy = v_com[1] - v_scatt[1] * mass_r;
  p2.vz = v_com[2] - v_scatt[2] * mass_r;

  // particle 2 emits photon in random direction
  // change of velocity in its own COM
  // not relativistic!
  // very small change so commented out:
  /*
  v = collisions[coll_i].E_min[react_i] * q_e /
  species[sp_2].mass / 2.99792458e8;
  cos_th = 1 - 2*frand();
  sin_th = sqrt(1 - cos_th*cos_th);
  phi = TWOPI*frand();
  p2.vx += v*sin_th*cos(phi);
  p2.vy += v*sin_th*sin(phi);
  p2.vz += v*cos_th;
  */
}

void ionisation_collision(int coll_i, int react_i, int cell_i, int sp_1,
                          int sp_2, int p_1, int p_2, REAL E) {
  // reaction_index = 3
  // e + A → e + (A+) + e

  // simple model with energy in COM equaly distributed between electrons
  // because m_e << m_A we assume that ion is stationary in COM at the end

  int sp_3;
  REAL m1, m2, v;
  REAL v_com[3], v_scatt[3]; ///< center of mass frame, scatter velocity
  Particle &p1 = cells[cell_i].part[sp_1][p_1];
  Particle &p2 = cells[cell_i].part[sp_2][p_2];
  Particle p3 = cells[cell_i].part[sp_2][p_2]; ///< position copy for new e/i

  m1 = species[sp_1].mass;
  m2 = species[sp_2].mass;

  COM_velocity(v_com, cell_i, sp_1, sp_2, p_1, p_2);
  // in COM reduce energy by ionization energy
  v = (E - collisions[coll_i].E_min[react_i]) * q_e; ///< E(2e)[J]
  v /= m1;
  v = sqrt(v);
  isotropic_scatter(v_scatt, v);

  // scatter first electron
  p1.vx = v_com[0] + v_scatt[0];
  p1.vy = v_com[1] + v_scatt[1];
  p1.vz = v_com[2] + v_scatt[2];

  // create and scatter new electron
  p3.vx = v_com[0] - v_scatt[0];
  p3.vy = v_com[1] - v_scatt[1];
  p3.vz = v_com[2] - v_scatt[2];
  cells[cell_i].part[sp_1].push_back(p3);

  // add new ion and destroy neutral
  p3.vx = v_com[0];
  p3.vy = v_com[1];
  p3.vz = v_com[2];
  sp_3 = collisions[coll_i].sp_out[react_i][1];
  cells[cell_i].part[sp_3].push_back(p3);
  species[sp_3].n += 1;
  remove_particle(cell_i, sp_2, p_1);
}

void recombination_collision(int coll_i, int react_i, int cell_i, int sp_1,
                             int sp_2, int p_1, int p_2, REAL E) {
  // reaction_index = 4
  // e + (A+) → A + γ
  std::cout << "Collision happened (not yet implemented)!\n";
}

void charge_exchange_collision(int coll_i, int react_i, int cell_i, int sp_1,
                               int sp_2, int p_1, int p_2, REAL E) {
  // reaction_index = 5
  // A + (B+) → (A+) + B
  Particle p1 = cells[cell_i].part[sp_1][p_1];
  Particle p2 = cells[cell_i].part[sp_2][p_2];
  std::vector<int> &sp_out = collisions[coll_i].sp_out[react_i];

  // change particle specie
  cells[cell_i].part[sp_out[0]].push_back(p1);
  cells[cell_i].part[sp_out[1]].push_back(p2);
  species[sp_out[0]].n += 1;
  species[sp_out[1]].n += 1;
  remove_particle(cell_i, sp_1, p_1);
  remove_particle(cell_i, sp_2, p_2);
}

/*
For anisotropic scattering (maybe) use eq.(3) set with parameter alpha from
D.Tskhakaya/et al., Contrib. Plasma Phys. 48, No. 1-3, 121 . 125 (2008)
*/
