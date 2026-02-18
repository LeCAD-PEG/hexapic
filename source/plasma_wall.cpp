/**
 * SPDX-License-Identifier: EUPL-1.2
 * SPDX-FileCopyrightText: 2025 Stefan Costea, Miha Radez, LeCAD-PEG
 *
 * @file plasma_wall.cpp
 * @brief Boundary conditions and plasma–wall interaction (sheath models,
 * secondary emission, reflection).
 */

/* HEXAPIC plasma wall interactions section */

#include "hexapic.hpp"

bool boundary_in_rank(int bound_i);
void wall_interaction(size_t cell_i, int sp, int p_i, int dir);
void reflect_particle(size_t cell_i, int sp, int p_i, int dir);
void random_order(std::vector<int> &order, int n);
void sigmd_E(REAL &sigmd, REAL &E, size_t cell_i, int sp, int p_i, int dir);
void secondary_emission(REAL &E, REAL sigmd, int sp, int sec_i, int dir,
                        REAL xy_along, REAL z);
REAL gamma(REAL E, REAL sigmd, REAL gamma0, REAL E0, REAL Ep, int sec_type);
REAL gamma_e(REAL E, REAL sigmd, REAL gamma0, REAL E0, REAL Ep);
REAL gamma_i(REAL E, REAL sigmd, REAL gamma0, REAL E0, REAL Ep, int amu);
void inject_secondary(int sp2, int inj_type, REAL &E, REAL T0, int dir,
                      REAL xy_along, REAL z);

void PSI_init() {
  // if boundary not in rank it does not get any additional data
  for (auto &bc : grid.boundaries) {
    bc.sp_sec = {};
    for (int i = 0; i < species.size(); i++)
      bc.sp_sec.push_back({});
  }

  for (int i = 0; i < PSI_init_data.size(); i++) {
    std::vector<std::string> psi = PSI_init_data[i];
    SecondaryEmission se_e;
    int n = psi.size();

    se_e.sp1 = specie_index(psi[1]);
    se_e.sp2 = {};
    for (int j = 0; j < n - 9; j++)
      se_e.sp2.push_back(specie_index(psi[2 + j]));
    se_e.sec_type = std::stoi(psi[n - 7]);
    se_e.inj_type = std::stoi(psi[n - 6]);
    se_e.Ew = std::stod(psi[n - 5]);
    se_e.E0 = std::stod(psi[n - 4]);
    se_e.gamma0 = std::stod(psi[n - 3]);
    se_e.T0 = std::stod(psi[n - 2]);
    se_e.Ep = std::stod(psi[n - 1]);

    if (psi[0] == "xmin")
      if (boundary_in_rank(0))
        grid.boundaries[0].sp_sec[se_e.sp1].push_back(se_e);
      else if (psi[0] == "xmax")
        if (boundary_in_rank(1))
          grid.boundaries[1].sp_sec[se_e.sp1].push_back(se_e);
        else if (psi[0] == "ymin")
          if (boundary_in_rank(2))
            grid.boundaries[2].sp_sec[se_e.sp1].push_back(se_e);
          else if (psi[0] == "ymax")
            if (boundary_in_rank(3))
              grid.boundaries[3].sp_sec[se_e.sp1].push_back(se_e);
            else
              std::cout << "Error: PSI input data has improper shape";
  }

  // Clear the data and release allocated memory
  for (auto &inner : PSI_init_data) {
    inner.clear();
    inner.shrink_to_fit();
  }
  PSI_init_data.clear();
  PSI_init_data.shrink_to_fit();
}

bool boundary_in_rank(int bound_i) {
  // assumes 4 boundaries in the shape of a rectangle
  Boundary bc = grid.boundaries[bound_i];

  return ((bc.j1 >= grid.N0x && bc.j1 <= grid.N0x + grid.Ncx) &&
          (bc.k1 >= grid.N0y && bc.k1 <= grid.N0y + grid.Ncy));
}

void particle_wall_interaction(size_t cell_i, int sp, int p_i, int lrbt) {
  /*
        |      |
     5  |  3   |  7
  ______|______|_____
            |      |
     0  | cell |  1
  ______|______|_____
            |      |
     4  |  2   |  6
        |      |

  TODO: for corner crossing, choose based on travel distance instead of rand()
  */

  if (lrbt < 4) { // just one interaction
    if (lrbt == 0)
      wall_interaction(cell_i, sp, p_i, 0);
    else if (lrbt == 1)
      wall_interaction(cell_i, sp, p_i, 1);
    else if (lrbt == 2)
      wall_interaction(cell_i, sp, p_i, 2);
    else
      wall_interaction(cell_i, sp, p_i, 3);
  }

  else {
    int r1, r2; ///< which wall interaction
    Cell &cell = cells[cell_i];

    if (lrbt == 4) {
      if (frand() < 0.5) { // hit x wall
        cell.part[sp][p_i].y = 0.0001;
        wall_interaction(cell_i, sp, p_i, 0);
      } else { // hit y wall
        cell.part[sp][p_i].x = 0.0001;
        wall_interaction(cell_i, sp, p_i, 2);
      }
    }

    else if (lrbt == 5) {
      if (frand() < 0.5) { // hit x wall
        cell.part[sp][p_i].y = REAL(grid.Ncy) - 0.0001;
        wall_interaction(cell_i, sp, p_i, 0);
      } else { // hit y wall
        cell.part[sp][p_i].x = 0.0001;
        wall_interaction(cell_i, sp, p_i, 3);
      }
    }

    else if (lrbt == 6) {
      if (frand() < 0.5) { // hit x wall
        cell.part[sp][p_i].y = 0.0001;
        wall_interaction(cell_i, sp, p_i, 1);
      } else { // hit y wall
        cell.part[sp][p_i].x = REAL(grid.Ncx) - 0.0001;
        wall_interaction(cell_i, sp, p_i, 2);
      }
    }

    else {
      if (frand() < 0.5) { // hit x wall
        cell.part[sp][p_i].y = REAL(grid.Ncy) - 0.0001;
        wall_interaction(cell_i, sp, p_i, 1);
      } else { // hit y wall
        cell.part[sp][p_i].x = REAL(grid.Ncx) - 0.0001;
        ;
        wall_interaction(cell_i, sp, p_i, 3);
      }
    }
  }
}

void wall_interaction(size_t cell_i, int sp, int p_i, int dir) {
  // direction of boundary dir=0, 1, 2, 3 = xmin, xmax, ymin, ymax
  Boundary &bound = grid.boundaries[dir];

  if (species[sp].charge == 0) { // Neutrals get reflected except if open bc
    if (!(bound.name == "Dielectric" && bound.reflection == 0.0)) // !open
      reflect_particle(cell_i, sp, p_i, dir);
  }

  else if (frand() < bound.reflection) { // Non-neutrals can also get reflected
    reflect_particle(cell_i, sp, p_i, dir);
  }

  else { // PSI
    std::vector<int> order;
    REAL E, sigmd, xy_along;
    Particle &p = cells[cell_i].part[sp][p_i];
    // primary particle E[eV], 1/cos(inc. angle), distance along the other dim.
    sigmd_E(sigmd, E, cell_i, sp, p_i, dir);
    random_order(order, bound.sp_sec[sp].size());

    if ((dir < 2 && p.y >= 0 && p.y < grid.Ncy) ||
        (dir >= 2 && p.x >= 0 && p.x < grid.Ncx)) {

      xy_along = (dir < 2) ? p.y : p.x;

      for (int i = 0; i < order.size(); i++)
        // E is a reference and is decreased iside the function
        secondary_emission(E, sigmd, sp, order[i], dir, xy_along, p.z);
    }

    else {
      std::cout << "Warning: PSI particle outside of grid: x=" << p.x
                << ", y=" << p.y << ", Nx=" << grid.Ncx << ", Ny=" << grid.Ncy
                << ", dir=" << dir << std::endl;
    }
  }
}

void reflect_particle(size_t cell_i, int sp, int p_i, int dir) {
  Particle p = cells[cell_i].part[sp][p_i];
  if (dir < 2) {
    p.vx = -p.vx;
    if (dir == 0)
      p.x = -p.x;
    else
      p.x = 2 * grid.Ncx - p.x;
  } else {
    p.vy = -p.vy;
    if (dir == 2)
      p.y = -p.y;
    else
      p.y = 2 * grid.Ncy - p.y;
  }

  // find proper cell
  cell_i = p.x;
  cell_i *= grid.Ncy;
  cell_i += p.y;
  cells[cell_i].part[sp].push_back(p);
  species[sp].n += 1;

  // particle outside the boundary is removed later
}

void random_order(std::vector<int> &order, int n) {
  if (n < 1)
    order = {};
  else if (n == 1)
    order = {0};
  else {
    order = {};
    for (int i = 0; i < n; i++)
      order.push_back(i);
    for (int i = n - 1; i >= 1; i--)
      std::swap(order[i], order[int(frand() * (i + 1))]);
  }
}

void sigmd_E(REAL &sigmd, REAL &E, size_t cell_i, int sp, int p_i, int dir) {
  // 1/cos(incident angle), energy of a particle
  Particle &p = cells[cell_i].part[sp][p_i];
  REAL v_sq = p.vx * p.vx + p.vy + p.vy + p.vz + p.vz;

  E = species[sp].mass * v_sq / (2 * q_e); ///< [eV]
  if (dir < 2)
    sigmd = sqrt(v_sq) / fabs(p.vx); ///< x direction
  else
    sigmd = sqrt(v_sq) / fabs(p.vy); ///< y direction
}

void secondary_emission(REAL &E, REAL sigmd, int sp, int sec_i, int dir,
                        REAL xy_along, REAL z) {
  // sp - specie of primary particle
  // sec_i - which of secondary emissions
  // direction of boundary 0, 1, 2, 3 = xmin, xmax, ymin, ymax

  SecondaryEmission &SE = grid.boundaries[dir].sp_sec[sp][sec_i];
  REAL gamma_sec; ///< yield

  if (SE.sec_type == 0) { // electron-induced SEE
    if (E <= 0)
      return;

    gamma_sec = gamma_e(E, sigmd, SE.gamma0, SE.E0, SE.Ep);

    // Optional physical cap: max secondaries limited by available energy scale
    // Ew
    int Ncap = INT_MAX;
    if (SE.Ew > 0) {
      Ncap = (int)floor(E / SE.Ew);
      if (Ncap < 0)
        Ncap = 0;
    }

    // Sample count: floor + Bernoulli(frac), fast and unbiased
    int n_emit = (int)gamma_sec;
    REAL frac = gamma_sec - (REAL)n_emit;
    if (frand() < frac)
      ++n_emit;
    if (n_emit > Ncap)
      n_emit = Ncap;
    if (n_emit <= 0)
      return;

    // Secondary energy: reuse Ew (~few eV). If you later add a sampler, seed
    // with T0.
    REAL Esec = (SE.Ew > 0) ? SE.Ew : (REAL)5.0;

    for (int k = 0; k < n_emit; ++k)
      inject_secondary(SE.sp2[0], SE.inj_type, Esec, SE.T0, // T0 kept as-is
                       dir, xy_along, z);
  }

  else if (SE.sec_type == 1) { // ion induced SEE
                               // Expected yield
    REAL delta = gamma_i(E, sigmd, SE.gamma0, SE.E0, SE.Ep,
                         int(species[1].mass / 1e-27));

    // Optional cap by available primary energy scale Ew (few-eV per secondary)
    int Ncap = INT_MAX;
    if (SE.Ew > 0) {
      Ncap = std::max(0, (int)floor(E / SE.Ew));
    }

    // floor + Bernoulli(frac) sampling
    int n_emit = (int)delta;
    REAL frac = delta - (REAL)n_emit;
    if (frand() < frac)
      ++n_emit;
    n_emit = std::min(n_emit, Ncap);
    if (n_emit <= 0)
      return;

    // Few-eV secondary energy
    REAL Esec = (SE.Ew > 0) ? SE.Ew : (REAL)5.0;

    for (int k = 0; k < n_emit; ++k) {
      inject_secondary(SE.sp2[0], SE.inj_type, Esec, SE.T0, dir, xy_along, z);
    }
  }

  else if (SE.sec_type == 3) { // H+ recycling
    gamma_sec = gamma(E, sigmd, SE.gamma0, SE.E0, SE.Ep, 3);
    if (gamma_sec > frand()) { // H injection
      // injection type is 1, SE.T0 replaced by hardcoded 0.5 for scaling
      inject_secondary(SE.sp2[0], 1, E, 0.5, dir, xy_along, z);

      E = 0.0; ///< Check if this is sensible.
    }
  }

  else
    std::cout << "Secondary emission of type " << SE.sec_type
              << " not yet implemented!\n";
}

REAL gamma(REAL E, REAL sigmd, REAL gamma0, REAL E0, REAL Ep, int sec_type) {
  // yield of new particles
  if (sec_type == 0) {
    gamma0 *= (sigmd > 2) ? sigmd : 2.0;
    REAL Eratio = E / E0;
    gamma0 *= Eratio;
    gamma0 *= exp(-2 * (1 - sqrt(Eratio)));
    return gamma0;
  }

  else if (sec_type == 1) {
    // vp saved in E an v0 saved in E0
    gamma0 *= (sigmd > 2) ? sigmd : 2.0;
    gamma0 *= (E - E0);
    return gamma0;
  }

  else if (sec_type == 3) {
    gamma0 = E / (sigmd * sigmd);
    gamma0 = 1.0 / (1 + 0.4 * powf(gamma0, 0.39));
    return gamma0;
  }

  else {
    std::cout << "Sec_type not yet implemented!";
    return 0.0;
  }
}

// One-hump SEE model (Sternglass/Vaughan-type)
// gamma0->delta_max, E0->Emax (auto-fix if too small), Ep->shape s (clamped),
// sigmd->surface factor (fallback 1.0). sec_type==0 => electron-induced SEE.
REAL gamma_e(REAL E, REAL sigmd, REAL gamma0, REAL E0, REAL Ep) {
  if (E <= 0 || gamma0 <= 0)
    return 0;

  // Map & sanitize inputs
  REAL delta_max = gamma0; // peak yield
  REAL Emax =
      (E0 > 50.0) ? E0 : 350.0; // tungsten-ish default if E0 looks wrong
  REAL s_raw = (Ep > 0) ? Ep : (REAL)1.35;
  REAL s = std::min((REAL)1.8, std::max((REAL)0.8, s_raw));
  REAL f_surf = (sigmd > 0) ? sigmd : (REAL)1.0;

  REAL x = E / Emax;
  if (x <= 0)
    return 0;

  // δ(E) = f_surf * δ_max * (E/Emax)^s * exp[s*(1 - E/Emax)]
  REAL delta = f_surf * delta_max * pow(x, s) * exp(s * (1.0 - x));
  return (delta > 0) ? delta : 0;
}

REAL gamma_i(REAL E, REAL sigmd, REAL gamma0, REAL E0, REAL Ep, int amu) {
  if (E <= 0 || gamma0 <= 0)
    return 0;

  // Surface multiplier (oxide/conditioning/roughness)
  const REAL f_surf = (sigmd > 0) ? sigmd : (REAL)1.0;

  // Shape exponent for ions (slightly narrower than electrons)
  const REAL s_raw = (Ep > 0) ? Ep : (REAL)1.2;
  const REAL s = std::min((REAL)1.6, std::max((REAL)0.8, s_raw));

  // Peak energy: use provided E0 if reasonable; otherwise per-amu fallback
  // Typical kinetic IISEE peaks: tens of keV/amu for H/D on metals.
  const REAL k_keV_per_amu = (REAL)30.0e3; // 30 keV/amu default
  const REAL Emax = (E0 > 5.0e3) ? E0 : (k_keV_per_amu * amu); // eV

  const REAL x = E / Emax;
  if (x <= 0)
    return 0;

  const REAL delta = f_surf * gamma0 * pow(x, s) * exp(s * (1.0 - x));
  return (delta > 0) ? delta : 0;
}

void inject_secondary(int sp2, int inj_type, REAL &E, REAL T0, int dir,
                      REAL xy_along, REAL z) {
  // inject new secondary particles of sp2
  // from boundary dir at other dimension xy_along
  Particle p;

  // different injection options
  if (inj_type == 0) { // Maxwell distributed with T0
    REAL v_th = sqrt(q_e * T0 / species[sp2].mass);
    if (dir < 2) {
      p.vx = normvel_vFv() * v_th;
      p.vy = normvel_MaxBol() * v_th;
      p.y = xy_along;
      p.x = std::min(0.9999, p.vx * dt * frand() + 0.0001);
      if (dir == 1) {
        p.vx *= -1;
        p.x = REAL(grid.Ncx) - p.x;
      }
    } else {
      p.vx = normvel_MaxBol() * v_th;
      p.vy = normvel_vFv() * v_th;
      p.x = xy_along;
      p.y = std::min(0.9999, p.vy * dt * frand() + 0.0001);
      if (dir == 3) {
        p.vy *= -1;
        p.y = REAL(grid.Ncy) - p.y;
      }
    }
    p.vz = normvel_MaxBol() * v_th;
    p.z = z;
  }

  else if (inj_type == 1) {
    REAL v_sec, r = frand(), theta = TWOPI * frand();
    v_sec = E * T0;
    v_sec /= (v_sec + 1.0);
    v_sec = E * powf(frand(), v_sec) * (2 * q_e) / species[sp2].mass;
    v_sec = sqrt(v_sec);

    if (dir < 2) {
      p.vx = sqrt(r) * v_sec;
      r = sqrt(1 - r) * v_sec;
      p.vy = r * sin(theta);
      p.vz = r * cos(theta);
      p.y = xy_along;
      p.x = std::min(0.9999, p.vx * dt * frand() + 0.0001);
      if (dir == 1) {
        p.vx *= -1;
        p.x = REAL(grid.Ncx) - p.x;
      }
    } else {
      p.vy = sqrt(r) * v_sec;
      r = sqrt(1 - r) * v_sec;
      p.vx = r * sin(theta);
      p.vz = r * cos(theta);
      p.x = xy_along;
      p.y = std::min(0.9999, p.vy * dt * frand() + 0.0001);
      if (dir == 3) {
        p.vy *= -1;
        p.y = REAL(grid.Ncy) - p.y;
      }
    }
    p.z = z;
  }

  // implement other types

  else {
    std::cout << "Injection type " << inj_type << " not yet implemented.\n";
    return;
  }

  size_t cell_i = p.x * grid.Ncy + p.y;
  // std::cout<<"DEBUG: sp2="<<sp2<<", inj_type="<<inj_type<<", dir="<<dir<<",
  // p.x="<<p.x<<", p.y="<<p.y<<", cell_i="<<cell_i<<std::endl;
  cells[cell_i].part[sp2].push_back(p);
  species[sp2].n += 1;
}
