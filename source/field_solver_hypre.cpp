/**
 * SPDX-License-Identifier: EUPL-1.2
 * SPDX-FileCopyrightText: 2025 Stefan Costea, LeCAD-PEG
 *
 * HYPRE implementation of the field solver
 * For now just basic version with fixed solver and preconditioner
 * It uses Struct as for our grid this gives us the most optimized solutions
 * Beacause Struct does not have ghost values direct access with API,
 * 	manual MPI send/recive of this values is done

*/

#include "hexapic.hpp"

constexpr int HYPRE_TAG_V = 101;
constexpr int HYPRE_TAG_RHO = 102;

void hypre_init_geo() {
  if (grid.boundaries[0].name == "Periodic")
    hypre.periodic[0] = grid.Nxg - 1;
  else
    hypre.periodic[0] = 0;
  if (grid.boundaries[2].name == "Periodic")
    hypre.periodic[1] = grid.Nyg - 1;
  else
    hypre.periodic[1] = 0;

  // if it continues with next rank-> edge nodes are solved in next domain
  hypre.ilower[0] = grid.N0x;
  hypre.ilower[1] = grid.N0y;

  if ((grid.neighbours.rxmax[0] == rmpi) && (!hypre.periodic[0])) // boundary
    hypre.iupper[0] = hypre.ilower[0] + grid.Ncx;
  else
    hypre.iupper[0] = hypre.ilower[0] + grid.Ncx - 1;

  if ((grid.neighbours.rymax[0] == rmpi) && (!hypre.periodic[1])) // boundary
    hypre.iupper[1] = hypre.ilower[1] + grid.Ncy;
  else
    hypre.iupper[1] = hypre.ilower[1] + grid.Ncy - 1;

  hypre.size_extend = (grid.Ncx + 1) * (grid.Ncy + 1);
  hypre.nx = hypre.iupper[0] - hypre.ilower[0] + 1;
  hypre.ny = hypre.iupper[1] - hypre.ilower[1] + 1;
  hypre.size = hypre.nx * hypre.ny;

  hypre.size_s.resize(nmpi);
  hypre.size_r.resize(nmpi);
  hypre.V_s_idx.resize(nmpi);
  hypre.V_r_idx.resize(nmpi);
  hypre.send_buff.resize(nmpi);
  hypre.recv_buff.resize(nmpi);

  hypre.bcs = {0, 0, 0, 0};
  std::vector<Boundary> &bcs = grid.boundaries;

  if (bcs[0].j1 == hypre.ilower[0]) {
    if (bcs[0].name == "Equipotential")
      hypre.bcs[0] = 1;
    else if (bcs[0].name == "Dielectric")
      hypre.bcs[0] = 2;
  }

  if (bcs[1].j1 == hypre.iupper[0]) {
    if (bcs[1].name == "Equipotential")
      hypre.bcs[1] = 1;
    else if (bcs[1].name == "Dielectric")
      hypre.bcs[1] = 2;
  }

  if (bcs[2].k1 == hypre.ilower[1]) {
    if (bcs[2].name == "Equipotential")
      hypre.bcs[2] = 1;
    else if (bcs[2].name == "Dielectric")
      hypre.bcs[2] = 2;
  }

  if (bcs[3].k1 == hypre.iupper[1]) {
    if (bcs[3].name == "Equipotential")
      hypre.bcs[3] = 1;
    else if (bcs[3].name == "Dielectric")
      hypre.bcs[3] = 2;
  }
}

void hypre_init_ghost() {
  Neighbours &neigh = grid.neighbours;
  int r, r_prev, r_c;

  // sending and receiving V

  // corner
  r_c = neigh.corners[0];
  if (r_c != rmpi)
    hypre.V_s_idx[r_c].push_back(0);

  // left boundary
  r_prev = r_c;
  for (int j = 0; j < neigh.rxmin.size(); j++) {
    r = neigh.rxmin[j];
    if (r == rmpi)
      continue;

    if (r != r_prev)
      hypre.V_s_idx[r].push_back(j);
    r_prev = r;

    if (j != neigh.rxmin.size() - 1)
      hypre.V_s_idx[r].push_back(j + 1);
    else if (grid.Ny == hypre.ny)
      hypre.V_s_idx[r].push_back(j + 1);
  }

  // bottom boundary
  r_prev = r_c;
  for (int i = 0; i < neigh.rymin.size(); i++) {
    r = neigh.rymin[i];
    if (r == rmpi)
      continue;

    if (r != r_prev)
      hypre.V_s_idx[r].push_back(i * grid.Ny);
    r_prev = r;

    if (i != neigh.rymin.size() - 1)
      hypre.V_s_idx[r].push_back((i + 1) * grid.Ny);
    else if (grid.Nx == hypre.nx)
      hypre.V_s_idx[r].push_back((i + 1) * grid.Ny);
  }

  // right boundary
  for (int j = 0; j < neigh.rxmax.size(); j++) {
    r = neigh.rxmax[j];
    if (r == rmpi)
      continue;

    hypre.V_r_idx[r].push_back(grid.Ny * (grid.Nx - 1) + j);
  }

  // top boundary
  for (int i = 0; i < neigh.rymax.size(); i++) {
    r = neigh.rymax[i];
    if (r == rmpi)
      continue;

    hypre.V_r_idx[r].push_back(grid.Ny * (i + 1) - 1);
  }

  // top corner
  r_c = neigh.corners[3];
  if (r_c != rmpi) {
    hypre.V_r_idx[r_c].push_back(grid.Ny * grid.Nx - 1);
  } else {
    r_c = neigh.rxmax[grid.Ncy - 1];
    if (r_c != rmpi) {
      hypre.V_r_idx[r_c].push_back(grid.Ny * grid.Nx - 1);
    } else {
      r_c = neigh.rymax[grid.Ncx - 1];
      if (r_c != rmpi) {
        hypre.V_r_idx[r_c].push_back(grid.Ny * grid.Nx - 1);
      }
    }
  }
  // if double boundary then no receiving...

  // rho sending and receiving is just in the opposite direction
  // also using opposite send recive buffers!

  // size, buffers
  for (int i = 0; i < nmpi; i++) {
    hypre.size_r[i] = hypre.V_r_idx[i].size();
    hypre.size_s[i] = hypre.V_s_idx[i].size();
    hypre.send_buff[i].resize(hypre.size_s[i]);
    hypre.recv_buff[i].resize(hypre.size_r[i]);
  }

  // --- Build persistent MPI requests for halo exchanges ---
  hypre_persist_v_reqs.clear();
  hypre_persist_rho_reqs.clear();

  // For V exchange: Recv into recv_buff (size_r), Send from send_buff (size_s)
  for (int r : grid.neighbours.all) {
    if (hypre.size_r[r] > 0) {
      MPI_Request req;
      MPI_Recv_init(hypre.recv_buff[r].data(), hypre.size_r[r], MPI_DOUBLE, r,
                    HYPRE_TAG_V, comm, &req);
      hypre_persist_v_reqs.push_back(req);
    }
  }
  for (int r : grid.neighbours.all) {
    if (hypre.size_s[r] > 0) {
      MPI_Request req;
      MPI_Send_init(hypre.send_buff[r].data(), hypre.size_s[r], MPI_DOUBLE, r,
                    HYPRE_TAG_V, comm, &req);
      hypre_persist_v_reqs.push_back(req);
    }
  }

  // For rho exchange (opposite buffers): Recv into send_buff (size_s), Send
  // from recv_buff (size_r)
  for (int r : grid.neighbours.all) {
    if (hypre.size_s[r] > 0) {
      MPI_Request req;
      MPI_Recv_init(hypre.send_buff[r].data(), hypre.size_s[r], MPI_DOUBLE, r,
                    HYPRE_TAG_RHO, comm, &req);
      hypre_persist_rho_reqs.push_back(req);
    }
  }
  for (int r : grid.neighbours.all) {
    if (hypre.size_r[r] > 0) {
      MPI_Request req;
      MPI_Send_init(hypre.recv_buff[r].data(), hypre.size_r[r], MPI_DOUBLE, r,
                    HYPRE_TAG_RHO, comm, &req);
      hypre_persist_rho_reqs.push_back(req);
    }
  }
}

void hypre_init_solver() {
  // grid
  HYPRE_StructGridCreate(comm, 2, &grid_h);
  HYPRE_StructGridSetExtents(grid_h, hypre.ilower, hypre.iupper);
  HYPRE_StructGridSetPeriodic(grid_h, hypre.periodic);
  HYPRE_StructGridAssemble(grid_h);

  // 5-point stencil
  HYPRE_StructStencilCreate(2, 5, &stencil_h);
  int offsets[5][2] = {{0, 0}, {-1, 0}, {1, 0}, {0, -1}, {0, 1}};
  for (int i = 0; i < 5; i++)
    HYPRE_StructStencilSetElement(stencil_h, i, offsets[i]);

  // matrix
  int nentries = 5;
  double values[5];
  int stencil_indices[5] = {0, 1, 2, 3, 4};
  HYPRE_StructMatrixCreate(comm, grid_h, stencil_h, &A_h);
  HYPRE_StructMatrixInitialize(A_h);

  dx2_inv = 1.0;
  dy2_inv = (grid.dx * grid.dx) / (grid.dy * grid.dy);

  std::vector<Boundary> &bcs = grid.boundaries;
  std::vector<int> bcs_ij;

  for (int j = hypre.ilower[1]; j <= hypre.iupper[1]; j++) {
    for (int i = hypre.ilower[0]; i <= hypre.iupper[0]; i++) {
      int index[2] = {i, j};
      values[0] = values[1] = values[2] = values[3] = values[4] = 0.0;

      // for periodic normal stencil
      // for dirichlet fixed point
      // for open the difference to outer point is 0

      bcs_ij = {0, 0, 0, 0};
      bool d = false;

      if ((i == bcs[0].j1 && i == bcs[0].j2)) { // xmin boundary
        if (bcs[0].name == "Equipotential") {
          bcs_ij[0] = 1;
          d = true;
        } else if (bcs[0].name == "Dielectric")
          bcs_ij[0] = 2;
      }

      if (i == bcs[1].j1 && i == bcs[1].j2) { // xmax boundary
        if (bcs[1].name == "Equipotential") {
          bcs_ij[1] = 1;
          d = true;
        } else if (bcs[1].name == "Dielectric")
          bcs_ij[1] = 2;
      }

      if (j == bcs[2].k1 && j == bcs[2].k2) { // ymin boundary
        if (bcs[2].name == "Equipotential") {
          bcs_ij[2] = 1;
          d = true;
        } else if (bcs[2].name == "Dielectric")
          bcs_ij[2] = 2;
      }

      if (j == bcs[3].k1 && j == bcs[3].k2) { // ymax boundary
        if (bcs[3].name == "Equipotential") {
          bcs_ij[3] = 1;
          d = true;
        } else if (bcs[3].name == "Dielectric")
          bcs_ij[3] = 2;
      }

      // if any dirichlet -> fixed boundary!
      if (d)
        values[0] = dx2_inv;
      else {
        if (bcs_ij[0] != 2) {
          values[0] -= dx2_inv;
          values[1] += dx2_inv;
        }
        if (bcs_ij[1] != 2) {
          values[0] -= dx2_inv;
          values[2] += dx2_inv;
        }
        if (bcs_ij[2] != 2) {
          values[0] -= dy2_inv;
          values[3] += dy2_inv;
        }
        if (bcs_ij[3] != 2) {
          values[0] -= dy2_inv;
          values[4] += dy2_inv;
        }
      }

      HYPRE_StructMatrixSetValues(A_h, index, nentries, stencil_indices,
                                  values);
    }
  }
  HYPRE_StructMatrixAssemble(A_h);

  // --- PFMG solver (geometric) ---
  // HYPRE_StructSolver pfmg;
  HYPRE_StructPFMGCreate(comm, &pfmg);
  HYPRE_StructPFMGSetTol(pfmg, 1e-5);
  HYPRE_StructPFMGSetMaxIter(pfmg, 200);
  HYPRE_StructPFMGSetRelChange(pfmg, 0);
  HYPRE_StructPFMGSetNumPreRelax(pfmg, 1);
  HYPRE_StructPFMGSetNumPostRelax(pfmg, 1);
  HYPRE_StructPFMGSetRelaxType(
      pfmg, 2); // Red-Black Gauss-Siedel, often better than Jacobi
  HYPRE_StructPFMGSetLogging(pfmg, 0);

  V_h = (double *)calloc(hypre.size, sizeof(double));
  rho_h = (double *)calloc(hypre.size, sizeof(double));
  V_local = (double *)calloc(hypre.size_extend, sizeof(double));

  HYPRE_StructVectorCreate(comm, grid_h, &b_h);
  HYPRE_StructVectorCreate(comm, grid_h, &x_h);
  HYPRE_StructVectorInitialize(b_h);
  HYPRE_StructVectorInitialize(x_h);

  HYPRE_StructVectorSetBoxValues(b_h, hypre.ilower, hypre.iupper, rho_h);
  HYPRE_StructVectorSetBoxValues(x_h, hypre.ilower, hypre.iupper, V_h);

  HYPRE_StructVectorAssemble(b_h);
  HYPRE_StructVectorAssemble(x_h);

  // One-time PFMG setup (outside timestep loop)
  HYPRE_StructPFMGSetup(pfmg, A_h, b_h, x_h);
}

void hypre_init() {
  hypre_init_geo();
  hypre_init_ghost();
  hypre_init_solver();
}

/* Split send and receive in order to overlap with computation */

// ---- V halo (pack + start) ----
void hypre_post_V_halo() {
  // pack
  for (int r : grid.neighbours.all)
    for (int j = 0; j < hypre.size_s[r]; ++j)
      hypre.send_buff[r][j] = V_local[hypre.V_s_idx[r][j]];

  if (!hypre_persist_v_reqs.empty())
    MPI_Startall((int)hypre_persist_v_reqs.size(), hypre_persist_v_reqs.data());
}

void hypre_finish_V_halo() {
  if (!hypre_persist_v_reqs.empty())
    MPI_Waitall((int)hypre_persist_v_reqs.size(), hypre_persist_v_reqs.data(),
                MPI_STATUSES_IGNORE);

  // unpack
  for (int r : grid.neighbours.all)
    for (int j = 0; j < hypre.size_r[r]; ++j)
      V_local[hypre.V_r_idx[r][j]] = hypre.recv_buff[r][j];
}

// ---- rho halo (pack + start) ----
void hypre_post_rho_halo() {
  // note: your rho packs into recv_buff in the opposite direction
  for (int r : grid.neighbours.all)
    for (int j = 0; j < hypre.size_r[r]; ++j)
      hypre.recv_buff[r][j] = rho[hypre.V_r_idx[r][j]];

  if (!hypre_persist_rho_reqs.empty())
    MPI_Startall((int)hypre_persist_rho_reqs.size(),
                 hypre_persist_rho_reqs.data());
}

void hypre_finish_rho_halo() {
  if (!hypre_persist_rho_reqs.empty())
    MPI_Waitall((int)hypre_persist_rho_reqs.size(),
                hypre_persist_rho_reqs.data(), MPI_STATUSES_IGNORE);

  // accumulate
  for (int r : grid.neighbours.all)
    for (int j = 0; j < hypre.size_s[r]; ++j)
      rho[hypre.V_s_idx[r][j]] += hypre.send_buff[r][j];
}

// void hypre_exchange_ghost_V() {
//     // Pack send buffers from V_local
//     for (int r : grid.neighbours.all) {
//         for (int j = 0; j < hypre.size_s[r]; ++j) {
//             hypre.send_buff[r][j] = V_local[hypre.V_s_idx[r][j]];
//         }
//     }

//     // Start all persistent recvs + sends (created in hypre_init_ghost)
//     if (!hypre_persist_v_reqs.empty()) {
//         MPI_Startall(static_cast<int>(hypre_persist_v_reqs.size()),
//         hypre_persist_v_reqs.data());
//         MPI_Waitall(static_cast<int>(hypre_persist_v_reqs.size()),
//         hypre_persist_v_reqs.data(), MPI_STATUSES_IGNORE);
//     }

//     // Unpack into V_local
//     for (int r : grid.neighbours.all) {
//         for (int j = 0; j < hypre.size_r[r]; ++j) {
//             V_local[hypre.V_r_idx[r][j]] = hypre.recv_buff[r][j];
//         }
//     }
// }

void hypre_update_Vlocal() {
  int i;
  for (int x_i = 0; x_i < hypre.nx; x_i++) {
    i = x_i * (grid.Ncy + 1);
    for (int y_i = 0; y_i < hypre.ny; y_i++) {
      V_local[i] = V_h[y_i * hypre.nx + x_i];
      i++;
    }
  }
}

void hypre_get_V() {
  HYPRE_StructVectorGetBoxValues(x_h, hypre.ilower, hypre.iupper, V_h);

  // first overwrite what you already have then send, receive and write
  hypre_update_Vlocal();

  // hypre_exchange_ghost_V();
  hypre_post_V_halo();
}

void hypre_field_solver() { HYPRE_StructPFMGSolve(pfmg, A_h, b_h, x_h); }

// void hypre_exchange_ghost_rho() {
//     // Pack from rho into recv_buff (note: opposite buffers by design)
//     for (int r : grid.neighbours.all) {
//         for (int j = 0; j < hypre.size_r[r]; ++j) {
//             hypre.recv_buff[r][j] = rho[hypre.V_r_idx[r][j]];
//         }
//     }

//     // Start all persistent recvs + sends for rho
//     if (!hypre_persist_rho_reqs.empty()) {
//         MPI_Startall(static_cast<int>(hypre_persist_rho_reqs.size()),
//         hypre_persist_rho_reqs.data());
//         MPI_Waitall(static_cast<int>(hypre_persist_rho_reqs.size()),
//         hypre_persist_rho_reqs.data(), MPI_STATUSES_IGNORE);
//     }

//     // Accumulate into rho from send_buff (received data placed there)
//     for (int r : grid.neighbours.all) {
//         for (int j = 0; j < hypre.size_s[r]; ++j) {
//             rho[hypre.V_s_idx[r][j]] += hypre.send_buff[r][j];
//         }
//     }
// }

void hypre_rho_update() {

  int i;
  for (int x_i = 0; x_i < hypre.nx; x_i++) {
    i = x_i * (grid.Ncy + 1);
    for (int y_i = 0; y_i < hypre.ny; y_i++) {
      rho_h[y_i * hypre.nx + x_i] = rho[i];
      i++;
    }
  }

  // correct boundary entries
  // just dirichlet, periodic and open are already okay!
  double b_value;
  if (hypre.bcs[0] == 1) { // Equipotential
    b_value = dx2_inv * grid.boundaries[0].C;
    for (int i = 0; i < hypre.ny; i++) {
      rho_h[i * hypre.nx] = b_value;
    }
  }
  if (hypre.bcs[1] == 1) { // Equipotential
    b_value = dx2_inv * grid.boundaries[1].C;
    for (int i = 0; i < hypre.ny; i++) {
      rho_h[(i + 1) * hypre.nx - 1] = b_value;
    }
  }
  if (hypre.bcs[2] == 1) { // Equipotential
    b_value = dx2_inv * grid.boundaries[2].C;
    for (int i = 0; i < hypre.nx; i++) {
      rho_h[i] = b_value;
    }
  }
  if (hypre.bcs[3] == 1) { // Equipotential
    b_value = dx2_inv * grid.boundaries[3].C;
    int j = hypre.nx * (hypre.ny - 1);
    for (int i = 0; i < hypre.nx; i++) {
      rho_h[j + i] = b_value;
    }
  }
}

void hypre_source_update() {
  // send, receive then write
  // hypre_exchange_ghost_rho();

  hypre_finish_rho_halo();

  hypre_rho_update();

  // finalize vector
  HYPRE_StructVectorSetBoxValues(b_h, hypre.ilower, hypre.iupper, rho_h);
  HYPRE_StructVectorAssemble(b_h);
}

void hypre_cleanup() {
  free(V_h);
  free(rho_h);
  free(V_local);
  // Free persistent MPI requests if created
  for (auto &req : hypre_persist_v_reqs) {
    MPI_Request_free(&req);
  }
  for (auto &req : hypre_persist_rho_reqs) {
    MPI_Request_free(&req);
  }
  hypre_persist_v_reqs.clear();
  hypre_persist_rho_reqs.clear();

  HYPRE_StructPFMGDestroy(pfmg);

  HYPRE_StructMatrixDestroy(A_h);
  HYPRE_StructVectorDestroy(b_h);
  HYPRE_StructVectorDestroy(x_h);
  HYPRE_StructGridDestroy(grid_h);
  HYPRE_StructStencilDestroy(stencil_h);
}
