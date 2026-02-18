/**
 * SPDX-License-Identifier: EUPL-1.2
 * SPDX-FileCopyrightText: 2025 Stefan Costea, LeCAD-PEG
 *
 * @file diagnostics.cpp
 * @brief Diagnostics and output: field and particle moments, probes, and file
 * writers.
 */

/* HEXAPIC diagnostics section */

#include "hexapic.hpp"

void create_output() {

  int i, j;
  REAL vth_x, vth_y, vth_z, xth_X, var;
  std::vector<REAL> vx_store, vy_store, vz_store;

  // simulation time
  t_store[0] = tstep * dt;

  // electric potential
  for (i = 0; i < grid.Nx * grid.Ny; i++)
    grid.V_store[i] = V_local[i];

  // particle density, fluid velocity and temperature
  for (j = 0; j < species.size(); j++) {
    Species &sp = species[j];
    std::fill(sp.n_store.begin(), sp.n_store.end(), 0);
    std::fill(sp.v_store.begin(), sp.v_store.end(), 0);
    std::fill(sp.T_store.begin(), sp.T_store.end(), 0);
    vx_store = sp.v_store;
    vy_store = sp.v_store;
    vz_store = sp.v_store;

    for (i = 0; i < cells.size(); i++)
      for (auto &p : cells[i].part[j]) {
        sp.n_store[i]++; ///< particle count
        vx_store[i] += p.vx;
        vy_store[i] += p.vy;
        vz_store[i] += p.vz;
      }

    for (i = 0; i < grid.Ncx * grid.Ncy; i++) {

      vx_store[i] /= sp.n_store[i]; ///< fluid velocity
      vy_store[i] /= sp.n_store[i];
      vz_store[i] /= sp.n_store[i];
      sp.v_store[i] =
          std::sqrt(vx_store[i] * vx_store[i] + vy_store[i] * vy_store[i] +
                    vz_store[i] * vz_store[i]);
    }

    for (i = 0; i < cells.size(); i++)
      for (auto &p : cells[i].part[j]) {
        vth_x = p.vx - vx_store[i];
        vth_y = p.vy - vy_store[i];
        vth_z = p.vz - vz_store[i];
        sp.T_store[i] += vth_x * vth_x;
        sp.T_store[i] += vth_y * vth_y;
        sp.T_store[i] += vth_z * vth_z;
      }

    for (i = 0; i < grid.Ncx * grid.Ncy; i++) {
      var = sp.T_store[i] / sp.n_store[i]; ///< variance 3d
      // T(eV) = variance * m / 3 * e
      // (dividing by 3 works as 'averaging' temperature)
      sp.T_store[i] = var * sp.mass / (3 * q_e); ///< temperature in eV
    }

    for (auto &c : sp.n_store)
      c *= grid.p2n; ///< particle density
  }
}

void save_output() {

  // TODO: assign pointers to output buffer instead of copying the values

  // Clear output buffer
  for (auto &vec : chunks_float)
    vec.clear();
  chunks_float.clear();

  int s, i, j, J, J0;

  Iteration it = series.writeIterations()[tstep];
  Datatype datatype_real = determineDatatype<REAL>();

  // save time
  MeshRecordComponent time_mesh =
      it.meshes["time"][MeshRecordComponent::SCALAR];
  Dataset dataset2 = Dataset(datatype_real, {1ul});
  time_mesh.resetDataset(dataset2);
  chunks_float.push_back(t_store);
  time_mesh.storeChunk(chunks_float.back(), {0ul}, {1ul});

  // save potential
  MeshRecordComponent V_mesh = it.meshes["V"][MeshRecordComponent::SCALAR];
  Extent global_extent1 = {1ul * grid.Nxg, 1ul * grid.Nyg};
  Dataset dataset1 = Dataset(datatype_real, global_extent1);
  Offset chunk_offset1 = {1ul * grid.N0x, 1ul * grid.N0y};
  Extent chunk_extent1 = {1ul * grid.Nx, 1ul * grid.Ny};
  V_mesh.resetDataset(dataset1);
  chunks_float.push_back(grid.V_store);
  V_mesh.storeChunk(chunks_float.back(), chunk_offset1, chunk_extent1);

  Extent global_extent = {1ul * (grid.Nxg - 1), 1ul * (grid.Nyg - 1)};
  Dataset dataset = Dataset(datatype_real, global_extent);
  Offset chunk_offset = {1ul * grid.N0x, 1ul * grid.N0y};
  Extent chunk_extent = {1ul * grid.Ncx, 1ul * grid.Ncy};

  // save density
  for (int i = 0; i < species.size(); i++) {
    MeshRecordComponent n_mesh =
        it.meshes["n" + std::to_string(i)][MeshRecordComponent::SCALAR];
    n_mesh.resetDataset(dataset);
    chunks_float.push_back(species[i].n_store);
    n_mesh.storeChunk(chunks_float.back(), chunk_offset, chunk_extent);
  }

  // save fluid velocity
  for (int i = 0; i < species.size(); i++) {
    MeshRecordComponent v_mesh =
        it.meshes["v" + std::to_string(i)][MeshRecordComponent::SCALAR];
    v_mesh.resetDataset(dataset);
    chunks_float.push_back(species[i].v_store);
    v_mesh.storeChunk(chunks_float.back(), chunk_offset, chunk_extent);
  }

  // save temperature
  for (int i = 0; i < species.size(); i++) {
    MeshRecordComponent T_mesh =
        it.meshes["T" + std::to_string(i)][MeshRecordComponent::SCALAR];
    T_mesh.resetDataset(dataset);
    chunks_float.push_back(species[i].T_store);
    T_mesh.storeChunk(chunks_float.back(), chunk_offset, chunk_extent);
  }

  // save phase-space
  if (ParticleDiagnosticFlag) {
    // collect info from all ranks
    int nsp = species.size();
    unsigned long np_sum[nsp], np_sp[nsp], np_r[nsp], np_sum_r[nsp];

    for (i = 0; i < nsp; i++) {
      np_sp[i] = species[i].n; ///< local extent
      np_sum_r[i] = 0;         ///< offset
    }

    // compute global extent
    // calculate the sum per index across all ranks
    MPI_Allreduce(np_sp, np_sum, nsp, MPI_UNSIGNED_LONG, MPI_SUM,
                  MPI_COMM_WORLD);

    // compute offset
    // Send data to higher-ranked processes
    for (i = rmpi + 1; i < nmpi; i++)
      MPI_Send(&np_sp, nsp, MPI_UNSIGNED_LONG, i, 1, MPI_COMM_WORLD);
    // Receive data from lower-ranked processes and accumulate the sum
    for (i = 0; i < rmpi; i++) {
      MPI_Recv(&np_r, nsp, MPI_UNSIGNED_LONG, i, 1, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      for (j = 0; j < nsp; j++)
        np_sum_r[j] += np_r[j];
    }

    // save particle data
    for (int isp = 0; isp < nsp; isp++) {

      std::vector<REAL> x_store, y_store, z_store, vx_store, vy_store, vz_store;

      ParticleSpecies specie = it.particles[std::to_string(isp)];

      // specie["charge"][RecordComponent::SCALAR]
      //         .resetDataset({Datatype::DOUBLE, {1}})
      //         .makeConstant(q[isp]);

      // specie["mass"][RecordComponent::SCALAR]
      //         .resetDataset({Datatype::DOUBLE, {1}})
      //         .makeConstant(m[isp]);

      // specie["weight"][RecordComponent::SCALAR]
      //         .resetDataset({Datatype::DOUBLE, {1}})
      //         .makeConstant(wt[isp]);

      for (Cell &cell : cells)
        for (auto &p : cell.part[isp]) {
          x_store.push_back(p.x + grid.N0x);
          y_store.push_back(p.y + grid.N0y);
          z_store.push_back(p.z);
          vx_store.push_back(p.vx);
          vy_store.push_back(p.vy);
          vz_store.push_back(p.vz);
        }

      RecordComponent rec_x = specie["position"]["x"];
      RecordComponent rec_y = specie["position"]["y"];
      RecordComponent rec_z = specie["position"]["z"];
      RecordComponent rec_vx = specie["velocity"]["x"];
      RecordComponent rec_vy = specie["velocity"]["y"];
      RecordComponent rec_vz = specie["velocity"]["z"];

      Extent gl_ext = {np_sum[isp]};
      Offset off = {np_sum_r[isp]};
      Extent ext = {np_sp[isp]};

      Dataset ds = Dataset(datatype_real, gl_ext);

      rec_x.resetDataset(ds);
      rec_y.resetDataset(ds);
      rec_z.resetDataset(ds);
      rec_vx.resetDataset(ds);
      rec_vy.resetDataset(ds);
      rec_vz.resetDataset(ds);

      if (x_store.size()) {
        chunks_float.push_back(x_store);
        x_store.clear();
        rec_x.storeChunk(chunks_float.back(), off, ext);

        chunks_float.push_back(y_store);
        y_store.clear();
        rec_y.storeChunk(chunks_float.back(), off, ext);

        chunks_float.push_back(z_store);
        z_store.clear();
        rec_z.storeChunk(chunks_float.back(), off, ext);

        chunks_float.push_back(vx_store);
        vx_store.clear();
        rec_vx.storeChunk(chunks_float.back(), off, ext);

        chunks_float.push_back(vy_store);
        vy_store.clear();
        rec_vy.storeChunk(chunks_float.back(), off, ext);

        chunks_float.push_back(vz_store);
        vz_store.clear();
        rec_vz.storeChunk(chunks_float.back(), off, ext);
      }
    }
  }
}
