/**
 * SPDX-License-Identifier: EUPL-1.2
 * SPDX-FileCopyrightText: 2025 Stefan Costea, LeCAD-PEG
 *
 * @file domain_decomposition.cpp
 * @brief Cartesian grid topology, ghost cells, and MPI domain decomposition
 * utilities.
 */

/* HEXAPIC domain decomposition */

#include "hexapic.hpp"
#include <iomanip>

std::vector<int> rank_domain(int Ncx, int Ncy, int n_mpi, int r_mpi) {
  // Decomposition by bisecting larger dimension in each step
  int n = n_mpi;
  // remaining ranks in current part of domain
  int r = r_mpi;
  // which of remaining n in current part of domain
  int x0 = 0, y0 = 0, nx = Ncx, ny = Ncy;
  // Recursivly split the domain
  // x0, y0 is bottom left corner of current domain
  // nx, ny is number of grid points in this part of domain per dimension
  int d;

  while (true) {
    // remaining part belongs to remaining rank
    if (n == 1)
      break;

    // we need at least 2 nodes per cell in each dimension
    else if (nx < 2 && ny < 2) {
      if (n > 1 && r > 0) {
        // some ranks are empty
        std::cout << "To many MPI processes to efficently decompose"
                  << " full domain by bisections\n";
        nx = 0;
        ny = 0;
      }
      break;
    }
    if (nx >= ny) {
      d = n / 2;
      if (r < d) {
        nx = std::ceil(double(nx) * d / n);
        n = d;
      } else {
        r = r - d;
        x0 = x0 + std::ceil(double(nx) * d / n);
        nx = nx - std::ceil(double(nx) * d / n);
        n = n - d;
      }
    } else {
      d = n / 2;
      if (r < d) {
        ny = std::ceil(double(ny) * d / n);
        n = d;
      } else {
        r = r - d;
        y0 = y0 + std::ceil(double(ny) * d / n);
        ny = ny - std::ceil(double(ny) * d / n);
        n = n - d;
      }
    }
  }
  std::vector<int> domain_part = {x0, y0, nx, ny};
  return domain_part;
}

void plot_decomposed_domain() {

  int Ncx = grid.Ncx;
  int Ncy = grid.Ncy;

  int x0, y0, nx, ny;
  std::vector<int> domain_i;
  // we plot cells (along dimension #cells = #nodes-1)
  int plot[Ncy][Ncy] = {};

  for (int i = 0; i < nmpi; i++) {
    domain_i = rank_domain(Ncx, Ncy, nmpi, i);
    x0 = domain_i[0];
    y0 = domain_i[1];
    nx = domain_i[2];
    ny = domain_i[3];
    for (int x = x0; x < x0 + nx; x++) {
      for (int y = y0; y < y0 + ny; y++) {
        plot[x][y] = i;
      }
    }
  }
  for (int y = Ncy - 1; y >= 0; y--) {
    for (int x = 0; x < Ncx; x++) {
      std::cout << std::setw(2) << plot[x][y] << " ";
    }
    std::cout << "\n";
  }
}

void decompose_domain() {
  int x0, y0, nx, ny, val;
  int i, x, y;
  std::vector<std::vector<int>> rank_x0y0 = {};

  int count = (grid.Ncx + 1) * 2 + (grid.Ncy + 1) * 2 + 8 + nmpi * 2;

  if (!rmpi) {
    std::vector<std::vector<int>> rank_map(grid.Ncx,
                                           std::vector<int>(grid.Ncy, -1));

    std::vector<std::vector<int>> rdom;
    rdom.resize(nmpi);
    for (i = 0; i < nmpi; i++)
      rdom[i].reserve(count);

    for (i = 0; i < nmpi; i++) {
      auto rd = rank_domain(grid.Ncx, grid.Ncy, nmpi, i);
      rdom[i].push_back(rd[0]);
      rdom[i].push_back(rd[1]);
      rdom[i].push_back(rd[2]);
      rdom[i].push_back(rd[3]);
      x0 = rdom[i][0];
      y0 = rdom[i][1];
      nx = rdom[i][2];
      ny = rdom[i][3];
      rank_x0y0.push_back({x0, y0});
      for (x = x0; x < x0 + nx; x++)
        for (y = y0; y < y0 + ny; y++)
          rank_map[x][y] = i;
    }

    int print_domain = 0;
    if (print_domain) {
      // print domain decomposition
      for (y = grid.Ncy - 1; y >= 0; y--) {
        for (x = 0; x < grid.Ncx; x++)
          std::cout << std::setw(2) << rank_map[x][y] << " ";
        std::cout << "\n";
      }
    }

    // populate neighbours
    for (i = 0; i < nmpi; i++) {
      x0 = rdom[i][0];
      y0 = rdom[i][1];
      nx = rdom[i][2];
      ny = rdom[i][3];
      // xmin neighbours
      for (y = y0; y < y0 + ny; y++) {
        if (x0 == 0)
          val = i;
        else
          val = rank_map[x0 - 1][y];
        rdom[i].push_back(val);
        if (!i)
          grid.neighbours.rxmin.push_back(val);
      }
      // xmax neighbours
      for (y = y0; y < y0 + ny; y++) {
        if (x0 + nx == grid.Ncx)
          val = i;
        else
          val = rank_map[x0 + nx + 1][y];
        rdom[i].push_back(val);
        if (!i)
          grid.neighbours.rxmax.push_back(val);
      }
      // ymin neighbours
      for (x = x0; x < x0 + nx; x++) {
        if (y0 == 0)
          val = i;
        else
          val = rank_map[x][y0 - 1];
        rdom[i].push_back(val);
        if (!i)
          grid.neighbours.rymin.push_back(val);
      }
      // ymax neighbours
      for (x = x0; x < x0 + nx; x++) {
        if (y0 + ny == grid.Ncy)
          val = i;
        else
          val = rank_map[x][y0 + ny + 1];
        rdom[i].push_back(val);
        if (!i)
          grid.neighbours.rymax.push_back(val);
      }

      // corners
      // bottom-left
      if (x0 == 0 || y0 == 0)
        val = i;
      else
        val = rank_map[x0 - 1][y0 - 1];
      rdom[i].push_back(val);
      if (!i)
        grid.neighbours.corners.push_back(val);
      // top-left
      if (x0 == 0 || y0 + ny == grid.Ncy)
        val = i;
      else
        val = rank_map[x0 - 1][y0 + ny];
      rdom[i].push_back(val);
      if (!i)
        grid.neighbours.corners.push_back(val);
      // bottom-right
      if (x0 + nx == grid.Ncx || y0 == 0)
        val = i;
      else
        val = rank_map[x0 + nx][y0 - 1];
      rdom[i].push_back(val);
      if (!i)
        grid.neighbours.corners.push_back(val);
      if (x0 + nx == grid.Ncx || y0 + ny == grid.Ncy)
        val = i;
      else
        val = rank_map[x0 + nx][y0 + ny];
      rdom[i].push_back(val);
      if (!i)
        grid.neighbours.corners.push_back(val);

      // x0, y0 for each rank to have shifts between neighbours
      for (int j = 0; j < nmpi; j++) {
        rdom[i].push_back(rdom[j][0]);
        rdom[i].push_back(rdom[j][1]);
      }

      if (i)
        MPI_Send(rdom[i].data(), rdom[i].size(), MPI_INT, i, 0, comm);
    }

    grid.Nxg = grid.Ncx + 1;
    grid.Nyg = grid.Ncy + 1;
    grid.N0x = rdom[0][0];
    grid.N0y = rdom[0][1];
    grid.Ncx = rdom[0][2];
    grid.Ncy = rdom[0][3];
    grid.Nx = rdom[0][2] + 1;
    grid.Ny = rdom[0][3] + 1;

  } else {
    int buf[count];
    MPI_Status *status;
    MPI_Recv(&buf, count, MPI_INT, 0, 0, comm, status);

    grid.Nxg = grid.Ncx + 1;
    grid.Nyg = grid.Ncy + 1;
    grid.N0x = buf[0];
    grid.N0y = buf[1];
    grid.Ncx = buf[2];
    grid.Ncy = buf[3];
    grid.Nx = buf[2] + 1;
    grid.Ny = buf[3] + 1;

    for (i = 4; i < 4 + (grid.Ncy); i++)
      grid.neighbours.rxmin.push_back(buf[i]);
    for (i = 4 + (grid.Ncy); i < 4 + 2 * (grid.Ncy); i++)
      grid.neighbours.rxmax.push_back(buf[i]);
    for (i = 4 + 2 * (grid.Ncy); i < 4 + 2 * (grid.Ncy) + (grid.Ncx); i++)
      grid.neighbours.rymin.push_back(buf[i]);
    for (i = 4 + 2 * (grid.Ncy) + (grid.Ncx);
         i < 4 + 2 * (grid.Ncy) + 2 * (grid.Ncx); i++)
      grid.neighbours.rymax.push_back(buf[i]);
    for (i = 4 + 2 * (grid.Ncy) + 2 * (grid.Ncx);
         i < 8 + 2 * (grid.Ncy) + 2 * (grid.Ncx); i++)
      grid.neighbours.corners.push_back(buf[i]);
    for (int j = 0; j < nmpi; j++) {
      i = 8 + 2 * (grid.Ncy) + 2 * (grid.Ncx) + j * 2;
      rank_x0y0.push_back({buf[i], buf[i + 1]});
    }
  }

  // add all possible connections
  for (int r : grid.neighbours.rxmax)
    grid.neighbours.all.insert(r);
  for (int r : grid.neighbours.rxmin)
    grid.neighbours.all.insert(r);
  for (int r : grid.neighbours.rymax)
    grid.neighbours.all.insert(r);
  for (int r : grid.neighbours.rymin)
    grid.neighbours.all.insert(r);
  for (int r : grid.neighbours.corners)
    grid.neighbours.all.insert(r);
  grid.neighbours.all.erase(rmpi);

  // add all possible neighbour shifts
  std::vector<std::vector<REAL>> empty(nmpi);
  grid.neighbours.shift_xy = empty;
  int shift_x, shift_y;
  for (int r : grid.neighbours.all) {
    shift_x = rank_x0y0[r][0] - grid.N0x;
    shift_y = rank_x0y0[r][1] - grid.N0y;
    grid.neighbours.shift_xy[r].push_back(shift_x);
    grid.neighbours.shift_xy[r].push_back(shift_y);
  }

  // update grid lengths
  grid.x1s += grid.N0x * grid.dx;
  grid.x1f = grid.x1s + (grid.Ncx) * grid.dx;
  grid.x2s += grid.N0y * grid.dy;
  grid.x2f = grid.x2s + (grid.Ncy) * grid.dy;

  // update particle load and injection locations
  for (auto &sp : species) {
    for (auto &ld : sp.load) {
      // if load location is outside the local domain
      if (ld.x1s > grid.x1f || ld.x2s > grid.x2f || ld.x1f < grid.x1s ||
          ld.x2f < grid.x2s)
        ld.density = 0; ///< no density
      // otherwise truncate if needed
      else {
        if (ld.x1s < grid.x1s)
          ld.x1s = grid.x1s;
        if (ld.x1f > grid.x1f)
          ld.x1f = grid.x1f;
        if (ld.x2s < grid.x2s)
          ld.x2s = grid.x2s;
        if (ld.x2f > grid.x2f)
          ld.x2f = grid.x2f;
      }
    }
    for (auto &inj : sp.inject) {
      // if inject location is outside the local domain
      if (inj.j1 > grid.N0x + (grid.Ncx) || inj.k1 > grid.N0y + (grid.Ncy) ||
          inj.j2 < grid.N0x || inj.k2 < grid.N0y)
        inj.I = 0; ///< no current

      else {
        // change current to current density
        if (inj.j2 - inj.j1)
          inj.I /= inj.j2 - inj.j1;
        if (inj.k2 - inj.k1)
          inj.I /= inj.k2 - inj.k1;
        // truncate if needed
        if (inj.j1 < grid.N0x)
          inj.j1 = grid.N0x;
        if (inj.j2 > grid.N0x + (grid.Ncx))
          inj.j2 = grid.N0x + (grid.Ncx);
        if (inj.k1 < grid.N0y)
          inj.k1 = grid.N0y;
        if (inj.k2 > grid.N0y + (grid.Ncy))
          inj.k2 = grid.N0y + (grid.Ncy);
      }
    }
  }

  // prepare send and receive buffer
  // (circumference) * (#part(sp) per cell) * (#species) * (size of particle)
  grid.neighbours.max_recv_size =
      (2 * grid.Nx + 2 * grid.Ny) * 1000 * species.size() * 6;
  std::vector<std::vector<REAL>> empty0(nmpi);
  std::vector<REAL> empty1(grid.neighbours.max_recv_size);
  grid.neighbours.send_data = empty0;
  grid.neighbours.recv_data = empty0;
  for (int r : grid.neighbours.all)
    grid.neighbours.recv_data[r] = empty1;
  std::vector<std::vector<std::vector<REAL>>> empty2(
      species.size(), std::vector<std::vector<REAL>>(nmpi));
  grid.neighbours.send_data_3d = empty2;

  int print_neighbours = 0;
  if (print_neighbours) {
    // print neighbours
    printf("\nRank %d: x0=%d, nx=%d, y0=%d, ny=%d\n", rmpi, grid.N0x, grid.Ncx,
           grid.N0y, grid.Ncy);
    printf("rxmin:");
    for (auto j : grid.neighbours.rxmin)
      printf("%d ", j);
    printf("\n");
    printf("rymin:");
    for (auto j : grid.neighbours.rymin)
      printf("%d ", j);
    printf("\n");
    printf("rxmax:");
    for (auto j : grid.neighbours.rxmax)
      printf("%d ", j);
    printf("\n");
    printf("rymax:");
    for (auto j : grid.neighbours.rymax)
      printf("%d ", j);
    printf("\n");
    // print corners
    printf("corners: b-l=%d, t-l=%d, b-r=%d, t-r=%d \n",
           grid.neighbours.corners[0], grid.neighbours.corners[1],
           grid.neighbours.corners[2], grid.neighbours.corners[3]);
    printf("all neighbours and shifts:\n");
    for (int r : grid.neighbours.all) {
      std::cout << rmpi << " to neighbour " << r
                << " shift_xy=" << grid.neighbours.shift_xy[r][0] << ", "
                << grid.neighbours.shift_xy[r][1] << "\n";
    }
    printf("\n");
  }
}
