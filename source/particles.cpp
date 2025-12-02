/**
 * SPDX-License-Identifier: EUPL-1.2
 * SPDX-FileCopyrightText: 2025 Stefan Costea, LeCAD-PEG
 * 
 * @file particles.cpp
 * @brief Particle pusher, weighting/interpolation, and particle management (injection, movement, deposition).
 */

/* HEXAPIC particles section */

#include "hexapic.hpp"

inline static bool is_edge(int i, int j, int Nx, int Ny) {
    return i==0 || j==0 || i==Nx-1 || j==Ny-1;
}

void particles_init() {
	REAL tb2;

    for(auto& sp : species) {
    	sp.chom = q_e*sp.charge/sp.mass;
		tb2 = 0.0;
		for(int j=0; j<3; j++)	{
			sp.tb[j] = hdt*sp.chom*Bf[j];
			tb2 += sp.tb[j]*sp.tb[j]; 
		}
		for(int j=0; j<3; j++)
			sp.sb[j] = 2*sp.tb[j]/(1+tb2);
	}
	srand(time(NULL)*(rmpi+1));
}

static void cell_init(int i, int j) {
	int index = i*grid.Ncy + j;

	// set position in local grid
	cells[index].n0x = i;
	cells[index].n0y = j;
	// empty particle vector
	cells[index].part = {};
	for (int k=0; k<species.size(); k++)
		cells[index].part.push_back({});

	// add neighbouring cells
	std::vector<int> neighb_cells = {index-grid.Ncy, index+grid.Ncy,
	index-1, index+1, index-grid.Ncy-1, index-grid.Ncy+1,
	index+grid.Ncy-1, index+grid.Ncy+1};
	cells[index].neighbours = neighb_cells;
}

static void cells_init_interior () {
	cells_interior.reserve((grid.Ncx-2)*(grid.Ncy-2));
	for (int i=1; i<grid.Ncx-1; i++) 
		for (int j=1; j<grid.Ncy-1; j++) {
			cell_init(i, j);
			int index = i*grid.Ncy + j;
			cells_interior.push_back(&cells[index]); // create pointer to interior cell
			}
	
}

static void cells_init_edge () {
	int i, j, k;
	cells_edge.reserve(2*grid.Ncx+2*grid.Ncy);
	// bottom edge
	for (j=0; j<grid.Ncy; j++) { 
		cell_init(0, j);
		cells_edge.push_back(&cells[j]); // create pointer to edge cell
		}
	// top edge
	for (j=0; j<grid.Ncy; j++) { 
		cell_init(grid.Ncx-1, j);
		cells_edge.push_back(&cells[(grid.Ncx-1)*grid.Ncy+j]); // create pointer to edge cell
		}
	// left edge
	for (i=0; i<grid.Ncx; i++) { 
		cell_init(i, 0);
		cells_edge.push_back(&cells[i*grid.Ncy]); // create pointer to edge cell
		}
	// right edge
	for (i=0; i<grid.Ncx; i++) { 
		cell_init(i, grid.Ncy-1);
		cells_edge.push_back(&cells[i*grid.Ncy+grid.Ncy-1]); // create pointer to edge cell
		}
	// add connections to other ranks as: -(r + 1)
	Neighbours& g_neigh = grid.neighbours;
	for (i=0; i<grid.Ncy; i++) {
		j = grid.Ncy * (grid.Ncx-1) + i;
		cells[i].neighbours[0] = -(g_neigh.rxmin[i]+1);
		cells[j].neighbours[1] = -(g_neigh.rxmax[i]+1);
		if (i > 0) {
			cells[i].neighbours[4] = -(g_neigh.rxmin[i-1]+1);
			cells[j].neighbours[6] = -(g_neigh.rxmax[i-1]+1);
		}
		else {
			cells[i].neighbours[4] = -(g_neigh.corners[0]+1);
			cells[j].neighbours[6] = -(g_neigh.corners[2]+1);
		}
		if (i < grid.Ncy-1) {
			cells[i].neighbours[5] = -(g_neigh.rxmin[i+1]+1);
			cells[j].neighbours[7] = -(g_neigh.rxmax[i+1]+1);
		}
		else {
			cells[i].neighbours[5] = -(g_neigh.corners[1]+1);
			cells[j].neighbours[7] = -(g_neigh.corners[3]+1);
		}
	}
	for (i=0; i<grid.Ncx; i++) {
		j = i*grid.Ncy;
		k = (i+1)*grid.Ncy-1;
		cells[j].neighbours[2] = -(g_neigh.rymin[i]+1);
		cells[k].neighbours[3] = -(g_neigh.rymax[i]+1);
		if (i > 0) {
			cells[j].neighbours[4] = -(g_neigh.rymin[i-1]+1);
			cells[k].neighbours[5] = -(g_neigh.rymax[i-1]+1);
		}
		if (i < grid.Ncx-1) {
			cells[j].neighbours[6] = -(g_neigh.rymin[i+1]+1);
			cells[k].neighbours[7] = -(g_neigh.rymax[i+1]+1);
		}
	}
}

void cells_init() {
	cells = {};
	cells.resize(grid.Ncx*grid.Ncy); // allocate once so pointers are stable
	cells_init_interior();
	cells_init_edge();	
}

void grid_init(int argc, char **args) {

	int i, j;
	grid.sort_boundaries();

	dV = grid.dx*grid.dy*1; ///< cell volume
	grid.p2n = np2c_global/dV;
	Sconst = -q_e*grid.p2n/EPS0; ///< particle to charge density
	Sconst *= grid.dx*grid.dx; ///< from multiplication of matrix A
	hdt = dt/2; ///< half time step
	// number of cells in each direction


	// allocate memory for meshes

	int n_corners = grid.Nx*grid.Ny;
	rho.reserve(n_corners); // computed at cell corners
	grid.V_store.reserve(n_corners);

	for(i=0; i<n_corners; i++) {
		rho.push_back(0);
		grid.V_store.push_back(0);
		}

	int n_cells = grid.Ncx*grid.Ncy;
	for(auto& sp : species) {
		sp.n_store.reserve(n_cells); ///< computed at cell center
		sp.T_store.reserve(n_cells);
		sp.v_store.reserve(n_cells);
		for(i=0; i<n_cells; i++) {
			sp.n_store.push_back(0);
			sp.T_store.push_back(0);
			sp.v_store.push_back(0);
			}
		}
	t_store.reserve(1); ///< time
	t_store.push_back(0);
}

void inject_particle(int sp, REAL x, REAL y) {
	Particle p;
	if ((!x)&&(!y)) {
		p.x = (grid.Nx-2)*frand()+1;
		p.y = (grid.Ny-2)*frand()+1;
		}
	else {
		p.x = x;
		p.y = y;
		}
	p.z = 0;
	p.vx = 0; ///< (10+frand())/20;
	p.vy = 0; ///< (10+frand())/20;
	p.vz = 0; ///< (10+frand())/20;
	int cell_i = int(p.y) + int(p.x) * grid.Ncy;
	species[sp].n += 1;
	cells[cell_i].part[sp].push_back(p);
}

// work on a single cell
static void part2grid_cell(Cell& cell) {
/* Map particles charge to grid 
  J,(i,j)_______px__________ J+grid.Ny,(i+1,j)
		|       |          |
		|  S1   |    S4    |
	  py|_______|__________|
		|       |          |
		|  S2   |    S3    |
		|       |          |
		|_______|__________|
 J+1,(i,j+1)		   J+1+grid.Ny,(i+1,j+1)
*/
	int i, j, k, J;
	REAL px, py, S1, S2, S3, S4;
    int q;

    i = cell.n0x;
	j = cell.n0y;

	J = i*grid.Ny+j;
	for (k=0; k<species.size(); k++) {
		Species& sp = species[k];
		q = sp.charge;
		if (q == 0) continue;
		for (auto& p : cell.part[k]) {
			px = p.x-i;
			py = p.y-j;
			S1 = px*py;
			S2 = px*(1-py);
			S3 = (1-px)*(1-py);
			S4 = (1-px)*py;
			// The charge to a grid point is proportional to opposite area
			rho[J]		+= q*S3;
			rho[J+1]	+= q*S4;
			rho[J+grid.Ny]	+= q*S2;
			rho[J+grid.Ny+1] += q*S1;
			}
		}
}

void part2grid_edge() {
	// must be called first, before part2grid_interior()
	size_t i, j, J;

	// zero the charge density array
	std::fill(rho.begin(), rho.end(), 0);

	for (Cell* p : cells_edge) part2grid_cell(*p);

	// J = i*grid.Ny+j;

	// update horizontal edges
	for(i=0; i<grid.Nx; i++) {
		J = i*grid.Ny;
		rho[J] *= Sconst;
		J += grid.Ny-1;
		rho[J] *= Sconst;
		}

	// update vertical edges
	for(j=0; j<grid.Ny; j++) {
		J = j;
		rho[J] *= Sconst;
		J += (grid.Nx-1)*grid.Ny;
		rho[J] *= Sconst;
		}

	// deploy for exchange
	hypre_post_rho_halo();
}

void part2grid_interior() {
	// must be called second, after part2grid_edge()
	size_t i, j, J;

	for (Cell* p : cells_interior) part2grid_cell(*p);

	// update interior cells
	for(i=1; i<grid.Nx-1; i++) 
		for(j=1; j<grid.Ny-1; j++) {
			J = i*grid.Ny+j;
			rho[J] *= Sconst;
			}
}


void grid2part(REAL Px, REAL Py, double* V, REAL *Ef) {
	/* Map electric field to particle position 
	            
  J,(i,j)_______px_________J+grid.Ny,(i+1,j)
		|       |          |
		|       |          |
	  py|_______|__________| 
		|       |          |
		|       |          |
		|       |          |
		|_______|__________|
J+1,(i,j+1)				J+1+grid.Ny,(i+1,j+1) */


	int i, j, J;
	REAL px, py; // position inside the cell
	i = Px; // cell index in X direction
	j = Py; // cell index in Y direction
	px = Px-i; // X position inside the cell
	py = Py-j; // Y position inside the cell	
	J = i*grid.Ny+j; // global index of cell in linearized mesh

	// The solution of the matrix is the exact potential
	Ef[0] = -((V[J+grid.Ny]-V[J])*(1-py) + (V[J+grid.Ny+1]-V[J+1])*py)/grid.dx;
	Ef[1] = -((V[J+1]-V[J])*(1-px) + (V[J+grid.Ny+1]-V[J+grid.Ny])*px)/grid.dy;
}

inline static void grid2cell(int J, double* V, REAL *V_cell) {
	/* Map electric field to cell nodes
	b-l, t-l, b-r, t-r
	*/
	V_cell[0] = V[J];
	V_cell[1] = V[J+1];
	V_cell[2] = V[J+grid.Ny];
	V_cell[3] = V[J+grid.Ny+1];
}

inline static void cell2part(const REAL px, const REAL py, const REAL* V_cell, REAL *Ef) {
	Ef[0] = -((V_cell[2]-V_cell[0])*(1-py) + (V_cell[3]-V_cell[1])*py)/grid.dx;
	Ef[1] = -((V_cell[1]-V_cell[0])*(1-px) + (V_cell[3]-V_cell[2])*px)/grid.dy;
}



inline static void particle_mover_boris_cell(Cell& cell) {
	int i, j, k, J;
    REAL V_cell[4], Ef[2], px, py, vx, vy, vz, chom, *sb, *tb;
	int sp, sc;

	i = cell.n0x;
	j = cell.n0y;
	J = j + grid.Ny*i;
	grid2cell(J, V_local, &V_cell[0]);

    for (k=0; k<species.size(); k++) {
    	Species& sp = species[k];
    	if(tstep % sp.subcycle != 0) continue;
		if(sp.charge!=0.0) for (auto& p : cell.part[k]) { // charged particles
			/* Calculate electric field at particle position */
			px = p.x - i;
			py = p.y - j;
		    cell2part(px, py, V_cell, &Ef[0]);

			/* First Half-Step Velocity Update with Electric Field */
		    p.vx += sp.chom*Ef[0]*hdt;
			p.vy += sp.chom*Ef[1]*hdt;

			/* Intermediate velocity update due to rotation in magnetic field */
			vx = p.vx + (p.vy*sp.tb[2]-p.vz*sp.tb[1]);
			vy = p.vy + (p.vz*sp.tb[0]-p.vx*sp.tb[2]);
			vz = p.vz + (p.vx*sp.tb[1]-p.vy*sp.tb[0]);

			/* Final velocity update due to rotation in magnetic field */
			p.vx += vy*sp.sb[2]-vz*sp.sb[1];
			p.vy += vz*sp.sb[0]-vx*sp.sb[2];
			p.vz += vx*sp.sb[1]-vy*sp.sb[0];

			/* Second Half-Step Velocity Update with Electric Field */
			p.vx += sp.chom*Ef[0]*hdt;
			p.vy += sp.chom*Ef[1]*hdt;

			/* Full-Step Position Update */
			p.x += sp.subcycle*p.vx*dt/grid.dx; ///< normalized to cell size
			p.y += sp.subcycle*p.vy*dt/grid.dy; ///< normalized to cell size
			p.z += sp.subcycle*p.vz*dt/grid.dy; ///< normalized to cell size

			// /* Save trajectory if test particle exists */
			// if(sp.part.size()==1) {
			// 	trajx.push_back(p.x);
	    	// 	trajy.push_back(p.y);
	    	// 	trajz.push_back(p.z);
	    	// 	}
			}
		else for (auto& p : cell.part[k]) { // neutrals
			/* Full-Step Position Update */
			p.x += sp.subcycle*p.vx*dt/grid.dx;
			p.y += sp.subcycle*p.vy*dt/grid.dy;
			p.z += sp.subcycle*p.vz*dt/grid.dy;
			}
		}
}

void particle_mover_boris(const std::vector<Cell*>& cells_x) {
	for (Cell* p : cells_x) 
		particle_mover_boris_cell(*p);
}

void particle_boundaries() {

	// erase particles out of boundary for open boundary
	// put them to the other side if periodic boundary
	// reflect if neumann boundary condition

	bool erase_particle;
	int n_sp, r, lrbt; ///< neighbour position
	std::vector<std::vector<std::vector<REAL>>>&
	send_data_3d = grid.neighbours.send_data_3d;


	for (int i=0; i<species.size(); i++) {
		Species& sp = species[i];
		if(tstep % sp.subcycle != 0) continue;

		for (int cell_i=0; cell_i<cells.size(); cell_i++) {
			Cell& cell = cells[cell_i];
			std::vector<Particle>& part = cell.part[i];
		    for (int p_i=0; p_i < part.size(); ) {
		    	Particle p = part[p_i];
		        erase_particle = false; 
		        
		        // x_min boundary
		        if (p.x <= cell.n0x) {
		        	// corner
		        	if (p.y <= cell.n0y) lrbt = 4;
		        	// corner
		        	else if (p.y >= cell.n0y+1) lrbt = 5;
		        	else lrbt = 0;
		            erase_particle = true;
		        }

		        // x_max boundary
		        else if (p.x >= cell.n0x+1) {
		        	// corner
					if (p.y <= cell.n0y) lrbt = 6;
		        	// corner
		        	else if (p.y >= cell.n0y+1) lrbt = 7;
		        	else lrbt = 1;
		            erase_particle = true;
		        }

		        // y_min boundary
		        else if (p.y <= cell.n0y) {
		            // corners already taken into account
		            lrbt = 2;
		            erase_particle = true;
		        }

		        // y_max boundary
		        else if (p.y >= cell.n0y+1) {
		            // corners already taken into account
		            lrbt = 3;
		            erase_particle = true;
		        }

		        if (erase_particle) {
		        	r = cell.neighbours[lrbt];
		        	if (r < 0) { // reserved for sending to other rank or PWI
		      			r = -r - 1;
		        		if (r!=rmpi)  // sending to other rank
		       				send_data_3d[i][r].insert(send_data_3d[i][r].end(),
		        			{p.x, p.y, p.z, p.vx, p.vy, p.vz});
		        		else  // particle reachead an actual boundary
		        			particle_wall_interaction(cell_i, i, p_i, lrbt);
		        		remove_particle(cell_i, i, p_i);
					}

		        	else { // sending to different cell inside the same rank
		        		cells[r].part[i].push_back(p);
		            		std::swap(part[p_i], part[part.size()-1]);
	        			part.pop_back();
	        		}
		        }
		        else p_i++;
		    }
		}
	}

	std::vector<std::vector<REAL>>& send_data = grid.neighbours.send_data;
	for (int j : grid.neighbours.all) {
		send_data[j].push_back(static_cast<double>(species.size()));
		for (int i=0; i<species.size(); i++) send_data[j].push_back(0.0);
	}
	
	for (int j : grid.neighbours.all)
		for (int i=0; i<species.size(); i++) {
	    	n_sp = (send_data_3d[i][j].size()) / 6;
	    	send_data[j][i+1] = n_sp;
	    	if (n_sp > 0) {
	    		send_data[j].insert(send_data[j].end(),
	    			send_data_3d[i][j].begin(), send_data_3d[i][j].end());
	    	}
	    	send_data_3d[i][j] = {};
	    }

	send_receive_particles();
	add_received_particles();
}


void add_received_particles() {
	for (int r : grid.neighbours.all) {

		std::vector<REAL>& recv_data = grid.neighbours.recv_data[r];
		std::vector<std::vector<std::vector<REAL>>>&
		send_data_3d = grid.neighbours.send_data_3d;
		size_t i, j, lrbt, cell_i, sp, n_sp;
		int r2;

		n_sp = static_cast<size_t>(recv_data[0]);
		i = 1 + n_sp;

		for (sp=0; sp<n_sp; sp++) {
			// add particles per each specie in proper cell
			species[sp].n += static_cast<int>(recv_data[sp+1]);
			for (j=0; j<static_cast<int>(recv_data[sp+1]); j++) {
				Particle p;
				p.x = recv_data[i] + grid.neighbours.shift_xy[r][0];
				p.y = recv_data[i+1] + grid.neighbours.shift_xy[r][1];
				p.z = recv_data[i+2];
				p.vx = recv_data[i+3];
				p.vy = recv_data[i+4];
				p.vz = recv_data[i+5];
				i += 6;

				// Check if particle inside new rank
				if (p.x >= 0 && p.x < grid.Ncx &&
					p.y >= 0 && p.y < grid.Ncy) [[likely]] {
					cell_i = size_t(p.y) + size_t(p.x) * grid.Ncy;
					cells[cell_i].part[sp].push_back(p);
				}
				else [[unlikely]] { // send it on in the next step or erase
					std::cout << "Warrning: Received particle of specie "
					<< species[sp].name << " jumped over multiple ranks!\n";
					}
				} 
			}
		}
}


void send_receive_particles() {
    std::vector<MPI_Request> data_requests;

    // Send and receive the actual data
    for (int r : grid.neighbours.all) {
        // Send data
	    std::vector<REAL>& data = grid.neighbours.send_data[r];
	    data_requests.push_back(MPI_Request());
	    MPI_Isend(data.data(), data.size(), MPI_DOUBLE, r, 1,
	    	MPI_COMM_WORLD, &data_requests.back());
        // Receive data
	    data_requests.push_back(MPI_Request());
	    MPI_Irecv(grid.neighbours.recv_data[r].data(),
	    	grid.neighbours.max_recv_size, MPI_DOUBLE, r, 1,
	    	MPI_COMM_WORLD, &data_requests.back());
    }
    MPI_Waitall(data_requests.size(), data_requests.data(),
    	MPI_STATUSES_IGNORE);

    // empty send buffer
    for (int r : grid.neighbours.all)
    	std::vector<REAL>().swap(grid.neighbours.send_data[r]);
}


void initial_particle_load() {
	int i, j, Npart, cell_i;
	for(i=0; i<species.size(); i++) {
		Species& sp = species[i];
		for(auto ld : sp.load) {
			if(grid.x1s > ld.x1f || grid.x1f < ld.x1s ||
			   grid.x2s > ld.x2f || grid.x2f < ld.x2s)
				continue; ///< if load area is outside the grid
			if(ld.method==0) {
				Npart = ld.density*(ld.x1f-ld.x1s)*(ld.x2f-ld.x2s)/np2c_global;
				for(j=0; j<Npart; j++) {
					Particle p;
					p.x = (frand()*(ld.x1f-ld.x1s)+ld.x1s)/grid.dx - grid.N0x;
					p.y = (frand()*(ld.x2f-ld.x2s)+ld.x2s)/grid.dy - grid.N0y;
					p.z = 0.0;
					p.vx = ld.v1drift + ld.v1thermal*normvel_MaxBol();
					p.vy = ld.v2drift + ld.v2thermal*normvel_MaxBol();
					p.vz = ld.v3drift + ld.v3thermal*normvel_MaxBol();
					cell_i = int(p.y) + int(p.x) * grid.Ncy;
					cells[cell_i].part[i].push_back(p);
					sp.n += 1;
					}
				}
			else if (!rmpi){
				std::cout << "LOAD: Initial load method " << ld.method << 
					" not implemented. Use method=0, i.e. uniform." << "\n";
				exit(1);
				}
			}
		}
}

void remove_particle(int cell_i, int sp, int p) {
	cells[cell_i].part[sp][p] = cells[cell_i].part[sp].back();
	cells[cell_i].part[sp].pop_back();
	species[sp].n -= 1;
}

void inject_particles() {
	int i, j, Npart, intx, inty, cell_i;
	REAL src;
	for(i=0; i<species.size(); i++) {
		Species& sp = species[i];
		for(auto& inj : sp.inject) {
			if(inj.I>0) {
				intx = inj.j2-inj.j1;
				inty = inj.k2-inj.k1;
				src = inj.I*dt/abs(sp.charge)/q_e/np2c_global;
				// just injection along one of dimensions
				if(intx)	src *= intx;
				if(inty)	src *= inty;
				inj.f += src;
				Npart = inj.f;
				for(j=0; j<Npart; j++) {
					Particle p;
					if(intx) {
						p.x = frand()*intx + inj.j1 - grid.N0x;
						p.vx = inj.v1thermal*normvel_MaxBol();
						}
					else {
						p.x = frand() + inj.j1 - grid.N0x;
						p.vx = inj.normal*inj.v1thermal*normvel_vFv();
						}

					if(inty) {
						p.y = frand()*inty + inj.k1 - grid.N0y;
						p.vy = inj.v2thermal*normvel_MaxBol();
						}
					else {
						p.y = frand() + inj.k1 - grid.N0y;
						p.vy = inj.normal*inj.v2thermal*normvel_vFv();
						}
					p.z = 0.0;
					p.vz = inj.v3thermal*normvel_MaxBol();

					cell_i = p.x;
					cell_i *= grid.Ncy;
					cell_i += p.y;

					cells[cell_i].part[i].push_back(p);
					sp.n += 1;
					}
				inj.f -= Npart;
				}
			}
	}
}


REAL normvel_MaxBol() {

	/* Box-Muller method
	 W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. Y. Vetterling.
     Numerical Recipes in C: The Art of Scientiﬁc Programming. Cambridge
     University Press, Cambridge, UK, 1988.	*/

	static std::vector<REAL> vel;
	REAL n1=1, n2=1, R2, val;

	// remove previously used value
	if(vel.size())	vel.pop_back();

	// return available value 		
	if(vel.size())	return(vel.back());		

	// else generate new values
	while (n1*n1+n2*n2 > 1) {
		n1 = frand()*2-1;
		n2 = frand()*2-1;
		}

	R2 = n1*n1+n2*n2;
	val = sqrt(-2*log(R2)/R2);

	vel.push_back(n1*val);
	vel.push_back(n2*val);

	return(vel.back());	
}


REAL normvel_vFv() {

	/* 
	 W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. Y. Vetterling.
     Numerical Recipes in C: The Art of Scientiﬁc Programming. Cambridge
     University Press, Cambridge, UK, 1988.	*/

	return sqrt(-log(frand()));
}
