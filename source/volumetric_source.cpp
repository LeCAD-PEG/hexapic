/**
 * SPDX-License-Identifier: EUPL-1.2
 * SPDX-FileCopyrightText: 2025 Stefan Costea, Miha Radez, LeCAD-PEG
 * 
 * @file volumetric_source.cpp
 * @brief Volumetric source terms (heating/ionization profiles, source sampling).
 */

/* HEXAPIC volumetric source and heating section */

#include "hexapic.hpp"

void order_sources();
void make_inject_basis();
void source_decomposition();
void source_part(int source_i);
void source_specie();
void source_inject(int sp, int sp_src, int src);
void heating_cell(int cell_i, int sp, int sp_src, int src);


void source_init() {
	order_sources();
	make_inject_basis();
	source_decomposition();
	source_specie();
}

void order_sources() {
	// order sources index in vector in case of different order in input
	int n = sources.size();
	std::vector<Source> new_sources;
	for (int i=0; i<n; i++) new_sources.push_back({});
	for (Source src : sources)
		new_sources[src.index] = src;
	sources = new_sources;
	// test
	for (int i=0; i < sources.size(); i++)
		if (sources[i].index != i)
			std::cout << "ERROR: Source order is not ok!\n";
}

void make_inject_basis() {
	// x' in B irection, other two perpendicular
	REAL basis[3][3] = {0.0}; 
	REAL norm = sqrt(Bf[0]*Bf[0] + Bf[1]*Bf[1] + Bf[2]*Bf[2]);
	if (norm == 0.0) {
		basis[0][0] = 1.0;
		basis[1][1] = 1.0;
		basis[2][2] = 1.0;
		for (int i=0; i<3; i++)
    	for (int j=0; j<3; j++)
    		Bb[i][j] = basis[i][j];
		return;
	}
	// direction of magnetic field
	basis[0][0] = Bf[0]/norm;
	basis[0][1] = Bf[1]/norm;
	basis[0][2] = Bf[2]/norm;

    REAL temp[3] = {0.0, 0.0, -1.0};
    if (fabs(basis[0][2]) > 0.9) {  // Avoid collinearity
        temp[1] = 1.0;
        temp[2] = 0.0;
    }
    basis[1][0] = basis[0][1] * temp[2] - basis[0][2] * temp[1];
    basis[1][1] = basis[0][2] * temp[0] - basis[0][0] * temp[2];
    basis[1][2] = basis[0][0] * temp[1] - basis[0][1] * temp[0];
    norm = sqrt(basis[1][0] * basis[1][0] +
            	basis[1][1] * basis[1][1] +
            	basis[1][2] * basis[1][2]);
    basis[1][0] /= norm;
    basis[1][1] /= norm;
    basis[1][2] /= norm;

    basis[2][0] = basis[0][1] * basis[1][2] - basis[0][2] * basis[1][1];
    basis[2][1] = basis[0][2] * basis[1][0] - basis[0][0] * basis[1][2];
    basis[2][2] = basis[0][0] * basis[1][1] - basis[0][1] * basis[1][0];

    for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
    	Bb[i][j] = basis[i][j];
}


void source_decomposition() {
	/*
	check active status
	translate to local grid (add cx, cy)
	calculate injection and heat per timestep
	*/
	for (int i=0; i<sources.size(); i++) {
		Source& src = sources[i];
		if (src.j1 >= grid.N0x+grid.Ncx || src.j2 <= grid.N0x ||
			src.k1 >= grid.N0y+grid.Ncy || src.k2 <= grid.N0y) {
			src.active = 0;
			continue;
		}
		src.active = 1;
		src.Lx = src.j2 - src.j1;
		src.Ly = src.k2 - src.k1;
		src.cx = double(src.j1) + double(src.Lx)/2;
		src.cy = double(src.k1) + double(src.Ly)/2;
		// translate
		src.cx -= grid.N0x;
		src.cy -= grid.N0y;
		src.j1 -= grid.N0x;
		src.j2 -= grid.N0x;
		src.k1 -= grid.N0y;
		src.k2 -= grid.N0y;
		// just inside the local grid
		src.j1 = std::max(src.j1, 0);
		src.j2 = std::min(src.j2, grid.Ncx);
		src.k1 = std::max(src.k1, 0);
		src.k2 = std::min(src.k2, grid.Ncy);
		src.lx = src.j2 - src.j1;
		src.ly = src.k2 - src.k1;

		src.n_heat = src.heat_freq * dt;
		// calculate injection portion of full source
		source_part(i);
	}
}

void source_part(int source_i) {
	/*
	calculates n_inj_step for this rank
	normalisation is done in a way that some n_inj means same number of
	injected particles globally no matter what the source shape is.
	Number of injected particles is therefore just t*n_cell*dV*n_inj 
	*/ 
	Source& src = sources[source_i];
	REAL F, F1, remain, stddev, j1_0, k1_0;
	src.n_inj_step = dV * dt * src.n_inj / np2c_global;

	// cumulative distributions along x
	if (src.shape_x == 0) // uniform distribution
		src.n_inj_step *= src.lx;

	else if (src.shape_x == 1) { // cos distribution
		j1_0 = src.cx - src.Lx / 2;
		src.n_inj_step *= src.Lx;
		src.Fx0 = (1 - cos((src.j1-j1_0)*PI/src.Lx))/2;
		src.Fx1 = (1 - cos((src.j2-j1_0)*PI/src.Lx))/2;
		src.n_inj_step *= src.Fx1 - src.Fx0;
	}
	else { // normal distribution
		/* 
		normal distribution with boundary at three sigma
		remain outside the distribution
		scale with 1 / (1-remain) at the end
		Fx and Fy should not be scaled
		*/
		stddev = src.Lx / 6;
		src.Fx0 = 1 + std::erf((src.j1 - src.cx) / (stddev * std::sqrt(2)));
		src.Fx0 *= 0.5;
		src.Fx1 = 1 + std::erf((src.j2 - src.cx) / (stddev * std::sqrt(2)));
		src.Fx1 *= 0.5;
		F = src.Fx1 - src.Fx0;
		src.n_inj_step *= F;
		src.n_inj_step *= src.Lx;
		remain = 2.699796063e-3;
		src.n_inj_step *= 1 / (1 - remain);

		src.Fxcell = {};
		for (int i=src.j1; i<src.j2; i++) {
			F1 = 1 + std::erf(((i+1) - src.cx) / (stddev * std::sqrt(2)));
			F1 *= 0.5;
			src.Fxcell.insert({(F1-src.Fx0)/F, i});
		}
	}

	// cumulative distributions along y
	if (src.shape_y == 0)
		src.n_inj_step *= src.ly;
	else if (src.shape_y == 1) {
		k1_0 = src.cy - src.Ly / 2;
		src.n_inj_step *= src.Ly;
		src.Fy0 = (1 - cos((src.k1-k1_0)*PI/src.Ly))/2;
		src.Fy1 = (1 - cos((src.k2-k1_0)*PI/src.Ly))/2;
		src.n_inj_step *= src.Fy1 - src.Fy0;
	}
	else {
		stddev = src.Ly / 6;
		src.Fy0 = 1 + std::erf((src.k1 - src.cy) / (stddev * std::sqrt(2)));
		src.Fy0 *= 0.5;
		src.Fy1 = 1 + std::erf((src.k2 - src.cy) / (stddev * std::sqrt(2)));
		src.Fy1 *= 0.5;
		F = src.Fy1 - src.Fy0;
		src.n_inj_step *= F;
		src.n_inj_step *= src.Ly;
		remain = 2.699796063e-3;
		src.n_inj_step *= 1 / (1 - remain);

		src.Fycell = {};
		for (int i=src.k1; i<src.k2; i++) {
			F1 = 1 + std::erf(((i+1) - src.cy) / (stddev * std::sqrt(2)));
			F1 *= 0.5;
			src.Fycell.insert({(F1-src.Fy0)/F, i});
		}
	}
}

void source_specie() {
	for (auto& sp : species) {
		for (auto& sp_src : sp.source) {
			Source src = sources[sp_src.index];
			sp_src.vtpar = sqrt(src.Tpar * q_e / sp.mass);
			sp_src.vtper = sqrt(src.Tper * q_e / sp.mass);
		}
	}
}

void source() {
	for (int sp=0; sp<species.size(); sp++) {
		for (int sp_src=0; sp_src<species[sp].source.size(); sp_src++) {
			source_inject(sp, sp_src, species[sp].source[sp_src].index);
		}
	}
}

void source_inject(int sp, int sp_src, int src) {
	Source& source = sources[src];
	SourceSpecie& source_sp = species[sp].source[sp_src];
	int cell_i, Ninj;

	if (!source.active) return;

	else if (source.t_dep)
		if (tstep*dt <= source.t0 || tstep*dt >= source.t3) return;
		else if (tstep*dt <= source.t1) // linear increase
			Ninj = source.n_inj_step * source_sp.source_frac *
				(tstep*dt - source.t0) / (source.t1 - source.t0) + frand();
		else if (tstep*dt <= source.t2) // fully active
			Ninj = source.n_inj_step * source_sp.source_frac + frand();
		else // linear decrease
			Ninj = source.n_inj_step * source_sp.source_frac *
				(source.t3 - tstep*dt) / (source.t3 - source.t2) + frand();

	else
		Ninj = source.n_inj_step * source_sp.source_frac + frand();


	Particle p;
	REAL v, r;
	REAL Fx = source.Fx1 - source.Fx0;
	REAL Fy = source.Fy1 - source.Fy0;
	p.z = 0.0;
	species[sp].n += Ninj;

	for (int i=0; i<Ninj; i++) {

		// x coordinate
		if (!source.shape_x)
			p.x = source.j1 + frand()*source.lx;
		else if (source.shape_x == 1) {
			r = source.Fx0 + frand()*Fx;
			r = asin(2*r - 1) * source.Lx / PI;
			p.x = source.cx + r;
		}
		else { // normal distribution (uniform inside cell)
			auto it_next = source.Fxcell.upper_bound({frand(), 1});
			p.x = it_next->second + frand();
		}

		// y coordinate
		if (!source.shape_y)
			p.y = source.k1 + frand()*source.ly;
		else if (source.shape_y == 1) {
			r = source.Fy0 + frand()*Fy;
			r = asin(2*r - 1) * source.Ly / PI;
			p.y = source.cy + r;
		}
		else { // normal distribution (uniform inside cell)
			auto it_next = source.Fycell.upper_bound({frand(), 1});
			p.y = it_next->second + frand();
		}


		// velocity
		if (source.Tpar == source.Tper) { // no need for basis transform
			p.vx = source_sp.vtpar*normvel_MaxBol();
			p.vy = source_sp.vtpar*normvel_MaxBol();
			p.vz = source_sp.vtpar*normvel_MaxBol();
		}
		else { // transform from B basis
			v = source_sp.vtpar*normvel_MaxBol();
			p.vx = Bb[0][0] * v;
			p.vy = Bb[0][1] * v;
			p.vz = Bb[0][2] * v;
			v = source_sp.vtper*normvel_MaxBol();
			p.vx += Bb[1][0] * v;
			p.vy += Bb[1][1] * v;
			p.vz += Bb[1][2] * v;
			v = source_sp.vtper*normvel_MaxBol();
			p.vx += Bb[2][0] * v;
			p.vy += Bb[2][1] * v;
			p.vz += Bb[2][2] * v;
		}

		// add the particle to a proper cell
		cell_i = int(p.y) + int(p.x) * grid.Ncy;
		cells[cell_i].part[sp].push_back(p);
	}
}

void source_heating() {
	int index, cell_i, cell_i_0;
	for (int sp=0; sp<species.size(); sp++) {
		for (int sp_src=0; sp_src<species[sp].source.size(); sp_src++) {
			index = species[sp].source[sp_src].index;
			Source& source = sources[index];
			if (!source.active) continue;

			cell_i_0 = source.k1 + source.j1 * grid.Ncy;
			for (int i=0; i<source.lx; i++) {
				cell_i = cell_i_0;
				for (int j=0; j<source.ly; j++) {
					heating_cell(cell_i, sp, sp_src, index);
					cell_i += 1;
				}
				cell_i_0 += grid.Ncy;
			}
		}
	}
}

void heating_cell(int cell_i, int sp, int sp_src, int src) {
	Source& source = sources[src];
	std::vector<Particle>& part = cells[cell_i].part[sp];
	SourceSpecie& source_sp = species[sp].source[sp_src];
	REAL v;
	int Nheat, part_size = part.size();

	if (source.t_dep)
		if (tstep*dt <= source.t0 || tstep*dt >= source.t3) return;
		else if (tstep*dt <= source.t1) // linear increase
			Nheat =  part_size * source.n_heat * source_sp.heat_frac *
				(tstep*dt - source.t0) / (source.t1 - source.t0) + frand();
		else if (tstep*dt <= source.t2) // fully active
			Nheat =  part_size * source.n_heat * source_sp.heat_frac + frand();
		else // linear decrease
			Nheat =  part_size * source.n_heat * source_sp.heat_frac *
				(source.t3 - tstep*dt) / (source.t3 - source.t2) + frand();
	else
		Nheat =  part_size * source.n_heat * source_sp.heat_frac + frand();

	for (int i=0; i<Nheat; i++) {
		Particle& p = part[int(part_size*frand())];
		// heating as changing velocity distribution into Maxwellian
		if (source.Tpar == source.Tper) { // no need for basis transform
			p.vx = source_sp.vtpar*normvel_MaxBol();
			p.vy = source_sp.vtpar*normvel_MaxBol();
			p.vz = source_sp.vtpar*normvel_MaxBol();
		}
		else { // transform from B basis
			v = source_sp.vtpar*normvel_MaxBol();
			p.vx = Bb[0][0] * v;
			p.vy = Bb[0][1] * v;
			p.vz = Bb[0][2] * v;
			v = source_sp.vtper*normvel_MaxBol();
			p.vx += Bb[1][0] * v;
			p.vy += Bb[1][1] * v;
			p.vz += Bb[1][2] * v;
			v = source_sp.vtper*normvel_MaxBol();
			p.vx += Bb[2][0] * v;
			p.vy += Bb[2][1] * v;
			p.vz += Bb[2][2] * v;
		}
	}
}
