/**
 * SPDX-License-Identifier: EUPL-1.2
 * SPDX-FileCopyrightText: 2025 Stefan Costea, LeCAD-PEG
 * 
 * @file main.cpp
 * @brief Program entry point: initializes MPI/PETSc, reads input, sets up the simulation, runs the time loop, and finalizes.
 */

/* HEXAPIC main file */

#define EXTERNALS
#include "hexapic.hpp"

int main(int argc, char *argv[]) {

	double start, end;
	unsigned long int Ntot, Nlocal=0;

	int steps=100; // total simulation steps
	dsteps = 1; 		// steps between diagnostics
	davg = 0;			// steps for averaging diagnostics

	MPI_Init(&argc, &argv); ///< initialization of MPI

    MPI_Comm_size(MPI_COMM_WORLD, &nmpi); ///< nmpi is the total number of CPUs
    MPI_Comm_rank(MPI_COMM_WORLD, &rmpi); ///< rmpi is the rank of current CPU

	read_input_file(argv);
	decompose_domain();
	//num_param_init(); ///< not compatible (yet) with domain decomposition
	grid_init(argc, argv);
	if(!rmpi)	std::cout << "HEXAPIC: Nodes initialised \n";
	cells_init();
	if(!rmpi)	std::cout << "HEXAPIC: Cells initialised \n";
	hypre_init();
	if(!rmpi)	std::cout << "HEXAPIC: Hypre solver initialised \n";
	particles_init();
	if(!rmpi)	std::cout << "HEXAPIC: Particles initialised \n";
	source_init();
	if(!rmpi)	std::cout << "HEXAPIC: Souces initialised \n";
	collisions_init();
	if(!rmpi)	std::cout << "HEXAPIC: MCC initialised \n";
	PSI_init();
	if(!rmpi)	std::cout << "HEXAPIC: PSI initialised \n";

	MPI_Barrier(MPI_COMM_WORLD); // otherwise starting phi is not always OK

	if(!rmpi)	std::cout << "HEXAPIC: Loading initial particles... \n";
	initial_particle_load();
	if(!rmpi)	std::cout << "HEXAPIC: Loading initial particles finished.\n";

	for (int i = 1; i < argc - 1; ++i) {
        std::string arg = argv[i];
        if (arg == "-steps") steps = std::atoi(argv[i + 1]);
        else if (arg == "-dsteps") dsteps = std::atoi(argv[i + 1]);
        else if (arg == "-davg") davg = std::atoi(argv[i + 1]);
    }
    if (davg > dsteps) davg = dsteps;

	// // open file for streaming
	// std::string file_stream = std::string(argv[1])+".sst";
    // series = Series(file_stream, Access::CREATE, MPI_COMM_WORLD,
    //         R"(
    //         {
    //           "iteration_encoding": "variable_based",
    //           "adios2": {
    //             "engine": {
    //               "parameters": {
    //                 "BufferGrowthFactor": "2.0",
    //                 "DataTransport": "WAN",
    //                 "RendezvousReaderCount": "0",
    //                 "QueueLimit": "2",
    //                 "QueueFullPolicy": "Discard",
    //                 "use_steps": "True",
    //                 "modifiable_attributes": "True"
    //                 }
    //             }
    //           }
    //         })");
    //     // set constant global attributes
    //     // series.setAttribute("var", var);

    // open file for writing to disk
    std::string file_write = std::string(argv[1])+".bp4";
    std::string const config = R"(
            {
              "iteration_encoding": "variable_based",
              "adios2": {
                "engine": {
                   "parameters": {
                        "use_steps2": "True",
                        "modifiable_attributes": "True"
                        }
                    }
                } 
             })";
    series = Series(file_write, Access::CREATE, MPI_COMM_WORLD, config);

    series.setAttribute("nsp", species.size());
    series.setAttribute("Nx", grid.Nxg);
    series.setAttribute("Ny", grid.Nyg);

	for(int i=0; i<10; i++)
		et[i] = 0.0; ///< profiling

	/* Main loop */

	for(tstep=0; tstep<steps; tstep++) {
		Nlocal = 0;
		for (auto sp : species)
			Nlocal += sp.n;
		MPI_Reduce(&Nlocal, &Ntot, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
		if(!rmpi)
			std::cout << "RUN: tstep=" << tstep << 
				", Nparticles=" << Ntot <<"\n";

		start = MPI_Wtime();
		part2grid_edge();
		part2grid_interior();
		end = MPI_Wtime();
		et[0] += end-start;

		start = MPI_Wtime();
		hypre_source_update();
		end = MPI_Wtime();
		et[1] += end-start;

		start = MPI_Wtime();
		hypre_field_solver();
		end = MPI_Wtime();
		et[2] += end-start;

		start = MPI_Wtime();
		hypre_get_V();
		end = MPI_Wtime();
		et[1] += end-start;

		start = MPI_Wtime();
		particle_mover_boris(cells_interior);
		end = MPI_Wtime();
		et[3] += end-start;

		start = MPI_Wtime();
		hypre_finish_V_halo();
		end = MPI_Wtime();
		et[1] += end-start;

		start = MPI_Wtime();
		particle_mover_boris(cells_edge);
		end = MPI_Wtime();
		et[3] += end-start;

		start = MPI_Wtime();
		inject_particles();
		source();	
		end = MPI_Wtime();
		et[5] += end-start;

		start = MPI_Wtime();
		particle_boundaries();
		end = MPI_Wtime();
		et[4] += end-start;

		start = MPI_Wtime();
		MCC();
		source_heating();
		end = MPI_Wtime();
		et[6] += end-start;

		if(tstep%dsteps==0) { // diagnostics
			MPI_Barrier(MPI_COMM_WORLD); ///< otherwise diagnostics can have problems
			start = MPI_Wtime();
			create_output();
			save_output();
			series.flush();
			end = MPI_Wtime();
			et[7] += end-start;
			}

		}
	/* end of main loop */

	MPI_Barrier(MPI_COMM_WORLD); ///< synchronization

	series.close(); 
	plot_profiling();

	hypre_cleanup();

	MPI_Barrier(MPI_COMM_WORLD); ///< synchronization
	if(!rmpi)	std::cout << "HEXAPIC: Finished. \n";
	MPI_Finalize();

	return 0;
}


void plot_profiling() {
	FILE *f1;
	if(rmpi)    return; ///< only rank 0 is allowed to create output file

    // Open file for writing
    
	f1 = fopen("profiling.txt","w");
	fprintf(f1, "nmpi\tpart2grid()\thypre_source_update()\t"); 
	fprintf(f1, "hypre_field_solver()\t"); 
	fprintf(f1, "particle_mover_boris()\tparticle_boundaries()\t");
	fprintf(f1, "inject_particles()\tMCC+heat\tdiagnostics\ttotal\n");


    for(int i=0; i<8; i++)
			et[8] += et[i]; ///< total time
    
    fprintf(f1, "%d\t", nmpi);
    for(int i=0;i<9;i++)
    	fprintf(f1, "%f\t", et[i]);
    fprintf(f1, "\n");

    fclose(f1);

}
