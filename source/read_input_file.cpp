/**
 * SPDX-License-Identifier: EUPL-1.2
 * SPDX-FileCopyrightText: 2025 Stefan Costea, LeCAD-PEG
 *
 * @file read_input_file.cpp
 * @brief Parser for the input file and runtime configuration (geometry,
 * species, sources, solver parameters).
 */

/* HEXAPIC read input file */

#include "hexapic.hpp"
#include <filesystem>

void process_section(std::string word, std::string ref) {
  static int word_found = FALSE, status = 0, specie_index = -1;
  static std::string name, subname;
  static Species_load *lv;   ///< pointer to load vector element
  static Species_inject *iv; ///< pointer to inject vector element

  if (!status && word == ref) {
    status = 1;
    name = word;
    // std::cout << "CHECK: name=" << name << "\n";
    return;
  }

  if (status == 1 && word == "{") {
    status = 2;
    if (name == "Species") {
      species.push_back({}); ///< zero-initialise
      species.back().source = {};
    }
    if (name == "Grid")
      grid = {}; ///< zero-initialise
    if (name == "Equipotential" || name == "Dielectric") {
      grid.boundaries.push_back({}); ///< zero-initialise
      grid.boundaries.back().name = name;
    }
    if (name == "MCC") {
      if (MCC_init_data.size() == 0)
        MCC_init_data = {};
      MCC_init_data.push_back({});
    }
    if (name == "PSI") {
      if (PSI_init_data.size() == 0)
        PSI_init_data = {};
      PSI_init_data.push_back({});
    }
    if (name == "Source") {
      sources.push_back({});
    }

    return;
  }

  if (status == 2) {
    if (name == "Species") {
      if (word == "name" || word == "m" || word == "q" || word == "subcycle" ||
          word.substr(0, 4) == "heat" || word.substr(0, 6) == "source") {
        word_found = TRUE;
        subname = word;
        return;
      }
      if (word_found && word != "=") {
        if (subname == "name")
          species.back().name = word;
        if (subname == "m")
          species.back().mass = std::stof(word);
        if (subname == "q")
          species.back().charge = int(std::stof(word) / q_e);
        if (subname == "subcycle")
          species.back().subcycle = std::stoi(word);
        if (subname.substr(0, 6) == "source") {
          SourceSpecie source;
          source.index = std::stoi(subname.substr(6));
          source.source_frac = std::stof(word);
          species.back().source.push_back(source);
        }
        if (subname.substr(0, 4) == "heat")
          species.back().source.back().heat_frac = std::stof(word);

        word_found = FALSE;
        subname = "";
      }
      if (word == "}") {
        name = "";
        status = 0;
        return;
      }
    }
    if (name == "Grid") {
      if (word == "J" || word == "x1s" || word == "x1f" || word == "K" ||
          word == "x2s" || word == "x2f" || word == "PeriodicFlagX1" ||
          word == "PeriodicFlagX2") {
        word_found = TRUE;
        subname = word;
        return;
      }
      if (word_found && word != "=") {
        // std::cout << subname << "=" << word << "\n";
        if (subname == "J")
          grid.Ncx = std::stoi(word);
        if (subname == "x1s")
          grid.x1s = std::stof(word);
        if (subname == "x1f")
          grid.x1f = std::stof(word);
        if (subname == "K")
          grid.Ncy = std::stoi(word);
        if (subname == "x2s")
          grid.x2s = std::stof(word);
        if (subname == "x2f")
          grid.x2f = std::stof(word);
        if (subname == "PeriodicFlagX1") {
          int PeriodicFlagX1 = std::stoi(word);
          if (PeriodicFlagX1 == 1) {
            grid.boundaries.push_back({});
            grid.boundaries.back().add_parameters("Periodic", 0, 0, 0,
                                                  grid.Ncy);
            grid.boundaries.push_back({});
            grid.boundaries.back().add_parameters("Periodic", grid.Ncx,
                                                  grid.Ncx, 0, grid.Ncy);
          }
        }
        if (subname == "PeriodicFlagX2") {
          int PeriodicFlagX2 = std::stoi(word);
          if (PeriodicFlagX2 == 1) {
            grid.boundaries.push_back({});
            grid.boundaries.back().add_parameters("Periodic", 0, grid.Ncx, 0,
                                                  0);
            grid.boundaries.push_back({});
            grid.boundaries.back().add_parameters("Periodic", 0, grid.Ncx,
                                                  grid.Ncy, grid.Ncy);
          }
        }

        word_found = FALSE;
        subname = "";
      }
      if (word == "}") {
        grid.dx = (grid.x1f - grid.x1s) / grid.Ncx;
        grid.dy = (grid.x2f - grid.x2s) / grid.Ncy;
        name = "";
        status = 0;
        return;
      }
    }
    if (name == "Control") {
      if (word == "dt" || word == "B01" || word == "B02" || word == "B03" ||
          word == "FieldSolverFlag" || word == "ParticleDiagnosticFlag") {
        word_found = TRUE;
        subname = word;
        return;
      }
      if (word_found && word != "=") {
        // std::cout << subname << "=" << word << "\n";
        if (subname == "dt")
          dt = std::stof(word);
        if (subname == "B01")
          Bf[0] = std::stof(word);
        if (subname == "B02")
          Bf[1] = std::stof(word);
        if (subname == "B03")
          Bf[2] = std::stof(word);
        if (subname == "FieldSolverFlag")
          FieldSolverFlag = std::stoi(word);
        if (subname == "ParticleDiagnosticFlag")
          ParticleDiagnosticFlag = std::stoi(word);
        word_found = FALSE;
        subname = "";
      }
      if (word == "}") {
        name = "";
        status = 0;
        return;
      }
    }
    if (name == "Equipotential" || name == "Dielectric") {
      if (word == "C" || word == "j1" || word == "j2" || word == "k1" ||
          word == "k2" || word == "reflection" || word == "QuseFlag") {
        word_found = TRUE;
        subname = word;
        return;
      }
      if (word_found && word != "=") {
        // std::cout << subname << "=" << word << "\n";
        if (subname == "C")
          grid.boundaries.back().C = std::stof(word);
        if (subname == "QuseFlag")
          grid.boundaries.back().QuseFlag = std::stoi(word);
        if (subname == "reflection")
          grid.boundaries.back().reflection = std::stof(word);
        if (subname == "j1")
          grid.boundaries.back().j1 = std::stoi(word);
        if (subname == "j2")
          grid.boundaries.back().j2 = std::stoi(word);
        if (subname == "k1")
          grid.boundaries.back().k1 = std::stoi(word);
        if (subname == "k2")
          grid.boundaries.back().k2 = std::stoi(word);
        word_found = FALSE;
        subname = "";
      }
      if (word == "}") {
        name = "";
        status = 0;
        return;
      }
    }
    if (name == "Load") {
      if (word == "speciesName" || word == "LoadMethodFlag" ||
          word == "x1MinMKS" || word == "x1MaxMKS" || word == "x2MinMKS" ||
          word == "x2MaxMKS" || word == "density" || word == "np2c" ||
          word == "v1thermal" || word == "v2thermal" || word == "v3thermal" ||
          word == "v1drift" || word == "v2drift" || word == "v3drift") {
        word_found = TRUE;
        subname = word;
        return;
      }
      if (word_found && word != "=") {
        // std::cout << subname << "=" << word << "\n";
        if (subname == "speciesName") { // must be first entry in section
          // find specie index to which loading applies
          for (int i = 0; i < species.size(); i++)
            if (species[i].name == word) {
              specie_index = i;
              break;
            }
          if (specie_index == -1) {
            std::cout << "INPUT Error: Specie <" << word
                      << "> in Load not found in Species"
                      << "\n";
            exit(1);
          }
          species[specie_index].load.push_back({}); ///< zero-initialise
          lv = &species[specie_index].load.back();  ///< assign pointer
        }
        if (subname == "x1MinMKS")
          lv->x1s = std::stof(word);
        if (subname == "x1MaxMKS")
          lv->x1f = std::stof(word);
        if (subname == "x2MinMKS")
          lv->x2s = std::stof(word);
        if (subname == "x2MaxMKS")
          lv->x2f = std::stof(word);
        if (subname == "density")
          lv->density = std::stof(word);
        if (subname == "np2c")
          lv->np2c = (long long)std::stof(word);
        if (subname == "v1thermal")
          lv->v1thermal = std::stof(word);
        if (subname == "v2thermal")
          lv->v2thermal = std::stof(word);
        if (subname == "v3thermal")
          lv->v3thermal = std::stof(word);
        if (subname == "v1drift")
          lv->v1drift = std::stof(word);
        if (subname == "v2drift")
          lv->v2drift = std::stof(word);
        if (subname == "v3drift")
          lv->v3drift = std::stof(word);
        if (subname == "LoadMethodFlag")
          lv->method = std::stoi(word);
        word_found = FALSE;
        subname = "";
      }
      if (word == "}") {
        name = "";
        status = 0;
        return;
      }
    }

    if (name == "EmitPort") {
      if (word == "speciesName" || word == "I" || word == "j1" ||
          word == "j2" || word == "k1" || word == "k2" || word == "np2c" ||
          word == "normal" || word == "v1thermal" || word == "v2thermal" ||
          word == "v3thermal" || word == "v1drift" || word == "v2drift" ||
          word == "v3drift") {
        word_found = TRUE;
        subname = word;
        return;
      }
      if (word_found && word != "=") {
        // std::cout << subname << "=" << word << "\n";
        if (subname == "speciesName") { // must be first entry in section
          // find specie index to which loading applies
          for (int i = 0; i < species.size(); i++)
            if (species[i].name == word) {
              specie_index = i;
              break;
            }
          if (specie_index == -1) {
            std::cout << "INPUT Error: Specie <" << word
                      << "> in EmitPort not found in Species"
                      << "\n";
            exit(1);
          }
          species[specie_index].inject.push_back({}); ///< initialise
          iv = &species[specie_index].inject.back();  ///< assign pointer
        }
        if (subname == "I")
          iv->I = fabs(std::stof(word));
        if (subname == "j1")
          iv->j1 = std::stoi(word);
        if (subname == "j2")
          iv->j2 = std::stoi(word);
        if (subname == "k1")
          iv->k1 = std::stoi(word);
        if (subname == "k2")
          iv->k2 = std::stoi(word);
        if (subname == "normal")
          iv->normal = std::stoi(word);
        if (subname == "np2c")
          iv->np2c = (long long)std::stof(word);
        if (subname == "v1thermal")
          iv->v1thermal = std::stof(word);
        if (subname == "v2thermal")
          iv->v2thermal = std::stof(word);
        if (subname == "v3thermal")
          iv->v3thermal = std::stof(word);
        if (subname == "v1drift")
          iv->v1drift = std::stof(word);
        if (subname == "v2drift")
          iv->v2drift = std::stof(word);
        if (subname == "v3drift")
          iv->v3drift = std::stof(word);
        word_found = FALSE;
        subname = "";
      }
      if (word == "}") {
        name = "";
        status = 0;
        return;
      }
    }
    if (name == "MCC") {
      if (word == "collName" || word == "specie1Name" ||
          word == "specie2Name" || word == "reactionIdx" ||
          word == "collXSectionFile") {
        word_found = TRUE;
        subname = word;
        return;
      }
      if (word_found && word != "=") {
        // std::cout << subname << "=" << word << "\n";
        if (subname == "collName" || subname == "specie1Name" ||
            subname == "specie2Name" || subname == "reactionIdx" ||
            subname == "collXSectionFile")
          MCC_init_data.back().push_back(word);
        word_found = FALSE;
        subname = "";
      }
      if (word == "}") {
        name = "";
        status = 0;
        return;
      }
    }

    if (name == "PSI") {
      if (word == "boundary" || word == "specie" || word == "secondary1" ||
          word == "secondary2" || word == "secondary3" || word == "sec_type" ||
          word == "inj_type" || word == "E0" || word == "Ew" ||
          word == "gamma0" || word == "T0" || word == "Ep") {
        word_found = TRUE;
        subname = word;
        return;
      }
      if (word_found && word != "=") {
        // std::cout << subname << "=" << word << "\n";
        if (subname == "boundary" || subname == "specie" ||
            subname == "secondary1" || subname == "secondary2" ||
            subname == "secondary3" || subname == "sec_type" ||
            subname == "inj_type" || subname == "E0" || subname == "Ew" ||
            subname == "gamma0" || subname == "T0" || subname == "Ep")
          PSI_init_data.back().push_back(word);
        word_found = FALSE;
        subname = "";
      }
      if (word == "}") {
        name = "";
        status = 0;
        return;
      }
    }

    if (name == "Source") {
      if (word == "index" || word == "n_inj" || word == "shape_x" ||
          word == "shape_y" || word == "j1" || word == "j2" || word == "k1" ||
          word == "k2" || word == "Tpar" || word == "Tper" ||
          word == "heat_freq" || word == "t_dep" || word == "t0" ||
          word == "t1" || word == "t2" || word == "t3") {
        word_found = TRUE;
        subname = word;
        return;
      }
      if (word_found && word != "=") {
        if (subname == "index")
          sources.back().index = std::stoi(word);
        if (subname == "n_inj")
          sources.back().n_inj = std::stof(word);
        if (subname == "shape_x")
          sources.back().shape_x = std::stoi(word);
        if (subname == "shape_y")
          sources.back().shape_y = std::stoi(word);
        if (subname == "j1")
          sources.back().j1 = std::stoi(word);
        if (subname == "j2")
          sources.back().j2 = std::stoi(word);
        if (subname == "k1")
          sources.back().k1 = std::stoi(word);
        if (subname == "k2")
          sources.back().k2 = std::stoi(word);
        if (subname == "Tpar")
          sources.back().Tpar = std::stof(word);
        if (subname == "Tper")
          sources.back().Tper = std::stof(word);
        if (subname == "heat_freq")
          sources.back().heat_freq = std::stof(word);
        if (subname == "t_dep")
          sources.back().t_dep = std::stoi(word);
        if (subname == "t0")
          sources.back().t0 = std::stof(word);
        if (subname == "t1")
          sources.back().t1 = std::stof(word);
        if (subname == "t2")
          sources.back().t2 = std::stof(word);
        if (subname == "t3")
          sources.back().t3 = std::stof(word);

        word_found = FALSE;
        subname = "";
      }
      if (word == "}") {
        name = "";
        status = 0;
        return;
      }
    }
  }
}

void read_input_file(char **argv) {
  FieldSolverFlag = 1;
  ParticleDiagnosticFlag = 0;
  QuietFlag = 0;
  // state descriptors: 0=not found, 1=found, 2=reading, 3=finished reading
  int Description = 0, Region = 0, Load = 0;

  // variables to hold read data
  std::vector<std::vector<REAL>> species_data, load_data;
  std::vector<REAL> grid_data, control_data;

  // filestream variable file
  std::fstream file;
  std::string word, filename, title, description;

  // name of the input file
  filename = std::string(argv[1]);

  // opening file
  file.open(filename.c_str());

  // extracting words from the file
  while (file >> word) {
    if (Description < 3) {
      if (!Description) {
        if (word == "{")
          Description = 2;
        else
          title.append(word);
        continue;
      }
      if (Description) {
        if (word == "}")
          Description = 3;
        else {
          description.append(word);
          description.append(" ");
        }
      }
    }

    if (Region < 3 && Description == 3) {
      if (!Region && word == "Region") {
        Region = 1;
        continue;
      }
      if (Region == 1 && word == "{") {
        Region = 2;
        continue;
      }
      if (Region == 2) {
        process_section(word, "Species");
        process_section(word, "Grid");
        process_section(word, "Control");
        process_section(word, "Equipotential");
        process_section(word, "Dielectric");
        process_section(word, "Load");
        process_section(word, "EmitPort");
        process_section(word, "MCC");
        process_section(word, "PSI");
        process_section(word, "Source");

      } // end of if(Region==2)

    } // end of if(Region<3)

  } // end of while

  file.close();

  // Update MCC file paths
  std::filesystem::path exe_dir =
      std::filesystem::canonical("/proc/self/exe").parent_path();
  std::string abs_path = exe_dir.string();
  for (auto &mcc : MCC_init_data) {
    mcc[4] = abs_path + "/../mcc/" + mcc[4];
    // std::cout << "Loaded MCC file: " << mcc[4] << std::endl;
  }

  /* ------- print read parameters ----------------------*/

  if (!rmpi) {
    // description
    std::cout << "INPUT: "
              << "title = " << title << "\n";
    std::cout << "INPUT: "
              << "description =" << description << "\n";

    // species
    for (auto sp : species) {
      printf("INPUT: Specie:\n");
      printf("INPUT: name=%s\n", sp.name.c_str());
      printf("INPUT: mass=%e\n", sp.mass);
      printf("INPUT: charge=%d\n", sp.charge);
      printf("INPUT: subcycle=%d\n", sp.subcycle);
      for (auto src : sp.source) {
        std::cout << "INPUT: source_i=" << src.index << "\n";
        std::cout << "INPUT: source_frac=" << src.source_frac << "\n";
        std::cout << "INPUT: heat_frac=" << src.heat_frac << "\n";
      }
    }

    // grid
    std::cout << "INPUT: Ncx = " << grid.Ncx << "\n";
    std::cout << "INPUT: x1s = " << grid.x1s << "\n";
    std::cout << "INPUT: x1f = " << grid.x1f << "\n";
    std::cout << "INPUT: dx = " << grid.dx << "\n";
    std::cout << "INPUT: Ncy = " << grid.Ncy << "\n";
    std::cout << "INPUT: x2s = " << grid.x2s << "\n";
    std::cout << "INPUT: x2f = " << grid.x2f << "\n";
    std::cout << "INPUT: dy = " << grid.dy << "\n";

    // control
    std::cout << "INPUT: dt = " << dt << "\n";
    std::cout << "INPUT: Bx = " << Bf[0] << "\n";
    std::cout << "INPUT: By = " << Bf[1] << "\n";
    std::cout << "INPUT: Bz = " << Bf[2] << "\n";

    // boundaries
    for (auto bc : grid.boundaries) {
      std::cout << "INPUT: " << bc.name << "\n";
      if (bc.name == "Equipotential") {
        std::cout << "INPUT: C = " << bc.C << "\n";
      }
      if (bc.name == "Dielectric") {
        std::cout << "INPUT: QuseFlag = " << bc.QuseFlag << "\n";
        std::cout << "INPUT: reflection = " << bc.reflection << "\n";
      }
      std::cout << "INPUT: j1 = " << bc.j1 << "\n";
      std::cout << "INPUT: j2 = " << bc.j2 << "\n";
      std::cout << "INPUT: k1 = " << bc.k1 << "\n";
      std::cout << "INPUT: k2 = " << bc.k2 << "\n";
    }

    // emit port
    for (auto inj : species[0].inject) {
      std::cout << "INPUT: EmitPort \n";
      std::cout << "INPUT: I = " << inj.I << "\n";
      std::cout << "INPUT: j1 = " << inj.j1 << "\n";
      std::cout << "INPUT: j2 = " << inj.j2 << "\n";
      std::cout << "INPUT: k1 = " << inj.k1 << "\n";
      std::cout << "INPUT: k2 = " << inj.k2 << "\n";
      std::cout << "INPUT: normal = " << inj.normal << "\n";
    }

    // Sources
    for (auto src : sources) {
      std::cout << "INPUT: Source \n";
      std::cout << "INPUT: index = " << src.index << "\n";
      std::cout << "INPUT: n_inj = " << src.n_inj << "\n";
      std::cout << "INPUT: shape_x = " << src.shape_x << "\n";
      std::cout << "INPUT: shape_y = " << src.shape_y << "\n";
      std::cout << "INPUT: j1 = " << src.j1 << "\n";
      std::cout << "INPUT: j2 = " << src.j2 << "\n";
      std::cout << "INPUT: k1 = " << src.k1 << "\n";
      std::cout << "INPUT: k2 = " << src.k2 << "\n";
      std::cout << "INPUT: Tpar = " << src.Tpar << "\n";
      std::cout << "INPUT: Tper = " << src.Tper << "\n";
      std::cout << "INPUT: heat_freq = " << src.heat_freq << "\n";
      std::cout << "INPUT: t_dep = " << src.t_dep << "\n";
      if (src.t_dep) {
        std::cout << "INPUT: t0 = " << src.t0 << "\n";
        std::cout << "INPUT: t1 = " << src.t1 << "\n";
        std::cout << "INPUT: t2 = " << src.t2 << "\n";
        std::cout << "INPUT: t3 = " << src.t3 << "\n";
      }
    }

    // MCC
    for (auto mcc : MCC_init_data) {
      std::cout << "INPUT: MCC \n";
      std::cout << "INPUT: collName = " << mcc[0] << "\n";
      std::cout << "INPUT: specie1Name = " << mcc[1] << "\n";
      std::cout << "INPUT: specie2Name = " << mcc[2] << "\n";
      std::cout << "INPUT: reactionIdx = " << mcc[3] << "\n";
      std::cout << "INPUT: collXSectionFile = " << mcc[4] << "\n";
    }

    // PSI
    for (auto psi : PSI_init_data) {
      int n = psi.size() - 9;
      std::cout << "INPUT: PSI \n";
      std::cout << "INPUT: boundary = " << psi[0] << "\n";
      std::cout << "INPUT: specie = " << psi[1] << "\n";
      for (int i = 1; i <= n; i++) {
        std::cout << "INPUT: secondary" << i << " = " << psi[1 + i] << "\n";
      }
      n = psi.size();
      std::cout << "INPUT: sec_type = " << psi[n - 7] << "\n";
      std::cout << "INPUT: inj_type = " << psi[n - 6] << "\n";
      std::cout << "INPUT: Ew = " << psi[n - 5] << "\n";
      std::cout << "INPUT: E0 = " << psi[n - 4] << "\n";
      std::cout << "INPUT: gamma0 = " << psi[n - 3] << "\n";
      std::cout << "INPUT: T0 = " << psi[n - 2] << "\n";
      std::cout << "INPUT: Ep = " << psi[n - 1] << "\n";
    }
  }

  /* ------ end of print read parameters --------------------*/
}
