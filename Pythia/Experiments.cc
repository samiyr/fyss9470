#ifndef EXPERIMENTS_H
#define EXPERIMENTS_H

#include "ExperimentDefs.cc"

int THREAD_COUNT = 16;

int main() {
	// run_pT_experiment();
	// pT_cross_section(100'000'000 / 5, true, true, "MPI/Bias subdivision/1e8/data.csv", -1, true);

	// pT_cross_section(100'000 / 5, true, true, "No MPI/Bias subdivision/1e5/data", false);
	// pT_cross_section(1'000'000 / 5, true, true, "No MPI/Bias subdivision/1e6/data.csv", false);

	// pp_mpi_run(1'000'000, Process::HardQCD, MPIStrategy::Disabled, 1.5, "Data/pp/delta_phi/No MPI/Hard15/1e6/");
	// pp_mpi_run(1'000'000, Process::SoftQCDNonDiffractive, MPIStrategy::Disabled, 1.5, "Data/pp/delta_phi/No MPI/Soft/1e6/");

	// pp_mpi_run(10'000'000, Process::HardQCD, MPIStrategy::Disabled, 1.0, "Data/pp/delta_phi/No MPI/Hard10/1e7/");
	// pp_mpi_run(10'000'000, Process::HardQCD, MPIStrategy::Disabled, 1.5, "Data/pp/delta_phi/No MPI/Hard15/1e7/");
	// pp_mpi_run(10'000'000, Process::HardQCD, MPIStrategy::Disabled, 2.0, "Data/pp/delta_phi/No MPI/Hard20/1e7/");
	// pp_mpi_run(10'000'000, Process::SoftQCDNonDiffractive, MPIStrategy::Disabled, 1.5, "Data/pp/delta_phi/No MPI/Soft/1e7/");

	// pp_mpi_run(100'000'000, Process::HardQCD, MPIStrategy::Disabled, 1.5, "Data/pp/delta_phi/No MPI/Hard15/1e8/");
	// pp_mpi_run(100'000'000, Process::SoftQCDNonDiffractive, MPIStrategy::Disabled, 1.5, "Data/pp/delta_phi/No MPI/Soft/1e8/");

	// pp_mpi_run(1'000'000'000, Process::HardQCD, MPIStrategy::Disabled, 1.0, "Data/pp/delta_phi/No MPI/Hard10/1e9/");
	// pp_mpi_run(1'000'000'000, Process::HardQCD, MPIStrategy::Disabled, 1.5, "Data/pp/delta_phi/No MPI/Hard15/1e9/");
	// pp_mpi_run(1'000'000'000, Process::HardQCD, MPIStrategy::Disabled, 2.0, "Data/pp/delta_phi/No MPI/Hard20/1e9/");
	// pp_mpi_run(1'000'000'000, Process::SoftQCDNonDiffractive, MPIStrategy::Disabled, 1.5, "Data/pp/delta_phi/No MPI/Soft/1e9/");

	// pp_dps_run(10'000'000, Process::HardQCD, MPIStrategy::DPS, 1.5, "Data/pp/delta_phi/DPS/Hard15/1e7/");
	// pp_dps_run(10'000'000, Process::SoftQCDNonDiffractive, MPIStrategy::DPS, 1.5, "Data/pp/delta_phi/DPS/Soft/1e7/");

	// pp_mpi_run(10'000'000, Process::HardQCD, MPIStrategy::PythiaMPI, 1.5, "Data/pp/delta_phi/MPI/Hard15/1e7/");
	// pp_mpi_run(10'000'000, Process::SoftQCDNonDiffractive, MPIStrategy::PythiaMPI, 1.5, "Data/pp/delta_phi/MPI/Soft/1e7/");

	// pp_dps_run(100'000'000, Process::HardQCD, MPIStrategy::DPS, 1.5, "Data/pp/delta_phi/DPS/Hard15/1e8/");
	// pp_dps_run(100'000'000, Process::SoftQCDNonDiffractive, MPIStrategy::DPS, 1.5, "Data/pp/delta_phi/DPS/Soft/1e8/");

	// pp_dps_run(1'000'000'000, Process::HardQCD, MPIStrategy::DPS, 1.5, "Data/pp/delta_phi/DPS/Hard15/1e9/");
	// pp_dps_run(1'000'000'000, Process::SoftQCDNonDiffractive, MPIStrategy::DPS, 1.5, "Data/pp/delta_phi/DPS/Soft/1e9/");

	// Al_run(100'000'000, Process::HardQCD, MPIStrategy::DPS, 1.5, "Data/pAl/delta_phi/DPS/Hard15/1e8/", false);
	// Al_run(100'000'000, Process::HardQCD, MPIStrategy::DPS, 1.5, "Data/pAl/delta_phi/DPS/Hard15 nPDF/1e8/", true);
	// Al_run(100'000'000, Process::SoftQCDNonDiffractive, MPIStrategy::DPS, 1.5, "Data/pAl/delta_phi/DPS/Soft/1e8/", false);

	// Au_run(100'000'000, Process::HardQCD, MPIStrategy::DPS, 1.5, "Data/pAu/delta_phi/DPS/Hard15/1e8/", false);
	// Au_run(100'000'000, Process::HardQCD, MPIStrategy::DPS, 1.5, "Data/pAu/delta_phi/DPS/Hard15 nPDF/1e8/", true);
	// Au_run(100'000'000, Process::SoftQCDNonDiffractive, MPIStrategy::DPS, 1.5, "Data/pAu/delta_phi/DPS/Soft/1e8/", false);

	// Al_run(1'000'000'000, Process::HardQCD, MPIStrategy::DPS, 1.5, "Data/pAl/delta_phi/DPS/Hard15/1e9/", false);
	// Al_run(1'000'000'000, Process::HardQCD, MPIStrategy::DPS, 1.5, "Data/pAl/delta_phi/DPS/Hard15 nPDF/1e9/", true);
	// Al_run(1'000'000'000, Process::SoftQCDNonDiffractive, MPIStrategy::DPS, 1.5, "Data/pAl/delta_phi/DPS/Soft/1e9/", false);

	// Au_run(1'000'000'000, Process::HardQCD, MPIStrategy::DPS, 1.5, "Data/pAu/delta_phi/DPS/Hard15/1e9/", false);
	// Au_run(1'000'000'000, Process::HardQCD, MPIStrategy::DPS, 1.5, "Data/pAu/delta_phi/DPS/Hard15 nPDF/1e9/", true);
	// Au_run(1'000'000'000, Process::SoftQCDNonDiffractive, MPIStrategy::DPS, 1.5, "Data/pAu/delta_phi/DPS/Soft/1e9/", false);

	// pp_dps_run(10'000'000'000, Process::SoftQCDNonDiffractive, MPIStrategy::DPS, 1.5, "Data/pp/delta_phi/DPS/Soft/1e10/");
	// Al_run(10'000'000'000, Process::SoftQCDNonDiffractive, MPIStrategy::DPS, 1.5, "Data/pAl/delta_phi/DPS/Soft/1e10/", false);
	// Au_run(10'000'000'000, Process::SoftQCDNonDiffractive, MPIStrategy::DPS, 1.5, "Data/pAu/delta_phi/DPS/Soft/1e10/", false);

	// pp_dps_run(50'000'000'000, Process::SoftQCDNonDiffractive, MPIStrategy::DPS, 1.5, "Data/pp/delta_phi/DPS/Soft/5e10/");
	// Al_run(50'000'000'000, Process::SoftQCDNonDiffractive, MPIStrategy::DPS, 1.5, "Data/pAl/delta_phi/DPS/Soft/5e10/", false);
	// Au_run(50'000'000'000, Process::SoftQCDNonDiffractive, MPIStrategy::DPS, 1.5, "Data/pAu/delta_phi/DPS/Soft/5e10/", false);

	// pp_dps_run(10'000'000, Process::SoftQCDNonDiffractive, MPIStrategy::DPS, 1.5, "Data/pp/delta_phi/DPS/Soft/1e7/STARC/", Normalization::STARC);
	// Al_run(10'000'000, Process::SoftQCDNonDiffractive, MPIStrategy::DPS, 1.5, "Data/pAl/delta_phi/DPS/Soft/1e7/STARC/", false, Normalization::STARC);
	// Au_run(10'000'000, Process::SoftQCDNonDiffractive, MPIStrategy::DPS, 1.5, "Data/pAu/delta_phi/DPS/Soft/1e7/STARC/", false, Normalization::STARC);

	// pp_dps_run(2'000'000'000, Process::HardQCD, MPIStrategy::DPS, 1.5, "Data/pp/delta_phi/DPS/Hard15/2e9/STARC/", Normalization::STARC);
	// pp_dps_run(2'000'000'000, Process::HardQCD, MPIStrategy::DPS, 2.5, "Data/pp/delta_phi/DPS/Hard25/2e9/STARC/", Normalization::STARC);

	// Al_run(2'000'000'000, Process::HardQCD, MPIStrategy::DPS, 1.5, "Data/pAl/delta_phi/DPS/Hard15/2e9/STARC/", false, Normalization::STARC);
	// Al_run(2'000'000'000, Process::HardQCD, MPIStrategy::DPS, 1.5, "Data/pAl/delta_phi/DPS/Hard15 nPDF/2e9/STARC/", true, Normalization::STARC);

	// Al_run(2'000'000'000, Process::HardQCD, MPIStrategy::DPS, 2.5, "Data/pAl/delta_phi/DPS/Hard25/2e9/STARC/", false, Normalization::STARC);
	// Al_run(2'000'000'000, Process::HardQCD, MPIStrategy::DPS, 2.5, "Data/pAl/delta_phi/DPS/Hard25 nPDF/2e9/STARC/", true, Normalization::STARC);

	// Au_run(2'000'000'000, Process::HardQCD, MPIStrategy::DPS, 1.5, "Data/pAu/delta_phi/DPS/Hard15/2e9/STARC/", false, Normalization::STARC);
	// Au_run(2'000'000'000, Process::HardQCD, MPIStrategy::DPS, 1.5, "Data/pAu/delta_phi/DPS/Hard15 nPDF/2e9/STARC/", true, Normalization::STARC);

	// Au_run(2'000'000'000, Process::HardQCD, MPIStrategy::DPS, 2.5, "Data/pAu/delta_phi/DPS/Hard25/2e9/STARC/", false, Normalization::STARC);
	// Au_run(2'000'000'000, Process::HardQCD, MPIStrategy::DPS, 2.5, "Data/pAu/delta_phi/DPS/Hard25 nPDF/2e9/STARC/", true, Normalization::STARC);

	// pp_dps_run(10'000'000, Process::HardQCD, MPIStrategy::DPS, 1.5, "Data/temp/", Normalization::STARC);
	// Al_run(10'000'000, Process::HardQCD, MPIStrategy::DPS, 1.5, "Data/temp/pAl/", true, Normalization::STARC);
	// Au_run(10'000'000, Process::HardQCD, MPIStrategy::DPS, 1.5, "Data/temp/pAu/", true, Normalization::STARC);

	// ppAlAu_run(10'000'000, Process::HardQCD, MPIStrategy::DPS, 1.5, "Data/temp/pp5/", "Data/temp/pAl5/", "Data/temp/pAu5/", Normalization::STARC);

	// ppAlAu_run(10'000'000'000, Process::HardQCD, MPIStrategy::DPS, 1.5, "Data/pp/delta_phi/DPS/Hard15/1e10/STARCN/", "Data/pAl/delta_phi/DPS/Hard15/1e10/STARCN/", "Data/pAu/delta_phi/DPS/Hard15/1e10/STARCN/", Normalization::STARC);
	// ppAlAu_run(10'000'000'000, Process::HardQCD, MPIStrategy::DPS, 2.5, "Data/pp/delta_phi/DPS/Hard25/1e10/STARCN/", "Data/pAl/delta_phi/DPS/Hard25/1e10/STARCN/", "Data/pAu/delta_phi/DPS/Hard25/1e10/STARCN/", Normalization::STARC);

	// pp_mpi_run(1'000'000, Process::SoftQCDNonDiffractive, MPIStrategy::PythiaMPI, 1.0, "Data/temp/pp6/", Normalization::STARC);
	// pp_ncoll_run(1'000'000, Process::SoftQCDNonDiffractive, 1.0, "Data/temp/pp7");
	// Al_ncoll_run(1'000'000, Process::SoftQCDNonDiffractive, 1.0, "Data/temp/pAl7", 2);
	// Au_ncoll_run(1'000'000, Process::SoftQCDNonDiffractive, 1.0, "Data/temp/pAu7", 0);

	pp_ncoll_run(100'000'000, Process::SoftQCDNonDiffractive, 1.0, "Data/pp/delta_phi/ncoll0_offset/Soft/1e8/", 0, 1000);
	// pp_ncoll_run(100'000'000, Process::SoftQCDNonDiffractive, 1.0, "Data/pp/delta_phi/ncoll2/Soft/1e8/", 2);
	// pp_ncoll_run(100'000'000, Process::SoftQCDNonDiffractive, 1.0, "Data/pp/delta_phi/ncoll5/Soft/1e8/", 5);

	Al_ncoll_run(100'000'000, Process::SoftQCDNonDiffractive, 1.0, "Data/pAl/delta_phi/ncoll0_offset/Soft/1e8/", 0, 1000);
	Al_ncoll_run(100'000'000, Process::SoftQCDNonDiffractive, 1.0, "Data/pAl/delta_phi/ncoll2_offset/Soft/1e8/", 2, 1000);
	Al_ncoll_run(100'000'000, Process::SoftQCDNonDiffractive, 1.0, "Data/pAl/delta_phi/ncoll5_offset/Soft/1e8/", 5, 1000);

	Au_ncoll_run(100'000'000, Process::SoftQCDNonDiffractive, 1.0, "Data/pAu/delta_phi/ncoll0_offset/Soft/1e8/", 0, 1000);
	Au_ncoll_run(100'000'000, Process::SoftQCDNonDiffractive, 1.0, "Data/pAu/delta_phi/ncoll2_offset/Soft/1e8/", 2, 1000);
	Au_ncoll_run(100'000'000, Process::SoftQCDNonDiffractive, 1.0, "Data/pAu/delta_phi/ncoll5_offset/Soft/1e8/", 5, 1000);

	return 0;
}

#endif // EXPERIMENTS_H