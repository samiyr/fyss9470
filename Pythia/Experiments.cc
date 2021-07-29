#ifndef EXPERIMENT_H
#define EXPERIMENT_H

#include "ExperimentDefs.cc"

int THREAD_COUNT = 16;

int main() {
	// run_pT_experiment();
	// run_pT_error_experiment();
	// run_dps_error_experiment();
	// run_dps_mpi_experiment(10'000'000);
	// run_hard_soft_npdf_experiment(100'000'000);
	// run_experimental_dps_error_experiment();
	run_nuclear_experiment(1'000'000);

	return 0;
}

#endif // EXPERIMENT_H