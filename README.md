# Azimuthal correlations with Pythia

## General

This repository includes all the code used in generating the data for my research training thesis on pion azimuthal correlations. The code was written during the summer of 2021 while I was doing a summer traineeship in the Department of Physics in the University of Jyväskylä.

## Installation

### Compiling

Using this code requires a `C++17` compiler with `OpenMP` support, such as a modern version of `GCC`, and Pythia. See [this](https://pythia.org) for installations instruction for Pythia. This code has been written with Pythia 8.306. Once Pythia is installed, the code can be compiled. On macOS, the code can be compiled like this:

```
g++ Experiments.cc -Xpreprocessor -fopenmp -o Experiments -O3 -I /path/to/pythia8306/include -std=c++17 -pedantic -W -Wall -Wshadow -fPIC -L /path/to/pythia8306/lib -Wl,-rpath -Wl,/path/to/pythia8306/lib -lpythia8 -ldl -lomp
```

On a RedHat system, the compilation call looks more like this:

```
g++ Experiments.cc -fopenmp -o Experiments -O3 -I /path/to/pythia8306/include -std=c++17 -pedantic -W -Wall -Wshadow -fPIC -L /path/to/pythia8306/lib -Wl,-rpath -Wl,/path/to/pythia8306/lib -lpythia8 -ldl -lgomp
```

In both cases, `/path/to/pythia8306` should of course be replaced with the actual path to the Pythia installation.

### nPDFs

The code includes support for nuclear modifications to the parton distribution functions. Pythia has built-in support for a few sets of nPDFs. For documentation, see the file _PDF Selection_ under _Beams_ in the Pythia online manual, and look for the section _Nuclear modifications of parton densities_. The grid files of these sets are not built into Pythia, however. Pythia's documentation includes links for obtaining these grid files. The grid files, which for EPPS16 have filenames `EPPS16NLOR_A`, where `A` is the mass number of the nucleus. The required grid files must be copied to the folder `pythia8306/share/Pythia8/pdfdata`. Recompiling Pythia is not required.

### Running

After compiling, all that's left is to run the executable created. Continuing the example above, this amounts to `./Experiments`. A useful thing to do is call `time ./Experiments`, which will tell how long the execution took. When the program is running, status updates are provided to `stdout` on the number of events generated thus far. A useful trick is to invoke `screen` before running the program. With `screen` enabled, it's now possible to exit the current screen by pressing <kbd>ctrl</kbd>+<kbd>a</kbd> and then <kbd>d</kbd>. The program is now running in the background, and the status can be checked by returning to the screen with `screen -r`. For more information on using `screen`, please consult Google.


## Program structure

### `Experiment`

At the top level, an `Experiment` class encapsulates the idea of running an experiment. The class contains a number of parameters, such as the com energy and event count, and the function `run()` for actually running the experiment. This is an abstract class, and must be subclassed. 

There are two concrete subclasses of `Experiment`: `TransverseMomentumExperiment` and `AzimuthalCorrelationExperiment`. The file `ExperimentDefs.cc` contains templates for running these experiments, while the actual runs are defined in the main file `Experiments.cc`. 

### `TransverseMomentumExperiment`

This class is intended for studying the transverse momenta of pions produced in proton-proton collisions. It generates a histogram of transverse momenta normalized to the cross section `E * d^3 σ / d^3 p`.

### `AzimuthalCorrelationExperiment`

This class is for the main area of study, azimuthal correlations in proton-proton and proton-nucleus collisions. It pairs generated pions, calculates the difference in their azimuthal correlations and fills a histogram. It supports using both Pythia's MPI model, or the DPS model described in the research training thesis. Additionally, it supports two sets of disjoint transverse momentum and rapidity filters.

Internally, `AzimuthalCorrelationExperiment` loops over (in parallel) a list of phase space cuts. For each cut, it creates an instance of `EventGenerator`, which is responsible for generating and analyzing the events. `EventGenerator` returns a list of results, one for each run specified in `AzimuthalCorrelationExperiment`. Thus, the total output is a 2D list, whose elements are lists of results and correspond to a specific phase space cut. Results across phase space cuts are then combined into a 1D list, representing the results of each run.

### `EventGenerator`

For the azimuthal correlation experiment, the `EventGenerator` class is used and it abstracts the event generation and analysis. This is used to parallelize the event generation by the use of phase space cuts. The `EventGenerator` class is initialized with general parameters (`GeneratorParameters`), a phase space cut as well as a list of run parameters (`Analyzer::Parameters`). The list of run parameters are used to construct a list of `Analyzer` objects. 

Internally, `EventGenerator` initializes an instance of `ParticleGenerator` and passes to it the generator parameters and the phase space cut. `ParticleGenerator` is responsible for the event generation loop, and simply outputs lists of particles. `EventGenerator` takes these lists and hands them over to every `Analyzer` instance. These instances keep track of the results of each run (a run corresponds with a single `Analyzer` instance). After the desired amount of events have been generated, the results are collected from each `Analyzer` instance. Then `EventGenerator` returns a list of `EventGenerator::Result` objects, which represent the results of individual runs.

### `Analyzer`

The `Analyzer` class is responsible for analyzing events generated by `EventGenerator` and keeping track of the results. A single `EventGenerator` instance has a list of `Analyzer` instances, created from a list of analyzer parameters (`Analyzer::Parameters`). All the analyzer instances receive the same events. Phase space cuts, for example, generate different events based on the cut. This means that for each phase space cut, the events must be generated anew. On the other hand, filtering generated particles based on their transverse momentum or rapidity can be done on the same set of events. Runs are intended for applying different analysis procedures to the same underlaying events. 

During event generation, `Analyzer` is fed events through the function `Analyzer::book()`, which takes in a pointer to a list of particles. This input corresponds to the pions found in a single event. This particle list is then looped over and the pion pairs are formed. The pairing procedure is described in detail in Section 4 of the thesis. For each pair, the difference of azimuthal angles is calculated and a histogram is populated. 

### `ParticleGenerator`

This class is the most basic class. It is responsible for the event loop and interfacing with Pythia. `ParticleGenerator` is initialized with a set of parameters (`GeneratorParameters`) and a phase space cut. `ParticleGenerator` then initalized an instance of `Pythia` and passes the relevant settings. The event loop is started by calling `ParticleGenerator::generate()`. During the event loop, `ParticleGenerator` extracts the pions from each event using a `ParticleFilter`, and then calls a lambda expression with these pions. The particles, originally of type `Pythia8::Particle`, are mapped to a custom type `ParticleContainer`. 

Alternatively, `ParticleGenerator` can collect all the particles, event-by-event, to a 2D list and return this.

### `ParticleFilter`

This class provides filtration of particles. A particle, given either as a type `Pythia8::Particle` or `ParticleContainer`, is checked against a set of filters. The default filters are an ID filter (to check whether a particle is a neutral pion), a decay check (whether the particle has decayed and whether decayed particles should be kept), and a transverse momentum and rapidity range filters. These filters are implemented as functions that take in a single input of type `ParticleContainer` and have return types of `bool`.

### `ParticleContainer`

This class is a simple container for particles. It serves two purposes. First, it stores only the necessary information, reducing memory footprint. When constructing a `ParticleContainer`, a `Pythia8::Particle` instance is given as input, and the constructor extracts the necessary properties and stores them. Second, it keeps track of the event weight associated with the particle. Strictly speaking, event weights, as the name suggests, are associated with events. However, it's easier to associate the weights with individual particles, even if all the particles from the same event have the same weight. See the section _Biasing_ for more.

### Biasing

The physics of biasing are explained in the thesis. From the programming point of view, each event comes with a specific weight. The weight is then stored in individual particles. The weight of a particle/event must be taken into account when filling histograms. However, the support for biasing is incomplete. As is explained in the section _Histograms_, the relative error of a histogram bin is given by the hit count of that bin. If the events are weighted, the simple hit count formula no longer works. However, for the purposes of the thesis, this functionality is not needed and thus it has not been implemented.

### Histograms

There are two types of histograms: `Histogram` and `ValueHistogram`. Both are template classes, although some functionality is limited to histograms with types `double`.

The `Histogram` class implements a non-standard histogram. This histograms stores the values it has been filled with, instead of just the number of values. This is used in calculating the cross section for the transverse momentum experiment, and is explained in the thesis. This class is not used anywhere else. A `Histogram` class can be converted to a `ValueHistogram` through the function `Histogram::export_to_values()`. 

`ValueHistogram` implements a more traditional histogram. It's comprised of a list of instances of `RangedContainer`. A `RangedContainer` has a lower and upper bound (bin edges), a value and a hit count. The histogram is filled through the function `ValueHistogram::fill()`. This function loops over the list of `RangedContainer` and tries to find a bin which contains the inputted value. If no such bin is found, it does nothing. If a bin is found, by default, the value is incremented by 1. If a weight is given, then the value is incremented by that weight. The hit count is always incremented by 1. 

The `ValueHistogram` class contains different normalization procedures, and it defines basic arithmetic operations. Additionally, the function `ValueHistogram::combine()` takes a list of histograms and combines them into a single histogram by adding the values of each bin across the list.

When outputting to file, the `ValueHistogram` is exported as list of comma-separated values in the following format:

```
[bin center],[value],[error],[hitcount]
```

### File output

Both histograms and run parameters can be exported to files. The `Experiment` class includes a parameter `working_directory`, which acts as a base directory for the data output. For `TransverseMomentumExperiment`, an additional parameter `filename` specifies the name of the produced file. Two files are produced: a histogram file, the format of which was specified in the section _Histograms_, and a text file containing the values of different run parameters. Both files have the same filename, but different file extensions. By default, the histogram file has extension `csv` and the parameter file `txt`, although these can be changed with `histogram_file_extension` and `run_data_file_extension`, defined in the `Experiment` class. Thus, the two files produced have the paths

```
[working_directory]/[filename].csv
[working_directory]/[filename].txt
```

The working directory is relative to the location of the executable.


### Parallelization

Parallelization uses `OpenMP` and works by parallelizing the creationg and execution of `EventGenerator` instances. The `Experiment` class has a parameter `pT_hat_bins` and both `TransverseMomentumExperiment` and `AzimuthalCorrelationExperiment` parallelize over this list. What this means is that for each element of `pT_hat_bins`, which defines a phase space cut, an `EventGenerator` (or just a `ParticleGenerator`) is created and the phase space cut is passed along. The generator produces output, which is then combined across the different generators. For `TransverseMomentumExperiment`, five phase space cuts are applied, which means that five generators are created and five threads are used via `OpenMP`. 

For `AzimuthalCorrelationExperiment`, only a single phase space cut is needed. Still, the same cut can be replicated a number of times. This provides a convenient way of parallelization, even if it technically is abusing the system. The parameter `THREAD_COUNT`, defined as an `extern int` in `ExperimentDefs.cc` and initialized in `Experiments.cc` controls the number of copies of the phase space cut. Thus, this parameter effectively controls the number of threads used. In fact, each `EventGenerator` instance is given a target event count is `event_count / THREAD_COUNT`, since the total event count should be multiplied by the number of independent `EventGenerator` instances.

### Miscellaneous

- The set of transverse momentum filters in `AzimuthalCorrelationExperiment` (`pT_small`, `pT_large`) must be disjoint, similarly for `y_small` and `y_large`.
- In `AzimuthalCorrelationExperiment`, the `pT_range` and `y_range` are applied before `pT_small`, `pT_large`, `y_small` and `y_large`. Thus, restricting `pT_range` and `y_range` too much might mean that no particles will pass the two sets of transverse momentum and rapidity filters.
- The property `variable_seed` of `Experiment` is related to the initialization of multiple `EventGenerator` instances in `AzimuthalCorrelationExperiment`. Since the phase space cuts are just copies, the events produced would, by default, be identical, defeating the point of parallelization. This is because the PRNG used by Pythia is seeded with a specific number by default. The `variable_seed` flag ensures that each generator gets a unique random seed, resulting in an actual increase in statistics.
- Default parameter values are defined in the file `Constants.cc`.
- There are practically no checks. If you enter invalid data, or more likely forget to set some parameter, anything might happen. When weird results happen or the program crashes, the first to check is whether all parameters are properly initalized in classes like `Experiment` or `ParticleGenerator`.
