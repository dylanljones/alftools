# alftools

alftools is an alternative to the [pyalf] package provided by the [ALF] team.


## Running simulations

A new sdimulation can be run in a few steps:
```python
from alftools import init_simulation, run_simulation, Parameters

directory = "test"

# Create the simulation directory
init_simulation(directory, start_dir="Start")

# Update simulation parameters
params = Parameters(directory)
params.set_qmc("nbin", 100)
params.set_qmc("nsweep", 50)
params.set_hubbard("ham_u", 1.0)
...
params.save()

# Run simulation
run_simulation(directory)
```


An existing simulation can be continued as follows:
```python
from alftools import out_to_in, run_simulation

directory = "test"
out_to_in(directory)
run_simulation(directory)
```

[ALF]: (https://git.physik.uni-wuerzburg.de/ALF/pyALF)
[pyalf]: (https://git.physik.uni-wuerzburg.de/ALF/pyALF)
