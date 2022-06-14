# alftools

alftools is an alternative to the [pyalf] package provided by the [ALF] team.

## 🔧 Installation

Install via `pip` from github:
```commandline
pip install git+https://github.com/dylanljones/alftools.git@VERSION
```

or download/clone the package, navigate to the root directory and install via
````commandline
pip install .
````

## 🚀 Quick-Start

### Running simulations

A new sdimulation can be run in a few steps:
```python
from alftools import Simulation

directory = "test"

# Create the simulation directory
sim = Simulation(directory)
sim.init(start_dir="Start")

# Update simulation parameters
params = sim.parameters
params.set_qmc("nbin", 100)
params.set_qmc("nsweep", 50)
params.set_hubbard("ham_u", 1.0)
...
params.save()

# Run simulation
sim.run()
```


An existing simulation can be continued as follows:
```python
sim.out_to_in()
sim.run()
```

[ALF]: (https://git.physik.uni-wuerzburg.de/ALF/pyALF)
[pyalf]: (https://git.physik.uni-wuerzburg.de/ALF/pyALF)
