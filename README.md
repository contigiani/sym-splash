# splash

Collection of scripts and Jupyter notebooks to predict the splashback radius in symmetron gravity.

## Fortran solver
* symsplash/fsolver.f95 is a fortran solver for the static symmetron equation of motion. To create a python wrapper of the routine, f2py is needed (see <https://docs.scipy.org/doc/numpy/f2py/>). To compile the .so module run the following code:

```
f2py -c fsolver.f95 -m fsolver
```

## Notebooks & files
The workflow includes the following scripts

* SSinit.ipnyb : where the self-similar profile solution is obtained through iteration.

* runSSs.py: an example of how to use the class `SSHalo` to get the change in splashback position.

* Plot.ipnyb : a handy notebook to create plots.

Everything is self-contained and, apart from running the f2py command described above, the final products can be generated in a few minutes.

## Documentation

You can access every method's docstring using the help() function in python. Note that the under the hood machinery is defined in `symsplash/collapse.py`.
