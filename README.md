# splash

Collection of scripts and Jupyter notebooks to predict the splashback radius in symmetron gravity.

## Fortran solver
* symsplash/fsolver.f95 is a fortran solver. To create a python wrapper of the routine, f2py is needed. See <https://docs.scipy.org/doc/numpy/f2py/>. To compile the .so module run the following:

`f2py -c fsolver.f95 -m fsolver`

## Notebooks & files
The workflow includes the following scripts.

* SSinit.ipnyb : where the self-similar profile solution is obtained through iteration.

* runSSs.py: an example of how to use the class `SSHalo` to get the change in splashback position.

* Plot.ipnyb : a handy notebook to create plots.


## Documentation

You can access every method's docstring by using the help() function in python. All of the under the hood machinery is defined in `symsplash/colapse.py`.
