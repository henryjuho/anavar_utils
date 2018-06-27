# anavar_utils 

## Introduction

anavar_utils is a python package designed to facilitate the easy creation of 
control files for anavar ([Barton and Zeng, 2018](https://academic.oup.com/mbe/article/35/6/1536/4960016)) 
as well as providing classes for the easy interpretation of anavar's results files.

Anavar is available [here](http://zeng-lab.group.shef.ac.uk/wordpress/?page_id=28), 
the manual describes in much more detail the models available and the parameters in 
the control files than in these docs.

To use the package to create a simple control file:

```python
from __future__ import print_function
import anavar_utils as an

# data needed for anavar, SFS, callable sites and sample size
sample_size = 10
site_frequencies = [20, 19, 18, 17, 16, 15, 2, 1]
n_sites = 10000

# anavar_utils takes the sfs data in a dictionary
sfs_dict = {'SNP': (site_frequencies, n_sites)}

# intiate control file instance
control_file = an.Snp1ControlFile()

# set the data
control_file.set_data(sfs_dict, sample_size)

# construct control file string
control_contents = control_file.construct()

# output to a file or to stdout
print(control_contents)
```

This gives:

```
[algorithm_commands]
search_algorithm: NLOPT_LD_LBFGS
maxeval: 100000
maxtime: 600
num_searches: 500
nnoimp: 1
maximp: 3
optional: false

[model_commands]
model: SNP_1
n: 17
m: 10000
folded: false
sfs: 20, 19, 18, 17, 16, 15, 2, 1
dfe: discrete
c: 1
theta_range: 1e-06, 0.1
gamma_range: -250, 10
e_range: 0.0, 0.5
constraint: none

```

## API

### Control Files

#### Site frequency data format

In its simplest form a control file can be created for a given model 
with only the site frequency spectrum, the number of callable sites 
(in some cases and if per site parameter estimates are not required 
this can be set to 0) and the sample size. The control file classes 
take this data (using the ```.set_data()``` method) in the form of a 
dictionary, with keys specific to each model, as show in the table 
below.

| model  | dictionary format |
|:-------|:------------------|
| SNP_1  |```{'SNP': (sfs, n_sites)}``` |
| INDEL_1 | ```{'INS': (sfs, n_sites), 'DEL': (sfs, n_sites)}``` |
| gBGC_GLEMIN_EXTENDED_M1* | ```{'neutral_SNPs': (sfs, n_sites), 'ws_SNPs': (sfs, n_sites), 'sw_SNPs': (sfs, n_sites)}``` |
| neutralINDEL_vs_selectedINDEL | ```{'neutral_INS': (sfs, n_sites), 'neutral_DEL': (sfs, n_sites), 'selected_INS': (sfs, n_sites), 'selected_DEL': (sfs, n_sites)}``` |
| neutralSNP_vs_selectedSNP | ```{'neutral_SNP': (sfs, n_sites), 'selected_SNP': (sfs, n_sites)}``` |

The site frequency data needs to provide in the form of frequency counts 
from low freq to high freq with frequencies with no variants entered as 0.
For example this site frequency spectrum (SFS):

![sfs](sfs_example.png)

Would be expressed as:

```python
sfs = [100, 50, 33, 25, 20, 17, 14, 12, 11, 10, 9, 8, 8, 7, 7, 6, 6, 6, 5]
```

#### Classes

##### ```class Snp1ControlFile```

Snp1ControlFile()

Class to be used for creating a control file for the SNP_1 model. Also the parent class 
for all other control file classes.


```construct(self)``` 

Creates the final control file string with all specified settings

```set_alg_opts(self, alg='NLOPT_LD_LBFGS', maxeval=100000, maxtime=600, search=500, nnoimp=1, maximp=3, optional=False, size=10000, key=3, epsabs=1e-50, epsrel=1e-10, rftol=1e-10, init=())```

Sets algorithm options in control file

<b>Parameters:<b>

* alg - <i>(str)<i> - default'NLOPT_LD_LBFGS': the algorithm to use, choose from: 'NLOPT_LD_SLSQP', 'NLOPT_LD_LBFGS', 'NLOPT_LD_VAR1', 'NLOPT_LD_VAR2',
'NLOPT_LN_NELDERMEAD' or 'NLOPT_LD_TNEWTON_PRECOND_RESTART'

* maxeval - <i>(int)<i> - default=100000: search will be aborted when the likelihood function has been evaluated this number of times but the algorithm is still deemed non-converged
 
* maxtime - <i>(int)<i> - default=600: search will be aborted when the algorithm has been running for the given amount of time (in seconds)

* search - <i>(int)<i> - default=500: number of independent searches to be carried out

* nnoimp - <i>(int)<i> - default=1: if is set to 1, then no improvement search is attempted, otherwise, improvement searches will be conducted until nnoimp number of consecutive improvement attempts lead to no significant improvement in the ln-likelihood 

* maximp - <i>(int)<i> - default=3: maximum number of improvement searches

* optional - <i>(bool)<i> - default=False: if True, then optional commands will be included (optional commands are set separately, but will not be included until optional is set to True) 
  
* size - <i>(int)<i> - default=10000: controls QAG adaptive integration algorithm in GSL see: <https://www.gnu.org/software/gsl/manual/html node/QAG-adaptive- integration.html>

* key - <i>(int)<i> - default=3: controls QAG adaptive integration algorithm in GSL see: <https://www.gnu.org/software/gsl/manual/html node/QAG-adaptive- integration.html>

* epsabs - <i>(float)<i> - default (anavar_utils)=1e-50 - default (anavar1.4)=1e-5: controls QAG adaptive integration algorithm in GSL see: <https://www.gnu.org/software/gsl/manual/html node/QAG-adaptive- integration.html>

* epsrel - <i>(float)<i> - default (anavar_utils)=1e-10 - default (anavar1.4)=1e-8: controls QAG adaptive integration algorithm in GSL see: <https://www.gnu.org/software/gsl/manual/html node/QAG-adaptive- integration.html>

* rftol - <i>(float)<i> - default (anavar_utils)=1e-10 - default (anavar1.4)=1e-8: when consecutive search steps differ by less than the given value then the algorithm is considered to have converged

* init - <i>(tuple)<i> - default=(): takes a tuple of starting values of all parameters (in order of appearance in results file) for the model, default gives random starting values within parameter ranges

```set_constraint(self, constraint)```

<!---
     |      sets constraint
     |      :param constraint: str
     |      :return: sets str
     |  
     |  set_data(self, sfs_m, n, snp_fold=False, dfe='discrete', c=1, theta_r=(1e-06, 0.1), gamma_r=(-250, 10), error_r=(0.0, 0.5), shape_r=(0.001, 200), scale_r=(0.1, 1000.0), r_r=(0.05, 5.0))
     |      sets model and dfe commands in control file
     |      :param sfs_m: dict{str: tuple(list(float), int), ..}
     |      :param n: int
     |      :param snp_fold: bool
     |      :param dfe: str
     |      :param c: int
     |      :param theta_r: tuple(float, float)
     |      :param gamma_r: tuple(float, float)
     |      :param error_r: tuple(float, float)
     |      :param shape_r: tuple(float, float)
     |      :param scale_r: tuple(float, float)
     |      :param r_r: tuple(float, float)
     |      :return: NA
     |  
     |  set_dfe_optional_opts(self, optional=False, fraction=0.005, degree=50, delta=1e-05)
     |  
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |  
     |  __dict__
     |      dictionary for instance variables (if defined)
     |  
     |  __weakref__
     |      list of weak references to the object (if defined)


##### ```class GbgcControlFile```
     |  Method resolution order:
     |      GbgcControlFile
     |      Snp1ControlFile
     |      __builtin__.object
     |  
     |  Methods defined here:
     |  
     |  __init__(self)
     |  
     |  set_data(self, sfs_m, n, snp_fold=False, dfe='discrete', c=1, theta_r=(1e-06, 0.1), gamma_r=(-250, 10), error_r=(0.0, 0.5), shape_r=(0.001, 200), scale_r=(0.1, 1000.0), r_r=(0.05, 5.0))
     |      sets model and dfe commands in control file
     |      :param sfs_m: dict{str: tuple(list(float), int), ..}
     |      :param n: int
     |      :param snp_fold: bool
     |      :param dfe: str
     |      :param c: int
     |      :param theta_r: tuple(float, float)
     |      :param gamma_r: tuple(float, float)
     |      :param error_r: tuple(float, float)
     |      :param shape_r: tuple(float, float)
     |      :param scale_r: tuple(float, float)
     |      :param r_r: tuple(float, float)
     |      :return: NA
     |  
     |  ----------------------------------------------------------------------
     |  Methods inherited from Snp1ControlFile:
     |  
     |  construct(self)
     |  
     |  set_alg_opts(self, alg='NLOPT_LD_LBFGS', maxeval=100000, maxtime=600, search=500, nnoimp=1, maximp=3, optional=False, size=10000, key=3, epsabs=1e-50, epsrel=1e-10, rftol=1e-10, init=())
     |      sets algorithm options in control file
     |      :param alg: str
     |      :param maxeval: int
     |      :param maxtime: int
     |      :param search: int
     |      :param nnoimp: int
     |      :param maximp: int
     |      :param optional: bool
     |      :param size: int
     |      :param key: int
     |      :param epsabs: float
     |      :param epsrel: float
     |      :param rftol: float
     |      :param init: tuple
     |      :return: sets algorithm string
     |  
     |  set_constraint(self, constraint)
     |      sets constraint
     |      :param constraint: str
     |      :return: sets str
     |  
     |  set_dfe_optional_opts(self, optional=False, fraction=0.005, degree=50, delta=1e-05)
     |  
     |  ----------------------------------------------------------------------
     |  Data descriptors inherited from Snp1ControlFile:
     |  
     |  __dict__
     |      dictionary for instance variables (if defined)
     |  
     |  __weakref__
     |      list of weak references to the object (if defined)
    
    class Indel1ControlFile(Snp1ControlFile)
     |  Method resolution order:
     |      Indel1ControlFile
     |      Snp1ControlFile
     |      __builtin__.object
     |  
     |  Methods defined here:
     |  
     |  __init__(self)
     |  
     |  set_data(self, sfs_m, n, snp_fold=False, dfe='discrete', c=1, theta_r=(1e-06, 0.1), gamma_r=(-250, 10), error_r=(0.0, 0.5), shape_r=(0.001, 200), scale_r=(0.1, 1000.0), r_r=(0.05, 5.0))
     |      sets model and dfe commands in control file
     |      :param sfs_m: dict{str: tuple(list(float), int), ..}
     |      :param n: int
     |      :param snp_fold: bool
     |      :param dfe: str
     |      :param c: int
     |      :param theta_r: tuple(float, float)
     |      :param gamma_r: tuple(float, float)
     |      :param error_r: tuple(float, float)
     |      :param shape_r: tuple(float, float)
     |      :param scale_r: tuple(float, float)
     |      :param r_r: tuple(float, float)
     |      :return: NA
     |  
     |  ----------------------------------------------------------------------
     |  Methods inherited from Snp1ControlFile:
     |  
     |  construct(self)
     |  
     |  set_alg_opts(self, alg='NLOPT_LD_LBFGS', maxeval=100000, maxtime=600, search=500, nnoimp=1, maximp=3, optional=False, size=10000, key=3, epsabs=1e-50, epsrel=1e-10, rftol=1e-10, init=())
     |      sets algorithm options in control file
     |      :param alg: str
     |      :param maxeval: int
     |      :param maxtime: int
     |      :param search: int
     |      :param nnoimp: int
     |      :param maximp: int
     |      :param optional: bool
     |      :param size: int
     |      :param key: int
     |      :param epsabs: float
     |      :param epsrel: float
     |      :param rftol: float
     |      :param init: tuple
     |      :return: sets algorithm string
     |  
     |  set_constraint(self, constraint)
     |      sets constraint
     |      :param constraint: str
     |      :return: sets str
     |  
     |  set_dfe_optional_opts(self, optional=False, fraction=0.005, degree=50, delta=1e-05)
     |  
     |  ----------------------------------------------------------------------
     |  Data descriptors inherited from Snp1ControlFile:
     |  
     |  __dict__
     |      dictionary for instance variables (if defined)
     |  
     |  __weakref__
     |      list of weak references to the object (if defined)
    
    class IndelNeuSelControlFile(Snp1ControlFile)
     |  Method resolution order:
     |      IndelNeuSelControlFile
     |      Snp1ControlFile
     |      __builtin__.object
     |  
     |  Methods defined here:
     |  
     |  __init__(self)
     |  
     |  set_data(self, sfs_m, n, snp_fold=False, dfe='discrete', c=1, theta_r=(1e-06, 0.1), gamma_r=(-250, 10), error_r=(0.0, 0.5), shape_r=(0.001, 200), scale_r=(0.1, 1000.0), r_r=(0.05, 5.0))
     |      sets model and dfe commands in control file
     |      :param sfs_m: dict{str: tuple(list(float), int), ..}
     |      :param n: int
     |      :param snp_fold: bool
     |      :param dfe: str
     |      :param c: int
     |      :param theta_r: tuple(float, float)
     |      :param gamma_r: tuple(float, float)
     |      :param error_r: tuple(float, float)
     |      :param shape_r: tuple(float, float)
     |      :param scale_r: tuple(float, float)
     |      :param r_r: tuple(float, float)
     |      :return: NA
     |  
     |  ----------------------------------------------------------------------
     |  Methods inherited from Snp1ControlFile:
     |  
     |  construct(self)
     |  
     |  set_alg_opts(self, alg='NLOPT_LD_LBFGS', maxeval=100000, maxtime=600, search=500, nnoimp=1, maximp=3, optional=False, size=10000, key=3, epsabs=1e-50, epsrel=1e-10, rftol=1e-10, init=())
     |      sets algorithm options in control file
     |      :param alg: str
     |      :param maxeval: int
     |      :param maxtime: int
     |      :param search: int
     |      :param nnoimp: int
     |      :param maximp: int
     |      :param optional: bool
     |      :param size: int
     |      :param key: int
     |      :param epsabs: float
     |      :param epsrel: float
     |      :param rftol: float
     |      :param init: tuple
     |      :return: sets algorithm string
     |  
     |  set_constraint(self, constraint)
     |      sets constraint
     |      :param constraint: str
     |      :return: sets str
     |  
     |  set_dfe_optional_opts(self, optional=False, fraction=0.005, degree=50, delta=1e-05)
     |  
     |  ----------------------------------------------------------------------
     |  Data descriptors inherited from Snp1ControlFile:
     |  
     |  __dict__
     |      dictionary for instance variables (if defined)
     |  
     |  __weakref__
     |      list of weak references to the object (if defined)
    
    class MultiResultsFile(ResultsFile)
     |  Method resolution order:
     |      MultiResultsFile
     |      ResultsFile
     |      __builtin__.object
     |  
     |  Methods defined here:
     |  
     |  __init__(self, file_list)
     |      takes a list of anavar results files and merges them into one,
     |      then calls ResultsFile instance on merged file before deleting
     |      merged file
     |      :param file_list: list
     |  
     |  ----------------------------------------------------------------------
     |  Methods inherited from ResultsFile:
     |  
     |  bounds_hit(self, theta_r=(1e-06, 0.1), gamma_r=(-250, 10), error_r=(0.0, 0.5), shape_r=(0.001, 200), scale_r=(0.1, 1000.0), r_r=(0.05, 5.0))
     |      lists parameters too close to estimate boundaries
     |      :param theta_r: tuple
     |      :param gamma_r: tuple
     |      :param error_r: tuple
     |      :param shape_r: tuple
     |      :param scale_r: tuple
     |      :param r_r: tuple
     |      :return: list
     |  
     |  control_file(self)
     |      returns control file path
     |      :return: str
     |  
     |  converged(self)
     |      assesses convergence for results file
     |      :return: bool
     |  
     |  data_type(self)
     |      assigns data as either indel, snp or snp_indel
     |      :return: str
     |  
     |  dfe(self)
     |      returns type of distribution
     |      :return: str
     |  
     |  estimates(self, as_string=False)
     |      returns generator of all estimates in dictionary form
     |      :param as_string: bool
     |      :return: generator
     |  
     |  free_parameters(self)
     |      returns all free parameters in results file
     |      :return: tuple
     |  
     |  get_alpha(self, dn, ds, var_type)
     |      calculates alpha as described in equation 19 Barton and Zeng 2018
     |      :param dn: float
     |      :param ds: float
     |      :param var_type: str
     |      :return: str
     |  
     |  header(self)
     |      returns the column headers line from the file
     |      :return: tuple
     |  
     |  ml_estimate(self, as_string=False)
     |      returns the maximum likelihood estimate
     |      :param as_string: bool
     |      :return: dict
     |  
     |  num_class(self)
     |      returns number of classes in results file
     |      :return: int
     |  
     |  num_runs(self)
     |      gives the number of runs listed in the results file
     |      :return: int
     |  
     |  seed(self)
     |      returns the seed anavar was run with
     |      :return: int
     |  
     |  ----------------------------------------------------------------------
     |  Data descriptors inherited from ResultsFile:
     |  
     |  __dict__
     |      dictionary for instance variables (if defined)
     |  
     |  __weakref__
     |      list of weak references to the object (if defined)
    
    class ResultsFile(__builtin__.object)
     |  Methods defined here:
     |  
     |  __init__(self, anavar_results_file)
     |      Creates an anavar results file object
     |      :param anavar_results_file: file
     |  
     |  bounds_hit(self, theta_r=(1e-06, 0.1), gamma_r=(-250, 10), error_r=(0.0, 0.5), shape_r=(0.001, 200), scale_r=(0.1, 1000.0), r_r=(0.05, 5.0))
     |      lists parameters too close to estimate boundaries
     |      :param theta_r: tuple
     |      :param gamma_r: tuple
     |      :param error_r: tuple
     |      :param shape_r: tuple
     |      :param scale_r: tuple
     |      :param r_r: tuple
     |      :return: list
     |  
     |  control_file(self)
     |      returns control file path
     |      :return: str
     |  
     |  converged(self)
     |      assesses convergence for results file
     |      :return: bool
     |  
     |  data_type(self)
     |      assigns data as either indel, snp or snp_indel
     |      :return: str
     |  
     |  dfe(self)
     |      returns type of distribution
     |      :return: str
     |  
     |  estimates(self, as_string=False)
     |      returns generator of all estimates in dictionary form
     |      :param as_string: bool
     |      :return: generator
     |  
     |  free_parameters(self)
     |      returns all free parameters in results file
     |      :return: tuple
     |  
     |  get_alpha(self, dn, ds, var_type)
     |      calculates alpha as described in equation 19 Barton and Zeng 2018
     |      :param dn: float
     |      :param ds: float
     |      :param var_type: str
     |      :return: str
     |  
     |  header(self)
     |      returns the column headers line from the file
     |      :return: tuple
     |  
     |  ml_estimate(self, as_string=False)
     |      returns the maximum likelihood estimate
     |      :param as_string: bool
     |      :return: dict
     |  
     |  num_class(self)
     |      returns number of classes in results file
     |      :return: int
     |  
     |  num_runs(self)
     |      gives the number of runs listed in the results file
     |      :return: int
     |  
     |  seed(self)
     |      returns the seed anavar was run with
     |      :return: int
     |  
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |  
     |  __dict__
     |      dictionary for instance variables (if defined)
     |  
     |  __weakref__
     |      list of weak references to the object (if defined)
    
    class SNPNeuSelControlFile(Snp1ControlFile)
     |  Method resolution order:
     |      SNPNeuSelControlFile
     |      Snp1ControlFile
     |      __builtin__.object
     |  
     |  Methods defined here:
     |  
     |  __init__(self)
     |  
     |  set_data(self, sfs_m, n, snp_fold=False, dfe='discrete', c=1, theta_r=(1e-06, 0.1), gamma_r=(-250, 10), error_r=(0.0, 0.5), shape_r=(0.001, 200), scale_r=(0.1, 1000.0), r_r=(0.05, 5.0))
     |      sets model and dfe commands in control file
     |      :param sfs_m: dict{str: tuple(list(float), int), ..}
     |      :param n: int
     |      :param snp_fold: bool
     |      :param dfe: str
     |      :param c: int
     |      :param theta_r: tuple(float, float)
     |      :param gamma_r: tuple(float, float)
     |      :param error_r: tuple(float, float)
     |      :param shape_r: tuple(float, float)
     |      :param scale_r: tuple(float, float)
     |      :param r_r: tuple(float, float)
     |      :return: NA
     |  
     |  ----------------------------------------------------------------------
     |  Methods inherited from Snp1ControlFile:
     |  
     |  construct(self)
     |  
     |  set_alg_opts(self, alg='NLOPT_LD_LBFGS', maxeval=100000, maxtime=600, search=500, nnoimp=1, maximp=3, optional=False, size=10000, key=3, epsabs=1e-50, epsrel=1e-10, rftol=1e-10, init=())
     |      sets algorithm options in control file
     |      :param alg: str
     |      :param maxeval: int
     |      :param maxtime: int
     |      :param search: int
     |      :param nnoimp: int
     |      :param maximp: int
     |      :param optional: bool
     |      :param size: int
     |      :param key: int
     |      :param epsabs: float
     |      :param epsrel: float
     |      :param rftol: float
     |      :param init: tuple
     |      :return: sets algorithm string
     |  
     |  set_constraint(self, constraint)
     |      sets constraint
     |      :param constraint: str
     |      :return: sets str
     |  
     |  set_dfe_optional_opts(self, optional=False, fraction=0.005, degree=50, delta=1e-05)
     |  
     |  ----------------------------------------------------------------------
     |  Data descriptors inherited from Snp1ControlFile:
     |  
     |  __dict__
     |      dictionary for instance variables (if defined)
     |  
     |  __weakref__
     |      list of weak references to the object (if defined)
    
    
FUNCTIONS
    merge_results(file_list, name)
        takes list of results files and merges them
        :param file_list: list
        :param name: str
        :return:

DATA
    print_function = _Feature((2, 6, 0, 'alpha', 2), (3, 0, 0, 'alpha', 0)...

-->


### Results Files