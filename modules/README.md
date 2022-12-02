# LICHEM Modules

This directory contains templates for both
[Lmod-](https://lmod.readthedocs.io/en/latest/)
and [environment-](http://modules.sourceforge.net/)
style modules.

Both modules will be updated by `make module`.

- `LICHEM_ROOT`/`lichem_root` is set to the installation path
- the version (date installed) and commit
(most recent commit in the history) are updated

## Changes Needed

Users will need to load any packages relevant to their system, including the
QM and/or MM modules they intend to use with LICHEM.

The modules will need to be added to the existing module path, or append
the folder to their module program's search list.

## Using Modules

Modules have to be added to the module search path.
You can check which style of modules you are using with `module --version`.
Most TACC resources use Lmod.

Example Lmod Output:
```
Modules based on Lua: Version 8.7.1  2022-05-02 13:33 -05:00
    by Robert McLay mclay@tacc.utexas.edu
```

Example Environment Output:
```
Modules Release 4.4.1 (2020-01-03)
```

### Module Use (Both)

Both Lmod and Environment modules use `module use` to add a path to the
module search path.
To always search this path, add the `module use` command to the `~/.bashrc`.

```
module use /path/to/modulefiles/with/lichem
```

### Using `use.own` (Only Environment Modules)

Environment Modules offer the `use.own` module that users can use to...
use their own modules.

```
$ cp lichem.tcl $HOME/privatemodules
$ module load use.own
```
