#%Module1.0
#
# LICHEM modulefile
#
# set local variables
set     name            "LICHEM"
set     version         "Version: MOD_FIX_VERS_DATE, commit MOD_FIX_VERS_DATE"
set     url             "https://github.com/CisnerosResearch/LICHEM"
set     LIB_ROOT        "/usr/lib64"
set     LICHEM_ROOT     "MOD_FIX_PATH"

proc ModulesHelp { } {
        global name
        global url
        puts stderr "More information about $name can be found at:"
        puts stderr "    $url\n"
}

module-whatis   "Name: $name\nVersion: $version\nURL: $url\n"

# FIXME: Load module(s) for QM engine
#module load gaussian
#module load psi4
#module load nwchem

# FIXME: Load module(s) for MM engine
#module load tinker
#module load lammps

# FIXME: Load modules required for MPI (if applicable)

# Set up overall environment when module is loaded
prepend-path    PATH              "${LICHEM_ROOT}"
prepend-path    PATH              "${LIB_ROOT}"

# Set scratch directory environment variable if using Gaussian
setenv          GAUSS_SCRDIR      "/tmp"

