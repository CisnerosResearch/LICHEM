-- This is a template Lmod module file for LICHEM

local help_message=[[
More information about LICHEM can be found at
https://github.com/CisnerosResearch/LICHEM
]]

help(help_message, "\n")

whatis("Name: LICHEM")
whatis("Version: MOD_FIX_VERS_DATE, commit MOD_FIX_VERS_DATE")
whatis("URL: https://github.com/CisnerosResearch/LICHEM")

-- FIXME: Load module(s) for QM engine
--load("psi4")
--load("gaussian")
--load("nwchem")

-- FIXME: Load module(s) for MM Engine
--load("tinker")
--load("lammps")

-- FIXME: Load modules required for MPI (if applicable)
--load("mvapich2")

-- Create environment variables
local lichem_root    = "MOD_FIX_PATH"
local lib_root     = "/usr/lib64"

-- Set up overall environment when module is loaded
prepend_path( "PATH",   lichem_root)
prepend_path( "PATH",   lib_root)

-- Set scratch directory environment variable if using Gaussian
--setenv( "GAUSS_SCRDIR",  "/tmp")

