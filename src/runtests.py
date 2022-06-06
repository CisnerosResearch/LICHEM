###################################################
#                                                 #
#   LICHEM: Layered Interacting CHEmical Models   #
#                                                 #
#        Symbiotic Computational Chemistry        #
#                                                 #
###################################################

# LICHEM semi-automated test suite

###
# usage: ./runtests [-h] [-v] [-n [NCPUS]] [-d] [-t TESTS [TESTS ...]] [-a]
#                   [-q [QM]] [-m [MM]]
#
# Run the LICHEM test suite.
#
# optional arguments:
#   -h, --help            show this help message and exit
#   -v, --verbose         Use developer options and print debugging
#                         information (default: False)
#
# Job Settings (optional):
#   Define conditions for running tests.
#
#   -n [NCPUS], --ncpus [NCPUS]
#                         Number of CPUs (default: 1)
#   -d, --dry             Perform a dry run of all tests (default: False)
#   -t TESTS [TESTS ...], --tests TESTS [TESTS ...]
#                         Explicit tests to run. Options: ['HF', 'PBE0',
#                         'CCSD', 'PM6', 'Frequencies', 'NEB_TS', 'QSM_TS',
#                         'TIP3P', 'AMOEBA/GK', 'PBE0/TIP3P', 'PBE0/AMOEBA',
#                         'DFP/Pseudobonds'] (default: all)
#
# Wrappers (optional):
#   Specify QM and MM Wrappers for tests. Use (-a|--all) to search for
#   multiple versions of the same program (ex., g09 and g16). Only 1
#   executable of each name will be tested. (So two foo2 executables
#   wouldn't be tested, but foo2 and foo3 would be.)
#
#   -a, --all             Auto-run all tests for each available wrapper.
#                         (default: False)
#   -q [QM], --qm [QM]    The QM wrapper to test. Options: ['gaussian',
#                         'g09', 'g16', 'nwchem', 'psi4'] Note: Gaussian and
#                         g09 are equivalent. (default: all)
#   -m [MM], --mm [MM]    The MM wrapper to test. Options: ['tinker',
#                         'tinker9', 'lammps'] (default: all)
###

# TODO:
# - Rework the if/if process for QM/MM engines tests (shorten code)
# - Figure out where QSM TS tinker.key copy attempt is...

# DevNote:
# > The comments are removed to make this an executable. Docstrings will be
#    printed, however. Because this is "private" code, PEP8 does not actually
#    require docstrings. For the best of both worlds, write function
#    info in the Numpy docstring style, and then comment each line out
#    with '# '.

"""
NOTE: See ../src/runtests.py for commented code!

Requires Python 3.6+ to run. Update the path as necessary!

Usage and help information printed with:

    user:$ ./runtests -h
      or
    user:$ ./runtests --h
"""

# Modules
import argparse
import os
import platform
import subprocess
import sys
import time

# Start timer immediately
startTime = time.time()

# Initialize globals
TTxtLen = 30  # Number of characters for the test name
passCt = 0    # Number of tests passed
failCt = 0    # Number of tests failed
skipCt = 0    # Number of tests skipped
testCt = 0    # Current test number

# Development settings -- can be forced to True with `-v` option!
# NB: Modified by the Makefile (make devtestexe)
updateResults = False  # Bool to print energies to update tests
forceAll = False       # Bool to force it to do tests even if they will fail
debugMode = False      # Bool to print LICHEM command used if test fails

# Note: All round statements expect eV!
har2eV = 27.21138505  # Convert eV to au/Hartree

# Assume g09 by default
useg16 = False

# List of regions files needing to be updated with program versions
regions_files = ["ccsdreg.inp", "freqreg.inp", "hfreg.inp", "mmreg.inp",
                 "nebreg.inp", "pbereg.inp", "pboptreg.inp", "pchrgreg.inp",
                 "pm6reg.inp", "polreg.inp",  "qsmreg.inp", "solvreg.inp"]

# List of possible QM wrappers
QM_wrappers = ['gaussian', 'g09', 'g16', 'nwchem', 'psi4']
# List of possible MM wrappers
MM_wrappers = ['tinker', 'tinker9', 'lammps']
# List of possible tests to run
test_opts = ["HF", "PBE0", "CCSD", "PM6", "Frequencies", "NEB_TS", "QSM_TS",
             "TIP3P", "AMOEBA/GK", "PBE0/TIP3P", "PBE0/AMOEBA",
             "DFP/Pseudobonds"]


# ---------------------- #
# --- Define Classes --- #
# ---------------------- #


class ClrSet:
    # """
    # Define the colors for messages printed to the console.
    # """
    # Unicode colors
    Norm = '\033[0m'
    Red = '\033[91m'
    Bold = '\033[1m'
    Green = '\033[92m'
    Blue = '\033[94m'
    Cyan = '\033[36m'
    # Set colors
    TFail = Bold+Red    # Highlight failed tests
    TPass = Bold+Green  # Highlight passed tests
    TSkip = Bold+Cyan   # Highlight skipped tests (in result)
    Reset = Norm        # Reset to defaults


# ------------------------ #
# --- Define functions --- #
# ------------------------ #

def get_args():
    # """
    # Read in provided arguments and define command-line usage.
    #
    # Returns
    # -------
    # parser.parse_args() : argparse.Namespace
    #     Program details specified by command-line arguments.
    # """
    global QM_wrappers
    global MM_wrappers
    global test_opts
    parser = argparse.ArgumentParser(
                    prog="./runtests",
                    description="Run the LICHEM test suite.",
                    # Print the default values
                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-v", "--verbose", action='store_true',
                        help=("Use developer options and print debugging "
                              "information"))
    job = parser.add_argument_group(
                title="Job Settings (optional)",
                description="Define conditions for running tests.")
    # nargs = ? saves only 1 value
    job.add_argument("-n", "--ncpus", nargs='?', default=1, type=int,
                     help="Number of CPUs")
    # store_true will only set to True if argument given!
    job.add_argument("-d", "--dry", action='store_true',
                     help="Perform a dry run of all tests")
    # nargs = + saves list of values
    job.add_argument("-t", "--tests", nargs="+", default="all",
                     type=str.lower,
                     help=f"Explicit tests to run. Options: {test_opts}")
    # Create a group of wrappers to lists together in help!
    wrappers = parser.add_argument_group(
                title="Wrappers (optional)",
                description=("Specify QM and MM Wrappers for tests. "
                             "Use (-a|--all) to search for multiple versions "
                             "of the same program (ex., g09 and g16). "
                             "Only 1 executable of each name will be tested. "
                             "(So two foo2 executables wouldn't be tested, "
                             "but foo2 and foo3 would be.)"))
    wrappers.add_argument("-a", "--all", action='store_true',
                          help=("Auto-run all tests for each available "
                                "wrapper."))
    wrappers.add_argument("-q", "--qm", nargs="+", default="all",
                          type=str.lower,
                          help=("The QM wrapper to test.\n"
                                f"Options: {QM_wrappers}\n"
                                "  Note: Gaussian and g09 are equivalent."))
    wrappers.add_argument("-m", "--mm", nargs="+", default="all",
                          type=str.lower,
                          help=("The MM wrapper to test.\n"
                                f"Options: {MM_wrappers}"))
    return parser.parse_args()


def LocateLICHEM():
    # """
    # Check that the LICHEM binary is in the user's PATH.
    #
    # Returns
    # -------
    # saved_lichem_bin : str
    #     The path to the LICHEM binary.
    # """
    global allTests
    cmd = "which lichem"
    # Find the binary and save the colored and uncolored output
    try:
        LICHEMbin = subprocess.check_output(cmd, shell=True)
        saved_lichem_bin = LICHEMbin.decode('utf-8').strip()  # Uncolored
        LICHEMbin = ClrSet.TPass+LICHEMbin.decode('utf-8').strip()+ClrSet.Reset
    # Exit if the binary isn't found
    except subprocess.CalledProcessError:
        LICHEMbin = ClrSet.TFail+"N/A"+ClrSet.Reset
        print(f"LICHEM binary: {LICHEMbin}\n")
        print("Error: LICHEM binary not found!\n")
        exit(0)
    # Print the colored binary path output
    else:
        if allTests is True:
            print(f"LICHEM binary: {LICHEMbin}\n")
    return saved_lichem_bin


def LocateProgram(executable):
    # """
    # Check that the program binary is in the user's PATH.
    #
    # Returns
    # -------
    # saved_bin : str
    #     The path to the program binary.
    # """
    global allTests
    global useg16
    cmd = "which {}".format(executable)
    # Find the binary and save the colored and uncolored output
    try:
        prog_bin = subprocess.check_output(cmd, shell=True)
        saved_bin = prog_bin.decode('utf-8').strip()  # Uncolored
        prog_bin = ClrSet.TPass+prog_bin.decode('utf-8').strip()+ClrSet.Reset
    except subprocess.CalledProcessError:
        # print("\n")  # Print newline after "which" failure
        # If g09 fails, try g16
        if executable == 'g09':
            print("g09 failed I guess...")
            cmd = "which g16"
            try:
                prog_bin = subprocess.check_output(cmd, shell=True)
                saved_bin = prog_bin.decode('utf-8').strip()  # Uncolored
                prog_bin = (ClrSet.TPass+prog_bin.decode('utf-8').strip() +
                            ClrSet.Reset)
                useg16 = True  # Use Gaussian 16 for tests
            except subprocess.CalledProcessError:
                prog_bin = ClrSet.TFail+"N/A"+ClrSet.Reset
                saved_bin = "N/A"
        else:
            prog_bin = ClrSet.TFail+"N/A"+ClrSet.Reset
            saved_bin = "N/A"
    if allTests is True:
        saved_bin = prog_bin  # Use colored version
    return saved_bin


def PrepRegions(keyword, d_val, u_val, file):
    # """
    # Replace the current keyword argument in a region file with another.
    #
    # Parameters
    # ----------
    # keyword : str
    #     The keyword to search for in the LICHEM regions file.
    # d_val : str
    #     The default value for the keyword.
    # u_val : str
    #     The updated value for the keywork.
    # file : str
    #     The name of the file to operate on.
    # """
    global revertMode
    # Account for difference in sed in-line between OSX and Linux
    # ex: sed -i '' '/^QM_type/s/Gaussian/g16'
    if platform.system() == "Darwin":
        cmd = ("sed -i '' '/^{}/s/{}/{}/' {}".format(keyword, d_val,
                                                     u_val, file))
        subprocess.call(cmd, shell=True)
    elif platform.system() == "Linux":
        cmd = ("sed -i '/^{}/s/{}/{}/' {}".format(keyword, d_val,
                                                  u_val, file))
        subprocess.call(cmd, shell=True)
    else:
        # Print warning instead of raising error so they can correct for
        #  future, since this prints due to an unknown platform.
        print("WARNING: your OS doesn't seem to be either OSX or Linux.\n"
              f"If you're using {d_val} executable, you'll need to modify\n"
              f"  the '{keyword}' block to *reg.inp files to '{u_val}'\n"
              "   manually.\n"
              f"This is because the program name calls '{d_val}'.\n"
              "These tests will likely fail!")
    if debugMode is True and revertMode is False:
        print("Modifying the regions files for this wrapper combination...")
        print("\n{' ':6}- {cmd}\n")
    elif debugMode is True and revertMode is True:
        print("Reverting the regions files for this wrapper combination...")
        print("\n{' ':6}- {cmd}\n")
    return


def RunLICHEM(xName, rName, cName):
    # """
    # Call a LICHEM instance for testing.
    #
    # Parameters
    # ----------
    # xName : str
    #     Name of the LICHEM XYZ file.
    # rName : str
    #     Name of the regions file.
    # cName : str
    #     Name of the connectivity file.
    #
    # Returns
    # -------
    # cmd_printed : str
    #     Formatted version of the attempted LICHEM command (for printing).
    # """
    global Ncpus
    # -l captures stdout, 2>&1 captures stderr
    cmd = ("lichem -n {n} -x {x} -r {r} -c {c} -o trash.xyz "
           "-l tests.out 2>&1").format(n=Ncpus, x=xName, r=rName, c=cName)
    subprocess.call(cmd, shell=True)  # Run calculations
    # Save a terminal-formatted version of the command
    #  \\\n escapes a printed backslash and then prints newline
    cmd_printed = ("lichem -n {n} -x {x} -r {r} \\\n{sp:9}-c {c} "
                   "-o trash.xyz -l tests.out 2>&1").format(
                        n=Ncpus, x=xName, r=rName, sp=' ', c=cName)
    return cmd_printed


def CleanFiles():
    # """
    # Clean out various files after tests are performed.
    # """
    # Delete junk files
    cleanCmd = "rm -f"
    # Remove LICHEM files
    cleanCmd += " BASIS tests.out trash.xyz"
    cleanCmd += " BeadStartStruct.xyz BurstStruct.xyz"
    # Remove TINKER files
    cleanCmd += " tinker.key"
    # Remove LAMMPS files
    cleanCmd += ""
    # Remove Gaussian files
    cleanCmd += " *.chk"
    # Remove PSI4 files
    cleanCmd += " timer.* psi.* *.32 *.180"
    # Remove NWChem files
    cleanCmd += " *.movecs"
    # Delete the files
    subprocess.call(cleanCmd, shell=True)
    # Wait 1 second for files to be deleted (avoid race conditions for I/O)
    try:
        time.sleep(1)
    except KeyboardInterrupt:
        print("Fine, I won't nap, gosh.")
    return


def RecoverEnergy(txtLabel, itemNum):
    # """
    # Recover the energy from the LICHEM output.
    #
    # Parameters
    # ----------
    # txtLabel : str
    #     Descriptive label in the LICHEM output used to identify energy line.
    # itemNum : int
    #     The split index in the identified line to pull the energy from.
    #
    # Returns
    # -------
    # finalEnergy : float
    #     The energy from the LICHEM output, rounded to 3 decimal places.
    #     If failed, the energy is set to 0.0.
    # savedResult : str
    #     Either the unrounded finalEnergy or a crashed message.
    # units : str
    #     The unit following the finalEnergy in the output.
    # """
    cmd = 'grep -e "{}" tests.out | tail -1'.format(txtLabel)
    savedResult = "Crashed..."
    try:
        # Safely check energy
        finalEnergy = subprocess.check_output(cmd, shell=True)  # Get results
        finalEnergy = finalEnergy.decode('utf-8').split()
        units = str(finalEnergy[itemNum+1])  # For developer mode
        finalEnergy = float(finalEnergy[itemNum])
        savedResult = "Energy: "+str(finalEnergy)  # Save it for later
        finalEnergy = round(finalEnergy, 3)
    # TODO: Figure out what except codes this will throw!
    except:
        # Calculation failed
        finalEnergy = 0.0
        units = ":("
        time.sleep(2)
    return finalEnergy, savedResult, units


def RecoverFreqs():
    # """
    # Recover a list of frequencies from the LICHEM output.
    #
    # Returns
    # -------
    # freqList : lst
    #     Frequencies from LICHEM output.
    # """
    global debugMode
    # sed commands to use for pulling frequencies
    cmd = ""
    cmd += "sed '/Usage Statistics/,$d' tests.out | "
    cmd += "sed -n '/Frequencies:/,$p' | "
    cmd += "sed '/Frequencies:/d'"
    # Set fake units for printing in developer mode
    units = " "
    try:
        # Safely check energy
        freqList = []
        tmpFreqs = subprocess.check_output(cmd, shell=True)  # Get results
        tmpFreqs = tmpFreqs.decode('utf-8').strip()
        tmpFreqs = tmpFreqs.split()
        for freqVal in tmpFreqs:
            freqList.append(float(freqVal))
    # TODO: Likely need this to be subprocess.CalledProcessError and exception
    #  raised by not being able to split output...
    except:
        # Calculation failed
        if debugMode is True:
            print("\nFrequencies not recovered.")
        freqList = []
    return freqList, units


def CompareEnergy(QMMMPack, QMMMEnergy, known, known_units=None):
    # """
    # Checks the calculated energy against a known value.
    #
    # Parameters
    # ----------
    # QMMMPack : str
    #     Name of the package being tested.
    # QMMMEnergy : float
    #     Value from the LICHEM log file.
    # known : float
    #     Expected value for the program.
    # known_units : str
    #     Units for the expected value.
    # """
    expected_energy = known
    if (QMMMEnergy == expected_energy):
        passEnergy = True
        return passEnergy
    else:
        saved_fail = SaveFailure(QMMMEnergy, expected_energy)
        return saved_fail


def SaveFailure(calc, expected):
    # """
    # Print the expected and test-calculated energies.
    #
    # Parameters
    # ----------
    # calc : float
    #     QM/MM/QMMM energy from the test LICHEM output.
    # expected : float
    #     Energy expected from output.
    #
    # Returns
    # -------
    # saved_val : str
    #     Expected and calculated values formatted for console printing.
    # """
    # The function input order matches the if statements!!!
    # Note: You don't need to add saved values because all tests of a pairing
    #       will run, so the if statements are independent.
    saved_val = "\n{:4}---> Expected: {}, Calculated: {}".format(
        " ", expected, calc)
    return saved_val


def PrintLICHEMDebug(cmd_printed):
    # """
    # Prints the LICHEM command that was run to the console.
    #
    # Parameters
    # -------
    # cmd_printed : str
    #     Formatted version of the attempted LICHEM command.
    # """
    global debugMode
    if debugMode is True:
        print(f"\n{' ':6}Tried to execute:\n{' ':7}{cmd_printed}")
    return


def PrintPassFail(tName, testPass, updateResults, enVal, units):
    # """
    # Write pass/fail message to the console with coloration and timing.
    #
    # Parameters
    # ----------
    # tName: str
    #     Name of the test being run.
    # testPass : bool
    #     True if test passed.
    # updateResults : bool
    #     True to print energies to update tests (for developers).
    # enVal : float
    #     The computed energy value.
    # units : str
    #     The units of the computed energy value.
    # """
    global passCt
    global failCt
    global TTxtLen
    # Get run time
    cmd = 'grep -e "Total wall time: " tests.out'
    try:
        runTime = subprocess.check_output(cmd, shell=True)  # Get run time
        runTime = runTime.decode('utf-8').split()
        runTime = "{:.4f} {}".format(float(runTime[3]), runTime[4])
    except (subprocess.CalledProcessError, KeyboardInterrupt):
        runTime = "N/A"
    if updateResults is True:
        if testPass is True:
            print(f"    {tName:<{TTxtLen-5}} " +
                  f"{ClrSet.TPass}Pass{ClrSet.Reset}, {runTime}, " +
                  f"{enVal} {units}")
            passCt += 1
        else:
            print(f"{' ':4}{tName:<{TTxtLen-5}} " +
                  f"{ClrSet.TFail}Fail{ClrSet.Reset}, {runTime}, " +
                  f"{enVal} {units}")
            failCt += 1
    else:
        if testPass is True:
            print(f"    {tName:<{TTxtLen-5}} " +
                  f"{ClrSet.TPass}Pass{ClrSet.Reset}, {runTime}")
            passCt += 1
        else:
            print(f"    {tName:<{TTxtLen-5}} " +
                  f"{ClrSet.TFail}Fail{ClrSet.Reset}, {runTime}")
            failCt += 1
    pf_printed = True
    return pf_printed


def SkipPass(tName, pf_printed):
    # """
    # Color the skip message printed to the console.
    #
    # Parameters
    # ----------
    # tName: str
    #     Name of the test being run.
    # pf_printed : bool
    #     True if pass/fail has already been printed. This is intended to help
    #      avoid race conditions with SkipSequence() after a KeyboardInterrupt,
    #      resulting in the skip being printed after a test has passed.
    # """
    global TTxtLen
    global skipCt
    if pf_printed is False:
        runTime = "N/A"
        # Print to max width - original 4+1 spaces (because less than)
        print(f"\n    {tName:<{TTxtLen-5}} " +
              f"{ClrSet.TSkip}Skip{ClrSet.Reset}, {runTime}")
        skipCt += 1
    return


def SkipSequence(tName, pf_printed):
    # """
    # Either skip or quit a test based off user input.
    #
    # Parameters
    # ----------
    # tName: str
    #     Name of the test being run.
    # pf_printed : bool
    #     True if pass/fail has already been printed. This is intended to help
    #      avoid race conditions with SkipSequence() after a KeyboardInterrupt,
    #      resulting in the skip being printed after a test has passed.
    # """
    print(f"\n{' ':6}Do you want to...")
    print(f"{' ':7}1. Skip this test only? (s|skip) [default]")
    print(f"{' ':7}2. Quit all? (q|quit|e|exit)")
    print(f"{' ':7}3. Quit all without clearing most recent files? "
          "(dev|debug)")
    try:
        response = input("      Selection: ")
    # Case where someone aggressively hits Cntrl+C
    except KeyboardInterrupt:
        print(f"\n{' ':7}Stopping all tests.\n")
        CleanFiles()
        exit(0)

    # Selected skip
    if response.lower() in ("s", "skip", "default"):
        print(f"\n{' ':7}Skipping this test.\n")
        SkipPass(tName, pf_printed)
        CleanFiles()
    # Default skip (where user pressed enter)
    elif response == "":
        print(f"\n{' ':7}Skipping this test.\n")
        SkipPass(tName, pf_printed)
        CleanFiles()
    # Normal quit
    elif response.lower() in ("q", "quit", "e", "exit"):
        print(f"\n{' ':7}Stopping all tests.\n")
        CleanFiles()
        exit(0)
    # Developer mode
    elif response.lower() in ("dev", "debug"):
        print(f"\n{' ':7}Stopping all tests. Recent output won't be "
              "cleared.\n")
        exit(0)
    # A cat crossed the keyboard...
    else:
        print(f"\n{' ':7}Response not understood. Skipping this test...\n")
        SkipPass(tName, pf_printed)
        CleanFiles()
    return


# --------------------------------- #
# --- Define the Test Functions --- #
# --------------------------------- #

def CopyRequired(ifile, ofile):
    # """
    # Copies files required for LICHEM job.
    #
    # Parameters
    # -------
    # ifile : str
    #     Name of the original file.
    # ofile : str
    #     Name of the copied file.
    # Returns
    # -------
    # cmd_printed : str
    #     Formatted version of the attempted copy command (for printing).
    # """
    cmd = "cp {} {}".format(ifile, ofile)
    cmd_printed = "\n{:6}- {}\n".format(' ', cmd)
    subprocess.call(cmd, shell=True)
    return cmd_printed


def QMMMWrapperTest(name, psi4_e, psi4_u, g09_e, g09_u, g16_e, g16_u,
                    nwchem_e, nwchem_u, energy_str, energy_loc, copy_command,
                    xName, rName, cName, QMMMPack, pf_printed):
    # """
    # Test for a QM wrapper on its own or with an MM wrapper.
    #
    # Parameters
    # ----------
    # name : str
    #     Name of the test being run.
    # psi4_e : float
    #     Expected PSI4 energy in the LICHEM log.
    # psi4_u : str
    #     Units for the expected PSI4 energy in the LICHEM log.
    # g09_e : float
    #     Expected Gaussian09 energy in the LICHEM log.
    # g09_u : str
    #     Units for the expected Gaussian09 energy in the LICHEM log.
    # g16_e : float
    #     Expected Gaussian16 energy in the LICHEM log.
    # g16_u : str
    #     Units for the expected Gaussian16 energy in the LICHEM log.
    # nwchem_e : float
    #     Expected NWChem energy in the LICHEM log.
    # nwchem_u : str
    #     Units for the expected NWChem energy in the LICHEM log.
    # energy_str : str
    #     Term to search for within the LICHEM log.
    # energy_loc : int
    #     Location within split energy_str containing the energy value.
    # copy_command : list
    #     A list of functions for copying required files.
    #      Format: [CopyRequired(ifile, ofile)]
    # xName : str
    #     Name of the LICHEM XYZ file.
    # rName : str
    #     Name of the regions file.
    # cName : str
    #     Name of the connectivity file.
    # QMMMPack : str
    #     Name of the package being tested.
    # pf_printed : bool
    #     True if pass/fail has already been printed. This is intended to help
    #      avoid race conditions with SkipSequence() after a KeyboardInterrupt,
    #      resulting in the skip being printed after a test has passed.
    #
    # Returns
    # -------
    # pf_printed : bool
    # """
    global testCt
    global debugMode
    # Skip irrelevant tests
    if ((psi4_e is None) and (QMMMPack == "PSI4")):
        return
    if ((g09_e is None) and (QMMMPack in ("Gaussian", "Gaussian09", "g09"))):
        return
    if ((g16_e is None) and (QMMMPack in ("Gaussian16", "g16"))):
        return
    if ((nwchem_e is None) and (QMMMPack == "NWChem")):
        return
    # Increment for relevant tests
    testCt += 1
    if testCt == 1:
        print(f"{' ':2}Test {testCt}: {name}")
    else:
        print(f"\n{' ':2}Test {testCt}: {name}")
    # Initialize energy as failing
    passEnergy = False
    # Copy necessary files
    if copy_command is not None:
        # Iterate through list
        for command in copy_command:
            cmd_printed = command
            if debugMode is True:
                print(cmd_printed)
    # Run LICHEM
    runC = RunLICHEM(xName, rName, cName)
    QMMMEnergy, savedEnergy, units = RecoverEnergy(energy_str, energy_loc)
    if QMMMPack == "PSI4":
        comparison = CompareEnergy(
                            QMMMPack=QMPack, QMMMEnergy=QMMMEnergy,
                            known=psi4_e, known_units=psi4_u)
    if QMMMPack in ("Gaussian", "Gaussian09", "g09"):
        comparison = CompareEnergy(
                            QMMMPack=QMPack, QMMMEnergy=QMMMEnergy,
                            known=g09_e, known_units=g09_u)
    if QMMMPack in ("Gaussian16", "g16"):
        comparison = CompareEnergy(
                            QMMMPack=QMPack, QMMMEnergy=QMMMEnergy,
                            known=g16_e, known_units=g16_u)
    if QMMMPack == "NWChem":
        comparison = CompareEnergy(
                            QMMMPack=QMPack, QMMMEnergy=QMMMEnergy,
                            known=nwchem_e, known_units=nwchem_u)
    # Check if comparsion == passEnergy == True, or if it's saved_fail
    if comparison is True:
        # CompareEnergy resaves passEnergy as comparison
        pf_printed = PrintPassFail(name+":", comparison,
                                   updateResults, savedEnergy, units)
    else:
        pf_printed = PrintPassFail(name+":", passEnergy,
                                   updateResults, savedEnergy, units)
        print(comparison)
        PrintLICHEMDebug(runC)
    CleanFiles()  # Clean up files
    return pf_printed


def QMMMFreqTest(name, psi4_e, psi4_u, g09_e, g09_u, g16_e, g16_u,
                 nwchem_e, nwchem_u, xName, rName, cName,
                 QMMMPack, pf_printed):
    # """
    # Frequency test for a QM wrapper.
    #
    # Parameters
    # ----------
    # name : str
    #     Name of the test being run.
    # psi4_e : float
    #     Expected PSI4 energy in the LICHEM log.
    # psi4_u : str
    #     Units for the expected PSI4 energy in the LICHEM log.
    # g09_e : float
    #     Expected Gaussian09 energy in the LICHEM log.
    # g09_u : str
    #     Units for the expected Gaussian09 energy in the LICHEM log.
    # g16_e : float
    #     Expected Gaussian16 energy in the LICHEM log.
    # g16_u : str
    #     Units for the expected Gaussian16 energy in the LICHEM log.
    # nwchem_e : float
    #     Expected NWChem energy in the LICHEM log.
    # nwchem_u : str
    #     Units for the expected NWChem energy in the LICHEM log.
    # xName : str
    #     Name of the LICHEM XYZ file.
    # rName : str
    #     Name of the regions file.
    # cName : str
    #     Name of the connectivity file.
    # QMMMPack : str
    #     Name of the package being tested.
    # pf_printed : bool
    #     True if pass/fail has already been printed. This is intended to help
    #      avoid race conditions with SkipSequence() after a KeyboardInterrupt,
    #      resulting in the skip being printed after a test has passed.
    #
    # Returns
    # -------
    # pf_printed : bool
    #
    # Notes
    # -----
    # A few more functions are required for selecting frequencies from the
    # log file, which is why this is a separate function from
    # QMMMWrapperTest().
    # """
    global testCt
    # Skip irrelevant tests
    if ((psi4_e is None) and (QMMMPack == "PSI4")):
        return
    if ((g09_e is None) and (QMMMPack in ("Gaussian", "Gaussian09", "g09"))):
        return
    if ((g16_e is None) and (QMMMPack in ("Gaussian16", "g16"))):
        return
    if ((nwchem_e is None) and (QMMMPack == "NWChem")):
        return
    # Increment for relevant tests
    testCt += 1
    if testCt == 1:
        print(f"{' ':2}Test {testCt}: {name}")
    else:
        print(f"\n{' ':2}Test {testCt}: {name}")
    # Initialize energy as failing
    passEnergy = False
    # Run LICHEM
    runC = RunLICHEM(xName, rName, cName)
    # Get Frequency Info
    QMMMFreqs, units = RecoverFreqs()
    # Find lowest frequency
    try:
        QMMMEnergy = min(QMMMFreqs)
        if QMMMEnergy > 1e100:
            savedEnergy = "Crashed..."
            time.sleep(2)
        else:
            # Save string for devMode printing
            savedEnergy = "Freq:{:3}{}".format(" ", QMMMEnergy)
            # Round for comparison
            QMMMEnergy = round(QMMMEnergy, 0)
    # If list is empty
    except ValueError:
        savedEnergy = "Crashed..."
        # Set to 0 for debugMode
        QMMMEnergy = 0.0
        time.sleep(2)
    if QMMMPack == "PSI4":
        comparison = CompareEnergy(
                            QMMMPack=QMPack, QMMMEnergy=QMMMEnergy,
                            known=psi4_e, known_units=psi4_u)
    if QMMMPack in ("Gaussian", "Gaussian09", "g09"):
        comparison = CompareEnergy(
                            QMMMPack=QMPack, QMMMEnergy=QMMMEnergy,
                            known=g09_e, known_units=g09_u)
    if QMMMPack in ("Gaussian16", "g16"):
        comparison = CompareEnergy(
                            QMMMPack=QMPack, QMMMEnergy=QMMMEnergy,
                            known=g16_e, known_units=g16_u)
    if QMMMPack == "NWChem":
        comparison = CompareEnergy(
                            QMMMPack=QMPack, QMMMEnergy=QMMMEnergy,
                            known=nwchem_e, known_units=nwchem_u)
    # Check if comparsion == passEnergy == True, or if it's saved_fail
    if comparison is True:
        # CompareEnergy resaves passEnergy as comparison
        pf_printed = PrintPassFail(name+":", comparison,
                                   updateResults, savedEnergy, units)
    else:
        pf_printed = PrintPassFail(name+":", passEnergy,
                                   updateResults, savedEnergy, units)
        print(comparison)
        PrintLICHEMDebug(runC)
    CleanFiles()  # Clean up files
    return pf_printed


def MMWrapperTest(name, tinker_e, tinker_u, lammps_e, lammps_u,
                  energy_str, energy_loc, tinkerKey,
                  xName, rName, cName, QMMMPack, pf_printed):
    # """
    # Test for an MM wrapper on its own.
    #
    # Parameters
    # ----------
    # name : str
    #     Name of the test being run.
    # tinker_e : float
    #     Expected PSI4 energy in the LICHEM log.
    # tinker_u : str
    #     Units for the expected PSI4 energy in the LICHEM log.
    # lammps_e : float
    #     Expected Gaussian09 energy in the LICHEM log.
    # lammps_u : str
    #     Units for the expected Gaussian09 energy in the LICHEM log.
    # energy_str : str
    #     Term to search for within the LICHEM log.
    # energy_loc : int
    #     Location within split energy_str containing the energy value.
    # tinkerKey : str
    #     Name of the TINKER key file to copy. (None for LAMMPS-only test.)
    # xName : str
    #     Name of the LICHEM XYZ file.
    # rName : str
    #     Name of the regions file.
    # cName : str
    #     Name of the connectivity file.
    # QMMMPack : str
    #     Name of the package being tested.
    # pf_printed : bool
    #     True if pass/fail has already been printed. This is intended to help
    #      avoid race conditions with SkipSequence() after a KeyboardInterrupt,
    #      resulting in the skip being printed after a test has passed.
    #
    # Returns
    # -------
    # pf_printed : bool
    # """
    global testCt
    global debugMode
    # Skip irrelevant tests
    if ((tinker_e is None) and (QMMMPack in ("TINKER", "TINKER9"))):
        return
    if ((lammps_e is None) and (QMMMPack == "LAMMPS")):
        return
    # Increment for relevant tests
    testCt += 1
    if testCt == 1:
        print(f"{' ':2}Test {testCt}: {name}")
    else:
        print(f"\n{' ':2}Test {testCt}: {name}")
    # Copy key file if it's TINKER
    if QMMMPack in ("TINKER", "TINKER9"):
        cmd_printed = CopyRequired(tinkerKey, "tinker.key")
        if debugMode is True:
            print(cmd_printed)
    # Initialize energy as failing
    passEnergy = False
    runC = RunLICHEM(xName, rName, cName)
    QMMMEnergy, savedEnergy, units = RecoverEnergy(energy_str, energy_loc)
    if QMMMPack in ("TINKER", "TINKER9"):
        comparison = CompareEnergy(
                            QMMMPack=QMPack, QMMMEnergy=QMMMEnergy,
                            known=tinker_e, known_units=tinker_u)
    if QMMMPack == "LAMMPS":
        comparison = CompareEnergy(
                            QMMMPack=QMPack, QMMMEnergy=QMMMEnergy,
                            known=lammps_e, known_units=lammps_u)
    # Check if comparsion == passEnergy == True, or if it's saved_fail
    if comparison is True:
        # CompareEnergy resaves passEnergy as comparison
        pf_printed = PrintPassFail(name+":", comparison,
                                   updateResults, savedEnergy, units)
    else:
        pf_printed = PrintPassFail(name+":", passEnergy,
                                   updateResults, savedEnergy, units)
        print(comparison)
        PrintLICHEMDebug(runC)
    CleanFiles()  # Clean up files
    return pf_printed


# --------------------------------------------------- #
# --- Prepare the Tests and Verify Binaries Exist --- #
# --------------------------------------------------- #

# Print title
print("\n"
      "***************************************************\n"
      "*                                                 *\n"
      "*   LICHEM: Layered Interacting CHEmical Models   *\n"
      "*                                                 *\n"
      "*        Symbiotic Computational Chemistry        *\n"
      "*                                                 *\n"
      "***************************************************\n"
      )

# Read in the arguments
args = get_args()

# Print defaults warning when no command-line arguments are specified
if not len(sys.argv) > 1:
    print("No command-line arguments specified. Running with defaults: \n"
          "  --ncpus [1]\n"
          "  --qm [all detected]\n"
          "  --mm [all detected]\n"
          "  --tests [all applicable]\n\n"
          "Usage info can be printed with './runtests -h' or "
          "'./runtests --help'.\n")

# Set for now, go back and use args.dry = False
dryRun = False    # Only check packages
allTests = False  # Run all tests at once

# Create a set of each QM and MM package to test to avoid duplicates
qm_test_packs = set()
mm_test_packs = set()

# Check the arguments!
# TODO: Go back and remove Makefile line!
if args.verbose:
    updateResults = True  # Bool to print energies to update tests
    forceAll = True       # Bool to force it to do tests even if they will fail
    debugMode = True      # Bool to print LICHEM command used if test fails

# Set the Number of CPUs (default is set as 1)
Ncpus = int(args.ncpus)

if args.dry:
    dryRun = True  # Only check packages

if args.all:
    allTests = True  # Run all tests at once
    if args.qm != 'all' and args.mm != 'all':
        print("Ignoring provided QM/MM wrappers and auto-running with all "
              "available options.\n")
    elif args.qm != 'all':
        print("Ignoring provided QM wrappers and auto-running with all "
              "available options.\n")
    elif args.mm != 'all':
        print("Ignoring provided MM wrappers and auto-running with all "
              "available options.\n")
else:
    # TODO: look for both g09 and g16, and use the one found
    if args.qm:
        # Verify the QM options
        # Check for PSI4
        if 'psi4' in args.qm or 'psi' in args.qm:
            qm_test_packs.add("psi4")
        # Check for Gaussian
        if 'gaussian' in args.qm:
            qm_test_packs.add("g09")
        elif 'g09' in args.qm:
            qm_test_packs.add("g09")
        if 'gaussian16' in args.qm:
            qm_test_packs.add("g16")
        elif 'g16' in args.qm:
            qm_test_packs.add("g16")
        # Check for NWChem
        if 'nwchem' in args.qm:
            qm_test_packs.add("nwchem")
        if 'all' in args.qm:
            qm_test_packs = {"g09", "g16", "nwchem", "psi4"}
        # Check for None assigned
        if len(qm_test_packs) == 0:
            print(f"\nError: QM package name {args.qm} not recognized.\n")
            print(f"Valid QM package options: {QM_wrappers}\n")
            exit(0)
    # TODO: look for both tinker and tinker9
    if args.mm:
        # Verify the MM options
        # Check for Tinker
        if 'tinker' in args.mm:
            mm_test_packs.add("tinker")
        if 'tinker9' in args.mm:
            mm_test_packs.add("tinker9")
        if 'lammps' in args.mm:
            mm_test_packs.add("lammps")
        if 'all' in args.mm:
            mm_test_packs = {"tinker", "tinker9", "lammps"}
        # Check for None assigned
        if len(mm_test_packs) == 0:
            print(f"\nError: MM package name {args.mm} not recognized.\n")
            print(f"Valid MM package options: {MM_wrappers}\n")
            exit(0)
    if (args.qm == 'all') or (args.mm == 'all'):
        print("Assuming -a|--all and using all available options.\n "
              "Either explicitly set -q|--qm and -m|--mm or remove unwanted "
              "executables\n from your path to change this behavior.\n")
        allTests = True

# Make test options lowercase to test against input
test_opts_lower = [x.lower() for x in test_opts]

# Verify that test input matches a test
if args.tests:
    # If not running all tests
    if args.tests != "all":
        # For index, value
        for i, test in enumerate(args.tests):
            if test not in test_opts_lower:
                # Replace common mistakes
                if test == "freq":
                    args.tests[i] = "frequencies"
                # Spaces register as separate tests, but QSM_TS/NEB_TS are
                #  distinct without "TS"
                elif test == "ts":
                    args.tests.remove('ts')
                elif test in ("nebts", "neb-ts", "neb ts", "neb"):
                    args.tests[i] = "neb_ts"
                elif test in ("qsmts", "qsm-ts", "qsm ts", "qsm"):
                    args.tests[i] = "qsm_ts"
                elif test in ("amoebagk", "amoeba-gk", "amoeba gk"):
                    args.tests[i] = "amoeba/gk"
                elif test in ("pbe0tip3p", "pbe0-tip3p", "pbe0 tip3p"):
                    args.tests[i] = "pbe0/tip3p"
                elif test in ("pbe0amoeba", "pbe0-amoeba", "pbe0 amoeba"):
                    args.tests[i] = "pbe0/amoeba"
                elif test in ("dfppseudobonds", "dfp-pseudobonds",
                              "dfp", "pseudobonds", "pseudo"):
                    args.tests[i] = "dfp/pseudobonds"
                else:
                    print(f"\nError: Test '{test}' not understood.\n\n"
                          f" Valid options: {test_opts}\n\n"
                          " Note: The test list is case-insensitive.\n")
                    exit(0)
        # Remove potential duplicates (e.g., from both dfp & pseudo)
        args.tests = set(args.tests)

# Look for LICHEM
LICHEMbin = LocateLICHEM()

# Set up dictionaries for key = QM wrapper, value = binary location
qm_pack_dict = {}
mm_pack_dict = {}

# Identify all available wrappers
if allTests is True:
    # Identify QM wrappers
    print("Checking available QM wrappers:")
    # Search for PSI4
    QMbin = LocateProgram(executable='psi4')
    qm_pack_dict["PSI4"] = QMbin
    print(f" PSI4: {QMbin}")
    # Search for Gaussian
    QMbin = LocateProgram(executable='g09')
    qm_pack_dict["Gaussian09"] = QMbin
    print(f" Gaussian09: {QMbin}")
    # Check if both g09 and g16 co-exist
    if useg16 is False:
        QMbin = LocateProgram(executable='g16')
        if QMbin != "N/A":
            useg16 = True
        qm_pack_dict["Gaussian16"] = QMbin
        print(f" Gaussian16: {QMbin}")
    # Search for NWChem
    QMbin = LocateProgram(executable='nwchem')
    qm_pack_dict["NWChem"] = QMbin
    print(f" NWChem: {QMbin}")
    # Identify MM wrappers
    print("\nChecking available MM wrappers:")
    # Search for TINKER
    MMbin = LocateProgram(executable='analyze')
    mm_pack_dict["TINKER"] = MMbin
    print(f" TINKER: {MMbin}")
    # Search for TINKER9
    MMbin = LocateProgram(executable='tinker9')
    mm_pack_dict["TINKER9"] = MMbin
    print(f" TINKER9: {MMbin}")
    # Search for LAMMPS
    MMbin = LocateProgram(executable='lammps')
    mm_pack_dict["LAMMPS"] = MMbin
    print(f" LAMMPS: {MMbin}\n")
# Search for requested wrappers
else:
    badQM = True
    # Ask if wrapper is in since it's part of a set!
    if "psi4" in qm_test_packs:
        QMPack = "PSI4"
        QMbin = LocateProgram(executable='psi4')
        if QMbin != "N/A":
            badQM = False
            qm_pack_dict[QMPack] = QMbin
    # Search for Gaussian
    if "g09" in qm_test_packs:
        QMPack = "Gaussian09"
        QMbin = LocateProgram(executable='g09')
        if QMbin != "N/A":
            badQM = False
            qm_pack_dict[QMPack] = QMbin
    if "g16" in qm_test_packs:
        QMPack = "Gaussian16"
        QMbin = LocateProgram(executable='g16')
        if QMbin != "N/A":
            useg16 = True
            badQM = False
            qm_pack_dict[QMPack] = QMbin
    if useg16 is True:
        QMPack = "Gaussian16"
    if "nwchem" in qm_test_packs:
        QMPack = "NWChem"
        QMbin = LocateProgram(executable='nwchem')
        if QMbin != "N/A":
            badQM = False
            qm_pack_dict[QMPack] = QMbin
    if badQM is True:
        # Quit with an error
        print("\nError: No QM executables located for any of the\n"
              f" following: {qm_test_packs}.\n")
        exit(0)
    # Check for MM
    badMM = True
    if "tinker" in mm_test_packs:
        MMPack = "TINKER"
        MMbin = LocateProgram(executable='analyze')
        if MMbin != "N/A":
            badMM = False
            mm_pack_dict[MMPack] = MMbin
    # Search for TINKER9
    if "tinker9" in mm_test_packs:
        MMPack = "TINKER9"
        MMbin = LocateProgram(executable='tinker9')
        if MMbin != "N/A":
            badMM = False
            mm_pack_dict[MMPack] = MMbin
    # Search for LAMMPS
    if "lammps" in mm_test_packs:
        MMPack = "LAMMPS"
        MMbin = LocateProgram(executable='lammps')
        if MMbin != "N/A":
            badMM = False
            mm_pack_dict[MMPack] = MMbin
    if badMM is True:
        # Quit with error
        print("\nError: No MM executables located for any of the"
              f" following: {mm_test_packs}.\n")
        exit(0)

if forceAll is False:
    # Remove N/A values from QM_Packs
    qm_del_list = [key for key, value in qm_pack_dict.items()
                   if value == ClrSet.TFail+"N/A"+ClrSet.Reset]
    for x in qm_del_list:
        del qm_pack_dict[x]

    # Remove N/A values from MM_Packs
    mm_del_list = [key for key, value in mm_pack_dict.items()
                   if value == ClrSet.TFail+"N/A"+ClrSet.Reset]
    for x in mm_del_list:
        del mm_pack_dict[x]

# Create lists from the keys of tests to run
QMTests = list(qm_pack_dict.keys())
MMTests = list(mm_pack_dict.keys())

# Find total number of packages to test
total_valid_qm = len(QMTests)
total_valid_mm = len(MMTests)

# Print test settings
print(
    "Settings:\n"
    f" CPUs: {Ncpus}")
# LICHEM binary information
print(f" LICHEM binary: {LICHEMbin}")
# QM wrapper binary information
# Loop through dictionary and start counter at 1
for i, (qm_key, qm_value) in enumerate(qm_pack_dict.items(), 1):
    if total_valid_qm > 1:
        print(
            f" QM package {i}: {qm_key}\n"
            f"  Binary: {qm_value}")
    else:
        print(
            f" QM package: {qm_key}\n"
            f"  Binary: {qm_value}")
# MM wrapper binary information
for i, (mm_key, mm_value) in enumerate(mm_pack_dict.items(), 1):
    if total_valid_mm > 1:
        print(
            f" MM package {i}: {mm_key}\n"
            f"  Binary: {mm_value}\n")
    else:
        print(
            f" MM package: {mm_key}\n"
            f"  Binary: {mm_value}\n")

# # Print test settings
# print(
#     "Settings:\n"
#     f" CPUs: {Ncpus}")
# # LICHEM binary information
# print(f" LICHEM binary: {LICHEMbin}")
# # QM wrapper binary information
# # num_valid = Total number of QM packages with binaries not set to "N/A"
# num_valid = sum(1 for x in qm_pack_dict.values()
#                 if x != ClrSet.TFail+"N/A"+ClrSet.Reset)
# # j = Counter for number of N/As skipped
# j = 0
# # Loop through dictionary and start counter at 1
# for i, (qm_key, qm_value) in enumerate(qm_pack_dict.items(), 1):
#     if qm_value == ClrSet.TFail+"N/A"+ClrSet.Reset:
#         j += 1  # Update N/A counter
#         continue  # Exit the for loop for N/A
#     elif num_valid > 1:
#         print(
#             f" QM package {i-j}: {qm_key}\n"
#             f"  Binary: {qm_value}")
#     else:
#         print(
#             f" QM package: {qm_key}\n"
#             f"  Binary: {qm_value}")
# # MM wrapper binary information
# num_valid = sum(1 for x in mm_pack_dict.values()
#                 if x != ClrSet.TFail+"N/A"+ClrSet.Reset)
# j = 0
# for i, (mm_key, mm_value) in enumerate(mm_pack_dict.items(), 1):
#     if mm_value == ClrSet.TFail+"N/A"+ClrSet.Reset:
#         j += 1
#         continue
#     elif num_valid > 1:
#         print(
#             f" MM package {i-j}: {mm_key}\n"
#             f"  Binary: {mm_value}\n")
#     else:
#         print(
#             f" MM package: {mm_key}\n"
#             f"  Binary: {mm_value}\n")

# Mode information
if allTests is True:
    if forceAll is True:
        print(" Mode: Development")
    else:
        print(" Mode: All tests")
if dryRun is True:
    print(" Mode: Dry run\n")

# Escape for dry runs
if dryRun is True:
    # Quit without an error
    print("Dry run completed.\n")
    exit(0)

# Escape if binaries not found when running specific tests
if (((QMbin == "N/A") or (MMbin == "N/A")) and (allTests is False)):
    # Quit with an error
    print("Error: Missing binaries.\n")
    exit(0)

print(
    "***************************************************\n\n"
    "Running LICHEM tests...\n")

# # TODO: FIXME
# # Make a list of tests
# QMTests = []
# MMTests = []
# if allTests is True:
#     # Safely add PSI4
#     cmd = "which psi4"
#     try:
#         # Run PSI4 tests
#         packBin = subprocess.check_output(cmd, shell=True)
#         QMTests.append("PSI4")
#     except subprocess.CalledProcessError:
#         # Skip tests that will fail
#         if forceAll is True:
#             QMTests.append("PSI4")
#     # Safely add Gaussian
#     cmd = "which g09"
#     try:
#         # Run Gaussian tests
#         packBin = subprocess.check_output(cmd, shell=True)
#         QMTests.append("Gaussian09")
#     except subprocess.CalledProcessError:
#         cmd = "which g16"
#         try:
#             packBin = subprocess.check_output(cmd, shell=True)
#             QMTests.append("g16")
#         except subprocess.CalledProcessError:
#             # Skip tests that will fail
#             if forceAll is True:
#                 # QMTests.append("Gaussian")
#                 # Instead try g16 if no other Gaussian found
#                 QMTests.append("g16")
#     # Safely add NWChem
#     cmd = "which nwchem"
#     try:
#         # Run NWChem tests
#         packBin = subprocess.check_output(cmd, shell=True)
#         QMTests.append("NWChem")
#     except subprocess.CalledProcessError:
#         # Skip tests that will fail
#         if forceAll is True:
#             QMTests.append("NWChem")
#     # Safely add TINKER
#     cmd = "which analyze"
#     try:
#         # Run TINKER tests
#         packBin = subprocess.check_output(cmd, shell=True)
#         MMTests.append("TINKER")
#     except subprocess.CalledProcessError:
#         # Skip tests that will fail
#         if forceAll is True:
#             MMTests.append("TINKER")
#     # Safely add lammps
#     cmd = "which lammps"
#     try:
#         # Run LAMMPS tests
#         packBin = subprocess.check_output(cmd, shell=True)
#         MMTests.append("LAMMPS")
#     except subprocess.CalledProcessError:
#         # Skip tests that will fail
#         if forceAll is True:
#             MMTests.append("LAMMPS")
# else:
#     # Add only the specified packages
#     QMTests = list(qm_pack_dict.keys())
#     MMTests = list(mm_pack_dict.keys())

# NB: Tests are in the following order:
#     1) HF energy
#     2) PBE0 energy
#     3) CCSD energy
#     4) PM6 energy
#     5) Frequencies
#     6) NEB TS energy
#     7) QSM TS energy
#     8) TIP3P energy
#     9) AMOEBA/GK energy
#    10) PBE0/TIP3P energy
#    11) PBE0/AMOEBA energy
#    12) DFP/Pseudobonds

# --------------------- #
# --- Run the Tests --- #
# --------------------- #

# Loop over tests
for QMPack in QMTests:
    for MMPack in MMTests:

        # Set path based on packages
        dirPath = ""
        if (QMPack == "PSI4"):
            dirPath += "PSI4_"
        if QMPack in ("Gaussian", "g09", "g16", "Gaussian09", "Gaussian16"):
            dirPath += "Gau_"
        # elif (QMPack == "g16"):
        #     dirPath += "G16_"
        if (QMPack == "NWChem"):
            dirPath += "NWChem_"
        # TODO: Fix for Tinker9 (and against 8.4+ for anglep)
        dirPath += MMPack
        dirPath += "/"

        # Change directory
        os.chdir(dirPath)

        # Fix the regions files (if the specific one exits) to match program
        #  version requested
        revertMode = False
        for filename in regions_files:
            # Check that the file exists for sed
            if os.path.isfile("./"+filename):
                # Gaussian
                if QMPack in ("g16", "Gaussian16"):
                    PrepRegions("QM_type", "Gaussian", "g16", filename)
                # Tinker version
                if MMPack == "TINKER9":
                    PrepRegions("MM_type", "TINKER", "TINKER9", filename)

        # Start printing results
        print(f"{QMPack}/{MMPack} results:")

        # Only run each test if requested
        # Note: "if ('hf' or 'all') in args.tests" does not work!
        if 'hf' in args.tests or 'all' in args.tests:
            # Start each try block by declaring that the Pass/Fail
            #  report for that test hasn't yet printed.
            pf_printed = False
            # Use try/except to control Cntrl+C behavior with
            #  SkipSequence().
            try:
                # Gaussian and Psi4 Only
                # Check HF energy
                pf_printed = QMMMWrapperTest(
                    name="HF Energy",
                    psi4_e=round(-4136.93039814/har2eV, 3), psi4_u="a.u.",
                    g09_e=round(-4136.9317704519/har2eV, 3), g09_u="a.u.",
                    g16_e=round(-4136.9317704519/har2eV, 3), g16_u="a.u.",
                    nwchem_e=None, nwchem_u=None,
                    energy_str="QM energy:", energy_loc=2,
                    copy_command=None,
                    xName="waterdimer.xyz", rName="hfreg.inp",
                    cName="watercon.inp",
                    QMMMPack=QMPack, pf_printed=pf_printed)
            except KeyboardInterrupt:
                SkipSequence("HF energy:", pf_printed)

        if 'pbe0' in args.tests or 'all' in args.tests:
            pf_printed = False
            try:
                # Check DFT energy
                pf_printed = QMMMWrapperTest(
                    name="PBE0 Energy",
                    psi4_e=round(-4154.16836599/har2eV, 3), psi4_u="a.u.",
                    g09_e=round(-4154.1676114324/har2eV, 3), g09_u="a.u.",
                    g16_e=round(-4154.1676114324/har2eV, 3), g16_u="a.u.",
                    nwchem_e=round(-4154.1683939169/har2eV, 3),
                    nwchem_u="a.u.",
                    energy_str="QM energy:", energy_loc=2,
                    copy_command=None,
                    xName="waterdimer.xyz", rName="pbereg.inp",
                    cName="watercon.inp",
                    QMMMPack=QMPack, pf_printed=pf_printed)
            except KeyboardInterrupt:
                SkipSequence("PBE0 energy:", pf_printed)

        if 'ccsd' in args.tests or 'all' in args.tests:
            pf_printed = False
            try:
                # PSI4-only test
                pf_printed = QMMMWrapperTest(
                    name="CCSD Energy",
                    psi4_e=round(-4147.7304837/har2eV, 3), psi4_u="a.u.",
                    g09_e=None, g09_u=None,
                    g16_e=None, g16_u=None,
                    nwchem_e=None, nwchem_u=None,
                    energy_str="QM energy:", energy_loc=2,
                    copy_command=None,
                    xName="waterdimer.xyz", rName="ccsdreg.inp",
                    cName="watercon.inp",
                    QMMMPack=QMPack, pf_printed=pf_printed)
            except KeyboardInterrupt:
                SkipSequence("CCSD energy:", pf_printed)

        if 'pm6' in args.tests or 'all' in args.tests:
            pf_printed = False
            try:
                # Gaussian-only test
                pf_printed = QMMMWrapperTest(
                    name="PM6 Energy",
                    psi4_e=None, psi4_u=None,
                    g09_e=round(-4.8623027634995/har2eV, 3), g09_u="a.u.",
                    g16_e=round(-4.8623027634995/har2eV, 3), g16_u="a.u.",
                    nwchem_e=None, nwchem_u=None,
                    energy_str="QM energy:", energy_loc=2,
                    copy_command=None,
                    xName="waterdimer.xyz", rName="pm6reg.inp",
                    cName="watercon.inp",
                    QMMMPack=QMPack, pf_printed=pf_printed)
            except KeyboardInterrupt:
                SkipSequence("PM6 energy:", pf_printed)

        if 'frequencies' in args.tests or 'all' in args.tests:
            pf_printed = False
            try:
                # Gaussian-only test
                pf_printed = QMMMFreqTest(
                    name="Frequencies",
                    psi4_e=round(-40.507339, 0), psi4_u="cm^-1",
                    g09_e=round(5.33078515, 0), g09_u="cm^-1",
                    g16_e=round(5.33078515, 0), g16_u="cm^-1",
                    nwchem_e=round(-31.769945, 0), nwchem_u="cm^-1",
                    xName="methfluor.xyz", rName="freqreg.inp",
                    cName="methflcon.inp",
                    QMMMPack=QMPack, pf_printed=pf_printed)
            except KeyboardInterrupt:
                SkipSequence("Frequencies", pf_printed)

        if 'neb_ts' in args.tests or 'all' in args.tests:
            pf_printed = False
            try:
                pf_printed = QMMMWrapperTest(
                    name="NEB TS Energy",
                    psi4_e=round(0.39581219003957, 3), psi4_u="eV",
                    g09_e=round(0.39582668467028, 3), g09_u="eV",
                    g16_e=round(0.39582668467028, 3), g16_u="eV",
                    nwchem_e=round(-0.39582668467028, 3), nwchem_u="eV",
                    energy_str="Forward barrier", energy_loc=3,
                    copy_command=[CopyRequired(ifile="methflbeads.xyz",
                                               ofile="BeadStartStruct.xyz")],
                    xName="methfluor.xyz", rName="nebreg.inp",
                    cName="methflcon.inp",
                    QMMMPack=QMPack, pf_printed=pf_printed)
            except KeyboardInterrupt:
                SkipSequence("NEB Forward barrier:", pf_printed)

        if 'qsm_ts' in args.tests or 'all' in args.tests:
            pf_printed = False
            try:
                pf_printed = QMMMWrapperTest(
                    name="QSM TS Energy",
                    psi4_e=round(-239.27686429, 3), psi4_u="eV",
                    g09_e=round(-239.27681732, 3), g09_u="eV",
                    g16_e=round(-239.27681732, 3), g16_u="eV",
                    nwchem_e=round(-239.27681732, 3), nwchem_u="eV",
                    energy_str="TS Energy", energy_loc=4,
                    copy_command=[CopyRequired(ifile="methflbeads.xyz",
                                               ofile="BeadStartStruct.xyz")],
                    xName="methfluor.xyz", rName="qsmreg.inp",
                    cName="methflcon.inp",
                    QMMMPack=QMPack, pf_printed=pf_printed)
            except KeyboardInterrupt:
                SkipSequence("QSM TS energy:", pf_printed)

        # TINKER-only, but require key to be copied
        # For LAMMPS-only, set tinkerKey argument to 'None'
        if 'tip3p' in args.tests or 'all' in args.tests:
            pf_printed = False
            try:
                pf_printed = MMWrapperTest(
                    name="TIP3P Energy",
                    tinker_e=round(-0.2596903536223/har2eV, 3),
                    tinker_u="a.u.",
                    lammps_e=None, lammps_u=None,
                    energy_str="MM energy:", energy_loc=2,
                    tinkerKey="pchrg.key",
                    xName="waterdimer.xyz", rName="mmreg.inp",
                    cName="watercon.inp",
                    QMMMPack=MMPack, pf_printed=pf_printed)
            except KeyboardInterrupt:
                SkipSequence("TIP3P energy:", pf_printed)

        if 'amoeba/gk' in args.tests or 'all' in args.tests:
            pf_printed = False
            try:
                pf_printed = MMWrapperTest(
                    name="AMOEBA/GK Energy",
                    tinker_e=round(-0.0095434448906, 3), tinker_u="a.u.",
                    lammps_e=None, lammps_u=None,
                    energy_str="MM energy:", energy_loc=2,
                    tinkerKey="pol.key",
                    xName="waterdimer.xyz", rName="mmreg.inp",
                    cName="watercon.inp",
                    QMMMPack=MMPack, pf_printed=pf_printed)
            except KeyboardInterrupt:
                SkipSequence("TIP3P energy:", pf_printed)

        if 'pbe0/tip3p' in args.tests or 'all' in args.tests:
            pf_printed = False
            try:
                pf_printed = QMMMWrapperTest(
                    name="PBE0/TIP3P Energy",
                    psi4_e=round(-2077.2021947277/har2eV, 3), psi4_u="a.u.",
                    g09_e=round(-76.34642184669, 3), g09_u="a.u.",
                    g16_e=round(-76.34642184669, 3), g16_u="a.u.",
                    nwchem_e=round(-2077.2022117306/har2eV, 3),
                    nwchem_u="a.u.",
                    energy_str="QMMM energy:", energy_loc=2,
                    copy_command=[CopyRequired(ifile="pchrg.key",
                                               ofile="tinker.key")],
                    xName="waterdimer.xyz", rName="pchrgreg.inp",
                    cName="watercon.inp",
                    QMMMPack=QMPack, pf_printed=pf_printed)
            except KeyboardInterrupt:
                SkipSequence("PBE0/TIP3P energy:", pf_printed)

        if 'pbe0/amoeba' in args.tests or 'all' in args.tests:
            pf_printed = False
            try:
                pf_printed = QMMMWrapperTest(
                    name="PBE0/AMOEBA Energy",
                    psi4_e=round(-2077.1114201829/har2eV, 3), psi4_u="a.u.",
                    g09_e=round(-2077.1090319595/har2eV, 3), g09_u="a.u.",
                    g16_e=round(-2077.1090319595/har2eV, 3), g16_u="a.u.",
                    nwchem_e=round(-2077.1094168459/har2eV, 3),
                    nwchem_u="a.u.",
                    energy_str="QMMM energy:", energy_loc=2,
                    copy_command=[CopyRequired(ifile="pol.key",
                                               ofile="tinker.key")],
                    xName="waterdimer.xyz", rName="polreg.inp",
                    cName="watercon.inp",
                    QMMMPack=QMPack, pf_printed=pf_printed)
            except KeyboardInterrupt:
                SkipSequence("PBE0/AMOEBA energy:", pf_printed)

        if 'dfp/pseudobonds' in args.tests or 'all' in args.tests:
            pf_printed = False
            try:
                # Gaussian and NWChem only
                pf_printed = QMMMWrapperTest(
                    name="DFP/Pseudobonds", QMMMPack=QMPack,
                    # Copy two files!
                    copy_command=[CopyRequired(ifile="pbopt.key",
                                               ofile="tinker.key"),
                                  CopyRequired(ifile="pbbasis.txt",
                                               ofile="BASIS")],
                    xName="alkyl.xyz", rName="pboptreg.inp",
                    cName="alkcon.inp",
                    energy_str="Opt. step: 2", energy_loc=6,
                    psi4_e=None, psi4_u=None,
                    g09_e=round(-3015.0548490566/har2eV, 3), g09_u="a.u.",
                    g16_e=round(-3015.0548490566/har2eV, 3), g16_u="a.u.",
                    nwchem_e=round(-3015.2278310975/har2eV, 3),
                    nwchem_u="a.u.",
                    pf_printed=pf_printed)
            except KeyboardInterrupt:
                SkipSequence("DFP/Pseudobonds:", pf_printed)

        # Revert the regions files
        revertMode = True
        for filename in regions_files:
            # Check that the file exists for sed
            if os.path.isfile("./"+filename):
                # Gaussian
                if QMPack in ("Gaussian16", "g16"):
                    PrepRegions("QM_type", "g16", "Gaussian", filename)
                # Tinker version
                if MMPack == "TINKER9":
                    PrepRegions("MM_type", "TINKER9", "TINKER", filename)

        # Print blank line and change directory
        print()
        os.chdir("../")
        # Reset test counter
        testCt = 0

# ----------------------------- #
# --- Print Test Statistics --- #
# ----------------------------- #

# Start printing the statistics
print("***************************************************\n\n"
      "Statistics:\n"
      f" Tests {ClrSet.TPass}Passed{ClrSet.Reset}:  {passCt}\n"
      f" Tests {ClrSet.TFail}Failed{ClrSet.Reset}:  {failCt}\n"
      f" Tests {ClrSet.TSkip}Skipped{ClrSet.Reset}: {skipCt}\n")

# Stop timer
endTime = time.time()
totalTime = (endTime-startTime)

# Find the correct units
# Maybe update to datetime or use divmod in the future?
timeUnits = "seconds"
if (totalTime > 60):
    totalTime /= 60.0
    timeUnits = "minutes"
    if (totalTime > 60):
        totalTime /= 60.0
        timeUnits = "hours"
        if (totalTime > 24):
            totalTime /= 24.0
            timeUnits = "days"

# Finish printing the statistics
print(f" Total run time: {round(totalTime, 2):2f} {timeUnits}\n\n"
      "***************************************************\n\n"
      "Done.\n\n")

# Quit
exit(0)
