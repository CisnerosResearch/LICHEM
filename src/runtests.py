###################################################
#                                                 #
#   LICHEM: Layered Interacting CHEmical Models   #
#                                                 #
#        Symbiotic Computational Chemistry        #
#                                                 #
###################################################

# LICHEM semi-automated test suite

###
#  Usage:
#
#    user:$ ./runtests Ncpus All
#      or
#    user:$ ./runtests Ncpus QMPackage MMPackage
#     or
#    user:$ ./runtests Ncpus QMPackage MMPackage dry
####

# TODO:
# - Support "named test" option (aka, functionalize tests)
# - Rework the if/if process for QM/MM engines tests (shorten code)
# - Figure out where QSM TS tinker.key copy attempt is...

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

# Development settings
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
test_opts = ["HF", "PBE0", "CCSD", "PM6", "Frequencies", "NEB_TS", "TIP3P",
             "AMOEBA/GK", "PBE0/TIP3P", "PBE0/AMOEBA", "DFP/Pseudobonds"]


# ---------------------- #
# --- Define Classes --- #
# ---------------------- #


class ClrSet:
    """
    Define the colors for messages printed to the console.
    """
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
    """
    Read in provided arguments and define command-line usage.

    Returns
    -------
    parser.parse_args() : argparse.Namespace
        Program details specified by command-line arguments.
    """
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
    wrappers.add_argument("-q", "--qm", nargs="?", default="all",
                          type=str.lower,
                          help=("The QM wrapper to test.\n"
                                f"Options: {QM_wrappers}\n"
                                "  Note: Gaussian and g09 are equivalent."))
    wrappers.add_argument("-m", "--mm", nargs="?", default="all",
                          type=str.lower,
                          help=("The MM wrapper to test.\n"
                                f"Options: {MM_wrappers}"))
    return parser.parse_args()


def LocateLICHEM():
    """
    Check that the LICHEM binary is in the user's PATH.

    Returns
    -------
    saved_lichem_bin : str
        The path to the LICHEM binary.
    """
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
    """
    Check that the program binary is in the user's PATH.

    Returns
    -------
    saved_bin : str
        The path to the program binary.
    """
    global allTests
    global useg16
    cmd = "which {}".format(executable)
    # Find the binary and save the colored and uncolored output
    try:
        prog_bin = subprocess.check_output(cmd, shell=True)
        saved_bin = prog_bin.decode('utf-8').strip()  # Uncolored
        prog_bin = ClrSet.TPass+prog_bin.decode('utf-8').strip()+ClrSet.Reset
    except subprocess.CalledProcessError:
        # If g09 fails, try g16
        if executable == 'g09':
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
    """
    Replace the current keyword argument in a region file with another.

    Parameters
    ----------
    keyword : str
        The keyword to search for in the LICHEM regions file.
    d_val : str
        The default value for the keyword.
    u_val : str
        The updated value for the keywork.
    file : str
        The name of the file to operate on.
    """
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
    return


def RunLICHEM(xName, rName, cName):
    """
    Call a LICHEM instance for testing.

    Parameters
    ----------
    xName : str
        Name of the LICHEM XYZ file.
    rName : str
        Name of the regions file.
    cName : str
        Name of the connectivity file.

    Returns
    -------
    cmd_printed : str
        Formatted version of the attempted LICHEM command (for printing).
    """
    global Ncpus
    # -l captures stdout, 2>&1 captures stderr
    cmd = ("lichem -n {n} -x {x} -r {r} -c {c} -o trash.xyz "
           "-l tests.out 2>&1").format(n=Ncpus, x=xName, r=rName, c=cName)
    subprocess.call(cmd, shell=True)  # Run calculations
    # Save a terminal-formatted version of the command
    #  \\\n escapes a printed backslash and then prints newline
    cmd_printed = ("lichem -n {n} -x {x} -r {r} \\\n{sp:9}-c {c} "
                   "-o trash.xyz -l tests.out 2>&1").format(n=Ncpus,
                   x=xName, r=rName, sp=' ', c=cName)
    return cmd_printed


def CleanFiles():
    """
    Clean out various files after tests are performed.
    """
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
    """
    Recover the energy from the LICHEM output.

    Parameters
    ----------
    txtLabel : str
        Descriptive label in the LICHEM output used to identify energy line.
    itemNum : int
        The split index in the identified line to pull the energy from.

    Returns
    -------
    finalEnergy : float
        The energy from the LICHEM output, rounded to 3 decimal places.
        If failed, the energy is set to 0.0.
    savedResult : str
        Either the unrounded finalEnergy or a crashed message.
    units : str
        The unit following the finalEnergy in the output.
    """
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
    """
    Recover a list of frequencies from the LICHEM output.

    Returns
    -------
    freqList : lst
        Frequencies from LICHEM output.
    """
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
    """
    Checks the calculated energy against a known value.

    Parameters
    ----------
    QMMMPack : str
        Name of the package being tested.
    QMMMEnergy : float
        Value from the LICHEM log file.
    known : float
        Expected value for the program.
    known_units : str
        Units for the expected value.
    """
    expected_energy = known
    if (QMMMEnergy == expected_energy):
        passEnergy = True
        return passEnergy
    else:
        saved_fail = SaveFailure(QMMMEnergy, expected_energy)
        return saved_fail


def SaveFailure(calc, expected):
    """
    Print the expected and test-calculated energies.

    Parameters
    ----------
    calc : float
        QM/MM/QMMM energy from the test LICHEM output.
    expected : float
        Energy expected from output.

    Returns
    -------
    saved_val : str
        Expected and calculated values formatted for console printing.
    """
    # The function input order matches the if statements!!!
    # Note: You don't need to add saved values because all tests of a pairing
    #       will run, so the if statements are independent.
    saved_val = "\n{:4}---> Expected: {}, Calculated: {}".format(
        " ", expected, calc)
    return saved_val


def PrintLICHEMDebug(cmd_printed):
    """
    Prints the LICHEM command that was run to the console.

    Parameters
    -------
    cmd_printed : str
        Formatted version of the attempted LICHEM command.
    """
    global debugMode
    if debugMode is True:
        print(f"\n{' ':6}Tried to execute:\n{' ':7}{cmd_printed}")
    return


def PrintPassFail(tName, testPass, updateResults, enVal, units):
    """
    Write pass/fail message to the console with coloration and timing.

    Parameters
    ----------
    tName: str
        Name of the test being run.
    testPass : bool
        True if test passed.
    updateResults : bool
        True to print energies to update tests (for developers).
    enVal : float
        The computed energy value.
    units : str
        The units of the computed energy value.
    """
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
                  f"{ClrSet.TPass}Pass{ClrSet.Reset}, {runTime}, "+
                  f"{enVal} {units}")
            passCt += 1
        else:
            print(f"{' ':4}{tName:<{TTxtLen-5}} " +
                  f"{ClrSet.TFail}Fail{ClrSet.Reset}, {runTime}, "+
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
    """
    Color the skip message printed to the console.

    Parameters
    ----------
    tName: str
        Name of the test being run.
    pf_printed : bool
        True if pass/fail has already been printed. This is intended to help
         avoid race conditions with SkipSequence() after a KeyboardInterrupt,
         resulting in the skip being printed after a test has passed.
    """
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
    """
    Either skip or quit a test based off user input.

    Parameters
    ----------
    tName: str
        Name of the test being run.
    pf_printed : bool
        True if pass/fail has already been printed. This is intended to help
         avoid race conditions with SkipSequence() after a KeyboardInterrupt,
         resulting in the skip being printed after a test has passed.
    """
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
    """
    Copies files required for LICHEM job.

    Parameters
    -------
    ifile : str
        Name of the original file.
    ofile : str
        Name of the copied file.
    Returns
    -------
    cmd_printed : str
        Formatted version of the attempted copy command (for printing).
    """
    cmd = "cp {} {}".format(ifile, ofile)
    cmd_printed = "\n{:6}- {}\n".format(' ', cmd)
    subprocess.call(cmd, shell=True)
    return cmd_printed


# Parameters
# ----------
# name : str
#     Name of the test being run.
# QMMMPack : str
#     Name of the package being tested.
# pf_printed : bool
#     True if pass/fail has already been printed. This is intended to help
#      avoid race conditions with SkipSequence() after a KeyboardInterrupt,
#      resulting in the skip being printed after a test has passed.
# Returns
# -------
# pf_printed: bool

# def Test1(name, QMMMPack, pf_printed):
#     """
#     Test for HF energy.
#     """
#     global testCt
#     testCt += 1
#     print(f"{' ':2}Test {testCt}: {name}")
#     # Initialize energy as failing
#     passEnergy = False
#     runC = RunLICHEM("waterdimer.xyz", "hfreg.inp", "watercon.inp")
#     QMMMEnergy, savedEnergy, units = RecoverEnergy("QM energy:", 2)
#     if QMMMPack == "PSI4":
#         comparison = CompareEnergy(
#                             QMMMPack=QMPack, QMMMEnergy=QMMMEnergy,
#                             known=round(-4136.93039814/har2eV, 3),
#                             known_units='a.u.')
#     if QMMMPack in ("Gaussian", "Gaussian09", "g09"):
#         comparison = CompareEnergy(
#                             QMMMPack=QMPack, QMMMEnergy=QMMMEnergy,
#                             known=round(-4136.9317704519/har2eV, 3),
#                             known_units='a.u.')
#     if QMMMPack in ("Gaussian16", "g16"):
#         comparison = CompareEnergy(
#                             QMMMPack=QMPack, QMMMEnergy=QMMMEnergy,
#                             known=round(-4136.9317704519/har2eV, 3),
#                             known_units='a.u.')
#     # Check if comparsion == passEnergy == True, or if it's saved_fail
#     if comparison is True:
#         # CompareEnergy resaves passEnergy as comparison
#         pf_printed = PrintPassFail(name+":", comparison,
#                                    updateResults, savedEnergy, units)
#     else:
#         pf_printed = PrintPassFail(name+":", passEnergy,
#                                    updateResults, savedEnergy, units)
#         print(comparison)
#         PrintLICHEMDebug(runC)
#     CleanFiles()  # Clean up files
#     return pf_printed
#
#
# def Test2(name, QMMMPack, pf_printed):
#     """
#     Test for PBE0 energy.
#     """
#     global testCt
#     testCt += 1
#     print(f"\n{' ':2}Test {testCt}: {name}")
#     # Initialize energy as failing
#     passEnergy = False
#     runC = RunLICHEM("waterdimer.xyz", "pbereg.inp", "watercon.inp")
#     QMMMEnergy, savedEnergy, units = RecoverEnergy("QM energy:", 2)
#     if QMMMPack == "PSI4":
#         comparison = CompareEnergy(
#                             QMMMPack=QMPack, QMMMEnergy=QMMMEnergy,
#                             known=round(-4154.16836599/har2eV, 3),
#                             known_units='a.u.')
#     if QMMMPack in ("Gaussian", "Gaussian09", "g09"):
#         comparison = CompareEnergy(
#                             QMMMPack=QMPack, QMMMEnergy=QMMMEnergy,
#                             known=round(-4154.1676114324/har2eV, 3),
#                             known_units='a.u.')
#     if QMMMPack in ("Gaussian16", "g16"):
#         comparison = CompareEnergy(
#                             QMMMPack=QMPack, QMMMEnergy=QMMMEnergy,
#                             known=round(-4154.1676114324/har2eV, 3),
#                             known_units='a.u.')
#     if QMMMPack == "NWChem":
#         comparison = CompareEnergy(
#                             QMMMPack=QMPack, QMMMEnergy=QMMMEnergy,
#                             known=round(-4154.1683939169/har2eV, 3),
#                             known_units='a.u.')
#     # Check if comparsion == passEnergy == True, or if it's saved_fail
#     if comparison is True:
#         # CompareEnergy resaves passEnergy as comparison
#         pf_printed = PrintPassFail(name+":", comparison,
#                                    updateResults, savedEnergy, units)
#     else:
#         pf_printed = PrintPassFail(name+":", passEnergy,
#                                    updateResults, savedEnergy, units)
#         print(comparison)
#         PrintLICHEMDebug(runC)
#     CleanFiles()  # Clean up files
#     return pf_printed

# pf_printed = QMMMWrapperTest(
#     name="DFP/Pseudobonds", QMMMPack=QMPack,
#     # Copy two files!
#     copy_command=[CopyRequired(ifile="pol.key",
#                            ofile="tinker.key"),
#                   CopyRequired(ifile="pbbasis.txt",
#                            ofile="BASIS")],
#     xName="alkyl.xyz", rName="pboptreg.inp",
#     cName="alkcon.inp",
#     energy_str="Opt. step: 2", energy_loc=6,
#     psi4_e=None, psi4_u=None,
#     g09_e=round(-3015.0548490566/har2eV, 3), g09_u="a.u.",
#     g16_e=round(-3015.0548490566/har2eV, 3), g16_u="a.u.",
#     nwchem_e=round(-3015.2278310975/har2eV, 3), nwchem_u="a.u.",
#     pf_printed=pf_printed)

# Parameters
# ----------
# name : str
#     Name of the test being run.
# QMMMPack : str
#     Name of the package being tested.
# pf_printed : bool
#     True if pass/fail has already been printed. This is intended to help
#      avoid race conditions with SkipSequence() after a KeyboardInterrupt,
#      resulting in the skip being printed after a test has passed.
# Returns
# -------
# pf_printed: bool

def QMMMWrapperTest(name, psi4_e, psi4_u, g09_e, g09_u, g16_e, g16_u,
                 nwchem_e, nwchem_u, energy_str, energy_loc, copy_command,
                 xName, rName, cName, QMMMPack, pf_printed):
    """
    Test for a QM wrapper on its own or with an MM wrapper.

    Parameters
    ----------
    name : str
        Name of the test being run.
    psi4_e : float
        Expected PSI4 energy in the LICHEM log.
    psi4_u : str
        Units for the expected PSI4 energy in the LICHEM log.
    g09_e : float
        Expected Gaussian09 energy in the LICHEM log.
    g09_u : str
        Units for the expected Gaussian09 energy in the LICHEM log.
    g16_e : float
        Expected Gaussian16 energy in the LICHEM log.
    g16_u : str
        Units for the expected Gaussian16 energy in the LICHEM log.
    nwchem_e : float
        Expected NWChem energy in the LICHEM log.
    nwchem_u : str
        Units for the expected NWChem energy in the LICHEM log.
    energy_str : str
        Term to search for within the LICHEM log.
    energy_loc : int
        Location within split energy_str containing the energy value.
    copy_command : list
        A list of functions for copying required files.
         Format: [CopyRequired(ifile, ofile)]
    xName : str
        Name of the LICHEM XYZ file.
    rName : str
        Name of the regions file.
    cName : str
        Name of the connectivity file.
    QMMMPack : str
        Name of the package being tested.
    pf_printed : bool
        True if pass/fail has already been printed. This is intended to help
         avoid race conditions with SkipSequence() after a KeyboardInterrupt,
         resulting in the skip being printed after a test has passed.

    Returns
    -------
    pf_printed : bool
    """
    global testCt
    global debugMode
    # Skip irrelevant tests
    if ((psi4_e == None) and (QMMMPack == "PSI4")):
        return
    if ((g09_e == None) and (QMMMPack in ("Gaussian", "Gaussian09", "g09"))):
        return
    if ((g16_e == None) and (QMMMPack in ("Gaussian16", "g16"))):
        return
    if ((nwchem_e == None) and (QMMMPack == "NWChem")):
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
    if copy_command != None:
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
    """
    Frequency test for a QM wrapper.

    Parameters
    ----------
    name : str
        Name of the test being run.
    psi4_e : float
        Expected PSI4 energy in the LICHEM log.
    psi4_u : str
        Units for the expected PSI4 energy in the LICHEM log.
    g09_e : float
        Expected Gaussian09 energy in the LICHEM log.
    g09_u : str
        Units for the expected Gaussian09 energy in the LICHEM log.
    g16_e : float
        Expected Gaussian16 energy in the LICHEM log.
    g16_u : str
        Units for the expected Gaussian16 energy in the LICHEM log.
    nwchem_e : float
        Expected NWChem energy in the LICHEM log.
    nwchem_u : str
        Units for the expected NWChem energy in the LICHEM log.
    xName : str
        Name of the LICHEM XYZ file.
    rName : str
        Name of the regions file.
    cName : str
        Name of the connectivity file.
    QMMMPack : str
        Name of the package being tested.
    pf_printed : bool
        True if pass/fail has already been printed. This is intended to help
         avoid race conditions with SkipSequence() after a KeyboardInterrupt,
         resulting in the skip being printed after a test has passed.

    Returns
    -------
    pf_printed : bool

    Notes
    -----
    A few more functions are required for selecting frequencies from the
    log file, which is why this is a separate function from QMMMWrapperTest(),
    """
    global testCt
    # Skip irrelevant tests
    if ((psi4_e == None) and (QMMMPack == "PSI4")):
        return
    if ((g09_e == None) and (QMMMPack in ("Gaussian", "Gaussian09", "g09"))):
        return
    if ((g16_e == None) and (QMMMPack in ("Gaussian16", "g16"))):
        return
    if ((nwchem_e == None) and (QMMMPack == "NWChem")):
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
    """
    Test for an MM wrapper on its own.

    Parameters
    ----------
    name : str
        Name of the test being run.
    tinker_e : float
        Expected PSI4 energy in the LICHEM log.
    tinker_u : str
        Units for the expected PSI4 energy in the LICHEM log.
    lammps_e : float
        Expected Gaussian09 energy in the LICHEM log.
    lammps_u : str
        Units for the expected Gaussian09 energy in the LICHEM log.
    energy_str : str
        Term to search for within the LICHEM log.
    energy_loc : int
        Location within split energy_str containing the energy value.
    tinkerKey : str
        Name of the TINKER key file to copy. (None for LAMMPS-only test.)
    xName : str
        Name of the LICHEM XYZ file.
    rName : str
        Name of the regions file.
    cName : str
        Name of the connectivity file.
    QMMMPack : str
        Name of the package being tested.
    pf_printed : bool
        True if pass/fail has already been printed. This is intended to help
         avoid race conditions with SkipSequence() after a KeyboardInterrupt,
         resulting in the skip being printed after a test has passed.

    Returns
    -------
    pf_printed : bool
    """
    global testCt
    global debugMode
    # Skip irrelevant tests
    if ((tinker_e == None) and (QMMMPack in ("TINKER", "TINKER9"))):
        return
    if ((lammps_e == None) and (QMMMPack == "LAMMPS")):
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
            # Check g09 not given twice
            if 'g09' in args.qm:
                qm_test_packs.add("g09")
            else:
                qm_test_packs.add("g09")
        elif 'g09' in args.qm:
            qm_test_packs.add("g09")
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

# Look for LICHEM
LICHEMbin = LocateLICHEM()

# Identify all available wrappers
if allTests is True:
    # Identify QM wrappers
    print("Available QM wrappers:")
    # Search for PSI4
    QMbin = LocateProgram(executable='psi4')
    print(f" PSI4: {QMbin}")
    # Search for Gaussian
    QMbin = LocateProgram(executable='g09')
    print(f" Gaussian09: {QMbin}")
    # Check if both g09 and g16 co-exist
    if useg16 is False:
        QMbin = LocateProgram(executable='g16')
        if QMbin != "N/A":
            useg16 = True
        print(f" Gaussian16: {QMbin}")
    # Search for NWChem
    QMbin = LocateProgram(executable='nwchem')
    print(f" NWChem: {QMbin}")
    # Identify MM wrappers
    print("\nAvailable MM wrappers:")
    # Search for TINKER
    MMbin = LocateProgram(executable='analyze')
    print(f" TINKER: {MMbin}")
    # Search for TINKER9
    MMbin = LocateProgram(executable='tinker9')
    print(f" TINKER9: {MMbin}")
    # Search for LAMMPS
    MMbin = LocateProgram(executable='lammps')
    print(f" LAMMPS: {MMbin}")
# Search for requested wrappers
else:
    badQM = True
    # Ask if wrapper is in since it's part of a set!
    if "psi4" in qm_test_packs:
        QMPack = "PSI4"
        QMbin = LocateProgram(executable='psi4')
        if QMbin != "N/A":
            badQM = False
    # Search for Gaussian
    if "g09" in qm_test_packs:
        QMPack = "Gaussian09"
        QMbin = LocateProgram(executable='g09')
        if QMbin != "N/A":
            badQM = False
    if "g16" in qm_test_packs:
        QMPack = "Gaussian16"
        QMbin = LocateProgram(executable='g16')
        if QMbin != "N/A":
            useg16 = True
            badQM = False
    if useg16 is True:
        QMPack = "Gaussian16"
    if "nwchem" in qm_test_packs:
        QMPack = "NWChem"
        QMbin = LocateProgram(executable='nwchem')
        if QMbin != "N/A":
            badQM = False
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
    # Search for TINKER9
    if "tinker9" in mm_test_packs:
        MMPack = "TINKER9"
        MMbin = LocateProgram(executable='tinker9')
        if MMbin != "N/A":
            badMM = False
    # Search for LAMMPS
    if "lammps" in mm_test_packs:
        MMPack = "LAMMPS"
        MMbin = LocateProgram(executable='lammps')
        if MMbin != "N/A":
            badMM = False
    if badMM is True:
        # Quit with error
        print("\nError: No MM executables located for any of the"
              f" following: {mm_test_packs}.\n")
        exit(0)

# Print test settings
print(
    "Settings:\n"
    f" Threads: {Ncpus}")
if allTests is False:
    print(
        f" LICHEM binary: {LICHEMbin}\n"
        f" QM package: {QMPack}\n"
        f" Binary: {QMbin}\n"
        f" MM package: {MMPack}\n"
        f" Binary: {MMbin}\n")
else:
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

# TODO: FIXME
# Make a list of tests
QMTests = []
MMTests = []
if allTests is True:
    # Safely add PSI4
    cmd = "which psi4"
    try:
        # Run PSI4 tests
        packBin = subprocess.check_output(cmd, shell=True)
        QMTests.append("PSI4")
    except:
        # Skip tests that will fail
        if forceAll is True:
            QMTests.append("PSI4")
    # Safely add Gaussian
    cmd = "which g09"
    try:
        # Run Gaussian tests
        packBin = subprocess.check_output(cmd, shell=True)
        QMTests.append("Gaussian")
    except subprocess.CalledProcessError:
        cmd = "which g16"
        try:
            packBin = subprocess.check_output(cmd, shell=True)
            QMTests.append("g16")
        except:
            # Skip tests that will fail
            if forceAll is True:
                # QMTests.append("Gaussian")
                # Instead try g16 if no other Gaussian found
                QMTests.append("g16")
    # Safely add NWChem
    cmd = "which nwchem"
    try:
        # Run NWChem tests
        packBin = subprocess.check_output(cmd, shell=True)
        QMTests.append("NWChem")
    except:
        # Skip tests that will fail
        if forceAll is True:
            QMTests.append("NWChem")
    # Safely add TINKER
    cmd = "which analyze"
    try:
        # Run TINKER tests
        packBin = subprocess.check_output(cmd, shell=True)
        MMTests.append("TINKER")
    except:
        # Skip tests that will fail
        if forceAll is True:
            MMTests.append("TINKER")
    # Safely add lammps
    cmd = "which lammps"
    try:
        # Run LAMMPS tests
        packBin = subprocess.check_output(cmd, shell=True)
        MMTests.append("LAMMPS")
    except:
        # Skip tests that will fail
        if forceAll is True:
            MMTests.append("LAMMPS")
else:
    # Add only the specified packages
    QMTests.append(QMPack)
    MMTests.append(MMPack)

# NB: Tests are in the following order:
#     1) HF energy
#     2) PBE0 energy
#     3) CCSD energy
#     4) PM6 energy
#     5) Frequencies
#     6) NEB TS energy
#     7) TIP3P energy
#     8) AMOEBA/GK energy
#     9) PBE0/TIP3P energy
#    10) PBE0/AMOEBA energy
#    11) DFP/Pseudobonds

# --------------------- #
# --- Run the Tests --- #
# --------------------- #

# Loop over tests
# for qmTest in QMTests:
#     for mmTest in MMTests:
#         # Set packages
#         QMPack = qmTest
#         MMPack = mmTest

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
        #
        # # Start each try block by declaring that the Pass/Fail report for that
        # #  test hasn't yet printed.
        # pf_printed = False
        # # Use try/except to control Cntrl+C behavior with SkipSequence().
        # try:
        #     # Gaussian and Psi4 Only, so update counter/define test in the loop
        #     # Check HF energy
        #     if ((QMPack == "PSI4") or (QMPack == "Gaussian") or
        #        (QMPack == "Gaussian16")):
        #         pf_printed = Test1(name="HF Energy",
        #                            QMMMPack=QMPack, pf_printed=pf_printed)
        # except KeyboardInterrupt:
        #     SkipSequence("HF energy:", pf_printed)
        #
        # pf_printed = False
        # try:
        #     pf_printed = Test2(name="PBE0 Energy",
        #                        QMMMPack=QMPack, pf_printed=pf_printed)
        # except KeyboardInterrupt:
        #     SkipSequence("HF energy:", pf_printed)

        pf_printed = False
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

        pf_printed = False
        try:
            # Check DFT energy
            pf_printed = QMMMWrapperTest(
                name="PBE0 Energy",
                psi4_e=round(-4154.16836599/har2eV, 3), psi4_u="a.u.",
                g09_e=round(-4154.1676114324/har2eV, 3), g09_u="a.u.",
                g16_e=round(-4154.1676114324/har2eV, 3), g16_u="a.u.",
                nwchem_e=round(-4154.1683939169/har2eV, 3), nwchem_u="a.u.",
                energy_str="QM energy:", energy_loc=2,
                copy_command=None,
                xName="waterdimer.xyz", rName="pbereg.inp",
                cName="watercon.inp",
                QMMMPack=QMPack, pf_printed=pf_printed)
        except KeyboardInterrupt:
            SkipSequence("PBE0 energy:", pf_printed)

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

        pf_printed = False
        try:
            # Gaussian-only test
            pf_printed = QMMMFreqTest(
                name="Frequencies",
                psi4_e=round(-40.507339, 0), psi4_u="wavenumbers",
                g09_e=round(-31.769945, 0), g09_u="wavenumbers",
                g16_e=round(5.32988942, 0), g16_u="wavenumbers",
                nwchem_e=round(-31.769945, 0), nwchem_u="wavenumbers",
                xName="methfluor.xyz", rName="freqreg.inp",
                cName="methflcon.inp",
                QMMMPack=QMPack, pf_printed=pf_printed)
        except KeyboardInterrupt:
            SkipSequence("PM6 energy:", pf_printed)

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
        pf_printed = False
        try:
            pf_printed = MMWrapperTest(
                name="TIP3P Energy",
                tinker_e=round(-0.2596903536223/har2eV, 3), tinker_u="a.u.",
                lammps_e=None, lammps_u=None,
                energy_str="MM energy:", energy_loc=2,
                tinkerKey="pchrg.key",
                xName="waterdimer.xyz", rName="mmreg.inp",
                cName="watercon.inp",
                QMMMPack=MMPack, pf_printed=pf_printed)
        except KeyboardInterrupt:
            SkipSequence("TIP3P energy:", pf_printed)

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

        pf_printed = False
        try:
            pf_printed = QMMMWrapperTest(
                name="PBE0/TIP3P Energy",
                psi4_e=round(-2077.2021947277/har2eV, 3), psi4_u="a.u.",
                g09_e=round(-76.34642184669, 3), g09_u="a.u.",
                g16_e=round(-76.34642184669, 3), g16_u="a.u.",
                nwchem_e=round(-2077.2022117306/har2eV, 3), nwchem_u="a.u.",
                energy_str="QMMM energy:", energy_loc=2,
                copy_command=[CopyRequired(ifile="pchrg.key",
                                       ofile="tinker.key")],
                xName="waterdimer.xyz", rName="pchrgreg.inp",
                cName="watercon.inp",
                QMMMPack=QMPack, pf_printed=pf_printed)
        except KeyboardInterrupt:
            SkipSequence("PBE0/TIP3P energy:", pf_printed)

        pf_printed = False
        try:
            pf_printed = QMMMWrapperTest(
                name="PBE0/AMOEBA Energy",
                psi4_e=round(-2077.1114201829/har2eV, 3), psi4_u="a.u.",
                g09_e=round(-2077.1090319595/har2eV, 3), g09_u="a.u.",
                g16_e=round(-2077.1090319595/har2eV, 3), g16_u="a.u.",
                nwchem_e=round(-2077.1094168459/har2eV, 3), nwchem_u="a.u.",
                energy_str="QMMM energy:", energy_loc=2,
                copy_command=[CopyRequired(ifile="pol.key",
                                       ofile="tinker.key")],
                xName="waterdimer.xyz", rName="polreg.inp",
                cName="watercon.inp",
                QMMMPack=QMPack, pf_printed=pf_printed)
        except KeyboardInterrupt:
            SkipSequence("PBE0/AMOEBA energy:", pf_printed)

        pf_printed = False
        try:
            # Gaussian and NWChem only
            pf_printed = QMMMWrapperTest(
                name="DFP/Pseudobonds", QMMMPack=QMPack,
                # Copy two files!
                copy_command=[CopyRequired(ifile="pol.key",
                                       ofile="tinker.key"),
                              CopyRequired(ifile="pbbasis.txt",
                                       ofile="BASIS")],
                xName="alkyl.xyz", rName="pboptreg.inp",
                cName="alkcon.inp",
                energy_str="Opt. step: 2", energy_loc=6,
                psi4_e=None, psi4_u=None,
                g09_e=round(-3015.0548490566/har2eV, 3), g09_u="a.u.",
                g16_e=round(-3015.0548490566/har2eV, 3), g16_u="a.u.",
                nwchem_e=round(-3015.2278310975/har2eV, 3), nwchem_u="a.u.",
                pf_printed=pf_printed)
        except KeyboardInterrupt:
            SkipSequence("DFP/Pseudobonds:", pf_printed)

        # # Use try/except to control Cntrl+C behavior with SkipSequence().
        # try:
        #     # Gaussian and Psi4 Only, so update counter/define test in the loop
        #     # Check HF energy
        #     if ((QMPack == "PSI4") or (QMPack == "Gaussian") or
        #        (QMPack == "g16")):
        #         testCt += 1
        #         print(f"{' ':2}Test {testCt}: HF Energy")
        #         # Initialize energy as failing
        #         passEnergy = False
        #         runC = RunLICHEM("waterdimer.xyz", "hfreg.inp", "watercon.inp")
        #         QMMMEnergy, savedEnergy, units = RecoverEnergy("QM energy:", 2)
        #         # Check result
        #         if (QMPack == "PSI4"):
        #             expected_energy = round(-4136.93039814/har2eV, 3)
        #             # Check against saved energy
        #             if (QMMMEnergy == expected_energy):
        #                 passEnergy = True
        #             else:
        #                 saved_fail = SaveFailure(QMMMEnergy, expected_energy)
        #         if ((QMPack == "Gaussian") or (QMPack == "g16")):
        #             expected_energy = round(-4136.9317704519/har2eV, 3)
        #             # Check against saved energy
        #             if (QMMMEnergy == expected_energy):
        #                 passEnergy = True
        #             else:
        #                 saved_fail = SaveFailure(QMMMEnergy, expected_energy)
        #         pf_printed = PrintPassFail("HF energy:", passEnergy,
        #                                    updateResults, savedEnergy, units)
        #         if passEnergy is False:
        #             print(saved_fail)
        #             PrintLICHEMDebug(runC)
        #         CleanFiles()  # Clean up files
        # except KeyboardInterrupt:
        #     SkipSequence("HF energy:", pf_printed)
        #
        # pf_printed = False
        # try:
        #     testCt += 1
        #     print(f"\n{' ':2}Test {testCt}: PBE0 Energy")
        #     # Check DFT energy
        #     # line = ""
        #     passEnergy = False
        #     runC = RunLICHEM("waterdimer.xyz", "pbereg.inp", "watercon.inp")
        #     QMMMEnergy, savedEnergy, units = RecoverEnergy("QM energy:", 2)
        #     # Check result
        #     if (QMPack == "PSI4"):
        #         expected_energy = round(-4154.16836599/har2eV, 3)
        #         # Check against saved energy
        #         if (QMMMEnergy == expected_energy):
        #             passEnergy = True
        #         else:
        #             saved_fail = SaveFailure(QMMMEnergy, expected_energy)
        #     if ((QMPack == "Gaussian") or (QMPack == "g16")):
        #         expected_energy = round(-4154.1676114324/har2eV, 3)
        #         # Check against saved energy
        #         if (QMMMEnergy == expected_energy):
        #             passEnergy = True
        #         else:
        #             saved_fail = SaveFailure(QMMMEnergy, expected_energy)
        #     if (QMPack == "NWChem"):
        #         expected_energy = round(-4154.1683939169/har2eV, 3)
        #         # Check against saved energy
        #         if (QMMMEnergy == expected_energy):
        #             passEnergy = True
        #         else:
        #             saved_fail = SaveFailure(QMMMEnergy, expected_energy)
        #     pf_printed = PrintPassFail("PBE0 energy:", passEnergy,
        #                                updateResults, savedEnergy, units)
        #     if passEnergy is False:
        #         print(saved_fail)
        #         PrintLICHEMDebug(runC)
        #     CleanFiles()  # Clean up files
        # except KeyboardInterrupt:
        #     SkipSequence("PBE0 energy:", pf_printed)
        #
        # pf_printed = False
        # try:
        #     # PSI4 Only, so update counter/define test in the loop
        #     # Check CCSD energy
        #     if (QMPack == "PSI4"):
        #         testCt += 1
        #         print(f"\n{' ':2}Test {testCt}: CCSD Energy")
        #         # line = ""
        #         passEnergy = False
        #         runC = RunLICHEM("waterdimer.xyz", "ccsdreg.inp",
        #                          "watercon.inp")
        #         QMMMEnergy, savedEnergy, units = RecoverEnergy("QM energy:", 2)
        #         # Check result
        #         if (QMPack == "PSI4"):
        #             expected_energy = round(-4147.7304837/har2eV, 3)
        #             # Check against saved energy
        #             if (QMMMEnergy == expected_energy):
        #                 passEnergy = True
        #             else:
        #                 saved_fail = SaveFailure(QMMMEnergy, expected_energy)
        #         pf_printed = PrintPassFail("CCSD energy:", passEnergy,
        #                                    updateResults, savedEnergy, units)
        #         if passEnergy is False:
        #             print(saved_fail)
        #             PrintLICHEMDebug(runC)
        #         CleanFiles()
        # except KeyboardInterrupt:
        #     SkipSequence("CCSD energy:", pf_printed)
        #
        # pf_printed = False
        # try:
        #     # Gaussian Only, so update counter/define test in the loop
        #     # Check PM6 energy
        #     if ((QMPack == "Gaussian") or (QMPack == "g16")):
        #         testCt += 1
        #         print(f"\n{' ':2}Test {testCt}: PM6 Energy")
        #         # line = ""
        #         passEnergy = False
        #         runC = RunLICHEM("waterdimer.xyz", "pm6reg.inp",
        #                          "watercon.inp")
        #         QMMMEnergy, savedEnergy, units = RecoverEnergy("QM energy:", 2)
        #         # Check result
        #         if ((QMPack == "Gaussian") or (QMPack == "g16")):
        #             expected_energy = round(-4.8623027634995/har2eV, 3)
        #             # Check against saved energy
        #             if (QMMMEnergy == expected_energy):
        #                 passEnergy = True
        #             else:
        #                 saved_fail = SaveFailure(QMMMEnergy, expected_energy)
        #         pf_printed = PrintPassFail("PM6 energy:", passEnergy,
        #                                    updateResults, savedEnergy, units)
        #         if passEnergy is False:
        #             print(saved_fail)
        #             PrintLICHEMDebug(runC)
        #         CleanFiles()
        # except KeyboardInterrupt:
        #     SkipSequence("PM6 energy:", pf_printed)
        #
        # pf_printed = False
        # try:
        #     testCt += 1
        #     print(f"\n{' ':2}Test {testCt}: Frequencies")
        #     # Check imaginary frequencies
        #     # line = ""
        #     passEnergy = False
        #     runC = RunLICHEM("methfluor.xyz", "freqreg.inp", "methflcon.inp")
        #     QMMMFreqs, units = RecoverFreqs()
        #     # Find lowest frequency
        #     try:
        #         QMMMEnergy = min(QMMMFreqs)
        #         if QMMMEnergy > 1e100:
        #           savedEnergy = "Crashed..."
        #         else:
        #             savedEnergy = "Freq:{:3}{}".format(" ", QMMMEnergy)
        #     # If list is empty
        #     except ValueError:
        #         savedEnergy = "Crashed..."
        #     # Check results
        #     if (QMPack == "PSI4"):
        #         expected_energy = round(-40.507339, 0)
        #         # Check against saved frequency
        #         if (round(QMMMEnergy, 0) == expected_energy):
        #             passEnergy = True
        #         else:
        #             saved_fail = SaveFailure(round(QMMMEnergy, 0),
        #                                      expected_energy)
        #     if (QMPack == "Gaussian"):
        #         # Gaussian 09 only!
        #         expected_energy = round(-31.769945, 0)
        #         # Check against saved frequency
        #         if (round(QMMMEnergy, 0) == expected_energy):
        #             passEnergy = True
        #         else:
        #             saved_fail = SaveFailure(round(QMMMEnergy, 0),
        #                                      expected_energy)
        #     elif (QMPack == "g16"):
        #         # Gaussian 16 will produce a different result!
        #         expected_energy = round(5.32988942, 0)
        #         if (round(QMMMEnergy, 0) == expected_energy):
        #             passEnergy = True
        #         else:
        #             saved_fail = SaveFailure(round(QMMMEnergy, 0),
        #                                      expected_energy)
        #     if (QMPack == "NWChem"):
        #         expected_energy = round(-31.769945, 0)
        #         # Check against saved frequency
        #         if (round(QMMMEnergy, 0) == expected_energy):
        #             passEnergy = True
        #         else:
        #             saved_fail = SaveFailure(round(QMMMEnergy, 0),
        #                                      expected_energy)
        #     pf_printed = PrintPassFail("Frequencies:", passEnergy,
        #                                updateResults, savedEnergy, units)
        #     if passEnergy is False:
        #         print(saved_fail)
        #         PrintLICHEMDebug(runC)
        #     CleanFiles()
        # except KeyboardInterrupt:
        #     SkipSequence("Frequencies:", pf_printed)
        #
        # pf_printed = False
        # try:
        #     testCt += 1
        #     print(f"\n{' ':2}Test {testCt}: NEB TS Energy")
        #     # Check NEB optimization
        #     passEnergy = False
        #     cmd = "cp methflbeads.xyz BeadStartStruct.xyz"
        #     PrintCopyDebug(cmd)
        #     subprocess.call(cmd, shell=True)  # Copy restart file
        #     runC = RunLICHEM("methfluor.xyz", "nebreg.inp", "methflcon.inp")
        #     # This is searching for matching line, which is in eV units.
        #     #    au and kcal/mol are reported on next line.
        #     QMMMEnergy, savedEnergy, units = RecoverEnergy("Forward barrier", 3)
        #     # Check result
        #     if (QMPack == "PSI4"):
        #         expected_energy = round(0.39581219003957, 3)
        #         # Check against saved energy
        #         if (QMMMEnergy == expected_energy):
        #             passEnergy = True
        #         else:
        #             saved_fail = SaveFailure(QMMMEnergy, expected_energy)
        #     if ((QMPack == "Gaussian") or (QMPack == "g16")):
        #         expected_energy = round(0.39582668467028, 3)
        #         # Check against saved energy
        #         if (QMMMEnergy == expected_energy):
        #             passEnergy = True
        #         else:
        #             saved_fail = SaveFailure(QMMMEnergy, expected_energy)
        #     if (QMPack == "NWChem"):
        #         expected_energy = round(-0.39582668467028, 3)
        #         # Check against saved energy
        #         if (QMMMEnergy == expected_energy):
        #             passEnergy = True
        #         else:
        #             saved_fail = SaveFailure(QMMMEnergy, expected_energy)
        #     pf_printed = PrintPassFail("NEB Forward barrier:", passEnergy,
        #                                updateResults, savedEnergy, units)
        #     if passEnergy is False:
        #         print(saved_fail)
        #         PrintLICHEMDebug(runC)
        #     CleanFiles()
        # except KeyboardInterrupt:
        #     SkipSequence("NEB Forward barrier:", pf_printed)
        #
        # pf_printed = False
        # try:
        #     testCt += 1
        #     print(f"\n{' ':2}Test {testCt}: QSM TS energy")
        #     # Check QSM optimization
        #     passEnergy = False
        #     cmd = "cp methflbeads.xyz BeadStartStruct.xyz"
        #     PrintCopyDebug(cmd)
        #     subprocess.call(cmd, shell=True)  # Copy restart file
        #     runC = RunLICHEM("methfluor.xyz", "qsmreg.inp", "methflcon.inp")
        #     QMMMEnergy, savedEnergy, units = RecoverEnergy("TS Energy", 4)
        #     # Check result
        #     if (QMPack == "PSI4"):
        #         expected_energy = round(-239.27686429, 3)
        #         # Check against saved energy
        #         if (QMMMEnergy == expected_energy):
        #             passEnergy = True
        #         else:
        #             saved_fail = SaveFailure(QMMMEnergy, expected_energy)
        #     if ((QMPack == "Gaussian") or (QMPack == "g16")):
        #         # Given in a.u.
        #         expected_energy = round(-239.27681732, 3)
        #         # Check against saved energy
        #         if (QMMMEnergy == expected_energy):
        #             passEnergy = True
        #         else:
        #             saved_fail = SaveFailure(QMMMEnergy, expected_energy)
        #     if (QMPack == "NWChem"):
        #         expected_energy = round(-239.27681732, 3)
        #         # Check against saved energy
        #         if (QMMMEnergy == expected_energy):
        #             passEnergy = True
        #         else:
        #             saved_fail = SaveFailure(QMMMEnergy, expected_energy)
        #     pf_printed = PrintPassFail("QSM TS energy:", passEnergy,
        #                                updateResults, savedEnergy, units)
        #     if passEnergy is False:
        #         print(saved_fail)
        #         PrintLICHEMDebug(runC)
        #     CleanFiles()
        # except KeyboardInterrupt:
        #     SkipSequence("QSM TS energy:", pf_printed)
        #
        # # TINKER-only, so update counter/define test in the loop
        # if (MMPack == "TINKER"):
        #     pf_printed = False
        #     try:
        #         testCt += 1
        #         print(f"\n{' ':2}Test {testCt}: TIP3P energy")
        #         # Check MM energy
        #         passEnergy = False
        #         cmd = "cp pchrg.key tinker.key"
        #         PrintCopyDebug(cmd)
        #         subprocess.call(cmd, shell=True)  # Copy key file
        #         runC = RunLICHEM("waterdimer.xyz", "mmreg.inp",
        #                          "watercon.inp")
        #         QMMMEnergy, savedEnergy, units = RecoverEnergy("MM energy:", 2)
        #         # Expected TINKER energy
        #         expected_energy = round(-0.2596903536223/har2eV, 3)
        #         # Check result
        #         if (QMMMEnergy == expected_energy):
        #             # Check against saved energy
        #             passEnergy = True
        #         else:
        #             saved_fail = SaveFailure(QMMMEnergy, expected_energy)
        #         pf_printed = PrintPassFail("TIP3P energy:", passEnergy,
        #                                    updateResults, savedEnergy, units)
        #         if passEnergy is False:
        #             print(saved_fail)
        #             PrintLICHEMDebug(runC)
        #         CleanFiles()
        #     except KeyboardInterrupt:
        #         SkipSequence("TIP3P energy:", pf_printed)
        #
        #     pf_printed = False
        #     try:
        #         # New test instance
        #         testCt += 1
        #         print(f"\n{' ':2}Test {testCt}: AMOEBA/GK Energy")
        #         # Check MM energy
        #         line = ""
        #         passEnergy = False
        #         cmd = "cp pol.key tinker.key"
        #         PrintCopyDebug(cmd)
        #         subprocess.call(cmd, shell=True)  # Copy key file
        #         runC = RunLICHEM("waterdimer.xyz", "solvreg.inp",
        #                          "watercon.inp")
        #         QMMMEnergy, savedEnergy, units = RecoverEnergy("MM energy:", 2)
        #         # Expected TINKER energy (a.u.)
        #         expected_energy = round(-0.0095434448906, 3)
        #         # Check result
        #         if (QMMMEnergy == expected_energy):
        #             # Check against saved energy
        #             passEnergy = True
        #         else:
        #             saved_fail = SaveFailure(QMMMEnergy, expected_energy)
        #         pf_printed = PrintPassFail("AMOEBA/GK energy:", passEnergy,
        #                                    updateResults, savedEnergy, units)
        #         if passEnergy is False:
        #             print(saved_fail)
        #             PrintLICHEMDebug(runC)
        #         CleanFiles()
        #     except KeyboardInterrupt:
        #         SkipSequence("AMOEBA/GK energy:", pf_printed)
        #
        # pf_printed = False
        # # Back to tests non-specific to TINKER
        # try:
        #     testCt += 1
        #     print(f"\n{' ':2}Test {testCt}: PBE0/TIP3P Energy")
        #     # Check QMMM point-charge energy results
        #     passEnergy = False
        #     cmd = "cp pchrg.key tinker.key"
        #     PrintCopyDebug(cmd)
        #     subprocess.call(cmd, shell=True)  # Copy key file
        #     runC = RunLICHEM("waterdimer.xyz", "pchrgreg.inp",
        #                      "watercon.inp")
        #     QMMMEnergy, savedEnergy, units = RecoverEnergy("QMMM energy:", 2)
        #     # Expected TINKER energy
        #     expected_energy = round(-2077.2021947277/har2eV, 3)
        #     # Check result
        #     if (QMPack == "PSI4"):
        #         # Check against saved energy
        #         if (QMMMEnergy == expected_energy):
        #             passEnergy = True
        #         else:
        #             saved_fail = SaveFailure(QMMMEnergy, expected_energy)
        #     if ((QMPack == "Gaussian") or (QMPack == "g16")):
        #         # Expected in a.u.
        #         expected_energy = round(-76.34642184669, 3)
        #         # Check against saved energy
        #         if (QMMMEnergy == expected_energy):
        #             passEnergy = True
        #         else:
        #             saved_fail = SaveFailure(QMMMEnergy, expected_energy)
        #     if (QMPack == "NWChem"):
        #         expected_energy = round(-2077.2022117306/har2eV, 3)
        #         # Check against saved energy
        #         if (QMMMEnergy == expected_energy):
        #             passEnergy = True
        #         else:
        #             saved_fail = SaveFailure(QMMMEnergy, expected_energy)
        #     pf_printed = PrintPassFail("PBE0/TIP3P energy:", passEnergy,
        #                                updateResults, savedEnergy, units)
        #     if passEnergy is False:
        #         print(saved_fail)
        #         PrintLICHEMDebug(runC)
        #     CleanFiles()
        # except KeyboardInterrupt:
        #     SkipSequence("PBE0/TIP3P energy:", pf_printed)
        #
        # pf_printed = False
        # try:
        #     testCt += 1
        #     print(f"\n{' ':2}Test {testCt}: PBE0/AMOEBA Energy")
        #     # Check QMMM polarizable energy results
        #     passEnergy = False
        #     cmd = "cp pol.key tinker.key"
        #     PrintCopyDebug(cmd)
        #     subprocess.call(cmd, shell=True)  # Copy key file
        #     runC = RunLICHEM("waterdimer.xyz", "polreg.inp", "watercon.inp")
        #     QMMMEnergy, savedEnergy, units = RecoverEnergy("QMMM energy:", 2)
        #     # Check result
        #     if (QMPack == "PSI4"):
        #         expected_energy = round(-2077.1114201829/har2eV, 3)
        #         # Check against saved energy
        #         if (QMMMEnergy == expected_energy):
        #             passEnergy = True
        #         else:
        #             saved_fail = SaveFailure(QMMMEnergy, expected_energy)
        #     if ((QMPack == "Gaussian") or (QMPack == "g16")):
        #         # In a.u.
        #         expected_energy = round(-0.0095434448906, 3)
        #         # Check against saved energy
        #         if (QMMMEnergy == expected_energy):
        #             passEnergy = True
        #         else:
        #             saved_fail = SaveFailure(QMMMEnergy, expected_energy)
        #     if (QMPack == "NWChem"):
        #         expected_energy = round(-2077.1094168459/har2eV, 3)
        #         # Check against saved energy
        #         if (QMMMEnergy == expected_energy):
        #             passEnergy = True
        #         else:
        #             saved_fail = SaveFailure(QMMMEnergy, expected_energy)
        #     pf_printed = PrintPassFail("PBE0/AMOEBA energy:", passEnergy,
        #                                updateResults, savedEnergy, units)
        #     if passEnergy is False:
        #         print(saved_fail)
        #         PrintLICHEMDebug(runC)
        #     CleanFiles()
        # except KeyboardInterrupt:
        #     SkipSequence("PBE0/AMOEBA energy:", pf_printed)
        #
        # pf_printed = False
        # try:
        #     # Gaussian and NWChem only, so update counter/def test in the loop
        #     # Check pseudobond optimizations
        #     if ((QMPack == "Gaussian") or (QMPack == "g16") or
        #        (QMPack == "NWChem")):
        #         testCt += 1
        #         print(f"\n{' ':2}Test {testCt}: DFP/Pseudobonds")
        #         # Carry out test
        #         passEnergy = False
        #         cmd = "cp pbopt.key tinker.key"
        #         PrintCopyDebug(cmd)
        #         subprocess.call(cmd, shell=True)  # Copy key file
        #         cmd = "cp pbbasis.txt BASIS"
        #         PrintCopyDebug(cmd)
        #         subprocess.call(cmd, shell=True)  # Copy BASIS set file
        #         runC = RunLICHEM("alkyl.xyz", "pboptreg.inp", "alkcon.inp")
        #         QMMMEnergy, savedEnergy, units = RecoverEnergy("Opt. step: 2",
        #                                                                     6)
        #         # Check result
        #         if ((QMPack == "Gaussian") or (QMPack == "g16")):
        #             expected_energy = round(-3015.0548490566/har2eV, 3)
        #             # Check against saved energy
        #             if (QMMMEnergy == expected_energy):
        #                 passEnergy = True
        #             else:
        #                 saved_fail = SaveFailure(QMMMEnergy, expected_energy)
        #         if (QMPack == "NWChem"):
        #             expected_energy = round(-3015.2278310975/har2eV, 3)
        #             # Check against saved energy
        #             if (QMMMEnergy == expected_energy):
        #                 passEnergy = True
        #             else:
        #                 saved_fail = SaveFailure(QMMMEnergy, expected_energy)
        #         pf_printed = PrintPassFail("DFP/Pseudobonds:", passEnergy,
        #                                    updateResults, savedEnergy, units)
        #         if passEnergy is False:
        #             print(saved_fail)
        #             PrintLICHEMDebug(runC)
        #         CleanFiles()
        # except KeyboardInterrupt:
        #     SkipSequence("DFP/Pseudobonds:", pf_printed)

        # Revert the regions files, but only if you're not debugging!
        if debugMode is False:
            for filename in regions_files:
                # Check that the file exists for sed
                if os.path.isfile("./"+filename):
                    # Gaussian
                    if QMPack == "g16":
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
