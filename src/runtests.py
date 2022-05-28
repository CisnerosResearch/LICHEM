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
import subprocess
import time
import sys
import os
import platform

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
        units = str(finalEnergy[itemNum+1]) # For developer mode
        finalEnergy = float(finalEnergy[itemNum])
        savedResult = "Energy: "+str(finalEnergy)  # Save it for later
        finalEnergy = round(finalEnergy, 3)
    except:
        # Calculation failed
        finalEnergy = 0.0
        units = ":("
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
    # Likely need this to be subprocess.CalledProcessError and exception
    #  raised by not being able to split output...
    except:
        # Calculation failed
        if debugMode is True:
            print("\nFrequencies not recovered.")
        freqList = []
    return freqList, units


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

def PrintCopyDebug(cmd_printed):
    """
    Prints the copy command that was run to the console.

    Parameters
    -------
    cmd_printed : str
        The attempted copy command.
    """
    global debugMode
    if debugMode is True:
        print(f"\n{' ':6}- {cmd_printed}\n")
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

# Read arguments
dryRun = False    # Only check packages
allTests = False  # Run all tests at once
if (len(sys.argv) == 3):
    if ((sys.argv[2]).lower() == "all"):
        # Automatically run all tests
        allTests = True
if (len(sys.argv) < 4):
    line = ""
    if allTests is False:
        # Print help if arguments are missing
        print(
              "Usage:\n"
              " user:$ ./runtests Ncpus All\n"
              " or \n"
              " user:$ ./runtests Ncpus QMPackage MMPackage\n"
              " or \n"
              " user:$ ./runtests Ncpus QMPackage MMPackage dry\n"
              )
    # Find LICHEM
    cmd = "which lichem"
    try:
        # Find path
        LICHEMbin = subprocess.check_output(cmd, shell=True)
        LICHEMbin = ClrSet.TPass+LICHEMbin.decode('utf-8').strip()+ClrSet.Reset
    except:
        LICHEMbin = ClrSet.TFail+"N/A"+ClrSet.Reset
    print(f"LICHEM binary: {LICHEMbin}\n")
    # Identify QM wrappers
    print("Available QM wrappers:")
    # Search for PSI4
    cmd = "which psi4"
    try:
        # Find path
        QMbin = subprocess.check_output(cmd, shell=True)
        QMbin = ClrSet.TPass+QMbin.decode('utf-8').strip()+ClrSet.Reset
    except:
        QMbin = ClrSet.TFail+"N/A"+ClrSet.Reset
    print(f" PSI4: {QMbin}")
    # Search for Gaussian
    cmd = "which g09"
    try:
        # Find path
        QMbin = subprocess.check_output(cmd, shell=True)
        QMbin = ClrSet.TPass+QMbin.decode('utf-8').strip()+ClrSet.Reset
    except subprocess.CalledProcessError:
        cmd = "which g16"
        try:
            # Find path
            QMbin = subprocess.check_output(cmd, shell=True)
            QMbin = ClrSet.TPass+QMbin.decode('utf-8').strip()+ClrSet.Reset
            useg16 = True  # Use Gaussian 16 for tests
        except:
            QMbin = ClrSet.TFail+"N/A"+ClrSet.Reset
    print(f" Gaussian: {QMbin}")
    # Search for NWChem
    cmd = "which nwchem"
    try:
        # Find path
        QMbin = subprocess.check_output(cmd, shell=True)
        QMbin = ClrSet.TPass+QMbin.decode('utf-8').strip()+ClrSet.Reset
    except:
        QMbin = ClrSet.TFail+"N/A"+ClrSet.Reset
    print(f" NWChem: {QMbin}\n")
    # Identify MM wrappers
    print("Available MM wrappers:")
    # Search for TINKER
    cmd = "which analyze"
    try:
        # Find path
        MMbin = subprocess.check_output(cmd, shell=True)
        MMbin = ClrSet.TPass+MMbin.decode('utf-8').strip()+ClrSet.Reset
    except:
        MMbin = ClrSet.TFail+"N/A"+ClrSet.Reset
    print(f" TINKER: {MMbin}")
    # Search for LAMMPS
    cmd = "which lammps"
    try:
        # Find path
        MMbin = subprocess.check_output(cmd, shell=True)
        MMbin = ClrSet.TPass+MMbin.decode('utf-8').strip()+ClrSet.Reset
    except:
        MMbin = ClrSet.TFail+"N/A"+ClrSet.Reset
    print(f" LAMMPS: {MMbin}\n")
    if allTests is False:
        # Quit
        exit(0)

# Parse command line options
Ncpus = int(sys.argv[1])  # Threads
if allTests is False:
    QMPack = sys.argv[2]  # QM wrapper for calculations
    QMPack = QMPack.lower()
    MMPack = sys.argv[3]  # MM wrapper for calculations
    MMPack = MMPack.lower()
    if (len(sys.argv) > 4):
        if ((sys.argv[4]).lower() == "dry"):
            # Quit early
            dryRun = True

# Initialize variables
LICHEMbin = ""
QMbin = ""
MMbin = ""

# Check packages and identify missing binaries
badLICHEM = False  # Assume executable can be found.
cmd = "which lichem"
try:
    # Find path
    LICHEMbin = subprocess.check_output(cmd, shell=True)
    LICHEMbin = LICHEMbin.decode('utf-8').strip()
except:
    LICHEMbin = "N/A"
    badLICHEM = True
if badLICHEM is True:
    # Quit with an error
    print("Error: LICHEM binary not found!\n")
    exit(0)
if allTests is False:
    # Check for QM binaries
    badQM = True
    if ((QMPack == "psi4") or (QMPack == "psi")):
        QMPack = "PSI4"
        cmd = "which psi4"
        try:
            # Find path
            QMbin = subprocess.check_output(cmd, shell=True)
            QMbin = QMbin.decode('utf-8').strip()
            badQM = False
        except:
            QMbin = "N/A"
    if ((QMPack == "gaussian") or (QMPack == "g09")):
        QMPack = "Gaussian"
        cmd = "which g09"
        try:
            # Find path
            QMbin = subprocess.check_output(cmd, shell=True)
            QMbin = QMbin.decode('utf-8').strip()
            badQM = False
        except:
            QMbin = "N/A"
    if ((QMPack == "g16") or (useg16 is True)):
        QMPack = "g16"
        cmd = "which g16"
        try:
            # Find path
            QMbin = subprocess.check_output(cmd, shell=True)
            QMbin = QMbin.decode('utf-8').strip()
            badQM = False
        except:
            QMbin = "N/A"
    if (QMPack == "nwchem"):
        QMPack = "NWChem"
        cmd = "which nwchem"
        try:
            # Find path
            QMbin = subprocess.check_output(cmd, shell=True)
            QMbin = QMbin.decode('utf-8').strip()
            badQM = False
        except:
            QMbin = "N/A"
    if badQM is True:
        # Quit with an error
        print(f"\nError: QM package name '{QMPack}' not recognized.\n")
        exit(0)
    # Check for MM
    badMM = True
    if (MMPack == "tinker"):
        MMPack = "TINKER"
        cmd = "which analyze"
        try:
            # Find path
            MMbin = subprocess.check_output(cmd, shell=True)
            MMbin = MMbin.decode('utf-8').strip()
            badMM = False
        except:
            MMbin = "N/A"
    elif (MMPack == "tinker9"):
        MMPack = "TINKER"
        cmd = "which tinker9"
        try:
            # Find path
            MMbin = subprocess.check_output(cmd, shell=True)
            MMbin = MMbin.decode('utf-8').strip()
            badMM = False
        except:
            MMbin = "N/A"
    if (MMPack == "lammps"):
        MMPack = "LAMMPS"
        cmd = "which lammps"
        try:
            # Find path
            MMbin = subprocess.check_output(cmd, shell=True)
            MMbin = MMbin.decode('utf-8').strip()
            badMM = False
        except:
            MMbin = "N/A"
    if badMM is True:
        # Quit with error
        print(f"\nError: MM package name '{MMPack}' not recognized.\n")
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
for qmTest in QMTests:
    for mmTest in MMTests:
        # Set packages
        QMPack = qmTest
        MMPack = mmTest

        # Set path based on packages
        dirPath = ""
        if (QMPack == "PSI4"):
            dirPath += "PSI4_"
        if QMPack in ("Gaussian", "g09", "g16"):
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
                if QMPack == "g16":
                    PrepRegions("QM_type", "Gaussian", "g16", filename)
                # Tinker version
                if MMPack == "TINKER9":
                    PrepRegions("MM_type", "TINKER", "TINKER9", filename)

        # Start printing results
        print(f"{QMPack}/{MMPack} results:")

        # Start each try block by declaring that the Pass/Fail report for that
        #  test hasn't yet printed.
        pf_printed = False
        # Use try/except to control Cntrl+C behavior with SkipSequence().
        try:
            # Gaussian and Psi4 Only, so update counter/define test in the loop
            # Check HF energy
            if ((QMPack == "PSI4") or (QMPack == "Gaussian") or
               (QMPack == "g16")):
                testCt += 1
                print(f"{' ':2}Test {testCt}: HF Energy")
                # Initialize energy as failing
                passEnergy = False
                runC = RunLICHEM("waterdimer.xyz", "hfreg.inp", "watercon.inp")
                QMMMEnergy, savedEnergy, units = RecoverEnergy("QM energy:", 2)
                # Check result
                if (QMPack == "PSI4"):
                    expected_energy = round(-4136.93039814/har2eV, 3)
                    # Check against saved energy
                    if (QMMMEnergy == expected_energy):
                        passEnergy = True
                    else:
                        saved_fail = SaveFailure(QMMMEnergy, expected_energy)
                if ((QMPack == "Gaussian") or (QMPack == "g16")):
                    expected_energy = round(-4136.9317704519/har2eV, 3)
                    # Check against saved energy
                    if (QMMMEnergy == expected_energy):
                        passEnergy = True
                    else:
                        saved_fail = SaveFailure(QMMMEnergy, expected_energy)
                pf_printed = PrintPassFail("HF energy:", passEnergy,
                                           updateResults, savedEnergy, units)
                if passEnergy is False:
                    print(saved_fail)
                    PrintLICHEMDebug(runC)
                CleanFiles()  # Clean up files
        except KeyboardInterrupt:
            SkipSequence("HF energy:", pf_printed)

        pf_printed = False
        try:
            testCt += 1
            print(f"\n{' ':2}Test {testCt}: PBE0 Energy")
            # Check DFT energy
            # line = ""
            passEnergy = False
            runC = RunLICHEM("waterdimer.xyz", "pbereg.inp", "watercon.inp")
            QMMMEnergy, savedEnergy, units = RecoverEnergy("QM energy:", 2)
            # Check result
            if (QMPack == "PSI4"):
                expected_energy = round(-4154.16836599/har2eV, 3)
                # Check against saved energy
                if (QMMMEnergy == expected_energy):
                    passEnergy = True
                else:
                    saved_fail = SaveFailure(QMMMEnergy, expected_energy)
            if ((QMPack == "Gaussian") or (QMPack == "g16")):
                expected_energy = round(-4154.1676114324/har2eV, 3)
                # Check against saved energy
                if (QMMMEnergy == expected_energy):
                    passEnergy = True
                else:
                    saved_fail = SaveFailure(QMMMEnergy, expected_energy)
            if (QMPack == "NWChem"):
                expected_energy = round(-4154.1683939169/har2eV, 3)
                # Check against saved energy
                if (QMMMEnergy == expected_energy):
                    passEnergy = True
                else:
                    saved_fail = SaveFailure(QMMMEnergy, expected_energy)
            pf_printed = PrintPassFail("PBE0 energy:", passEnergy,
                                       updateResults, savedEnergy, units)
            if passEnergy is False:
                print(saved_fail)
                PrintLICHEMDebug(runC)
            CleanFiles()  # Clean up files
        except KeyboardInterrupt:
            SkipSequence("PBE0 energy:", pf_printed)

        pf_printed = False
        try:
            # PSI4 Only, so update counter/define test in the loop
            # Check CCSD energy
            if (QMPack == "PSI4"):
                testCt += 1
                print(f"\n{' ':2}Test {testCt}: CCSD Energy")
                # line = ""
                passEnergy = False
                runC = RunLICHEM("waterdimer.xyz", "ccsdreg.inp",
                                 "watercon.inp")
                QMMMEnergy, savedEnergy, units = RecoverEnergy("QM energy:", 2)
                # Check result
                if (QMPack == "PSI4"):
                    expected_energy = round(-4147.7304837/har2eV, 3)
                    # Check against saved energy
                    if (QMMMEnergy == expected_energy):
                        passEnergy = True
                    else:
                        saved_fail = SaveFailure(QMMMEnergy, expected_energy)
                pf_printed = PrintPassFail("CCSD energy:", passEnergy,
                                           updateResults, savedEnergy, units)
                if passEnergy is False:
                    print(saved_fail)
                    PrintLICHEMDebug(runC)
                CleanFiles()
        except KeyboardInterrupt:
            SkipSequence("CCSD energy:", pf_printed)

        pf_printed = False
        try:
            # Gaussian Only, so update counter/define test in the loop
            # Check PM6 energy
            if ((QMPack == "Gaussian") or (QMPack == "g16")):
                testCt += 1
                print(f"\n{' ':2}Test {testCt}: PM6 Energy")
                # line = ""
                passEnergy = False
                runC = RunLICHEM("waterdimer.xyz", "pm6reg.inp",
                                 "watercon.inp")
                QMMMEnergy, savedEnergy, units = RecoverEnergy("QM energy:", 2)
                # Check result
                if ((QMPack == "Gaussian") or (QMPack == "g16")):
                    expected_energy = round(-4.8623027634995/har2eV, 3)
                    # Check against saved energy
                    if (QMMMEnergy == expected_energy):
                        passEnergy = True
                    else:
                        saved_fail = SaveFailure(QMMMEnergy, expected_energy)
                pf_printed = PrintPassFail("PM6 energy:", passEnergy,
                                           updateResults, savedEnergy, units)
                if passEnergy is False:
                    print(saved_fail)
                    PrintLICHEMDebug(runC)
                CleanFiles()
        except KeyboardInterrupt:
            SkipSequence("PM6 energy:", pf_printed)

        pf_printed = False
        try:
            testCt += 1
            print(f"\n{' ':2}Test {testCt}: Frequencies")
            # Check imaginary frequencies
            # line = ""
            passEnergy = False
            runC = RunLICHEM("methfluor.xyz", "freqreg.inp", "methflcon.inp")
            QMMMFreqs, units = RecoverFreqs()
            # Find lowest frequency
            try:
                QMMMEnergy = min(QMMMFreqs)
                if QMMMEnergy > 1e100:
                  savedEnergy = "Crashed..."
                else:
                    savedEnergy = "Freq:{:3}{}".format(" ", QMMMEnergy)
            # If list is empty
            except ValueError:
                savedEnergy = "Crashed..."
            # Check results
            if (QMPack == "PSI4"):
                expected_energy = round(-40.507339, 0)
                # Check against saved frequency
                if (round(QMMMEnergy, 0) == expected_energy):
                    passEnergy = True
                else:
                    saved_fail = SaveFailure(round(QMMMEnergy, 0),
                                             expected_energy)
            if (QMPack == "Gaussian"):
                # Gaussian 09 only!
                expected_energy = round(-31.769945, 0)
                # Check against saved frequency
                if (round(QMMMEnergy, 0) == expected_energy):
                    passEnergy = True
                else:
                    saved_fail = SaveFailure(round(QMMMEnergy, 0),
                                             expected_energy)
            elif (QMPack == "g16"):
                # Gaussian 16 will produce a different result!
                expected_energy = round(5.32988942, 0)
                if (round(QMMMEnergy, 0) == expected_energy):
                    passEnergy = True
                else:
                    saved_fail = SaveFailure(round(QMMMEnergy, 0),
                                             expected_energy)
            if (QMPack == "NWChem"):
                expected_energy = round(-31.769945, 0)
                # Check against saved frequency
                if (round(QMMMEnergy, 0) == expected_energy):
                    passEnergy = True
                else:
                    saved_fail = SaveFailure(round(QMMMEnergy, 0),
                                             expected_energy)
            pf_printed = PrintPassFail("Frequencies:", passEnergy,
                                       updateResults, savedEnergy, units)
            if passEnergy is False:
                print(saved_fail)
                PrintLICHEMDebug(runC)
            CleanFiles()
        except KeyboardInterrupt:
            SkipSequence("Frequencies:", pf_printed)

        pf_printed = False
        try:
            testCt += 1
            print(f"\n{' ':2}Test {testCt}: NEB TS Energy")
            # Check NEB optimization
            passEnergy = False
            cmd = "cp methflbeads.xyz BeadStartStruct.xyz"
            PrintCopyDebug(cmd)
            subprocess.call(cmd, shell=True)  # Copy restart file
            runC = RunLICHEM("methfluor.xyz", "nebreg.inp", "methflcon.inp")
            # This is searching for matching line, which is in eV units.
            #    au and kcal/mol are reported on next line.
            QMMMEnergy, savedEnergy, units = RecoverEnergy("Forward barrier", 3)
            # Check result
            if (QMPack == "PSI4"):
                expected_energy = round(0.39581219003957, 3)
                # Check against saved energy
                if (QMMMEnergy == expected_energy):
                    passEnergy = True
                else:
                    saved_fail = SaveFailure(QMMMEnergy, expected_energy)
            if ((QMPack == "Gaussian") or (QMPack == "g16")):
                expected_energy = round(0.39582668467028, 3)
                # Check against saved energy
                if (QMMMEnergy == expected_energy):
                    passEnergy = True
                else:
                    saved_fail = SaveFailure(QMMMEnergy, expected_energy)
            if (QMPack == "NWChem"):
                expected_energy = round(-0.39582668467028, 3)
                # Check against saved energy
                if (QMMMEnergy == expected_energy):
                    passEnergy = True
                else:
                    saved_fail = SaveFailure(QMMMEnergy, expected_energy)
            pf_printed = PrintPassFail("NEB Forward barrier:", passEnergy,
                                       updateResults, savedEnergy, units)
            if passEnergy is False:
                print(saved_fail)
                PrintLICHEMDebug(runC)
            CleanFiles()
        except KeyboardInterrupt:
            SkipSequence("NEB Forward barrier:", pf_printed)

        pf_printed = False
        try:
            testCt += 1
            print(f"\n{' ':2}Test {testCt}: QSM TS energy")
            # Check QSM optimization
            passEnergy = False
            cmd = "cp methflbeads.xyz BeadStartStruct.xyz"
            PrintCopyDebug(cmd)
            subprocess.call(cmd, shell=True)  # Copy restart file
            runC = RunLICHEM("methfluor.xyz", "qsmreg.inp", "methflcon.inp")
            QMMMEnergy, savedEnergy, units = RecoverEnergy("TS Energy", 4)
            # Check result
            if (QMPack == "PSI4"):
                expected_energy = round(-239.27686429, 3)
                # Check against saved energy
                if (QMMMEnergy == expected_energy):
                    passEnergy = True
                else:
                    saved_fail = SaveFailure(QMMMEnergy, expected_energy)
            if ((QMPack == "Gaussian") or (QMPack == "g16")):
                expected_energy = round(-239.27681732, 3)
                # Check against saved energy
                if (QMMMEnergy == expected_energy):
                    passEnergy = True
                else:
                    saved_fail = SaveFailure(QMMMEnergy, expected_energy)
            if (QMPack == "NWChem"):
                expected_energy = round(-239.27681732, 3)
                # Check against saved energy
                if (QMMMEnergy == expected_energy):
                    passEnergy = True
                else:
                    saved_fail = SaveFailure(QMMMEnergy, expected_energy)
            pf_printed = PrintPassFail("QSM TS energy:", passEnergy,
                                       updateResults, savedEnergy, units)
            if passEnergy is False:
                print(saved_fail)
                PrintLICHEMDebug(runC)
            CleanFiles()
        except KeyboardInterrupt:
            SkipSequence("QSM TS energy:", pf_printed)

        # TINKER-only, so update counter/define test in the loop
        if (MMPack == "TINKER"):
            pf_printed = False
            try:
                testCt += 1
                print(f"\n{' ':2}Test {testCt}: TIP3P energy")
                # Check MM energy
                passEnergy = False
                cmd = "cp pchrg.key tinker.key"
                PrintCopyDebug(cmd)
                subprocess.call(cmd, shell=True)  # Copy key file
                runC = RunLICHEM("waterdimer.xyz", "mmreg.inp",
                                 "watercon.inp")
                QMMMEnergy, savedEnergy, units = RecoverEnergy("MM energy:", 2)
                # Expected TINKER energy
                expected_energy = round(-0.2596903536223/har2eV, 3)
                # Check result
                if (QMMMEnergy == expected_energy):
                    # Check against saved energy
                    passEnergy = True
                else:
                    saved_fail = SaveFailure(QMMMEnergy, expected_energy)
                pf_printed = PrintPassFail("TIP3P energy:", passEnergy,
                                           updateResults, savedEnergy, units)
                if passEnergy is False:
                    print(saved_fail)
                    PrintLICHEMDebug(runC)
                CleanFiles()
            except KeyboardInterrupt:
                SkipSequence("TIP3P energy:", pf_printed)

            pf_printed = False
            try:
                # New test instance
                testCt += 1
                print(f"\n{' ':2}Test {testCt}: AMOEBA/GK Energy")
                # Check MM energy
                line = ""
                passEnergy = False
                cmd = "cp pol.key tinker.key"
                PrintCopyDebug(cmd)
                subprocess.call(cmd, shell=True)  # Copy key file
                runC = RunLICHEM("waterdimer.xyz", "solvreg.inp",
                                 "watercon.inp")
                QMMMEnergy, savedEnergy, units = RecoverEnergy("MM energy:", 2)
                # Expected TINKER energy
                expected_energy = round(-1.2549403662026/har2eV, 3)
                # Check result
                if (QMMMEnergy == expected_energy):
                    # Check against saved energy
                    passEnergy = True
                else:
                    saved_fail = SaveFailure(QMMMEnergy, expected_energy)
                pf_printed = PrintPassFail("AMOEBA/GK energy:", passEnergy,
                                           updateResults, savedEnergy, units)
                if passEnergy is False:
                    print(saved_fail)
                    PrintLICHEMDebug(runC)
                CleanFiles()
            except KeyboardInterrupt:
                SkipSequence("AMOEBA/GK energy:", pf_printed)

        pf_printed = False
        # Back to tests non-specific to TINKER
        try:
            testCt += 1
            print(f"\n{' ':2}Test {testCt}: PBE0/TIP3P Energy")
            # Check QMMM point-charge energy results
            passEnergy = False
            cmd = "cp pchrg.key tinker.key"
            PrintCopyDebug(cmd)
            subprocess.call(cmd, shell=True)  # Copy key file
            runC = RunLICHEM("waterdimer.xyz", "pchrgreg.inp",
                             "watercon.inp")
            QMMMEnergy, savedEnergy, units = RecoverEnergy("QMMM energy:", 2)
            # Expected TINKER energy
            expected_energy = round(-2077.2021947277/har2eV, 3)
            # Check result
            if (QMPack == "PSI4"):
                # Check against saved energy
                if (QMMMEnergy == expected_energy):
                    passEnergy = True
                else:
                    saved_fail = SaveFailure(QMMMEnergy, expected_energy)
            if ((QMPack == "Gaussian") or (QMPack == "g16")):
                expected_energy = round(-2077.2018207808/har2eV, 3)
                # Check against saved energy
                if (QMMMEnergy == expected_energy):
                    passEnergy = True
                else:
                    saved_fail = SaveFailure(QMMMEnergy, expected_energy)
            if (QMPack == "NWChem"):
                expected_energy = round(-2077.2022117306/har2eV, 3)
                # Check against saved energy
                if (QMMMEnergy == expected_energy):
                    passEnergy = True
                else:
                    saved_fail = SaveFailure(QMMMEnergy, expected_energy)
            pf_printed = PrintPassFail("PBE0/TIP3P energy:", passEnergy,
                                       updateResults, savedEnergy, units)
            if passEnergy is False:
                print(saved_fail)
                PrintLICHEMDebug(runC)
            CleanFiles()
        except KeyboardInterrupt:
            SkipSequence("PBE0/TIP3P energy:", pf_printed)

        pf_printed = False
        try:
            testCt += 1
            print(f"\n{' ':2}Test {testCt}: PBE0/AMOEBA Energy")
            # Check QMMM polarizable energy results
            passEnergy = False
            cmd = "cp pol.key tinker.key"
            PrintCopyDebug(cmd)
            subprocess.call(cmd, shell=True)  # Copy key file
            runC = RunLICHEM("waterdimer.xyz", "polreg.inp", "watercon.inp")
            QMMMEnergy, savedEnergy, units = RecoverEnergy("QMMM energy:", 2)
            # Check result
            if (QMPack == "PSI4"):
                expected_energy = round(-2077.1114201829/har2eV, 3)
                # Check against saved energy
                if (QMMMEnergy == expected_energy):
                    passEnergy = True
                else:
                    saved_fail = SaveFailure(QMMMEnergy, expected_energy)
            if ((QMPack == "Gaussian") or (QMPack == "g16")):
                expected_energy = round(-2077.1090319595/har2eV, 3)
                # Check against saved energy
                if (QMMMEnergy == expected_energy):
                    passEnergy = True
                else:
                    saved_fail = SaveFailure(QMMMEnergy, expected_energy)
            if (QMPack == "NWChem"):
                expected_energy = round(-2077.1094168459/har2eV, 3)
                # Check against saved energy
                if (QMMMEnergy == expected_energy):
                    passEnergy = True
                else:
                    saved_fail = SaveFailure(QMMMEnergy, expected_energy)
            pf_printed = PrintPassFail("PBE0/AMOEBA energy:", passEnergy,
                                       updateResults, savedEnergy, units)
            if passEnergy is False:
                print(saved_fail)
                PrintLICHEMDebug(runC)
            CleanFiles()
        except KeyboardInterrupt:
            SkipSequence("PBE0/AMOEBA energy:", pf_printed)

        pf_printed = False
        try:
            # Gaussian and NWChem only, so update counter/def test in the loop
            # Check pseudobond optimizations
            if ((QMPack == "Gaussian") or (QMPack == "g16") or
               (QMPack == "NWChem")):
                testCt += 1
                print(f"\n{' ':2}Test {testCt}: DFP/Pseudobonds")
                # Carry out test
                passEnergy = False
                cmd = "cp pbopt.key tinker.key"
                PrintCopyDebug(cmd)
                subprocess.call(cmd, shell=True)  # Copy key file
                cmd = "cp pbbasis.txt BASIS"
                PrintCopyDebug(cmd)
                subprocess.call(cmd, shell=True)  # Copy BASIS set file
                runC = RunLICHEM("alkyl.xyz", "pboptreg.inp", "alkcon.inp")
                QMMMEnergy, savedEnergy, units = RecoverEnergy("Opt. step: 2",
                                                                            6)
                # Check result
                if ((QMPack == "Gaussian") or (QMPack == "g16")):
                    expected_energy = round(-3015.0548490566/har2eV, 3)
                    # Check against saved energy
                    if (QMMMEnergy == expected_energy):
                        passEnergy = True
                    else:
                        saved_fail = SaveFailure(QMMMEnergy, expected_energy)
                if (QMPack == "NWChem"):
                    expected_energy = round(-3015.2278310975/har2eV, 3)
                    # Check against saved energy
                    if (QMMMEnergy == expected_energy):
                        passEnergy = True
                    else:
                        saved_fail = SaveFailure(QMMMEnergy, expected_energy)
                pf_printed = PrintPassFail("DFP/Pseudobonds:", passEnergy,
                                           updateResults, savedEnergy, units)
                if passEnergy is False:
                    print(saved_fail)
                    PrintLICHEMDebug(runC)
                CleanFiles()
        except KeyboardInterrupt:
            SkipSequence("DFP/Pseudobonds:", pf_printed)

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
