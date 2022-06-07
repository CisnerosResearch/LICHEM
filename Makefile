###################################################
#                                                 #
#   LICHEM: Layered Interacting CHEmical Models   #
#                                                 #
#        Symbiotic Computational Chemistry        #
#                                                 #
###################################################

### Standard compiler settings ###

CXX= g++
CXXFLAGS= -static -O3 -fopenmp

UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
CXXFLAGS= -O3 -fopenmp
endif

### Install directory ### 
INSTALLBIN=/tmp/LICHEM1.1/bin

### Version Information ### 
MAKEDATE=2222-22-22
COMMIT=0a1b2c3

#############################
### LICHEM Makefile Usage ###
#############################
### General Usage:
### 1. make install        General LICHEM installation
### 2. make manual         Compile the LaTeX-based documentation
### 3. make tests          Clear pre-existing test output, remake test source
### 4. make module         Put path, commit, & date info in ./modules/*
### 5. make clean          Clear bin and manual
###
### Developer Usage:
###     ** These set additional flags in the process! **
### 1. make Dev            CPU-based LICHEM installation
### 2. make GPUDev         GPU-based LICHEM installation
### 3. make checksyntax    Check for syntax errors, print stats
### 4. make mantest        Only compile the manual
### 5. make mandel         Only remove the manual
### 6. make gitclean       Clear bin, tests, & manual, generalize INSTALLBIN
###
###############################################################################

######################################################
###  Main Environment Variables; Change as Needed  ###
######################################################

## LDFLAGS       Linker flags (for including additional libraries)
## PYPATH        Path to Python 3 executable
## SEDI          In-place flag for sed (GNU: -i, OSX: -i "")
## TEX           TeX engine (doc has been tested with pdflatex only!)
## BIB           BibTeX engine
## DEVFLAGS      Additional compiler flags for CPU-based install
## GPUFLAGS      Additional compiler flags for GPU-based install

# The local copy of Eigen is located in ./Eigen3/
LDFLAGS=-I./Eigen3/

PYPATH=$(shell which python3)
#PYPATH=/usr/bin/python

SEDI=-i

### LaTeX settings ###
# For OSX these can be replaced with dummy calls to cat
TEX=pdflatex
BIB=bibtex

### Advanced compiler settings for developers ###
DEVFLAGS=-g -Wall -std=c++14
GPUFLAGS=-fopenacc

############################################################
### You Should Not Need to Change Things Below This Line ###
############################################################

################################################
###  Compile Rules for Users and Developers  ###
################################################

# NB: By definition, these are written with a tab after the colon

install:	title binary testexe compdone

Dev:	title devbin testexe manual stats compdone

GPUDev:	title gpubin testexe manual stats compdone

tests:	title deltests testexe

clean:	title delbin compclean

gitclean:	title delbin deltests gitmk

##################################
### Combine Settings Variables ###
##################################

# NB: Do not modify this section!

FLAGSBIN=$(CXXFLAGS) $(LDFLAGS) -I./src/ -I./include/
FLAGSBINMPI=$(CXXFLAGS) $(LDFLAGS) -I./src/ -I./src/mpi/ -I./include/
FLAGSDEV=$(CXXFLAGS) $(DEVFLAGS) $(LDFLAGS) -I./src/ -I./include/
FLAGSGPU=$(CXXFLAGS) $(DEVFLAGS) $(GPUFLAGS) $(LDFLAGS) -I./src/ -I./include/

#####################################################
### Rules for building various parts of the code ###
#####################################################

devbin:
	@echo ""; \
	echo "### Compiling the LICHEM development binary ###"; \
	mkdir -p $(INSTALLBIN)
	$(CXX) ./src/LICHEM.cpp -o $(INSTALLBIN)/lichem $(FLAGSDEV)

gpubin:
	@echo ""; \
	echo "### Compiling the LICHEM GPU binary ###"; \
	mkdir -p $(INSTALLBIN)
	$(CXX) ./src/LICHEM.cpp -o $(INSTALLBIN)/lichem $(FLAGSGPU)

testexe:
	@echo ""; \
	echo "### Creating test suite executable ###"
	@echo 'echo "#!$(PYPATH)" > ./tests/runtests'; \
	echo "!!$(PYPATH)" > ./tests/runtests
	cat ./src/runtests.py >> ./tests/runtests
	@sed $(SEDI) 's/\#.*//g' ./tests/runtests; \
	sed $(SEDI) 's/\s*$$//g' ./tests/runtests; \
	sed $(SEDI) '/^$$/d' ./tests/runtests; \
	sed $(SEDI) 's/\!\!/\#\!/g' ./tests/runtests; \
	chmod a+x ./tests/runtests
	@echo ""

checksyntax:	title
	@echo ""; \
	echo "### Checking for warnings and syntax errors ###"
	$(CXX) -fsyntax-only ./src/LICHEM.cpp -o lichem $(FLAGSDEV)
	@echo ""; \
	echo "### Source code statistics ###"; \
	echo "Number of LICHEM source code files:"; \
	ls -al include/* src/* | wc -l; \
	echo "Total length of LICHEM (lines):"; \
	grep '' -aR ./src ./include | wc -l; \
	echo ""

manual:
	@echo ""; \
	echo "### Compiling the documentation ###"; \
	cd src/; \
	echo "$(TEX) manual"; \
	$(TEX) manual > doclog.txt; \
	$(BIB) manual > doclog.txt; \
	$(TEX) manual > doclog.txt; \
	$(TEX) manual > doclog.txt; \
	$(TEX) manual > doclog.txt; \
	$(TEX) manual > doclog.txt; \
	echo "$(BIB) manual"; \
	$(BIB) manual > doclog.txt; \
	$(TEX) manual > doclog.txt; \
	$(TEX) manual > doclog.txt; \
	$(TEX) manual > doclog.txt; \
	mv manual.pdf ../doc/LICHEM_manual.pdf; \
	rm -f manual.aux manual.bbl manual.blg; \
	rm -f manual.log manual.out manual.toc; \
	rm -f manual.lof manual.lot; \
	rm -f doclog.txt manual.synctex.gz acs-manual.bib

mantest:
	@echo ""; \
	echo "### Compiling the documentation ###"; \
	cd src/; \
	echo "$(TEX) manual"; \
	$(TEX) manual > doclog.txt; \
	$(BIB) manual > doclog.txt; \
	$(TEX) manual > doclog.txt; \
	$(TEX) manual > doclog.txt; \
	$(TEX) manual > doclog.txt; \
	$(TEX) manual > doclog.txt; \
	echo "$(BIB) manual"; \
	$(BIB) manual > doclog.txt; \
	$(TEX) manual > doclog.txt; \
	$(TEX) manual > doclog.txt; \
	$(TEX) manual > doclog.txt; \
	mv manual.pdf ../doc/LICHEM_manual.pdf

mandel:
	@echo ""; \
	echo "### Removing the LICHEM documentation ###"; \
	echo ""; \
	rm -rf ./src/manual.pdf ./doc/LICHEM_manual.pdf; \
	rm -rf ./src/manual.aux ./src/manual.bbl ./src/manual.blg; \
	rm -rf ./src/manual.log ./src/manual.out ./src/manual.toc; \
	rm -rf ./src/manual.lof ./src/manual.lot; \
	rm -rf ./src/doclog.txt ./src/manual.synctex.gz ./src/acs-manual.bib

module:
	@echo ""; \
	echo "### Adding path to module templates ###"; \
	sed $(SEDI) 's:MOD_FIX_VERS_DATE:$(MAKEDATE):g' ./modules/lichem.lua; \
	sed $(SEDI) 's:MOD_FIX_COMMIT:$(COMMIT):g' ./modules/lichem.lua; \
	sed $(SEDI) 's:MOD_FIX_PATH:$(INSTALLBIN):g' ./modules/lichem.lua; \
	sed $(SEDI) 's:MOD_FIX_VERS_DATE:$(MAKEDATE):g' ./modules/lichem.tcl; \
	sed $(SEDI) 's:MOD_FIX_COMMIT:$(COMMIT):g' ./modules/lichem.tcl; \
	sed $(SEDI) 's:MOD_FIX_PATH:$(INSTALLBIN):g' ./modules/lichem.tcl; \
	echo ""; \
	echo "Remember to add or copy the module to your module path!"; \
	echo "  Use '.lua' for LMOD modules and '.tcl' for environment modules."; \
	echo "";\
	echo "And don't forget to modify the file to load the QM and MM packages!";\
	echo ""

title:
	@echo ""; \
	echo "###################################################"; \
	echo "#                                                 #"; \
	echo "#   LICHEM: Layered Interacting CHEmical Models   #"; \
	echo "#                                                 #"; \
	echo "#        Symbiotic Computational Chemistry        #"; \
	echo "#                                                 #"; \
	echo "###################################################"

stats:
	@echo ""; \
        echo "### Source code statistics ###"; \
	echo "Number of LICHEM source code files:"; \
	ls -al include/* src/* | wc -l; \
	echo "Total length of LICHEM (lines):"; \
	grep '' -aR ./src ./include | wc -l; \
	echo "";

compdone:
	@echo ""; \
	echo "Done."; \
	echo ""
	@echo ""
	@echo "Installation complete"
	@echo "Please add lichem executable directory to your path:"
	@echo ""
	@echo "     export PATH="$(INSTALLBIN)":\$$PATH"
	@echo ""

compclean:
	@echo ""; \
	echo "Done."; \
	echo ""
	@echo ""
	@echo "LICHEM has been uninstalled."
	@echo "Please remove the deleted executable directory from your path:"
	@echo ""
	@echo "     "$(INSTALLBIN)""
	@echo ""

vroom:
	@echo ""; \
	if grep -q "JOKES = 1" include/LICHEM_options.h; then \
	echo '     ___'; \
	echo '    |_  |'; \
	echo '      \ \'; \
	echo '      |\ \'; \
	echo '      | \ \'; \
	echo '      \  \ \'; \
	echo '       \  \ \'; \
	echo '        \  \ \       <wrrr vroooom wrrr> '; \
	echo '         \__\ \________'; \
	echo '             |_________\'; \
	echo '             |__________|  ..,  ,.,. .,.,, ,..'; \
	echo ""; \
 	fi; \
        echo "";

delbin:	vroom
	@echo "Removing binary and manual..."; \
	rm -rf lichem ./doc/LICHEM_manual.pdf ./tests/runtests $(INSTALLBIN)

deltests: vroom
	@echo ""; \
	echo "Removing any output from previous runtests."; \
	rm -rf ./tests/*/LICH* ./tests/*/trash.xyz ./tests/*/tests.out; \
	rm -rf ./tests/*_TINKER/BeadStartStruct.xyz; \
	rm -rf ./tests/*_TINKER/BurstStruct.xyz; \
	rm -rf ./tests/*_TINKER/tinker.key; \

# Arbitrarily assign MAKEDATE and COMMIT in Makefile so you see actual
# updates in commit history.
gitmk:
	@echo "";\
	echo "Preparing the directory for a git commit!"; \
	sed $(SEDI) 's:$(INSTALLBIN):/tmp/LICHEM1.1/bin:g' ./Makefile; \
	sed $(SEDI) 's:$(MAKEDATE):2222-22-22:g' ./Makefile; \
	sed $(SEDI) 's:$(COMMIT):0a1b2c3:g' ./Makefile; \
	sed $(SEDI) 's/QM_type: g16/QM_type: Gaussian/g' \
		./tests/Gau_TINKER/*reg.inp; \
	sed $(SEDI) 's:$(INSTALLBIN):MOD_FIX_PATH:g' ./modules/lichem.lua; \
	sed $(SEDI) 's:$(MAKEDATE):MOD_FIX_VERS_DATE:g' ./modules/lichem.lua; \
	sed $(SEDI) 's:$(COMMIT):MOD_FIX_VERS_DATE:g' ./modules/lichem.lua; \
	sed $(SEDI) 's:$(INSTALLBIN):MOD_FIX_PATH:g' ./modules/lichem.tcl; \
	sed $(SEDI) 's:$(MAKEDATE):MOD_FIX_VERS_DATE:g' ./modules/lichem.tcl; \
	sed $(SEDI) 's:$(COMMIT):MOD_FIX_VERS_DATE:g' ./modules/lichem.tcl; \
	echo ""

# NB: binary MUST be defined last because ./configure appends
#     compiler-specific info!
binary:
	@echo ""; \
	echo "### Compiling the LICHEM binary ###"; \
	mkdir -p $(INSTALLBIN)
	$(CXX) ./src/LICHEM.cpp -o $(INSTALLBIN)/lichem $(FLAGSBIN)
	@strip $(INSTALLBIN)/lichem
