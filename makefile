ejec:=mod_
srcs:=$(wildcard *.cpp *.c *.cc *.cxx *.htp *.obj)
rscript:='source(file.path("plotFigures.R"))'

modelo := $(foreach program,$(ejec),$(program)$(modl))
datafile := $(foreach program,$(ejec),$(program)$(modl).dat)

all: $(modelo) cleanBuild

$(modelo): $(modelo).tpl
	admb $(modelo)

run: $(modelo)
	 ./$(modelo)

runmcmc: $(modelo)
	     ./$(modelo) -mcmc -mcsave 100
	     ./$(modelo) -mceval

cleanBuild: 
	@echo 'removing build files'
	@rm -rf *.c *.cc *.cxx *.htp *.obj

cleanRun: 
	@echo 'removing running files'
	@rm -rf $(modelo).[bcehlmprs]* sims variance admodel.* fmin.log

cleanExe: 
	@echo 'removing executable'
	@rm -rf $(modelo)


plots: $(modelo)
	-@echo $(rscript) | R --vanilla --slave

