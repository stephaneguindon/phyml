# EXTRA_DIST = doc bin examples


DEFS = $(REVISION)

if WANT_PHYTIME
bin_PROGRAMS = phytime
PROG = PHYTIME
else
if WANT_PHYCONT
bin_PROGRAMS = phycont
PROG = PHYCONT
else
if WANT_PART
bin_PROGRAMS = part
PROG = PART
else
if WANT_RWRAP
PROG = RWRAP
else
if WANT_TIPORDER
bin_PROGRAMS = tiporder
PROG = TIPORDER
else 
if WANT_M4
bin_PROGRAMS = m4
PROG = M4
else 
if WANT_RF
bin_PROGRAMS = rf
PROG = RF
else 
if WANT_MPI
bin_PROGRAMS = phyml-mpi
PROG = PHYML
else
if WANT_TEST
bin_PROGRAMS = test
PROG = TEST
else
if WANT_INVITEE
bin_PROGRAMS = invitee
PROG = INVITEE
else
if WANT_GEO
bin_PROGRAMS = phylogeo
PROG = GEO
else
if WANT_EVOLVE
bin_PROGRAMS = evolve
PROG = EVOLVE
else
if WANT_CHECKPOINT
bin_PROGRAMS = checkpoint
PROG = CHECKPOINT
else
if WANT_BEAGLE
bin_PROGRAMS = phyml-beagle
PROG = PHYML
else
if WANT_MIGREP
bin_PROGRAMS = migrep
PROG = MIGREP
else
bin_PROGRAMS = phyml
PROG = PHYML
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif

if WANT_PHYTIME
phytime_SOURCES = main.c \
utilities.c utilities.h\
optimiz.c optimiz.h\
lk.c lk.h\
bionj.c bionj.h\
models.c models.h\
free.c free.h\
help.c help.h\
simu.c simu.h\
eigen.c eigen.h\
pars.c pars.h\
alrt.c alrt.h\
interface.c interface.h\
cl.c cl.h\
spr.c spr.h\
times.c times.h\
m4.c m4.h\
draw.c draw.h\
rates.c rates.h\
mcmc.c mcmc.h\
stats.c stats.h\
mg.c mg.h\
tiporder.c tiporder.h\
io.c io.h\
make.c make.h\
mixt.c mixt.h\
init.c init.h\
nexus.c nexus.h\
xml.c xml.h\
migrep.c migrep.h
phytime_LDADD = -lm
else 
if WANT_PHYCONT
phycont_SOURCES = main.c \
utilities.c utilities.h\
optimiz.c optimiz.h\
lk.c lk.h\
bionj.c bionj.h\
models.c models.h\
free.c free.h\
help.c help.h\
simu.c simu.h\
eigen.c eigen.h\
pars.c pars.h\
alrt.c alrt.h\
interface.c interface.h\
cl.c cl.h\
spr.c spr.h\
times.c times.h\
m4.c m4.h\
draw.c draw.h\
rates.c rates.h\
mcmc.c mcmc.h\
stats.c stats.h\
mg.c mg.h\
tiporder.c tiporder.h
# continuous.c continuous.h
phycont_LDADD = -lm
else 
if WANT_RWRAP
lib_LTLIBRARIES = librwrap.la
librwrap_la_SOURCES = main.c \
utilities.c  utilities.h\
optimiz.c  optimiz.h\
lk.c  lk.h\
bionj.c bionj.h\
models.c  models.h\
free.c  free.h\
help.c  help.h\
simu.c  simu.h\
eigen.c  eigen.h\
pars.c  pars.h\
alrt.c  alrt.h\
interface.c  interface.h\
cl.c  cl.h\
spr.c spr.h\
draw.c  draw.h\
stats.c stats.h\
rates.c rates.h\
mcmc.c mcmc.h\
times.c times.h\
tiporder.c tiporder.h\
mg.c mg.h\
m4.c m4.h\
io.c io.h\
make.c make.h\
nexus.c nexus.h\
init.c init.h\
xml.c xml.h\
mixt.c mixt.h\
rwrapper.c rwrapper.h
librwrap_la_LIBADD = -lm
librwrap_la_LDFLAGS = -I/usr/share/R/include -shared -module -flat_namespace
librwrap_la_CFLAGS=-std=gnu99 -fPIC -Wl,-z,defs
else 
if WANT_PART
part_SOURCES = main.c \
utilities.c utilities.h\
optimiz.c optimiz.h\
lk.c lk.h\
bionj.c bionj.h\
models.c models.h\
free.c free.h\
help.c help.h\
simu.c simu.h\
eigen.c eigen.h\
pars.c pars.h\
alrt.c alrt.h\
interface.c interface.h\
cl.c cl.h\
mg.c mg.h\
spr.c spr.h\
m4.c m4.h\
draw.c draw.h\
stats.c stats.h\
tiporder.c tiporder.h
part_LDADD = -lm
else
if WANT_TIPORDER
tiporder_SOURCES = main.c \
utilities.c utilities.h\
optimiz.c optimiz.h\
lk.c lk.h\
bionj.c bionj.h\
models.c models.h\
free.c free.h\
help.c help.h\
simu.c simu.h\
eigen.c eigen.h\
pars.c pars.h\
alrt.c alrt.h\
interface.c interface.h\
cl.c cl.h\
mg.c mg.h\
times.c times.h\
mcmc.c mcmc.h\
rates.c rates.h\
spr.c spr.h\
m4.c m4.h\
draw.c draw.h\
stats.c stats.h\
tiporder.c tiporder.h
tiporder_LDADD = -lm
else
if WANT_M4
m4_SOURCES = main.c \
utilities.c utilities.h\
optimiz.c optimiz.h\
lk.c lk.h\
bionj.c bionj.h\
models.c models.h\
free.c free.h\
help.c help.h\
simu.c simu.h\
eigen.c eigen.h\
pars.c pars.h\
alrt.c alrt.h\
interface.c interface.h\
cl.c cl.h\
mg.c mg.h\
times.c times.h\
mcmc.c mcmc.h\
rates.c rates.h\
spr.c spr.h\
draw.c draw.h\
stats.c stats.h\
tiporder.c tiporder.h\
m4.c m4.h
m4_LDADD = -lm
else 
if WANT_RF
rf_SOURCES = main.c \
utilities.c utilities.h\
optimiz.c optimiz.h\
lk.c lk.h\
bionj.c bionj.h\
models.c models.h\
free.c free.h\
help.c help.h\
simu.c simu.h\
eigen.c eigen.h\
pars.c pars.h\
alrt.c alrt.h\
interface.c interface.h\
cl.c cl.h\
mg.c mg.h\
times.c times.h\
mcmc.c mcmc.h\
rates.c rates.h\
spr.c spr.h\
draw.c draw.h\
stats.c stats.h\
tiporder.c tiporder.h\
m4.c m4.h
rf_LDADD = -lm
else
if WANT_MPI
phyml_mpi_SOURCES = main.c \
utilities.c  utilities.h\
optimiz.c  optimiz.h\
lk.c  lk.h\
bionj.c bionj.h\
models.c  models.h\
free.c  free.h\
help.c  help.h\
simu.c  simu.h\
eigen.c  eigen.h\
pars.c  pars.h\
alrt.c  alrt.h\
interface.c  interface.h\
cl.c  cl.h\
spr.c spr.h\
draw.c  draw.h\
stats.c stats.h\
rates.c rates.h\
mcmc.c mcmc.h\
times.c times.h\
tiporder.c tiporder.h\
mg.c mg.h\
m4.c m4.h\
io.c io.h\
make.c make.h\
nexus.c nexus.h\
init.c init.h\
xml.c xml.h\
mixt.c mixt.h\
mpi_boot.c mpi_boot.h
phyml_mpi_LDADD = -lm
else
if WANT_TEST
test_SOURCES = main.c \
utilities.c  utilities.h\
optimiz.c  optimiz.h\
lk.c  lk.h\
bionj.c bionj.h\
models.c  models.h\
free.c  free.h\
help.c  help.h\
simu.c  simu.h\
eigen.c  eigen.h\
pars.c  pars.h\
alrt.c  alrt.h\
interface.c  interface.h\
cl.c  cl.h\
spr.c spr.h\
draw.c  draw.h\
stats.c stats.h\
rates.c rates.h\
mcmc.c mcmc.h\
times.c times.h\
tiporder.c tiporder.h\
mg.c mg.h\
m4.c m4.h\
io.c io.h\
make.c make.h\
nexus.c nexus.h\
init.c init.h\
xml.c xml.h\
mixt.c mixt.h
test_LDADD = -lm
else
if WANT_INVITEE
invitee_SOURCES = main.c \
utilities.c  utilities.h\
optimiz.c  optimiz.h\
lk.c  lk.h\
bionj.c bionj.h\
models.c  models.h\
free.c  free.h\
help.c  help.h\
simu.c  simu.h\
eigen.c  eigen.h\
pars.c  pars.h\
alrt.c  alrt.h\
interface.c  interface.h\
cl.c  cl.h\
spr.c spr.h\
draw.c  draw.h\
stats.c stats.h\
rates.c rates.h\
mcmc.c mcmc.h\
times.c times.h\
tiporder.c tiporder.h\
mg.c mg.h\
m4.c m4.h\
io.c io.h\
make.c make.h\
nexus.c nexus.h\
init.c init.h\
xml.c xml.h\
invitee.c invitee.h\
mixt.c mixt.h
invitee_LDADD = -lm 
else
if WANT_EVOLVE
evolve_SOURCES = main.c\
utilities.c  utilities.h\
optimiz.c  optimiz.h\
lk.c  lk.h\
bionj.c bionj.h\
models.c  models.h\
free.c  free.h\
help.c  help.h\
simu.c  simu.h\
eigen.c  eigen.h\
pars.c  pars.h\
alrt.c  alrt.h\
interface.c  interface.h\
cl.c  cl.h\
spr.c spr.h\
draw.c  draw.h\
stats.c stats.h\
rates.c rates.h\
mcmc.c mcmc.h\
times.c times.h\
tiporder.c tiporder.h\
mg.c mg.h\
m4.c m4.h\
io.c io.h\
make.c make.h\
nexus.c nexus.h\
init.c init.h\
xml.c xml.h\
mixt.c mixt.h
evolve_LDADD = -lm 
else
if WANT_GEO
phylogeo_SOURCES = main.c\
utilities.c  utilities.h\
optimiz.c  optimiz.h\
lk.c  lk.h\
bionj.c bionj.h\
models.c  models.h\
free.c  free.h\
help.c  help.h\
simu.c  simu.h\
eigen.c  eigen.h\
pars.c  pars.h\
alrt.c  alrt.h\
interface.c  interface.h\
cl.c  cl.h\
spr.c spr.h\
draw.c  draw.h\
stats.c stats.h\
rates.c rates.h\
mcmc.c mcmc.h\
times.c times.h\
tiporder.c tiporder.h\
mg.c mg.h\
m4.c m4.h\
io.c io.h\
make.c make.h\
nexus.c nexus.h\
init.c init.h\
xml.c xml.h\
mixt.c mixt.h\
geo.c geo.h
phylogeo_LDADD = -lm 
else
if WANT_CHECKPOINT
checkpoint_SOURCES = main.c\
utilities.c  utilities.h\
optimiz.c  optimiz.h\
lk.c  lk.h\
bionj.c bionj.h\
models.c  models.h\
free.c  free.h\
help.c  help.h\
simu.c  simu.h\
eigen.c  eigen.h\
pars.c  pars.h\
alrt.c  alrt.h\
interface.c  interface.h\
cl.c  cl.h\
spr.c spr.h\
draw.c  draw.h\
stats.c stats.h\
rates.c rates.h\
mcmc.c mcmc.h\
times.c times.h\
tiporder.c tiporder.h\
mg.c mg.h\
m4.c m4.h\
io.c io.h\
make.c make.h\
nexus.c nexus.h\
init.c init.h\
xml.c xml.h\
mixt.c mixt.h\
checkpoint.c checkpoint.h
checkpoint_LDADD = -lm
else
if WANT_MIGREP
migrep_SOURCES = main.c\
utilities.c  utilities.h\
optimiz.c  optimiz.h\
lk.c  lk.h\
bionj.c bionj.h\
models.c  models.h\
free.c  free.h\
help.c  help.h\
simu.c  simu.h\
eigen.c  eigen.h\
pars.c  pars.h\
alrt.c  alrt.h\
interface.c  interface.h\
cl.c  cl.h\
spr.c spr.h\
draw.c  draw.h\
stats.c stats.h\
rates.c rates.h\
mcmc.c mcmc.h\
times.c times.h\
tiporder.c tiporder.h\
mg.c mg.h\
m4.c m4.h\
io.c io.h\
make.c make.h\
nexus.c nexus.h\
init.c init.h\
xml.c xml.h\
mixt.c mixt.h\
migrep.c migrep.h
else
if WANT_BEAGLE
phyml_beagle_SOURCES = main.c \
utilities.c  utilities.h\
optimiz.c  optimiz.h\
lk.c  lk.h\
bionj.c bionj.h\
models.c  models.h\
free.c  free.h\
help.c  help.h\
simu.c  simu.h\
eigen.c  eigen.h\
pars.c  pars.h\
alrt.c  alrt.h\
interface.c  interface.h\
cl.c  cl.h\
spr.c spr.h\
draw.c  draw.h\
stats.c stats.h\
rates.c rates.h\
mcmc.c mcmc.h\
times.c times.h\
tiporder.c tiporder.h\
mg.c mg.h\
m4.c m4.h\
io.c io.h\
make.c make.h\
nexus.c nexus.h\
init.c init.h\
xml.c xml.h\
mixt.c mixt.h\
beagle_utils.c beagle_utils.h
phyml_beagle_LDADD = -lm -lhmsbeagle
else
phyml_SOURCES = main.c \
utilities.c  utilities.h\
optimiz.c  optimiz.h\
lk.c  lk.h\
bionj.c bionj.h\
models.c  models.h\
free.c  free.h\
help.c  help.h\
simu.c  simu.h\
eigen.c  eigen.h\
pars.c  pars.h\
alrt.c  alrt.h\
interface.c  interface.h\
cl.c  cl.h\
spr.c spr.h\
draw.c  draw.h\
stats.c stats.h\
rates.c rates.h\
mcmc.c mcmc.h\
times.c times.h\
tiporder.c tiporder.h\
mg.c mg.h\
m4.c m4.h\
io.c io.h\
make.c make.h\
nexus.c nexus.h\
init.c init.h\
xml.c xml.h\
mixt.c mixt.h
phyml_LDADD = -lm
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif
endif


all-am:	intro $(bin_PROGRAMS)
	@echo ""
	@echo "Done."

intro:	
	@echo ""
	@echo ""
	@echo ".: Building [$(bin_PROGRAMS)]. Version $(VERSION) :."
	@echo ""
	@echo ""

