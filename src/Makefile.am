# EXTRA_DIST = doc bin examples


DEFS = $(REVISION)
bin_PROGRAMS = 
if WANT_PHYML
bin_PROGRAMS += phyml
endif
if WANT_PHYTIME
bin_PROGRAMS += phytime
endif
if WANT_RF
bin_PROGRAMS += rf
endif
if WANT_MPI
bin_PROGRAMS += phyml-mpi
endif
if WANT_WIN
bin_PROGRAMS += phyml-windows
endif
if WANT_TEST
bin_PROGRAMS += test
endif
if WANT_EVOLVE
bin_PROGRAMS += evolve
endif
if WANT_PHYREX
bin_PROGRAMS += phyrex
endif
if WANT_PHYREXSIM
bin_PROGRAMS += phyrexsim
endif
if WANT_DATE
bin_PROGRAMS += date
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
draw.c draw.h\
rates.c rates.h\
mcmc.c mcmc.h\
stats.c stats.h\
mg.c mg.h\
io.c io.h\
make.c make.h\
mixt.c mixt.h\
init.c init.h\
nexus.c nexus.h\
date.c date.h\
tbe.c tbe.h\
sse.c sse.h\
avx.c avx.h\
ancestral.c ancestral.h\
xml.c xml.h
phytime_CFLAGS = -DPHYTIME $(AM_CFLAGS)
endif
# if WANT_RWRAP
# lib_LTLIBRARIES = librwrap.la
# librwrap_la_SOURCES = main.c \
# utilities.c  utilities.h\
# optimiz.c  optimiz.h\
# lk.c  lk.h\
# bionj.c bionj.h\
# models.c  models.h\
# free.c  free.h\
# help.c  help.h\
# simu.c  simu.h\
# eigen.c  eigen.h\
# pars.c  pars.h\
# alrt.c  alrt.h\
# interface.c  interface.h\
# cl.c  cl.h\
# spr.c spr.h\
# draw.c  draw.h\
# stats.c stats.h\
# rates.c rates.h\
# mcmc.c mcmc.h\
# times.c times.h\
# tiporder.c tiporder.h\
# mg.c mg.h\
# io.c io.h\
# make.c make.h\
# nexus.c nexus.h\
# init.c init.h\
# xml.c xml.h\
# mixt.c mixt.h\
# rwrapper.c rwrapper.h
# librwrap_la_LDFLAGS = -I/usr/share/R/include -shared -module -flat_namespace
# librwrap_la_CFLAGS=-std=gnu99 -fPIC -Wl,-z,defs
# else 
if WANT_RF
rf_SOURCES = main.c \
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
mg.c mg.h\
io.c io.h\
make.c make.h\
nexus.c nexus.h\
init.c init.h\
xml.c xml.h\
mixt.c mixt.h\
tbe.c tbe.h \
ancestral.c ancestral.h\
sse.c sse.h\
avx.c avx.h
rf_CFLAGS = -DRF $(AM_CFLAGS)
endif
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
io.c io.h\
make.c make.h\
nexus.c nexus.h\
init.c init.h\
xml.c xml.h\
mixt.c mixt.h\
tbe.c tbe.h \
sse.c sse.h\
avx.c avx.h\
ancestral.c ancestral.h\
mpi_boot.c mpi_boot.h
phyml_mpi_CFLAGS = -DPHYML $(AM_CFLAGS)
endif
if WANT_WIN
phyml_windows_SOURCES = main.c \
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
tbe.c tbe.h \
stats.c stats.h\
rates.c rates.h\
mcmc.c mcmc.h\
times.c times.h\
mg.c mg.h\
io.c io.h\
make.c make.h\
nexus.c nexus.h\
init.c init.h\
xml.c xml.h\
mixt.c mixt.h\
sse.c sse.h\
avx.c avx.h\
ancestral.c ancestral.h\
date.c date.h
endif
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
tbe.c tbe.h \
stats.c stats.h\
rates.c rates.h\
mcmc.c mcmc.h\
times.c times.h\
mg.c mg.h\
io.c io.h\
make.c make.h\
nexus.c nexus.h\
init.c init.h\
xml.c xml.h\
mixt.c mixt.h\
sse.c sse.h\
avx.c avx.h\
ancestral.c ancestral.h\
date.c date.h
test_CFLAGS = -DTEST $(AM_CFLAGS)
endif
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
tbe.c tbe.h \
stats.c stats.h\
rates.c rates.h\
mcmc.c mcmc.h\
times.c times.h\
mg.c mg.h\
io.c io.h\
make.c make.h\
nexus.c nexus.h\
init.c init.h\
xml.c xml.h\
mixt.c mixt.h
evolve_CFLAGS = -DEVOLVE $(AM_CFLAGS)
endif
if WANT_PHYREX
phyrex_SOURCES = main.c\
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
io.c io.h\
make.c make.h\
nexus.c nexus.h\
init.c init.h\
xml.c xml.h\
mixt.c mixt.h\
date.c date.h\
tbe.c tbe.h \
sse.c sse.h\
avx.c avx.h\
ancestral.c ancestral.h\
phyrex.c phyrex.h\
location.c location.h\
slfv.c slfv.h\
rw.c rw.h\
rrw.c rrw.h\
ibm.c ibm.h
phyrex_CFLAGS = -DPHYREX $(AM_CFLAGS)
endif
if WANT_PHYREXSIM
phyrexsim_SOURCES = main.c\
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
mg.c mg.h\
io.c io.h\
make.c make.h\
nexus.c nexus.h\
init.c init.h\
xml.c xml.h\
mixt.c mixt.h\
date.c date.h\
tbe.c tbe.h \
sse.c sse.h\
avx.c avx.h\
ancestral.c ancestral.h\
phyrex.c phyrex.h
phyrexsim_CFLAGS = -DPHYREXSIM $(AM_CFLAGS)
endif
if WANT_DATE
date_SOURCES = main.c \
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
mg.c mg.h\
io.c io.h\
make.c make.h\
nexus.c nexus.h\
init.c init.h\
xml.c xml.h\
tbe.c tbe.h \
mixt.c mixt.h\
date.c date.h
date_CFLAGS = -DDATE $(AM_CFLAGS)
endif
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
io.c io.h\
make.c make.h\
nexus.c nexus.h\
init.c init.h\
xml.c xml.h\
sse.c sse.h\
avx.c avx.h\
tbe.c tbe.h \
ancestral.c ancestral.h \
mixt.c mixt.h
phyml_CFLAGS = -DPHYML $(AM_CFLAGS)


all-am:	intro $(bin_PROGRAMS)
	@echo ""
	@echo "Done."

intro:	
	@echo ""
	@echo ""
	@echo ".:  Building [$(bin_PROGRAMS)]. Version $(VERSION) :."
	@echo ""
	@echo ""
