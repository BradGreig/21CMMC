# C compiler and flags
CPPFLAGS = -I/usr/local/include
LDFLAGS = -lgsl -lgslcblas -lfftw3f_omp -lfftw3f -lm
CC      = gcc -fopenmp -Ofast

# Directory hosting the non-project specific code
COSMO_DIR = ../Cosmo_c_files

# Directory containing parameter files, definint relevant constants
PARAMETER_DIR = ../Parameter_files

# non-project dependencies
COSMO_FILES = 	${PARAMETER_DIR}/COSMOLOGY.H \
	${COSMO_DIR}/ps.c \
	${COSMO_DIR}/cosmo_progs.c \
	${COSMO_DIR}/misc.c \
	${PARAMETER_DIR}/INIT_PARAMS.H \
	${PARAMETER_DIR}/ANAL_PARAMS.H \
	${PARAMETER_DIR}/HEAT_PARAMS.H \

# object files
OBJ_FILES = drive_21cmMC_streamlined
#########################################################################

drive_21cmMC_streamlined: drive_21cmMC_streamlined.c \
	${COSMO_FILES} \

	${CC} ${CPPFLAGS} -o drive_21cmMC_streamlined drive_21cmMC_streamlined.c ${LDFLAGS}

clean:
	rm ${OBJ_FILES} *~
