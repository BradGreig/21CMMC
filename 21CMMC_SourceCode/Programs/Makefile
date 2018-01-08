# C compiler and flags
CPPFLAGS = -I/opt/local/include
LDFLAGS = -L/opt/local/lib -lgsl -lgslcblas -lfftw3f_omp -lfftw3f -lm
CC      = /usr/local/bin/gcc -fopenmp -Ofast

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
OBJ_FILES = perturb_field \
	drive_21cmMC_streamlined\
	SplitMockObservation \
	CreateFcollTable\
	CreateSmoothedDensityBoxes \
	ComputingTau_e \
	create_dens_boxes_for_LC \
	Createfcoll_ionisation_LC

#########################################################################

init:	init.c \
	filter.c \
	${OBJ_FILES} \
	${COSMO_FILES}

	${CC} ${CPPFLAGS} -o init init.c ${LDFLAGS}

perturb_field:	perturb_field.c \
	${COSMO_FILES}

	${CC} ${CPPFLAGS} -o perturb_field perturb_field.c ${LDFLAGS}

drive_21cmMC_streamlined: drive_21cmMC_streamlined.c \
	${COSMO_FILES} \

	${CC} ${CPPFLAGS} -o drive_21cmMC_streamlined drive_21cmMC_streamlined.c ${LDFLAGS}

CreateFcollTable: CreateFcollTable.c \
	${COSMO_FILES} \

	${CC} ${CPPFLAGS} -o CreateFcollTable CreateFcollTable.c ${LDFLAGS}

ComputingTau_e: ComputingTau_e.c \
	${COSMO_FILES} \

	${CC} ${CPPFLAGS} -o ComputingTau_e ComputingTau_e.c ${LDFLAGS}

CreateSmoothedDensityBoxes: CreateSmoothedDensityBoxes.c \
	${COSMO_FILES} \

	${CC} ${CPPFLAGS} -o CreateSmoothedDensityBoxes CreateSmoothedDensityBoxes.c ${LDFLAGS}

SplitMockObservation: SplitMockObservation.c \
	${COSMO_FILES} \

	${CC} ${CPPFLAGS} -o SplitMockObservation SplitMockObservation.c ${LDFLAGS}

SplitMockObservation_Different: SplitMockObservation_Different.c \
	${COSMO_FILES} \

	${CC} ${CPPFLAGS} -o SplitMockObservation_Different SplitMockObservation_Different.c ${LDFLAGS}

create_dens_boxes_for_LC: create_dens_boxes_for_LC.c \
	${COSMO_FILES} \

	${CC} ${CPPFLAGS} -o create_dens_boxes_for_LC create_dens_boxes_for_LC.c ${LDFLAGS}

Createfcoll_ionisation_LC: Createfcoll_ionisation_LC.c \
	${COSMO_FILES} \

	${CC} ${CPPFLAGS} -o Createfcoll_ionisation_LC Createfcoll_ionisation_LC.c ${LDFLAGS}

clean:
	rm ${OBJ_FILES} *~
