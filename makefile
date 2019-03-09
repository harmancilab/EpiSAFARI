all: EpiSAFARI 

CC = g++
comp_flags = -c -O2 -Wall
gsl_flags = -lgsl -lgslcblas -lz
exec_name = ./bin/EpiSAFARI
LIB_DIR = ./src

# Define pattern rule for building object files.
%.o: %.cpp
	@echo Compiling $@
	@${CC} ${comp_flags} $< -o $@

%.o: %.c
	@echo Compiling $@
	@${CC} ${comp_flags} $< -o $@

objs = \
${LIB_DIR}/epsfr_main.o \
${LIB_DIR}/epsfr_episafari_utils.o \
${LIB_DIR}/epsfr_gsl_polyfit_utils.o \
${LIB_DIR}/epsfr_signal_track_tools.o \
${LIB_DIR}/epsfr_filter_utils.o \
${LIB_DIR}/epsfr_mapped_read_tools.o \
${LIB_DIR}/epsfr_min_max_utils.o \
${LIB_DIR}/epsfr_ansi_string.o \
${LIB_DIR}/epsfr_nucleotide.o \
${LIB_DIR}/epsfr_genomics_coords.o \
${LIB_DIR}/epsfr_nomenclature.o\
${LIB_DIR}/epsfr_annot_region_tools.o \
${LIB_DIR}/epsfr_rng.o \
${LIB_DIR}/epsfr_seed_manager.o \
${LIB_DIR}/epsfr_utils.o \
${LIB_DIR}/epsfr_xlog_math.o 

EpiSAFARI: ${objs}
	${CC} -O2 ${gsl_flags} -o ${exec_name} ${objs}

clean:
	rm -f ${objs} ${exec_name} 

