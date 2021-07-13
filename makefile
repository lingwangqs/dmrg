.SUFFIXES: .f90 .f .c .cpp .o

FC=mpiifort
F90=mpiifort
CC=mpiicc
CXX=mpiicpc
FFLAGS=-openmp -g -O3
F90FLAGS=-openmp -g -O3
CFLAGS=-openmp -g -O3
CCFLAGS=-openmp -g -O3
CXXFLAGS=-openmp -g -O3
lib=~/program/SU2_real_mpi


main_mpi: main_mpi.o tensor.o tensor_su2.o su2bond.o su2struct.o ran.o sort2.o dmrg_su2.o dmrg_su2_mpi.o lanczos_su2.o permutation_operator.o
	${F90} ${FFLAGS} -o $@ $< tensor.o tensor_su2.o su2bond.o su2struct.o ran.o sort2.o dmrg_su2.o dmrg_su2_mpi.o lanczos_su2.o permutation_operator.o -L/home/intel/impi/5.0.1.035/intel64/lib/  -lmpicxx -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lstdc++ -nofor_main 

main_mpi_aklt: main_mpi.o tensor.o tensor_su2.o su2bond.o su2struct.o ran.o sort2.o dmrg_su2_aklt.o dmrg_su2_mpi.o lanczos_su2.o
	${F90} ${FFLAGS} -o main_mpi_aklt main_mpi.o tensor.o tensor_su2.o su2bond.o su2struct.o ran.o sort2.o dmrg_su2_aklt.o dmrg_su2_mpi.o lanczos_su2.o -L/home/intel/impi/5.0.1.035/intel64/lib/  -lmpicxx -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lstdc++ -nofor_main

main_mpi_mea_stgm: main_mpi_mea_stgm.o tensor.o tensor_su2.o su2bond.o su2struct.o ran.o sort2.o dmrg_su2_mea_stgm.o dmrg_su2_mpi_mea_stgm.o lanczos_su2.o
	${F90} ${FFLAGS} -o $@ $< tensor.o tensor_su2.o su2bond.o su2struct.o ran.o sort2.o dmrg_su2_mea_stgm.o dmrg_su2_mpi_mea_stgm.o lanczos_su2.o -L/home/intel/impi/5.0.1.035/intel64/lib/  -lmpicxx -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lstdc++ -nofor_main 

main_mpi_mea_enr: main_mpi_mea_enr.o tensor.o tensor_su2.o su2bond.o su2struct.o ran.o sort2.o dmrg_su2_mea_enr.o dmrg_su2_mpi_mea_enr.o lanczos_su2.o
	${F90} ${FFLAGS} -o $@ $< tensor.o tensor_su2.o su2bond.o su2struct.o ran.o sort2.o dmrg_su2_mea_enr.o dmrg_su2_mpi_mea_enr.o lanczos_su2.o -L/home/intel/impi/5.0.1.035/intel64/lib/  -lmpicxx -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lstdc++ -nofor_main 

main_mpi_mea_permutation: main_mpi_mea_permutation.o tensor.o tensor_su2.o su2bond.o su2struct.o ran.o sort2.o dmrg_su2.o dmrg_su2_mpi_mea_permutation.o lanczos_su2.o permutation_operator.o
	${F90} ${FFLAGS} -o $@ $< tensor.o tensor_su2.o su2bond.o su2struct.o ran.o sort2.o dmrg_su2.o dmrg_su2_mpi_mea_permutation.o lanczos_su2.o permutation_operator.o -L/home/intel/impi/5.0.1.035/intel64/lib/  -lmpicxx -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lstdc++ -nofor_main 

main_mpi_mea_plqplq: main_mpi_mea_plqplq.o tensor.o tensor_su2.o su2bond.o su2struct.o ran.o sort2.o dmrg_su2_mea_plqplq.o dmrg_su2_mpi_mea_plqplq.o lanczos_su2.o permutation_operator.o
	${F90} ${FFLAGS} -o $@ $< tensor.o tensor_su2.o su2bond.o su2struct.o ran.o sort2.o dmrg_su2_mea_plqplq.o dmrg_su2_mpi_mea_plqplq.o lanczos_su2.o permutation_operator.o -L/home/intel/impi/5.0.1.035/intel64/lib/  -lmpicxx -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lstdc++ -nofor_main 

#main_mpi_mea_plqplq: main_mpi_mea_plqplq.o tensor.o tensor_su2.o su2bond.o su2struct.o ran.o sort2.o dmrg_su2_mea_plqplq.o dmrg_su2_mpi_mea_plqplq.o lanczos_su2.o permutation_operator.o
#	${F90} ${FFLAGS} -o $@ $< tensor.o tensor_su2.o su2bond.o su2struct.o ran.o sort2.o dmrg_su2_mea_plqplq.o dmrg_su2_mpi_mea_plqplq.o lanczos_su2.o permutation_operator.o -L/home/intel/impi/5.0.1.035/intel64/lib/  -lmpicxx -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lstdc++ -nofor_main 

main_mpi_mea_bindercumulant: main_mpi_mea_bindercumulant.o tensor.o tensor_su2.o su2bond.o su2struct.o ran.o sort2.o dmrg_su2_mea_bindercumulant.o dmrg_su2_mpi_mea_bindercumulant.o lanczos_su2.o permutation_operator.o
	${F90} ${FFLAGS} -o $@ $< tensor.o tensor_su2.o su2bond.o su2struct.o ran.o sort2.o dmrg_su2_mea_bindercumulant.o dmrg_su2_mpi_mea_bindercumulant.o lanczos_su2.o permutation_operator.o -L/home/intel/impi/5.0.1.035/intel64/lib/  -lmpicxx -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lstdc++ -nofor_main 

main_nonmpi_mea_bindercumulant: main_nonmpi_mea_bindercumulant.o tensor.o tensor_su2.o su2bond.o su2struct.o ran.o sort2.o dmrg_su2_mea_bindercumulant.o dmrg_su2_nonmpi_mea_bindercumulant.o lanczos_su2.o permutation_operator.o
	${F90} ${FFLAGS} -o $@ $< tensor.o tensor_su2.o su2bond.o su2struct.o ran.o sort2.o dmrg_su2_mea_bindercumulant.o dmrg_su2_nonmpi_mea_bindercumulant.o lanczos_su2.o permutation_operator.o -L/home/intel/impi/5.0.1.035/intel64/lib/  -lmpicxx -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lstdc++ -nofor_main 

main_mpi_mea_translation: main_mpi_mea_plqplq.o tensor.o tensor_su2.o su2bond.o su2struct.o ran.o sort2.o dmrg_su2_mea_plqplq.o dmrg_su2_mpi_mea_translation.o lanczos_su2.o permutation_operator.o
	${F90} ${FFLAGS} -o $@ main_mpi_mea_plqplq.o tensor.o tensor_su2.o su2bond.o su2struct.o ran.o sort2.o dmrg_su2_mea_plqplq.o dmrg_su2_mpi_mea_translation.o lanczos_su2.o permutation_operator.o -L/home/intel/impi/5.0.1.035/intel64/lib/  -lmpicxx -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lstdc++ -nofor_main 

main_mpi_mea_bond_enr: main_mpi_mea_bond_enr.o tensor.o tensor_su2.o su2bond.o su2struct.o ran.o sort2.o dmrg_su2_mea_bond_enr.o dmrg_su2_mpi_mea_bond_enr.o lanczos_su2.o
	${F90} ${FFLAGS} -o $@ $< tensor.o tensor_su2.o su2bond.o su2struct.o ran.o sort2.o dmrg_su2_mea_bond_enr.o dmrg_su2_mpi_mea_bond_enr.o lanczos_su2.o -L/home/intel/impi/5.0.1.035/intel64/lib/  -lmpicxx -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lstdc++ -nofor_main 

main_mpi_mea_glide: main_mpi_mea_glide.o tensor.o tensor_su2.o su2bond.o su2struct.o ran.o sort2.o dmrg_su2_mea_glide.o dmrg_su2_mpi_mea_glide.o lanczos_su2.o permutation_operator.o
	${F90} ${FFLAGS} -o $@ main_mpi_mea_glide.o tensor.o tensor_su2.o su2bond.o su2struct.o ran.o sort2.o dmrg_su2_mea_glide.o dmrg_su2_mpi_mea_glide.o lanczos_su2.o permutation_operator.o -L/home/intel/impi/5.0.1.035/intel64/lib/  -lmpicxx -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lstdc++ -nofor_main 

main_mpi_mea_reflectionx: main_mpi_mea_glide.o tensor.o tensor_su2.o su2bond.o su2struct.o ran.o sort2.o dmrg_su2_mea_glide.o dmrg_su2_mpi_mea_reflectionx.o lanczos_su2.o permutation_operator.o
	${F90} ${FFLAGS} -o $@ main_mpi_mea_glide.o tensor.o tensor_su2.o su2bond.o su2struct.o ran.o sort2.o dmrg_su2_mea_glide.o dmrg_su2_mpi_mea_reflectionx.o lanczos_su2.o permutation_operator.o -L/home/intel/impi/5.0.1.035/intel64/lib/  -lmpicxx -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lstdc++ -nofor_main 

main_mpi_mea_reflectiony: main_mpi_mea_reflectiony.o tensor.o tensor_su2.o su2bond.o su2struct.o ran.o sort2.o dmrg_su2_mea_reflectiony.o dmrg_su2_mpi_mea_reflectiony.o lanczos_su2.o permutation_operator.o
	${F90} ${FFLAGS} -o $@ main_mpi_mea_reflectiony.o tensor.o tensor_su2.o su2bond.o su2struct.o ran.o sort2.o dmrg_su2_mea_reflectiony.o dmrg_su2_mpi_mea_reflectiony.o lanczos_su2.o permutation_operator.o -L/home/intel/impi/5.0.1.035/intel64/lib/  -lmpicxx -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lstdc++ -nofor_main 


clean:
	/bin/rm -f *.o 

.f.o:
	${FC} ${F77FLAGS} -c $<  ${FCINCLUDE}
.f90.o:
	${F90} ${F90FLAGS} -c $< ${FCINCLUDE}
.c.o:
	${CC} ${CCFLAGS} -c $< ${CXXINCLUDE}
.cpp.o:
	${CXX} ${CXXFLAGS} -c $< ${CXXINCLUDE}
