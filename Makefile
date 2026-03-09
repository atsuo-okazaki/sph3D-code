##############################################################
###                                                        ###
###           Make File for the GRAPE SPH Code             ###
###     originally written by M. Bate and R. Klessen       ###
###                      March 1996                        ###
###     and later modified for SR11000 on 12 Feb, 2006.    ###
###                                                        ###
###     This compiles the code using the BINARY TREE       ###
###                                                        ###
### (NOTE: As well as the different Makefile, `igrape'     ###
### must be changed when switching between GRAPE and TREE) ###
###                                                        ###
##############################################################

#.KEEP_STATE:

#FORTRAN = mpif90_r -64 -i,L -i,EU -loglist -cpp -DSR16K
#FORTRAN = mpif90_r -64 -loglist -cpp -DSR16K
#FORTRAN = mpiifort
FORTRAN = ifort
#FORTRAN = gfortran
#FORTRAN = nagfor

##############################################################
### Prepare for debugging                                  ###
##############################################################
PROGRAM = sph3d_smp_2025

### ifort
ifeq ($(FORTRAN), ifort)
   FC = ifort
   LD = $(FC)
   FFLAGS = -O3 -xCORE-AVX512 -fp-model consistent -qopenmp -g -traceback -fpp -r8 -mcmodel=medium -convert big_endian
   LDFLAGS = $(FFLAGS)
endif

### gfortran
ifeq ($(FORTRAN), gfortran)
   FC = gfortran
   LD = $(FC)
   FFLAGS =  -O3 -fconvert=big-endian -fopenmp
   LDFLAGS =  $(FFLAGS)
endif

### NAG fortran
ifeq ($(FORTRAN), nagfor)
   FC = nagfor
   LD = $(FC)
   FFLAGS =  -O4 -convert=BIG_ENDIAN -openmp -PIC
   LDFLAGS =  $(FFLAGS)
endif

SOURCES = mod_idims.f90 mod_mpi.f90 mod_accnum.f90 \
	mod_accurpt.f90 mod_actio.f90 mod_active.f90 mod_angm.f90 \
	mod_artvb.f90 mod_avail.f90 mod_binary.f90 \
	mod_bodys.f90 mod_btree.f90 mod_call.f90 mod_capt.f90 \
	mod_carac.f90 mod_cgas.f90 mod_cooldata.f90 mod_crpart.f90 \
	mod_curlist.f90 mod_current.f90 mod_debpt.f90 mod_debug.f90 \
	mod_debugit.f90 mod_debugpt.f90 mod_delay.f90 mod_densi.f90 \
	mod_diskbd.f90 mod_dissi.f90 mod_divve.f90 mod_dum.f90 \
	mod_ener.f90 mod_eosq.f90 mod_expan.f90 mod_files.f90 \
	mod_flag.f90 mod_force.f90 mod_fracg.f90 mod_g3monit.f90 \
	mod_ghost.f90 mod_glrho.f90 mod_gravi.f90 mod_gtime.f90 \
	mod_hagain.f90 mod_ian.f90 mod_indexx.f90 mod_infor.f90 \
	mod_init.f90 mod_integ.f90 mod_isnpt.f90 mod_kerne.f90 \
	mod_latti.f90 mod_logun.f90 mod_maslos.f90 \
	mod_maspres.f90 mod_misali.f90 mod_neighbor_P.f90 \
	mod_new.f90 mod_nextmpt.f90 mod_numpa.f90 mod_out.f90 \
	mod_outneigh.f90 mod_part.f90 mod_perform.f90 mod_phase.f90 \
	mod_constants.f90 mod_physeos.f90 mod_polyk2.f90 \
	mod_pres.f90 mod_ptbin.f90 mod_ptdump.f90 mod_ptmass.f90 \
	mod_ptsoft.f90 mod_rbnd.f90 mod_recor.f90 mod_rotat.f90 \
	mod_savernd.f90 mod_secret.f90 mod_setbin.f90 mod_soft.f90 \
	mod_split.f90 mod_sphcom.f90 mod_stepopt.f90 mod_stop.f90 \
	mod_sync.f90 mod_table.f90 mod_task.f90 mod_timei.f90 \
	mod_tlist.f90 mod_tming.f90 mod_torq.f90 mod_trans.f90 \
	mod_treecom_P.f90 mod_typef.f90 mod_units.f90 mod_useles.f90 \
	mod_varet.f90 mod_vargam.f90 mod_vbound.f90 mod_visc.f90 \
	mod_unidis.f90 mod_winds.f90 mod_xforce.f90 mod_xtorq.f90 \
	mod_zzhp.f90 ran1.f90 \
	sph.f90 accrete.f90 addump.f90 angmom.f90 boundry.f90 \
	cartdis.f90 cartmas.f90 cartpres.f90 cartvel.f90 chanref.f90 \
	chekopt.f90 condense.f90 constan.f90 coriol.f90 cyldis.f90 \
	cylmas.f90 cylpres.f90 cylvel.f90 densityi_P.f90 derivi_P.f90 \
	endrun.f90 energ.f90 eospg.f90 error.f90 evol.f90 \
	externf.f90 extract.f90 file.f90 forcei_P.f90 getneigh.f90 \
	gforsa_P.f90 gforsn_P.f90 gforspt_P.f90 ghostp1.f90 \
	ghostp2.f90 ghostp3.f90 hcalc.f90 hdot_P.f90 header.f90 \
	homexp.f90 indexx2.f90 inform.f90 \
	inopts.f90 insulate_TREE_P.f90 integs.f90 ktable.f90 \
	labrun.f90 mainop.f90 mesop.f90 \
	modif.f90 mtree_P.f90 newrun.f90 options.f90 phoenix.f90 \
	place.f90 preset.f90 prout.f90 psplit.f90 quit.f90 \
	rdump.f90 reduce.f90 revtree_P.f90 save.f90 \
        scaling.f90 secmes.f90 setpart.f90 \
	sphdis.f90 sphmas.f90 sphpres.f90 sphvel.f90 step_P.f90 \
	toten.f90 treef_P.f90 unifdis.f90 unit.f90 \
	radcool.f90 velweight.f90 wdump.f90 wrinsph.f90

LIBS =
PMLIBS =

OBJECTS = $(SOURCES:.f90=.o)

$(PROGRAM): $(OBJECTS)
	$(LD) $(FFLAGS) -o $@ $(OBJECTS) $(LIBS) $(PMLIBS)

.SUFFIXES: .o .f90
.f90.o : 
	 $(FC) $(FFLAGS) -c $<

clean:
	rm -f *.o *.log *.L *.mod
