#OBJDIR = ~/myobj/TDPPGEOOBJ
#PROGDIR = ~/mybins

OBJDIR = .
PROGDIR = .

PROG =	$(addprefix $(PROGDIR)/, tdppgeo)

SRCS =	algeb.f alli.f alpha.f area_elem.f ass_covinv.f ass_eval.f \
	ass_ext_drift.f ass_integ_val.f ass_real_val.f \
	assemb_crossterm_into_band_no_sym.f \
	assemb_crossterm_into_shuffled_band_no_sym.f \
	assemb_crossterm_into_shuffled_sparse.f \
	assemb_crossterm_into_sparse.f assemb_lhs.f assemb_rhs.f \
	assemb_rhs_coupled.f assemble.f \
	assemble_derivative_into_band_no_sym.f \
	assemble_derivative_into_shuffled_band_no_sym.f \
	assemble_derivative_into_shuffled_sparse.f \
	assemble_derivative_into_sparse.f assemble_full_into_band_no_sym.f \
	assemble_full_into_band_sym.f \
	assemble_full_into_shuffled_band_no_sym.f \
	assemble_full_into_shuffled_sparse.f assemble_full_into_sparse.f \
	assemble_into_band_no_sym.f assemble_into_band_sym.f \
	assemble_into_shuffled_band_no_sym.f assemble_into_shuffled_sparse.f \
	assemble_into_sparse.f assemble_lhs_coupled.f assemble_shuffled.f \
	assemble_sym_into_band_no_sym.f assemble_sym_into_band_sym.f \
	assemble_sym_into_shuffled_band_no_sym.f \
	assemble_sym_into_shuffled_sparse.f assemble_sym_into_sparse.f \
	assemble_sym_nd_into_shuffled_sparse.f \
	assemble_sym_no_diag_into_band_no_sym.f \
	assemble_sym_no_diag_into_band_sym.f \
	assemble_sym_no_diag_into_shuffled_band_no_sym.f \
	assemble_sym_no_diag_into_sparse.f assemble_vec_elem_into_band.f \
	assemble_vec_elem_into_shuffled_band.f \
	assemble_vector_into_shuffled_sparse.f \
	assemble_vector_nod_into_shuffled_band.f \
	assemble_vector_node_into_band.f assemble_vector_node_into_sparse.f \
	assign.f auxxx.f averages_var.f balance.f balance_fl.f balance_tr.f \
	balance_write_aux.f balance_write_fl.f \
	balance_write_tr.f basisfunc_obs.f buweight_obs.f \
	calc_dens.f calc_flowdens.f calc_mean_poldrift.f calc_visc.f \
	calc_watvol.f ch_geo_4n.f check_cons.f check_convergence.f \
	check_divergence.f check_obs.f check_par.f check_par_angle.f \
	check_stop.f choose_linearizt_method.f comat_bc.f comflow.f \
	comflow_conc.f comflow_presc_conc.f comflow_presc_head.f comp_aflu.f \
	comp_atra.f comp_atra_rec_fod.f comp_bflu.f comp_bflu_buoyancy.f \
	comp_btra.f comp_buoyancy_cf.f comp_cflu.f \
	comp_der_atra.f comp_der_btra.f comp_der_dtra.f comp_der_flow.f \
	comp_dflu.f comp_div_q.f comp_dtra.f comp_flowgrav.f comp_grad.f \
	comp_grad_loc.f comp_hess.f \
	comp_ind_change.f comp_obs.f comp_obs_aux.f comp_param_flow.f \
	comp_param_tra.f comp_phi_quadr.f comp_time_constants.f \
	comp_trans_mat.f comp_vel_buoyancy.f comp_vel_grav.f comp_vol_nod.f \
	compare_init.f compute_resid_max.f comvel.f \
	comvelelement.f cond_est_sim.f conditional_simulation.f contribdev.f \
	correct_ddtra_rhs.f correct_ddtra_rhs_i.f correct_ddtra_rhs_l.f \
	correct_parz.f coupled_flow_transport.f cova_krig_matrix.f \
	covpar_matrix.f der_atra_fp_vd.f der_atra_vd.f der_leak.f der_param.f \
	der_presc_flow.f der_presc_head.f der_vd.f der_vdtra.f \
	der_vdtra_buoyancy.f der_watvol_gen.f der_watvol_h.f der_wtv_stg.f \
	derarr.f derarr_tpt.f derclk.f dercoe.f dercoe_arr_flu.f \
	dercoe_nod_flu.f dercrd.f derdfm.f derdsl.f derdst.f derfod.f \
	derpor.f derpor_flu.f derq_arr.f derq_gen.f derq_nud.f \
	derq_por.f derq_stg.f derq_tra.f derq_tra_buoyancy.f derstg.f \
	dertra.f dertra_buoyancy.f dertra_tensor.f determinant.f difmat_dmt.f \
	difmat_mtdz.f difmat_mz.f endat.f endatinicond.f endatinicond_aux.f \
	endatobs.f entdat_dens.f entdat_funnoli_grav.f \
	entdat_groups_zones.f entdatix_coor.f entdatix_zon.f entdatlx_elem.f \
	entdatlx_zon.f entdatnz.f entflow.f entvel.f error.f estimcov.f \
	ext_atr_assign.f fe_iadm_d.f fileex.f flow.f fun_par.f funnoli.f \
	funnoli_in.f generstat.f geo_calc.f get_zone_group.f \
	getparname.f getvalue.f gradient.f grav_proj.f gslib_routines.f \
	hessiano.f histogram.f increment_iter.f init_extrap.f \
	init_geoest.f init_integ.f init_inv_prob_info.f init_var.f \
	input_mass.f interp_crit.f interp_eznum.f interp_nznum.f inver.f \
	jac_c.f jac_coupled.f jac_h.f jvj.f kriging.f ldimen.f lec_cfe.f \
	lec_cfn.f lec_ft.f lecdim.f leel.f likelihood.f \
	make_adjacency_arrays.f marquardt.f max_conect.f max_par.f \
	measure_obs.f min_process_info.f minimization.f mod_rapson_criteria.f \
	modelselec.f modif_jac_c_rdch.f modif_maxiter.f \
	names.f obj_var.f \
	open_one_file.f opt_groups_zones.f order.f param_value.f \
	position.f presc_bc.f presc_bc_coupled.f presc_leak_bc.f \
	principal.f prod_mat_vec.f prodat.f \
	provisional.f random_pipo.f read_extdrift_zones.f read_lin_comb.f \
	read_locations.f read_par.f read_search_drift.f read_tra.f \
	read_variograms.f residuals.f rhs_in_chp.f \
	rhs_in_chp_trlum.f rhs_node_flux_inv.f \
	rhs_noli_in.f rotmat_univkrig_covmax.f \
	sim_jac.f so_obj.f solution.f solve.f source_rdch.f \
	spat_weight_obs.f src_ncard.f stat_outp.f \
	state_variable_init.f store_var_scratch.f \
	temp_weight_obs.f tempcoeff.f tracehjvj.f transdens.f \
	transport.f unshuffle_lhs.f unshuffle_rhs.f update_files_scratch.f \
	update_parz.f update_sol_time.f update_state_variable.f update_time.f \
	update_time_incr.f update_vectors_aux.f verify_bw.f watsolv_new.f \
	watvol_ini.f watvol_update.f weightresid.f wri_deriv.f \
	wri_param_history.f wri_part.f wri_zone.f writ.f \
	writ_device.f write_arrayn.f write_debug_info.f write_matrix.f \
	write_nonlin_info.f write_stats.f write_temp_state_vars.f \
	write_var_gs.f write_velocity_vmsh.f \
	write_visualmeshplot.f zero_row.f zone_geometry.f comp_obs_geo.f \
	Write_VTK_Unstructured_Grid.f GetVTKCellType.f calc_gimli.f90 forpy_mod.F90 

COMMONS = MAIN_COMM.FOR COMMON.FOR COMMON_DMT.FOR COMMON_MTDZ.FOR \
          EQUIV_MTDZ.FOR

OBJS =	$(addprefix $(OBJDIR)/, algeb.o alli.o alpha.o area_elem.o \
	ass_covinv.o ass_eval.o \
	ass_ext_drift.o ass_integ_val.o ass_real_val.o \
	assemb_crossterm_into_band_no_sym.o \
	assemb_crossterm_into_shuffled_band_no_sym.o \
	assemb_crossterm_into_shuffled_sparse.o \
	assemb_crossterm_into_sparse.o assemb_lhs.o assemb_rhs.o \
	assemb_rhs_coupled.o assemble.o \
	assemble_derivative_into_band_no_sym.o \
	assemble_derivative_into_shuffled_band_no_sym.o \
	assemble_derivative_into_shuffled_sparse.o \
	assemble_derivative_into_sparse.o assemble_full_into_band_no_sym.o \
	assemble_full_into_band_sym.o \
	assemble_full_into_shuffled_band_no_sym.o \
	assemble_full_into_shuffled_sparse.o assemble_full_into_sparse.o \
	assemble_into_band_no_sym.o assemble_into_band_sym.o \
	assemble_into_shuffled_band_no_sym.o assemble_into_shuffled_sparse.o \
	assemble_into_sparse.o assemble_lhs_coupled.o assemble_shuffled.o \
	assemble_sym_into_band_no_sym.o assemble_sym_into_band_sym.o \
	assemble_sym_into_shuffled_band_no_sym.o \
	assemble_sym_into_shuffled_sparse.o assemble_sym_into_sparse.o \
	assemble_sym_nd_into_shuffled_sparse.o \
	assemble_sym_no_diag_into_band_no_sym.o \
	assemble_sym_no_diag_into_band_sym.o \
	assemble_sym_no_diag_into_shuffled_band_no_sym.o \
	assemble_sym_no_diag_into_sparse.o assemble_vec_elem_into_band.o \
	assemble_vec_elem_into_shuffled_band.o \
	assemble_vector_into_shuffled_sparse.o \
	assemble_vector_nod_into_shuffled_band.o \
	assemble_vector_node_into_band.o assemble_vector_node_into_sparse.o \
	assign.o auxxx.o averages_var.o balance.o balance_fl.o balance_tr.o \
	balance_write_aux.o balance_write_fl.o \
	balance_write_tr.o basisfunc_obs.o buweight_obs.o \
	calc_dens.o calc_flowdens.o calc_mean_poldrift.o calc_visc.o \
	calc_watvol.o ch_geo_4n.o check_cons.o check_convergence.o \
	check_divergence.o check_obs.o check_par.o check_par_angle.o \
	check_stop.o choose_linearizt_method.o comat_bc.o comflow.o \
	comflow_conc.o comflow_presc_conc.o comflow_presc_head.o comp_aflu.o \
	comp_atra.o comp_atra_rec_fod.o comp_bflu.o comp_bflu_buoyancy.o \
	comp_btra.o comp_buoyancy_cf.o comp_cflu.o \
	comp_der_atra.o comp_der_btra.o comp_der_dtra.o comp_der_flow.o \
	comp_dflu.o comp_div_q.o comp_dtra.o comp_flowgrav.o comp_grad.o \
	comp_grad_loc.o comp_hess.o \
	comp_ind_change.o comp_obs.o comp_obs_aux.o comp_param_flow.o \
	comp_param_tra.o comp_phi_quadr.o comp_time_constants.o \
	comp_trans_mat.o comp_vel_buoyancy.o comp_vel_grav.o comp_vol_nod.o \
	compare_init.o compute_resid_max.o comvel.o \
	comvelelement.o cond_est_sim.o conditional_simulation.o contribdev.o \
	correct_ddtra_rhs.o correct_ddtra_rhs_i.o correct_ddtra_rhs_l.o \
	correct_parz.o coupled_flow_transport.o cova_krig_matrix.o \
	covpar_matrix.o der_atra_fp_vd.o der_atra_vd.o der_leak.o der_param.o \
	der_presc_flow.o der_presc_head.o der_vd.o der_vdtra.o \
	der_vdtra_buoyancy.o der_watvol_gen.o der_watvol_h.o der_wtv_stg.o \
	derarr.o derarr_tpt.o derclk.o dercoe.o dercoe_arr_flu.o \
	dercoe_nod_flu.o dercrd.o derdfm.o derdsl.o derdst.o derfod.o \
	derpor.o derpor_flu.o derq_arr.o derq_gen.o derq_nud.o \
	derq_por.o derq_stg.o derq_tra.o derq_tra_buoyancy.o derstg.o \
	dertra.o dertra_buoyancy.o dertra_tensor.o determinant.o difmat_dmt.o \
	difmat_mtdz.o difmat_mz.o endat.o endatinicond.o endatinicond_aux.o \
	endatobs.o entdat_dens.o entdat_funnoli_grav.o \
	entdat_groups_zones.o entdatix_coor.o entdatix_zon.o entdatlx_elem.o \
	entdatlx_zon.o entdatnz.o entflow.o entvel.o error.o estimcov.o \
	ext_atr_assign.o fe_iadm_d.o fileex.o flow.o fun_par.o funnoli.o \
	funnoli_in.o generstat.o geo_calc.o get_zone_group.o \
	getparname.o getvalue.o gradient.o grav_proj.o gslib_routines.o \
	hessiano.o histogram.o increment_iter.o init_extrap.o \
	init_geoest.o init_integ.o init_inv_prob_info.o init_var.o \
	input_mass.o interp_crit.o interp_eznum.o interp_nznum.o inver.o \
	jac_c.o jac_coupled.o jac_h.o jvj.o kriging.o ldimen.o lec_cfe.o \
	lec_cfn.o lec_ft.o lecdim.o leel.o likelihood.o \
	make_adjacency_arrays.o marquardt.o max_conect.o max_par.o \
	measure_obs.o min_process_info.o minimization.o mod_rapson_criteria.o \
	modelselec.o modif_jac_c_rdch.o modif_maxiter.o \
	names.o obj_var.o \
	open_one_file.o opt_groups_zones.o order.o param_value.o \
	position.o presc_bc.o presc_bc_coupled.o presc_leak_bc.o \
	principal.o prod_mat_vec.o prodat.o \
	provisional.o random_pipo.o read_extdrift_zones.o read_lin_comb.o \
	read_locations.o read_par.o read_search_drift.o read_tra.o \
	read_variograms.o residuals.o rhs_in_chp.o \
	rhs_in_chp_trlum.o rhs_node_flux_inv.o \
	rhs_noli_in.o rotmat_univkrig_covmax.o \
	sim_jac.o so_obj.o solution.o solve.o source_rdch.o \
	spat_weight_obs.o src_ncard.o stat_outp.o \
	state_variable_init.o store_var_scratch.o \
	temp_weight_obs.o tempcoeff.o tracehjvj.o transdens.o \
	transport.o unshuffle_lhs.o unshuffle_rhs.o update_files_scratch.o \
	update_parz.o update_sol_time.o update_state_variable.o update_time.o \
	update_time_incr.o update_vectors_aux.o verify_bw.o watsolv_new.o \
	watvol_ini.o watvol_update.o weightresid.o wri_deriv.o \
	wri_param_history.o wri_part.o wri_zone.o writ.o \
	writ_device.o write_arrayn.o write_debug_info.o write_matrix.o \
	write_nonlin_info.o write_stats.o write_temp_state_vars.o \
	write_var_gs.o write_velocity_vmsh.o \
	write_visualmeshplot.o zero_row.o zone_geometry.o comp_obs_geo.o \
	Write_VTK_Unstructured_Grid.o GetVTKCellType.o calc_gimli.o forpy_mod.o)

MODS =  forpy_mod.mod

LIBS =	

F90 = gfortran
#F90 = x86_64-conda_cos6-linux-gnu-gfortran 
#Debug
#F90FLAGS = -ggdb -fbounds-check -Wall -Wtabs -g3
#F90FLAGS = -ggdb -Wall -Wtabs -g3 
#Release
F90FLAGS = -O3 -finit-local-zero 
LDFLAGS = 
PYFLAGS =  -fno-lto `python3-config --ldflags` -lpython3.6
PYFLAGS =  `python3-config --ldflags` 
.SUFFIXES : .o .f .f90 .FOR .F90 .mod

all: $(PROG)

$(PROG): $(OBJS) $(MODS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS) $(PYFLAGS)

clean:
	rm -f $(PROG) $(OBJS)  $(MODS)

cleanobj:
	rm -f $(OBJS)

# Regla para los commons (no hace nada)
COMMON.FOR:
COMMON_DMT.FOR:
COMMON_MTDZ.FOR:
EQUIV_MTDZ.FOR:
MAIN_COMM.FOR:


# Regla implícita para pasar .f a .o (forma nueva "patern rule").
$(OBJDIR)/%.o : %.f
	$(F90) $(F90FLAGS) -c $< -o $@

$(OBJDIR)/%.o : %.f90
	$(F90) $(F90FLAGS) -c $< -o $@

$(OBJDIR)/calc_gimli.o: forpy_mod.mod

forpy_mod.mod: forpy_mod.F90
	$(F90)  -c  forpy_mod.F90

$(OBJDIR)/algeb.o: COMMON.FOR
$(OBJDIR)/alpha.o: tracehjvj.f
$(OBJDIR)/basisfunc_obs.o: COMMON.FOR
$(OBJDIR)/buweight_obs.o: leel.f
$(OBJDIR)/contribdev.o: position.f
$(OBJDIR)/covpar_matrix.o: leel.f
$(OBJDIR)/difmat_dmt.o: MAIN_COMM.FOR COMMON_DMT.FOR
$(OBJDIR)/difmat_mtdz.o: COMMON_MTDZ.FOR EQUIV_MTDZ.FOR
$(OBJDIR)/difmat_mz.o: COMMON_MTDZ.FOR
$(OBJDIR)/endatinicond_aux.o: leel.f
$(OBJDIR)/endatobs.o: leel.f
$(OBJDIR)/entdat.o: ldimen.o
$(OBJDIR)/entdat_dens.o: leel.f
$(OBJDIR)/entdat_funnoli_grav.o: leel.f
$(OBJDIR)/entdat_groups_zones.o: leel.f
$(OBJDIR)/entdatix_coor.o: leel.f
$(OBJDIR)/entdatix_zon.o: leel.f
$(OBJDIR)/entdatlx_elem.o: ldimen.f leel.f
$(OBJDIR)/entdatlx_zon.o: ldimen.f leel.f
$(OBJDIR)/entvel.o: leel.f
$(OBJDIR)/getvalue.o: position.f
$(OBJDIR)/init_extrap.o: tempcoeff.f
$(OBJDIR)/init_inv_prob_info.o: tempcoeff.f
$(OBJDIR)/interp_eznum.o: ldimen.o
$(OBJDIR)/jvj.o: getvalue.f
$(OBJDIR)/lec_cfe.o: leel.f
$(OBJDIR)/lec_cfn.o: leel.f
$(OBJDIR)/lec_ft.o: leel.f
$(OBJDIR)/lecdim.o: fileex.f leel.f
$(OBJDIR)/modif_maxiter.o: fileex.f
$(OBJDIR)/open_one_file.o: fileex.f
$(OBJDIR)/opt_groups_zones.o: leel.f
$(OBJDIR)/provisional.o: COMMON.FOR
$(OBJDIR)/read_extdrift_zones.o: leel.f
$(OBJDIR)/read_lin_comb.o: leel.f
$(OBJDIR)/read_locations.o: leel.f
$(OBJDIR)/read_par.o: leel.f
$(OBJDIR)/read_search_drift.o: leel.f
$(OBJDIR)/read_tra.o: leel.f
$(OBJDIR)/read_variograms.o: leel.f
$(OBJDIR)/residuals.o: weightresid.f histogram.f
$(OBJDIR)/sim_jac.o: MAIN_COMM.FOR COMMON_DMT.FOR COMMON_MTDZ.FOR
$(OBJDIR)/so_obj.o: contribdev.f histogram.f
$(OBJDIR)/stat_outp.o: histogram.f
$(OBJDIR)/tracehjvj.o: getvalue.f
$(OBJDIR)/transdens.o: COMMON.FOR MAIN_COMM.FOR
$(OBJDIR)/weightresid.o: getvalue.f
$(OBJDIR)/write_stats.o: getparname.f getvalue.f histogram.f
