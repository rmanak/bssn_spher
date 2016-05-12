.IGNORE:

AR = xiar
SHELL = /bin/sh
LIBS = -lm -lbbhutil -lvutil -lifcore 
LDFLAGS = -O3 -L/usr/loca/lib -L.
CC = icc
CFLAGS = -O3
CPPFLAGS = -I. -I/usr/local/include
CC_COMP = $(CC) -c $(CFLAGS) $(CPPFLAGS)
CC_LOAD = $(CC) $(LDFLAGS)
F77 = ifort
F77FLAGS = -fast -w -w90 -w95 -cm -Vaxlib
F77_COMP = $(F77) -c $(F77FLAGS)
F77_LOAD = $(F77) $(F77FLAGS) $(F77_LDFLAGS) $(LDFLAGS)


.f.o:
	$(F77_COMP) $*.f
.c.o:
	$(CC_COMP) $*.c

OBJECTS = ire_dph.o  ire_dtk.o  ire_dgx.o  init_a.o  init_b.o  init_tk.o  ire_daatheta.o  ire_daax.o  ire_dlamx.o  ire_eqipib2.o  ire_eqipsi.o  ire_eqrpib2.o  ire_eqrpsi.o  ire_dgtheta.o  compute_c1.o  compute_c5.o  compute_ddlx.o  compute_dipsi.o  compute_drpsi.o  compute_emph.o  compute_ptheta.o  compute_px.o  compute_rrtheta.o  compute_rrx.o  compute_utheta.o  compute_ux.o  compute_a1.o  compute_a2.o  compute_b1.o  compute_ipi.o  compute_jsthth.o  compute_jsxx.o  compute_metx.o  compute_phi.o  compute_phi2.o  compute_pi2.o  compute_psi.o  compute_rho.o  compute_rpi.o  compute_sx.o  compute_ts.o  compute_u.o  evolve_daatheta.o  evolve_daax.o  evolve_dlamx.o  evolve_dph.o  evolve_dtk.o  evolve_eqipib2.o  evolve_eqipsi.o  evolve_eqrpib2.o  evolve_eqrpsi.o  evolve_dgtheta.o  evolve_dgx.o  init_athth.o  init_axx.o  init_lamx.o  init_iphi.o  init_ipi.o  init_ipib2.o  init_ipsi.o  init_ipsidot.o  init_rphi.o  init_rpi.o  init_rpib2.o  init_rpsi.o  init_rpsidot.o  resid_daatheta.o  resid_daax.o  resid_dlamx.o  resid_dph.o  resid_dtk.o  resid_eqipib2.o  resid_eqipsi.o  resid_eqrpib2.o  resid_eqrpsi.o  resid_hs.o  resid_dgtheta.o  resid_dgx.o  solve_hs.o  compute_ddltheta.o  compute_mettheta.o init_alpha.o init_beta.o init_ctfm.o init_ctfmp.o init_octfmp.o ire_restx.o ire_resxx.o ire_restt.o ire_rpsi_direct.o ire_val_lamx.o ire_hs.o ire_resthth.o resid_hs3.o solve_hs3.o compute_c5psi.o evol_alphaopl.o resid_alphaopl.o evol_alphakd.o resid_alphakd.o compute_th.o evol_alphadyopl.o resid_alphadyopl.o evol_betagd.o resid_betagd.o  resid_betagd2.o evol_betagd2.o init_ctfm2.o init_ctfmp2.o init_octfmp2.o ire_mk.o irev_dph.o irev_dtk.o irev_dgx.o irev_daatheta.o irev_daax.o irev_dlamx.o irev_dgtheta.o compute_db.o evol_betakd.o resid_betakd.o eval_dlam.o evol_beta_hgd_adv.o evol_BB_hgd_adv.o evol_alpha_hgd_adv.o resid_BB_hgd_adv.o resid_alpha_hgd_adv.o resid_beta_hgd_adv.o 


all: libbssn.a install

libbssn.a: $(OBJECTS)
	$(AR) r libbssn.a $(OBJECTS)

install: libbssn.a
	/bin/cp libbssn.a ../lib/

clean:
	/bin/rm -rf $(OBJECTS) libbssn.a
