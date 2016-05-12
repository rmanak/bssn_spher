#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <bbhutil.h>
#include <sdf_read.h>
#include <utils.h>

#include <compute_tt.h>
#include <eval_dlam.h>
#include <compute_db.h>
#include <compute_dmdr.h>
#include <rel_diff.h>
#include <ire_mk.h>
#include <init_ctfm2.h>
#include <init_octfmp2.h>
#include <init_ctfmp2.h>

#include <evol_betakd.h>
#include <resid_betakd.h>
#include <evol_betagd2.h>
#include <resid_betagd2.h>
#include <compute_th.h>
#include <evol_alphakd.h>
#include <resid_alphakd.h>
#include <evol_alphaopl.h>
#include <resid_alphaopl.h>
#include <resid_hs3.h>
#include <solve_hs3.h>
#include <compute_c5psi.h>
#include <ire_restx.h>
#include <ire_resxx.h>
#include <ire_restt.h>
#include <ire_resthth.h>
#include <ire_val_lamx.h>
#include <ire_hs.h>
#include <ire_rpsi_direct.h>
#include <compute_b.h>
#include <compute_athth.h>
#include <init_ctfm.h>
#include <init_octfmp.h>
#include <init_ctfmp.h>
#include <init_a.h>
#include <init_b.h>
#include <init_tk.h>
#include <ire_dph.h>
#include <ire_dtk.h>
#include <ire_dgx.h>
#include <ire_daatheta.h>
#include <ire_daax.h>
#include <ire_dlamx.h>
#include <ire_eqipib2.h>
#include <ire_eqrpib2.h>
#include <ire_dgtheta.h>

#include <irev_dph.h>
#include <irev_dtk.h>
#include <irev_dgx.h>
#include <irev_daatheta.h>
#include <irev_daax.h>
#include <irev_dlamx.h>
#include <irev_dgtheta.h>

#include <compute_c1.h>
#include <compute_c5.h>
#include <compute_ddlx.h>
#include <compute_dipsi.h>
#include <compute_drpsi.h>
#include <compute_emph.h>
#include <compute_ptheta.h>
#include <compute_px.h>
#include <compute_rrtheta.h>
#include <compute_rrx.h>
#include <compute_utheta.h>
#include <compute_ux.h>
#include <compute_a1.h>
#include <compute_a2.h>
#include <compute_b1.h>
#include <compute_ipi.h>
#include <compute_jsthth.h>
#include <compute_jsxx.h>
#include <compute_metx.h>
#include <compute_phi.h>
#include <compute_phi2.h>
#include <compute_pi2.h>
#include <compute_psi.h>
#include <compute_rho.h>
#include <compute_rpi.h>
#include <compute_sx.h>
#include <compute_ts.h>
#include <compute_u.h>
#include <evolve_daatheta.h>
#include <evolve_daax.h>
#include <evolve_dlamx.h>
#include <evolve_dph.h>
#include <evolve_dtk.h>
#include <evolve_eqipib2.h>
#include <evolve_eqrpib2.h>
#include <evolve_dgtheta.h>
#include <evolve_dgx.h>
#include <init_athth.h>
#include <init_axx.h>
#include <init_lamx.h>
#include <init_iphi.h>
#include <init_ipib2.h>
#include <init_ipsi.h>
#include <init_rphi.h>
#include <init_rpib2.h>
#include <init_rpsi.h>
#include <resid_daatheta.h>
#include <resid_daax.h>
#include <resid_dlamx.h>
#include <resid_dph.h>
#include <resid_dtk.h>
#include <resid_eqipib2.h>
#include <resid_eqrpib2.h>
#include <resid_hs.h>
#include <resid_dgtheta.h>
#include <resid_dgx.h>
#include <solve_hs.h>
#include <compute_ddltheta.h>
#include <compute_mettheta.h>
#include <interpolator.h>
#include <compute_mass.h>

typedef int bool;
#define true 1
#define false 0

struct timeval start_time, stop_time, elapsed_time, estart_time,estop_time,eelapsed_time,
               rstart_time,rstop_time,relapsed_time,wstart_time,wstop_time,welapsed_time;

int check_bh = 1;
int check_dis = 1;
double gma,gmc,gmd;
double ck, epsal;
int slicing_c;
double mus, eta;

/* Shapes: */
int intplt_rc;


int Nx;
int shape[1];
int dim;
int level;
/* Coordinates: */
double *x;
/* Grid Functions: */

double *ctfm;
double *ctfmp;
double *octfmp;

double *n_DLamx;
double *np1_DLamx;

double *n_BB;
double *np1_BB;

double *n_divbeta;
double *np1_divbeta;

double *n_ire_daath;
double *n_ire_daax;
double *n_ire_dlamx;
double *n_ire_dph;
double *n_ire_dtk;
double *n_ire_dgth;
double *n_ire_dgx;

double *n_mk; 
double *n_mass;
double *n_mass2;
double *n_tmr;
double *nph_th;

double *n_restx;
double *n_restt;
double *n_resxx;
double *n_resthth;
double *n_reshs;

double *n_U;
double *np1_U;

double *n_psi;
double *no_psi;
double *np1_psi;

double *n_PHI2;
double *n_PI2;
double *n_iPIb2;
double *n_rPIb2;

double *np1_PHI2;
double *np1_PI2;
double *np1_iPIb2;
double *np1_rPIb2;

double *n_C1s;
double *n_C5psi;
double *n_C5s;
double *n_rpsi;
double *n_ipsi;
double *n_rPI;
double *n_iPI;
double *n_rPHI;
double *n_iPHI;

double *np1_rpsi;
double *np1_ipsi;
double *np1_rPI;
double *np1_iPI;
double *np1_rPHI;
double *np1_iPHI;

double *n_a2;
double *n_a1;
double *n_b1;

double *np1_a2;
double *np1_a1;
double *np1_b1;

double *nm1_a1;
double *nm1_b1;

double *n_A;
double *n_Athth;
double *n_Axx;
double *n_B;
double *n_DDLxx;
double *n_JSxx;
double *n_K;
double *n_Lamx;
double *n_Pxx;
double *n_RRxx;
double *n_Sx;
double *n_TS;
double *n_Uxx;
double *n_beta;
double *n_metxx;
double *n_DDLthetatheta;
double *n_JSthth;
double *n_Pthetatheta;
double *n_RRthetatheta;
double *n_Uthetatheta;
double *n_em4phi;
double *n_metthetatheta;
double *n_phi;
double *n_rho;
double *n_alpha;

double *nph_A;
double *nph_B;
double *nph_Axx;
double *nph_Athth;
double *nph_phi;
double *nph_Lamx;
double *nph_K;


double *nph_ipsi;
double *nph_rpsi;
double *nph_iPIb2;
double *nph_rPIb2;


double *np1_A;
double *np1_Athth;
double *np1_Axx;
double *np1_B;
double *np1_DDLxx;
double *np1_JSxx;
double *np1_K;
double *np1_Lamx;
double *np1_Pxx;
double *np1_RRxx;
double *np1_Sx;
double *np1_TS;
double *np1_Uxx;
double *np1_beta;
double *np1_metxx;
double *np1_DDLthetatheta;
double *np1_JSthth;
double *np1_Pthetatheta;
double *np1_RRthetatheta;
double *np1_Uthetatheta;
double *np1_em4phi;
double *np1_metthetatheta;
double *np1_phi;
double *np1_rho;
double *np1_alpha;


double *nm1_ipsi;
double *nm1_rpsi;
double *nm1_iPIb2;
double *nm1_rPIb2;

double *nm1_A;
double *nm1_Athth;
double *nm1_Axx;
double *nm1_B;
double *nm1_DDLxx;
double *nm1_JSxx;
double *nm1_K;
double *nm1_Lamx;
double *nm1_Pxx;
double *nm1_RRxx;
double *nm1_Sx;
double *nm1_TS;
double *nm1_Uxx;
double *nm1_beta;
double *nm1_metxx;
double *nm1_DDLthetatheta;
double *nm1_JSthth;
double *nm1_Pthetatheta;
double *nm1_RRthetatheta;
double *nm1_Uthetatheta;
double *nm1_em4phi;
double *nm1_metthetatheta;
double *nm1_phi;
double *nm1_rho;
double *nm1_alpha;

double *n_ipsidot;
double *n_rpsidot;



double maxtmr;

double res_psi;

/* Parameters: */
double init_tol, evol_tol;
double lv, rv;
double omega;
double mass;
double ht;
double myzero;
double vee;
double deltx, x0, amp;
/* Coordinate Parameters: */
double zepsdis;
double x_max;
double x_min;
double hx;
double bbox[2];
int phys_bdy[2];
/* Time Evolution Parameters: */
int steps;
int output_freq;
int output_freq_phi;
double lambda;
double time;
double tres;

double mabs(double val) {
    if (val > 0.0) { return val; }
    else { return -val;}
}
void set_zero(double *f, int N) {
    int i;
    for(i=0;i<N;i++) { f[i] = 0.0; }
}

void find_min(double *f, int N,double *min_th) {
    int i;
    *min_th = f[1];
    for(i=2;i<N-1;i++) { 
        if ( x[i] < x_min + (x_max-x_min)*0.9  ) {
            if ( f[i] < *min_th ) { *min_th = f[i]; }
        }
    }
}


void set_one(double *f, int N) {
    int i;
    for(i=0;i<N;i++) { f[i] = 1.0; }
}


void swap_levels(double **a, double **b) {
    double *t;
    t = *a;
    *a = *b;
    *b = t;
}



void read_params(char *p_file,double *deltx, double *x0, double *amp, double *mass, int *level,double *lambda,int *output_freq, int *output_freq_phi, double *myzero,double *x_max,int *Nx,double *x_min,int *steps,int *read_lo_lev, int *compute_ire, double *omega, double *init_tol, double *evol_tol, int *flat_spacetime, int *freeze_matter, double *lv, double *rv, int *nc,double *zepsdis, int *slicing_c, int *check_bh, int *check_dis, double* bh_threshold, double* dis_threshold, double *gma, double *gmb, double *ck, double *epsal, double *trap_threshold, double *mus, double *eta, double *gmc, double *gmd, double *vee, int *min_evol_iter, double *advc)
{
    *advc = 0.0;
    *min_evol_iter = 10;
    *vee = 1.0;
    *mus = 0.0;
    *eta = 0.0;
    *trap_threshold = 0.0;
    *ck = 0.0;
    *epsal= 0.0;
    *gma = 1.0;
    *gmb = 1.0;
    *gmc = 1.0;
    *gmd = 1.0;
    *bh_threshold = 0.99;
    *dis_threshold = 0.01;
    *check_bh = 1;
    *check_dis = 1;
    *slicing_c = 0;
    *zepsdis = 0.0;
    *nc = 0;
    *lv = 0.0;
    *rv = 1.0;
    *level = 0;
    *lambda = 0.1;
    *deltx = 0.1;
    *x0 = 0.5;
    *amp = 1.0;
    *output_freq = 1;
    *output_freq_phi = 1;
    *myzero = 0.0;
    *x_max = 1;
    *x_min = 0;
    *Nx = 128;
    *steps = 0;
    *mass = 1.0;
    *omega = 0.0;
    *init_tol = 1.0e-8;
    *evol_tol = 1.0e-8;
    *read_lo_lev = 0;
    *compute_ire = 0;
    *flat_spacetime = 0;
    *freeze_matter = 0;
    get_param(p_file,"level","long",1,level);
    get_param(p_file,"lambda","double",1,lambda);
    get_param(p_file,"deltx","double",1,deltx);
    get_param(p_file,"x0","double",1,x0);
    get_param(p_file,"amp","double",1,amp);
    get_param(p_file,"lambda","double",1,lambda);
    get_param(p_file,"output_freq","long",1,output_freq);
    get_param(p_file,"output_freq_phi","long",1,output_freq_phi);
    get_param(p_file,"myzero","double",1,myzero);
    get_param(p_file,"x_max","double",1,x_max);
    get_param(p_file,"Nx","long",1,Nx);
    get_param(p_file,"x_min","double",1,x_min);
    get_param(p_file,"steps","long",1,steps);
    get_param(p_file,"flat_spacetime","long",1,flat_spacetime);
    get_param(p_file,"freeze_matter","long",1,freeze_matter);
    get_param(p_file,"mass","double",1,mass);
    get_param(p_file,"omega","double",1,omega);
    get_param(p_file,"init_tol","double",1,init_tol);
    get_param(p_file,"evol_tol","double",1,evol_tol);
    printf("reading low lev\n");
    get_param(p_file,"read_lo_lev","long",1,read_lo_lev);
    printf("read low lev in function = %d \n",*read_lo_lev);
    get_param(p_file,"compute_ire","long",1,compute_ire);
    get_param(p_file,"lv","double",1,lv);
    get_param(p_file,"rv","double",1,rv);
    get_param(p_file,"epsdis","double",1,zepsdis);
    get_param(p_file,"compact_type","long",1,nc);
    get_param(p_file,"slicing_c","long",1,slicing_c);
    get_param(p_file,"check_bh","long",1,check_bh);
    get_param(p_file,"check_dis","long",1,check_dis);
    get_param(p_file,"bh_threshold","double",1,bh_threshold);
    get_param(p_file,"dis_threshold","double",1,dis_threshold);
    get_param(p_file,"gma","double",1,gma);
    get_param(p_file,"gmb","double",1,gmb);
    get_param(p_file,"gmc","double",1,gmc);
    get_param(p_file,"gmd","double",1,gmd);
    get_param(p_file,"ck","double",1,ck);
    get_param(p_file,"epsal","double",1,epsal);
    get_param(p_file,"trap_threshold","double",1,trap_threshold);
    get_param(p_file,"mus","double",1,mus);
    get_param(p_file,"eta","double",1,eta);
    get_param(p_file,"vee","double",1,vee);
    get_param(p_file,"advc","double",1,advc);
    get_param(p_file,"min_evol_iter","long",1,min_evol_iter);
}

void find_inner_radious(double *mass, int Nx, double *in_rd){
    int i;
    *in_rd = 0.0;
    i = 0;
    while (mass[i] < 0.01*mass[Nx-1]) {
        i++;
    }
    *in_rd = x[i];
}

int main(int argc, char **argv) {

    double bhmass;
    double tau;
    double trap_threshold;
    double in_rd;
    double min_theta;
    double bh_threshold,dis_threshold;
    int adaptive_time = 0;
    double dtphi_center;
    int first_adt_warn = 1;
    int i = 0;
    int flat_spacetime = 0;
    int freeze_matter = 0;
    int compute_ire = 0;
    char pfile[64];
    int read_lo_lev = 0;
    int nc;
    FILE* bh_h;
    FILE* dis_h;
    FILE* mf;
    FILE *phi_0_h;
    double gmb;
    double advc;
    int min_evol_iter;
    phys_bdy[0] = 1;
    phys_bdy[1] = 1;
    strcpy(pfile,argv[1]);

    /* Initialization of Coordinate: */
    dim =1;
    read_params(pfile,&deltx, &x0, &amp,&mass, &level,&lambda,&output_freq,&output_freq_phi,&myzero,&x_max,&Nx,&x_min,&steps,&read_lo_lev,&compute_ire,&omega,&init_tol,&evol_tol,&flat_spacetime,&freeze_matter,&lv,&rv,&nc,&zepsdis,&slicing_c,&check_bh,&check_dis,&bh_threshold,&dis_threshold,&gma,&gmb,&ck,&epsal,&trap_threshold,&mus,&eta,&gmc,&gmd,&vee,&min_evol_iter,&advc);

    if (level == 0) { read_lo_lev = 0; }

    Nx = Nx*(int)pow(2.0,(double)level)+1;
    steps = steps*(int)pow(2.0,(double)level);
    /* Allocating Memory to Grid Functions: */

    double *res_rpsi_direct;
    double *res_lamx_val;
    double *dmdr;
    double *np1_tt;
    double max_tt = 0.0;

    np1_tt = vec_alloc(Nx);

    n_DLamx = vec_alloc(Nx);
    np1_DLamx = vec_alloc(Nx);

    n_BB = vec_alloc(Nx);
    np1_BB = vec_alloc(Nx);

    n_divbeta = vec_alloc(Nx);
    np1_divbeta = vec_alloc(Nx);
    dmdr = vec_alloc(Nx);
    ctfm = vec_alloc(Nx);
    ctfmp = vec_alloc(Nx);
    octfmp = vec_alloc(Nx);
    nph_th = vec_alloc(Nx);
    n_mass = vec_alloc(Nx);
    n_mk = vec_alloc(Nx);
    n_mass2 = vec_alloc(Nx);
    n_tmr = vec_alloc(Nx);
    res_lamx_val = vec_alloc(Nx);
    res_rpsi_direct = vec_alloc(Nx);
    n_ipsidot= vec_alloc(Nx);
    n_rpsidot= vec_alloc(Nx);

    n_ire_daath = vec_alloc(Nx);
    n_ire_daax = vec_alloc(Nx);
    n_ire_dlamx = vec_alloc(Nx);
    n_ire_dph = vec_alloc(Nx);
    n_ire_dtk = vec_alloc(Nx);
    n_ire_dgx = vec_alloc(Nx);
    n_ire_dgth = vec_alloc(Nx);

    n_restx = vec_alloc(Nx);
    n_restt = vec_alloc(Nx);
    n_resxx = vec_alloc(Nx);
    n_resthth = vec_alloc(Nx);
    n_reshs = vec_alloc(Nx);
    nm1_b1 = vec_alloc(Nx);
    nm1_a1 = vec_alloc(Nx);

    x = vec_alloc(Nx);
    n_A = vec_alloc(1*Nx);

    n_psi = vec_alloc(1*Nx);
    no_psi = vec_alloc(1*Nx);
    np1_psi = vec_alloc(1*Nx);

    n_C1s = vec_alloc(1*Nx);
    n_C5s = vec_alloc(1*Nx);
    n_C5psi = vec_alloc(1*Nx);
    n_PHI2 = vec_alloc(1*Nx);
    n_PI2 = vec_alloc(1*Nx);
    n_U = vec_alloc(1*Nx);
    np1_U = vec_alloc(1*Nx);

    n_rpsi = vec_alloc(1*Nx);
    n_ipsi = vec_alloc(1*Nx);
    n_rPI = vec_alloc(1*Nx);
    n_iPI = vec_alloc(1*Nx);
    n_rPHI = vec_alloc(1*Nx);
    n_iPHI = vec_alloc(1*Nx);
    n_iPIb2 = vec_alloc(1*Nx);
    n_rPIb2 = vec_alloc(1*Nx);

    np1_rpsi = vec_alloc(1*Nx);
    np1_ipsi = vec_alloc(1*Nx);
    np1_rPI = vec_alloc(1*Nx);
    np1_iPI = vec_alloc(1*Nx);
    np1_rPHI = vec_alloc(1*Nx);
    np1_iPHI = vec_alloc(1*Nx);
    np1_iPIb2 = vec_alloc(1*Nx);
    np1_rPIb2 = vec_alloc(1*Nx);

    np1_PHI2 = vec_alloc(1*Nx);
    np1_PI2 = vec_alloc(1*Nx);


    n_a2 = vec_alloc(1*Nx);
    np1_a2 =vec_alloc(1*Nx);

    n_a1 = vec_alloc(1*Nx);
    np1_a1 = vec_alloc(1*Nx);

    n_b1 = vec_alloc(1*Nx);
    np1_b1 =vec_alloc(1*Nx);


    n_Athth = vec_alloc(1*Nx);
    n_Axx = vec_alloc(1*Nx);
    n_B = vec_alloc(1*Nx);
    n_DDLxx = vec_alloc(1*Nx);
    n_JSxx = vec_alloc(1*Nx);
    n_K = vec_alloc(1*Nx);
    n_Lamx = vec_alloc(1*Nx);
    n_Pxx = vec_alloc(1*Nx);
    n_RRxx = vec_alloc(1*Nx);
    n_Sx = vec_alloc(1*Nx);
    n_TS = vec_alloc(1*Nx);
    n_Uxx = vec_alloc(1*Nx);
    n_alpha = vec_alloc(1*Nx);
    n_beta = vec_alloc(1*Nx);
    n_metxx = vec_alloc(1*Nx);
    n_phi = vec_alloc(1*Nx);
    n_rho = vec_alloc(1*Nx);
    n_DDLthetatheta = vec_alloc(1*Nx);
    n_JSthth = vec_alloc(1*Nx);
    n_Pthetatheta = vec_alloc(1*Nx);
    n_RRthetatheta = vec_alloc(1*Nx);
    n_Uthetatheta = vec_alloc(1*Nx);
    n_em4phi = vec_alloc(1*Nx);
    n_metthetatheta = vec_alloc(1*Nx);



    nph_ipsi = vec_alloc(1*Nx);
    nph_rpsi = vec_alloc(1*Nx);
    nph_iPIb2 = vec_alloc(1*Nx);
    nph_rPIb2 = vec_alloc(1*Nx);



    nph_A = vec_alloc(1*Nx);
    nph_B = vec_alloc(1*Nx);
    nph_Axx = vec_alloc(1*Nx);
    nph_Athth = vec_alloc(1*Nx);
    nph_K = vec_alloc(1*Nx);
    nph_phi = vec_alloc(1*Nx);
    nph_Lamx = vec_alloc(1*Nx);

    np1_A = vec_alloc(1*Nx);
    np1_Athth = vec_alloc(1*Nx);
    np1_Axx = vec_alloc(1*Nx);
    np1_B = vec_alloc(1*Nx);
    np1_DDLxx = vec_alloc(1*Nx);
    np1_JSxx = vec_alloc(1*Nx);
    np1_K = vec_alloc(1*Nx);
    np1_Lamx = vec_alloc(1*Nx);
    np1_Pxx = vec_alloc(1*Nx);
    np1_RRxx = vec_alloc(1*Nx);
    np1_Sx = vec_alloc(1*Nx);
    np1_TS = vec_alloc(1*Nx);
    np1_Uxx = vec_alloc(1*Nx);
    np1_alpha = vec_alloc(1*Nx);
    np1_beta = vec_alloc(1*Nx);
    np1_metxx = vec_alloc(1*Nx);
    np1_phi = vec_alloc(1*Nx);
    np1_rho = vec_alloc(1*Nx);
    np1_DDLthetatheta = vec_alloc(1*Nx);
    np1_JSthth = vec_alloc(1*Nx);
    np1_Pthetatheta = vec_alloc(1*Nx);
    np1_RRthetatheta = vec_alloc(1*Nx);
    np1_Uthetatheta = vec_alloc(1*Nx);
    np1_em4phi = vec_alloc(1*Nx);
    np1_metthetatheta = vec_alloc(1*Nx);

    nm1_A = vec_alloc(1*Nx);
    nm1_Athth = vec_alloc(1*Nx);
    nm1_Axx = vec_alloc(1*Nx);
    nm1_B = vec_alloc(1*Nx);
    nm1_DDLxx = vec_alloc(1*Nx);
    nm1_JSxx = vec_alloc(1*Nx);
    nm1_K = vec_alloc(1*Nx);
    nm1_Lamx = vec_alloc(1*Nx);
    nm1_Pxx = vec_alloc(1*Nx);
    nm1_RRxx = vec_alloc(1*Nx);
    nm1_Sx = vec_alloc(1*Nx);
    nm1_TS = vec_alloc(1*Nx);
    nm1_Uxx = vec_alloc(1*Nx);
    nm1_alpha = vec_alloc(1*Nx);
    nm1_beta = vec_alloc(1*Nx);
    nm1_metxx = vec_alloc(1*Nx);
    nm1_phi = vec_alloc(1*Nx);
    nm1_rho = vec_alloc(1*Nx);
    nm1_DDLthetatheta = vec_alloc(1*Nx);
    nm1_JSthth = vec_alloc(1*Nx);
    nm1_Pthetatheta = vec_alloc(1*Nx);
    nm1_RRthetatheta = vec_alloc(1*Nx);
    nm1_Uthetatheta = vec_alloc(1*Nx);
    nm1_em4phi = vec_alloc(1*Nx);
    nm1_metthetatheta = vec_alloc(1*Nx);

    nm1_ipsi=vec_alloc(Nx);
    nm1_rpsi=vec_alloc(Nx);
    nm1_iPIb2=vec_alloc(Nx);
    nm1_rPIb2=vec_alloc(Nx);


    hx = (x_max-x_min)/(Nx-1);
    //dvumsh(x,Nx,x_min,x_max);
    for ( i=0; i < Nx; i++) {
        x[i] = x_min + ((double)(i))/((double)((Nx-1))) * (x_max-x_min);
    }
    x[0] = x_min;
    x[Nx-1] = x_max;
    ht = lambda*sqrt( 0.0+ hx*hx);
    double ht0 = ht;
    shape[0]=Nx;
    bbox[0]=x_min;
    bbox[1]=x_max;
    time=0.0;

    // Initialization 

    // initializing spacetime to flat static 

    set_zero(octfmp,Nx);
    set_one(ctfmp,Nx);
    rdvcpy(ctfm,x,Nx);


    // Initializing coordinate:

    if ( nc == 0 ) {
        printf("Initializing to non compact\n");
        printf("Nx = %d\n",Nx); 
        for(i=0;i<Nx;i++) { 
            ctfm[i] = x[i];
            ctfmp[i] = 1.0;
            octfmp[i] = 0.0;
        }
    } else {
        if (nc ==1) {
            printf("Initializing to compact type 1\n");
            init_ctfm_(x,&Nx,&lv,&rv,ctfm);
            init_ctfmp_(x,&Nx,&rv,ctfmp);
            init_octfmp_(x,&Nx,&rv,octfmp);
        }
        if (nc == 2) {
            printf("Initializing to compact type 2\n");
            init_ctfm2_(x,&Nx,&lv,&rv,ctfm);
            init_ctfmp2_(x,&Nx,&lv,&rv,ctfmp);
            init_octfmp2_(x,&Nx,&lv,&rv,octfmp);
            // init_ctfm2_(x,&Nx,&rv,ctfm);
            //init_ctfmp2_(x,&Nx,&rv,ctfmp);
            //init_octfmp2_(x,&Nx,&rv,octfmp);
        }
    }

    printf("X_max = %1.14e\n",x[Nx-1]);
    printf("X_min = %1.14e\n",x[0]);

    printf("R_max = %1.14e\n",ctfm[Nx-1]);
    printf("R_min = %1.14e\n",ctfm[0]);
    printf("R_med = %1.14e\n",ctfm[Nx/2]);

    gft_out_bbox("ctfm.sdf",time,shape,dim,bbox,ctfm);
    gft_out_bbox("octfmp.sdf",time,shape,dim,bbox,octfmp);
    gft_out_bbox("ctfmp.sdf",time,shape,dim,bbox,ctfmp);


    set_zero(n_restx,Nx);
    set_zero(n_resxx,Nx);
    set_zero(n_restt,Nx);
    set_zero(n_resthth,Nx);

    set_zero(n_reshs,Nx);

    /* Setting geometry to flat: */

    init_a_(x,&Nx,&myzero,n_A);

    init_b_(x,&Nx,&myzero,n_B);

    init_axx_(x,&Nx,&myzero,n_Axx);

    init_athth_(x,&Nx,&myzero,n_Athth);

    //init_alpha_(x,&Nx,&myzero,n_alpha);
    set_one(n_alpha,Nx);

    init_beta_(x,&Nx,&myzero,n_beta);

    init_lamx_(x,&Nx,&myzero,n_Lamx);

    init_tk_(x,&Nx,&myzero,n_K);

    set_one(n_psi,Nx);
    set_one(no_psi,Nx);

    compute_b1_(n_B,n_psi,&Nx,n_b1);


    compute_c1_(ctfm,ctfmp,n_A,n_B,n_Lamx,octfmp,x,&Nx,&hx,&myzero,phys_bdy,n_C1s);

    /* Gaussian rpsi, zero ipsi */
    init_rpsi_(ctfm,&Nx,&amp,&deltx,&x0,n_rpsi);
    init_ipsi_(x,&Nx,&myzero,n_ipsi);

    init_iphi_(ctfmp,x,&Nx,&myzero,n_iPHI);
    init_rphi_(ctfm,ctfmp,&Nx,&amp,&deltx,&hx,&x0,n_rPHI);

    init_ipsidot_(ctfm,&Nx,&amp,&deltx,&omega,&x0,n_ipsidot);
    init_rpsidot_(x,&Nx,&myzero,n_rpsidot);

    /* Steps to compute rho using geometry: */
    compute_a2_(n_A,n_psi,&Nx,n_a2);
    compute_a1_(n_a2,&Nx,n_a1);

    init_ipi_(n_a1,n_alpha,n_beta,n_iPHI,n_ipsidot,&Nx,n_iPI);
    init_rpi_(n_a1,n_alpha,n_beta,n_rPHI,n_rpsidot,&Nx,n_rPI);

    compute_pi2_(n_iPI,n_rPI,&Nx,n_PI2);
    compute_phi2_(n_iPHI,n_rPHI,&Nx,n_PHI2);

    compute_u_(n_ipsi,n_rpsi,&Nx,&mass,n_U);

    compute_rho_(n_PHI2,n_PI2,n_U,n_a2,&Nx,n_rho);

    compute_c5_(n_A,n_Athth,n_Axx,n_B,n_K,n_rho,x,&Nx,&myzero,phys_bdy,n_C5s);

    compute_c5psi_(n_C5s,n_psi,&Nx,n_C5psi);

    /* End of steps to compute rho and c's */

    time = 0;
    tau = 0.0;

    gft_out_bbox("ctfm.sdf",time,shape,dim,bbox,ctfm);
    gft_out_bbox("ctfmp.sdf",time,shape,dim,bbox,ctfmp);
    gft_out_bbox("octfmp.sdf",time,shape,dim,bbox,octfmp);

    /* For pure gauge evolution: */

    int pure_gauge = 0;
    if (pure_gauge == 1) {
        double ampa;
        ampa = -0.05;
        init_rpsi_(ctfm,&Nx,&ampa,&deltx,&x0,n_alpha);
        for(i=0;i<Nx;i++) { 
            n_alpha[i] = n_alpha[i]+1.0;
        }
    }


    int etime,rtime,wtime;
    double rd=1.0;

    printf("Starting Initial data solver \n");
    printf("read low lev = %d\n",read_lo_lev);
    if (flat_spacetime == 0) {
        if (read_lo_lev == 1) {
            printf("Reading psilm1.sdf...\n");
            SDF psilm1;
            int rc;
            rc = sdf_read(&psilm1,"./init/psilm1.sdf",1,&time);
            if (rc == 0) {
                printf("Unable to read psilm1.sdf file\n");
                exit(1);
            }

            if ( (Nx-1) / (psilm1.Nx - 1)  == 2 ) {
                int method = 2;intplt_(psilm1.x,psilm1.data,&psilm1.Nx,x,n_psi,&Nx,&method,&intplt_rc);
                if (intplt_rc != 1) {
                    printf("Interpolation failed\n"); 
                    exit(1);
                } else {
                    /* recalculating the rho and c's */
                    compute_a2_(n_A,n_psi,&Nx,n_a2);
                    compute_a1_(n_a2,&Nx,n_a1);

                    /* alpha might need update here */
                    init_ipi_(n_a1,n_alpha,n_beta,n_iPHI,n_ipsidot,&Nx,n_iPI);
                    init_rpi_(n_a1,n_alpha,n_beta,n_rPHI,n_rpsidot,&Nx,n_rPI);

                    compute_pi2_(n_iPI,n_rPI,&Nx,n_PI2);
                    compute_phi2_(n_iPHI,n_rPHI,&Nx,n_PHI2);

                    compute_u_(n_ipsi,n_rpsi,&Nx,&mass,n_U);

                    compute_rho_(n_PHI2,n_PI2,n_U,n_a2,&Nx,n_rho);

                    compute_c5_(n_A,n_Athth,n_Axx,n_B,n_K,n_rho,x,&Nx,&myzero,phys_bdy,n_C5s);

                    compute_c5psi_(n_C5s,n_psi,&Nx,n_C5psi);
                    rdvcpy(no_psi,n_psi,Nx);
                    printf("interpolation completed\n");
                }

            } else {
                printf("incompatible psilm1\n");
                exit(1);
            }
        }
        resid_hs3_(ctfm,ctfmp,n_A,n_B,n_C1s,n_C5psi,n_psi,octfmp,&Nx,&hx,&myzero,phys_bdy,&res_psi);
        printf("HC solver begins, resid = %1.14e\n",res_psi);
        i = 0;
        while( rd > init_tol )  {
            solve_hs3_(ctfm,ctfmp,n_A,n_B,n_C1s,n_C5psi,n_psi,octfmp,&Nx,&hx,&myzero,phys_bdy,no_psi);
            /* recalculating the rho and c's */
            compute_a2_(n_A,no_psi,&Nx,n_a2);
            compute_a1_(n_a2,&Nx,n_a1);

            /* alpha might need update here */
            init_ipi_(n_a1,n_alpha,n_beta,n_iPHI,n_ipsidot,&Nx,n_iPI);
            init_rpi_(n_a1,n_alpha,n_beta,n_rPHI,n_rpsidot,&Nx,n_rPI);

            compute_pi2_(n_iPI,n_rPI,&Nx,n_PI2);
            compute_phi2_(n_iPHI,n_rPHI,&Nx,n_PHI2);

            compute_u_(n_ipsi,n_rpsi,&Nx,&mass,n_U);

            compute_rho_(n_PHI2,n_PI2,n_U,n_a2,&Nx,n_rho);

            compute_c5_(n_A,n_Athth,n_Axx,n_B,n_K,n_rho,x,&Nx,&myzero,phys_bdy,n_C5s);

            compute_c5psi_(n_C5s,no_psi,&Nx,n_C5psi);

            resid_hs3_(ctfm,ctfmp,n_A,n_B,n_C1s,n_C5psi,n_psi,octfmp,&Nx,&hx,&myzero,phys_bdy,&res_psi);
            rel_diff_(n_psi,no_psi,&Nx,&rd);
            rdvcpy(n_psi,no_psi,Nx);
            if (i % (100000) == 0) { 
                printf("Iter %d of %d complete. resid = %1.14e, rel_dff = %1.14e\n",i/1000,Nx,res_psi,rd);
            }

            i++;
        }
        printf("Iteration complete after: %d iteration\n",i);
    } 
    else {
        set_one(n_psi,Nx);
    }





    /******** END OF INITIAL DATA SOLVER **************/

    /* Computing everything using new psi */
    compute_phi_(n_psi,&Nx,n_phi);
    compute_a2_(n_A,n_psi,&Nx,n_a2);
    compute_a1_(n_a2,&Nx,n_a1);
    compute_b1_(n_B,n_psi,&Nx,n_b1);

    init_ipi_(n_a1,n_alpha,n_beta,n_iPHI,n_ipsidot,&Nx,n_iPI);
    init_rpi_(n_a1,n_alpha,n_beta,n_rPHI,n_rpsidot,&Nx,n_rPI);

    compute_pi2_(n_iPI,n_rPI,&Nx,n_PI2);
    compute_phi2_(n_iPHI,n_rPHI,&Nx,n_PHI2);

    /* End of steps to compute rho and c's */

    init_ipib2_(n_b1,n_iPI,&Nx,n_iPIb2);
    init_rpib2_(n_b1,n_rPI,&Nx,n_rPIb2);


    compute_db_(ctfm,ctfmp,n_A,n_B,n_beta,x,&Nx,&hx,&myzero,phys_bdy,n_divbeta);

    compute_u_(n_ipsi,n_rpsi,&Nx,&mass,n_U);

    compute_rho_(n_PHI2,n_PI2,n_U,n_a2,&Nx,n_rho);

    compute_jsxx_(n_PHI2,n_PI2,n_U,n_a2,&Nx,n_JSxx);

    compute_jsthth_(n_PHI2,n_PI2,n_U,n_a2,n_b1,&Nx,n_JSthth);

    compute_ts_(n_PHI2,n_PI2,n_U,n_a2,&Nx,n_TS);

    compute_sx_(n_a1,n_iPHI,n_iPI,n_rPHI,n_rPI,&Nx,n_Sx);

    compute_emph_(n_phi,x,&Nx,&myzero,phys_bdy,n_em4phi);

    compute_ux_(n_A,n_Axx,n_K,x,&Nx,&myzero,phys_bdy,n_Uxx);



    compute_utheta_(n_Athth,n_B,n_K,x,&Nx,&myzero,phys_bdy,n_Uthetatheta);

    compute_px_(ctfmp,n_Axx,n_beta,n_divbeta,x,&Nx,&hx,&myzero,&vee,phys_bdy,n_Pxx);

    compute_ptheta_(ctfm,ctfmp,n_Athth,n_beta,n_divbeta,x,&Nx,&hx,&myzero,&vee,phys_bdy,n_Pthetatheta);


    compute_metx_(n_A,n_phi,x,&Nx,&myzero,phys_bdy,n_metxx);

    compute_mettheta_(n_B,n_phi,x,&Nx,&myzero,phys_bdy,n_metthetatheta);


    compute_ddlx_(ctfmp,n_alpha,n_metxx,octfmp,x,&Nx,&hx,&myzero,phys_bdy,n_DDLxx);

    compute_ddltheta_(ctfm,ctfmp,n_alpha,n_metthetatheta,n_metxx,x,&Nx,&hx,&myzero,phys_bdy,n_DDLthetatheta);

    compute_rrx_(ctfm,ctfmp,n_A,n_B,n_Lamx,n_phi,octfmp,x,&Nx,&hx,&myzero,phys_bdy,n_RRxx);

    compute_rrtheta_(ctfm,ctfmp,n_A,n_B,n_Lamx,n_phi,octfmp,x,&Nx,&hx,&myzero,phys_bdy,n_RRthetatheta);

    gft_out_bbox("C1.sdf",time,shape,dim,bbox,n_C1s);
    gft_out_bbox("C5.sdf",time,shape,dim,bbox,n_C5s);
    gft_out_bbox("C5psi.sdf",time,shape,dim,bbox,n_C5psi);

    gft_out_bbox("psi.sdf",time,shape,dim,bbox,n_psi);
    gft_out_bbox("rho.sdf",time,shape,dim,bbox,n_rho);


    gft_out_bbox("A.sdf",time,shape,dim,bbox,n_A);
    gft_out_bbox("B.sdf",time,shape,dim,bbox,n_B);

    gft_out_bbox("phi.sdf",time,shape,dim,bbox,n_phi);
    gft_out_bbox("K.sdf",time,shape,dim,bbox,n_K);

    gft_out_bbox("Axx.sdf",time,shape,dim,bbox,n_Axx);
    gft_out_bbox("Athth.sdf",time,shape,dim,bbox,n_Athth);
    gft_out_bbox("alpha.sdf",time,shape,dim,bbox,n_alpha);
    gft_out_bbox("beta.sdf",time,shape,dim,bbox,n_beta);
    gft_out_bbox("a1.sdf",time,shape,dim,bbox,n_a1);
    gft_out_bbox("b1.sdf",time,shape,dim,bbox,n_b1);
    gft_out_bbox("rpsi.sdf",time,shape,dim,bbox,n_rpsi);
    gft_out_bbox("Lamx.sdf",time,shape,dim,bbox,n_Lamx);



    gft_out_bbox("JSxx.sdf",time,shape,dim,bbox,n_JSxx);
    gft_out_bbox("JSthth.sdf",time,shape,dim,bbox,n_JSthth);
    gft_out_bbox("Sx.sdf",time,shape,dim,bbox,n_Sx);
    gft_out_bbox("TS.sdf",time,shape,dim,bbox,n_TS);

    gft_out_bbox("em4phi.sdf",time,shape,dim,bbox,n_em4phi);
    gft_out_bbox("metxx.sdf",time,shape,dim,bbox,n_metxx);
    gft_out_bbox("metthetatheta.sdf",time,shape,dim,bbox,n_metthetatheta);
    gft_out_bbox("RRxx.sdf",time,shape,dim,bbox,n_RRxx);
    gft_out_bbox("RRthetatheta.sdf",time,shape,dim,bbox,n_RRthetatheta);

    gft_out_bbox("Pxx.sdf",time,shape,dim,bbox,n_Pxx);
    gft_out_bbox("Pthetatheta.sdf",time,shape,dim,bbox,n_Pthetatheta);
    gft_out_bbox("Uxx.sdf",time,shape,dim,bbox,n_Uxx);
    gft_out_bbox("Uthetatheta.sdf",time,shape,dim,bbox,n_Uthetatheta);
    gft_out_bbox("DDLxx.sdf",time,shape,dim,bbox,n_DDLxx);
    gft_out_bbox("DDLthetatheta.sdf",time,shape,dim,bbox,n_DDLthetatheta);
    gft_out_bbox("U.sdf",time,shape,dim,bbox,n_U);


    gft_out_bbox("iPIb2.sdf",time,shape,dim,bbox,n_iPIb2);
    gft_out_bbox("rPIb2.sdf",time,shape,dim,bbox,n_rPIb2);
    gft_out_bbox("PI2.sdf",time,shape,dim,bbox,n_PI2);
    gft_out_bbox("PHI2.sdf",time,shape,dim,bbox,n_PHI2);
    gft_out_bbox("iPHI.sdf",time,shape,dim,bbox,n_iPHI);
    gft_out_bbox("rPHI.sdf",time,shape,dim,bbox,n_rPHI);
    gft_out_bbox("iPI.sdf",time,shape,dim,bbox,n_iPI);
    gft_out_bbox("rPI.sdf",time,shape,dim,bbox,n_rPI);

    gft_out_bbox("ipsi.sdf",time,shape,dim,bbox,n_ipsi);
    gft_out_bbox("rpsi.sdf",time,shape,dim,bbox,n_rpsi);

    gft_out_bbox("ipsidot.sdf",time,shape,dim,bbox,n_ipsidot);
    gft_out_bbox("rpsidot.sdf",time,shape,dim,bbox,n_rpsidot);


    gft_out_bbox("restx.sdf",time,shape,dim,bbox,n_restx);
    gft_out_bbox("resxx.sdf",time,shape,dim,bbox,n_resxx);
    gft_out_bbox("restt.sdf",time,shape,dim,bbox,n_restt);
    gft_out_bbox("resthth.sdf",time,shape,dim,bbox,n_resthth);


    compute_mass_(n_A,n_B,n_psi, ctfm,ctfmp,n_rho,n_b1,n_a1,n_mass,n_mass2,n_tmr,&maxtmr,&bhmass,&hx,&Nx);
    mf=fopen("mass","w");
    fprintf(mf,"%1.14e %1.14e\n",time, n_mass[Nx-2]);

    double i_in_rd;
    find_inner_radious(n_mass2,Nx,&i_in_rd);

    compute_dmdr_(ctfm,n_rho,&Nx,dmdr);

    gft_out_bbox("dmdr.sdf",time,shape,dim,bbox,dmdr);

    gft_out_bbox("mass.sdf",time,shape,dim,bbox,n_mass);
    gft_out_bbox("tmr.sdf",time,shape,dim,bbox,n_tmr);

    ire_mk_(ctfm,ctfmp,n_A,n_Athth,n_Axx,n_B,n_K,n_Sx,n_psi,x,&Nx,&hx,&myzero,phys_bdy,n_mk);
    ire_val_lamx_(ctfm,ctfmp,n_A,n_B,n_Lamx,&Nx,&hx,res_lamx_val);
    ire_hs_(ctfm,ctfmp,n_A,n_Athth,n_Axx,n_B,n_K,n_Lamx,n_psi,n_rho,octfmp,&Nx,&hx,n_reshs);


    eval_dlam_(ctfm,ctfmp,n_A,n_Athth,n_Axx,n_B,n_K,n_Lamx,n_Sx,n_alpha,n_beta,n_divbeta,n_phi,octfmp,x,&Nx,&hx,&myzero,&vee,phys_bdy,n_DLamx);
    set_zero(n_BB,Nx);



    gft_out_bbox("reshs.sdf",time,shape,dim,bbox,n_reshs);
    gft_out_bbox("resmk.sdf",time,shape,dim,bbox,n_mk);
    gft_out_bbox("reslamx.sdf",time,shape,dim,bbox,res_lamx_val);



    wtime = 0;
    etime = 0;
    rtime = 0;

    double tot_res = 1.0;
    int j = 0;
    double tot_ire = 1.0;
    double norm_phi,norm_A,norm_B,norm_K,norm_Axx,norm_Athth,norm_Lamx,norm_ipsi,norm_rpsi,norm_rPIb2,norm_iPIb2;
    double norm_alpha,norm_beta;
    bool wave_reached_org = false;
    double norm_BB;

    gettimeofday(&start_time,NULL);


    rdvcpy(np1_phi,n_phi,Nx);
    rdvcpy(np1_A,n_A,Nx);
    rdvcpy(np1_B,n_B,Nx);
    rdvcpy(np1_Axx,n_Axx,Nx);
    rdvcpy(np1_Athth,n_Athth,Nx);
    rdvcpy(np1_Lamx,n_Lamx,Nx);
    rdvcpy(np1_K,n_K,Nx);

    rdvcpy(np1_alpha,n_alpha,Nx);
    rdvcpy(np1_beta,n_beta,Nx);


    rdvcpy(np1_psi,n_psi,Nx);
    rdvcpy(np1_a2,n_a2,Nx);
    rdvcpy(np1_a1,n_a1,Nx);
    rdvcpy(np1_b1,n_b1,Nx);

    rdvcpy(np1_rpsi,n_rpsi,Nx);
    rdvcpy(np1_ipsi,n_ipsi,Nx);
    rdvcpy(np1_rPHI,n_rPHI,Nx);
    rdvcpy(np1_iPHI,n_iPHI,Nx);
    rdvcpy(np1_iPI,n_iPI,Nx);
    rdvcpy(np1_rPI,n_rPI,Nx);
    rdvcpy(np1_PHI2,n_PHI2,Nx);
    rdvcpy(np1_PI2,n_PI2,Nx);
    rdvcpy(np1_iPIb2,n_iPIb2,Nx);
    rdvcpy(np1_rPIb2,n_rPIb2,Nx);

    rdvcpy(np1_U,n_U,Nx);
    rdvcpy(np1_rho,n_rho,Nx);
    rdvcpy(np1_JSxx,n_JSxx,Nx);
    rdvcpy(np1_Sx,n_Sx,Nx);
    rdvcpy(np1_TS,n_TS,Nx);
    rdvcpy(np1_JSthth,n_JSthth,Nx);

    rdvcpy(np1_divbeta,n_divbeta,Nx);

    rdvcpy(np1_em4phi,n_em4phi,Nx);

    rdvcpy(np1_metxx,n_metxx,Nx);
    rdvcpy(np1_metthetatheta,n_metthetatheta,Nx);

    rdvcpy(np1_DDLxx,n_DDLxx,Nx);
    rdvcpy(np1_DDLthetatheta,n_DDLthetatheta,Nx);

    rdvcpy(np1_RRxx,n_RRxx,Nx);
    rdvcpy(np1_RRthetatheta,n_RRthetatheta,Nx);

    rdvcpy(np1_Uxx,n_Uxx,Nx);
    rdvcpy(np1_Uthetatheta,n_Uthetatheta,Nx);

    rdvcpy(np1_Pxx,n_Pxx,Nx);
    rdvcpy(np1_Pthetatheta,n_Pthetatheta,Nx);

    rdvcpy(np1_DLamx,n_DLamx,Nx);
    rdvcpy(np1_BB,n_BB,Nx);




    phi_0_h=fopen("phi_0_t","w");

    fprintf(phi_0_h,"%1.14e %1.14e %1.14e\n",tau,n_rpsi[0],0.0E0);

    printf("Starting Evolution...\n");
    for (i=0; i<steps; i++) {
        j = 0;
        tot_res = 1.0;
        norm_phi = l2norm(Nx,n_phi);
        norm_A = l2norm(Nx,n_A);
        norm_B = l2norm(Nx,n_B);
        norm_K = l2norm(Nx,n_K);
        norm_Axx = l2norm(Nx,n_Axx);
        norm_Athth = l2norm(Nx,n_Athth);
        norm_Lamx = l2norm(Nx,n_Lamx);
        norm_ipsi = l2norm(Nx,n_ipsi);
        norm_rpsi = l2norm(Nx,n_rpsi);
        norm_rPIb2 = l2norm(Nx,n_rPIb2);
        norm_iPIb2 = l2norm(Nx,n_iPIb2);
        if (slicing_c == 1) {
            norm_alpha = l2norm(Nx,n_alpha);
        }
        if (slicing_c == 2) {
            norm_alpha = l2norm(Nx,n_alpha);
        }
        if (slicing_c == 3) {
            norm_alpha = l2norm(Nx,n_alpha);
        }
        if (slicing_c == 4) {
            norm_alpha = l2norm(Nx,n_alpha);
            norm_beta = l2norm(Nx,n_beta);
        }

        if (slicing_c == 5) {
            norm_alpha = l2norm(Nx,n_alpha);
            norm_beta = l2norm(Nx,n_beta);
        }


        if (slicing_c == 6) {
            norm_alpha = l2norm(Nx,n_alpha);
            norm_beta = l2norm(Nx,n_beta);
        }

        if (slicing_c == 7) {
            norm_alpha = l2norm(Nx,n_alpha);
            norm_beta = l2norm(Nx,n_beta);
        }

        if (slicing_c == 8) {
            norm_alpha = l2norm(Nx,n_alpha);
            norm_beta = l2norm(Nx,n_beta);
            norm_BB = l2norm(Nx,n_BB);
        }

        while ( (tot_res > evol_tol ) ||  ( j < min_evol_iter ) ) {
            j++;

            gettimeofday(&estart_time,NULL);

            if (flat_spacetime == 0) {
                evolve_dph_(ctfmp,n_K,n_alpha,n_beta,n_divbeta,n_phi,np1_K,np1_alpha,np1_beta,np1_divbeta,np1_phi,x,&Nx,&ht,&hx,&myzero,&zepsdis,phys_bdy,np1_phi);

                evolve_dgx_(ctfmp,n_A,n_Axx,n_alpha,n_beta,n_divbeta,np1_A,np1_Axx,np1_alpha,np1_beta,np1_divbeta,x,&Nx,&ht,&hx,&myzero,&vee,&zepsdis,phys_bdy,np1_A);

                compute_b_(np1_A,&Nx,np1_B);
                //evolve_dgtheta_(ctfm,ctfmp,n_Athth,n_B,n_alpha,n_beta,n_divbeta,np1_Athth,np1_B,np1_alpha,np1_beta,np1_divbeta,x,&Nx,&ht,&hx,&myzero,&vee,&zepsdis,phys_bdy,np1_B);

                evolve_dtk_(ctfmp,n_A,n_Athth,n_Axx,n_B,n_DDLthetatheta,n_DDLxx,n_K,n_TS,n_alpha,n_beta,n_metthetatheta,n_metxx,n_rho,np1_A,np1_Athth,np1_Axx,np1_B,np1_DDLthetatheta,np1_DDLxx,np1_K,np1_TS,np1_alpha,np1_beta,np1_metthetatheta,np1_metxx,np1_rho,x,&Nx,&ht,&hx,&myzero,&zepsdis,phys_bdy,np1_K);

                evolve_daax_(n_Axx,n_DDLthetatheta,n_DDLxx,n_JSthth,n_JSxx,n_Pxx,n_RRthetatheta,n_RRxx,n_Uxx,n_alpha,n_em4phi,n_metthetatheta,n_metxx,np1_Axx,np1_DDLthetatheta,np1_DDLxx,np1_JSthth,np1_JSxx,np1_Pxx,np1_RRthetatheta,np1_RRxx,np1_Uxx,np1_alpha,np1_em4phi,np1_metthetatheta,np1_metxx,x,&Nx,&ht,&hx,&myzero,&zepsdis,phys_bdy,np1_Axx);


                compute_athth_(np1_Axx,np1_A,np1_B,&Nx,np1_Athth);
                //evolve_daatheta_(n_Athth,n_DDLthetatheta,n_DDLxx,n_JSthth,n_JSxx,n_Pthetatheta,n_RRthetatheta,n_RRxx,n_Uthetatheta,n_alpha,n_em4phi,n_metthetatheta,n_metxx,np1_Athth,np1_DDLthetatheta,np1_DDLxx,np1_JSthth,np1_JSxx,np1_Pthetatheta,np1_RRthetatheta,np1_RRxx,np1_Uthetatheta,np1_alpha,np1_em4phi,np1_metthetatheta,np1_metxx,x,&Nx,&ht,&hx,&myzero,&zepsdis,phys_bdy,np1_Athth);


                evolve_dlamx_(ctfm,ctfmp,n_A,n_Athth,n_Axx,n_B,n_K,n_Lamx,n_Sx,n_alpha,n_beta,n_divbeta,n_phi,np1_A,np1_Athth,np1_Axx,np1_B,np1_K,np1_Lamx,np1_Sx,np1_alpha,np1_beta,np1_divbeta,np1_phi,octfmp,x,&Nx,&ht,&hx,&myzero,&vee,&zepsdis,phys_bdy,np1_Lamx);


                /*
                   rdvcpy(np1_A,nph_A,Nx);
                   rdvcpy(np1_B,nph_B,Nx);
                   rdvcpy(np1_Axx,nph_Axx,Nx);
                   rdvcpy(np1_Athth,nph_Athth,Nx);
                   rdvcpy(np1_Lamx,nph_Lamx,Nx);
                   rdvcpy(np1_phi,nph_phi,Nx);
                   rdvcpy(np1_K,nph_K,Nx);
                   */
            }

            gettimeofday(&estop_time,NULL);

            timersub(&estop_time,&estart_time,&eelapsed_time);
            etime += eelapsed_time.tv_sec*1000000+eelapsed_time.tv_usec;

            gettimeofday(&wstart_time,NULL);

            compute_psi_(np1_phi,&Nx,np1_psi);

            if (slicing_c == 1) { 
                // 1 + log  zero shift
                evol_alphadyopl_(n_K,n_alpha,np1_K,np1_alpha,x,&Nx,&ht,&myzero,phys_bdy,np1_alpha);
            }

            if (slicing_c == 2) {
                // 1 + log dynamical zero shift
                evol_alphaopl_(n_K,n_alpha,np1_K,np1_alpha,x,&Nx,&gma,&gmb,&ht,&myzero,phys_bdy,np1_alpha);
            }

            if (slicing_c == 3) {
                // K drive
                evol_alphakd_(n_K,n_alpha,np1_K,np1_alpha,x,&Nx,&ck,&epsal,&ht,&myzero,phys_bdy,np1_alpha);
            }

            if (slicing_c == 4) {
                evol_alphadyopl_(n_K,n_alpha,np1_K,np1_alpha,x,&Nx,&ht,&myzero,phys_bdy,np1_alpha);
                evol_betagd_(n_Lamx,n_beta,np1_Lamx,np1_beta,x,&Nx,&eta,&ht,&mus,&myzero,phys_bdy,np1_beta);
            }
            if (slicing_c == 5) {

                evol_alphaopl_(n_K,n_alpha,np1_K,np1_alpha,x,&Nx,&gma,&gmb,&ht,&myzero,phys_bdy,np1_alpha);

                evol_betagd_(n_Lamx,n_beta,np1_Lamx,np1_beta,x,&Nx,&eta,&ht,&mus,&myzero,phys_bdy,np1_beta);

            }

            if (slicing_c == 6) {

                evol_alphadyopl_(n_K,n_alpha,np1_K,np1_alpha,x,&Nx,&ht,&myzero,phys_bdy,np1_alpha);

                evol_betakd_(n_K,n_beta,np1_K,np1_beta,x,&Nx,&eta,&ht,&hx,&mus,&myzero,phys_bdy,np1_beta);


            }

            if (slicing_c == 7) {
                evol_alphaopl_(n_K,n_alpha,np1_K,np1_alpha,x,&Nx,&gma,&gmb,&ht,&myzero,phys_bdy,np1_alpha);

                evol_betagd2_(n_Lamx,n_beta,np1_Lamx,np1_beta,x,&Nx,&ck,&epsal,&ht,&myzero,phys_bdy,np1_beta);

            }
            if (slicing_c == 8) {

                evol_beta_hgd_adv_(ctfmp,n_BB,n_beta,np1_BB,np1_beta,x,&Nx,&advc,&ht,&hx,&mus,&myzero,phys_bdy,np1_beta);

                evol_bb_hgd_adv_(ctfmp,n_BB,n_DLamx,n_beta,np1_BB,np1_DLamx,np1_beta,x,&Nx,&advc,&eta,&ht,&hx,&myzero,phys_bdy,np1_BB);

                evol_alpha_hgd_adv_(ctfmp,n_K,n_alpha,n_beta,np1_K,np1_alpha,np1_beta,x,&Nx,&advc,&ht,&hx,&myzero,phys_bdy,np1_alpha);

                eval_dlam_(ctfm,ctfmp,np1_A,np1_Athth,np1_Axx,np1_B,np1_K,np1_Lamx,np1_Sx,np1_alpha,np1_beta,np1_divbeta,np1_phi,octfmp,x,&Nx,&hx,&myzero,&vee,phys_bdy,np1_DLamx);

            }

            compute_db_(ctfm,ctfmp,np1_A,np1_B,np1_beta,x,&Nx,&hx,&myzero,phys_bdy,np1_divbeta);

            compute_emph_(np1_phi,x,&Nx,&myzero,phys_bdy,np1_em4phi);

            compute_ux_(np1_A,np1_Axx,np1_K,x,&Nx,&myzero,phys_bdy,np1_Uxx);

            compute_utheta_(np1_Athth,np1_B,np1_K,x,&Nx,&myzero,phys_bdy,np1_Uthetatheta);

            compute_px_(ctfmp,np1_Axx,np1_beta,np1_divbeta,x,&Nx,&hx,&myzero,&vee,phys_bdy,np1_Pxx);

            compute_ptheta_(ctfm,ctfmp,np1_Athth,np1_beta,np1_divbeta,x,&Nx,&hx,&myzero,&vee,phys_bdy,np1_Pthetatheta);

            compute_metx_(np1_A,np1_phi,x,&Nx,&myzero,phys_bdy,np1_metxx);

            compute_mettheta_(np1_B,np1_phi,x,&Nx,&myzero,phys_bdy,np1_metthetatheta);

            compute_ddlx_(ctfmp,np1_alpha,np1_metxx,octfmp,x,&Nx,&hx,&myzero,phys_bdy,np1_DDLxx);

            compute_ddltheta_(ctfm,ctfmp,np1_alpha,np1_metthetatheta,np1_metxx,x,&Nx,&hx,&myzero,phys_bdy,np1_DDLthetatheta);

            compute_rrx_(ctfm,ctfmp,np1_A,np1_B,np1_Lamx,np1_phi,octfmp,x,&Nx,&hx,&myzero,phys_bdy,np1_RRxx);

            compute_rrtheta_(ctfm,ctfmp,np1_A,np1_B,np1_Lamx,np1_phi,octfmp,x,&Nx,&hx,&myzero,phys_bdy,np1_RRthetatheta);


            compute_a2_(np1_A,np1_psi,&Nx,np1_a2);

            compute_a1_(np1_a2,&Nx,np1_a1);

            compute_b1_(np1_B,np1_psi,&Nx,np1_b1);


            gettimeofday(&wstop_time,NULL);

            timersub(&wstop_time,&wstart_time,&welapsed_time);
            wtime += welapsed_time.tv_sec*1000000+welapsed_time.tv_usec;


            gettimeofday(&estart_time,NULL);

            if (freeze_matter == 0 ) {

                //	evolve_eqipsi_(n_a1,n_alpha,n_beta,n_iPHI,n_iPI,n_ipsi,np1_a1,np1_alpha,np1_beta,np1_iPHI,np1_iPI,np1_ipsi,x,&Nx,&ht,&myzero,&zepsdis,phys_bdy,np1_ipsi);

                evolve_eqrpsi_(n_a1,n_alpha,n_beta,n_rPHI,n_rPI,n_rpsi,np1_a1,np1_alpha,np1_beta,np1_rPHI,np1_rPI,np1_rpsi,x,&Nx,&ht,&myzero,&zepsdis,phys_bdy,np1_rpsi);

                //	evolve_eqipib2_(ctfm,ctfmp,n_a1,n_alpha,n_b1,n_beta,n_iPIb2,n_ipsi,np1_a1,np1_alpha,np1_b1,np1_beta,np1_iPIb2,np1_ipsi,octfmp,x,&Nx,&ht,&hx,&mass,&myzero,&zepsdis,phys_bdy,np1_iPIb2);

                evolve_eqrpib2_(ctfm,ctfmp,n_a1,n_alpha,n_b1,n_beta,n_rPIb2,n_rpsi,np1_a1,np1_alpha,np1_b1,np1_beta,np1_rPIb2,np1_rpsi,octfmp,x,&Nx,&ht,&hx,&mass,&myzero,&zepsdis,phys_bdy,np1_rPIb2);

                /*
                   rdvcpy(np1_ipsi,nph_ipsi,Nx);
                   rdvcpy(np1_rpsi,nph_rpsi,Nx);
                   rdvcpy(np1_iPIb2,nph_iPIb2,Nx);
                   rdvcpy(np1_rPIb2,nph_rPIb2,Nx);
                   */

            }


            gettimeofday(&estop_time,NULL);

            timersub(&estop_time,&estart_time,&eelapsed_time);
            etime += eelapsed_time.tv_sec*1000000+eelapsed_time.tv_usec;


            gettimeofday(&wstart_time,NULL);

            compute_ipi_(np1_b1,np1_iPIb2,&Nx,np1_iPI);

            compute_rpi_(np1_b1,np1_rPIb2,&Nx,np1_rPI);


            compute_dipsi_(ctfm,ctfmp,np1_ipsi,x,&Nx,&hx,&myzero,phys_bdy,np1_iPHI);

            compute_drpsi_(ctfm,ctfmp,np1_rpsi,x,&Nx,&hx,&myzero,phys_bdy,np1_rPHI);

            compute_phi2_(np1_iPHI,np1_rPHI,&Nx,np1_PHI2);

            compute_pi2_(np1_iPI,np1_rPI,&Nx,np1_PI2);


            compute_u_(np1_ipsi,np1_rpsi,&Nx,&mass,np1_U);


            // compute new energy-momentum here


            compute_rho_(np1_PHI2,np1_PI2,np1_U,np1_a2,&Nx,np1_rho);



            compute_jsxx_(np1_PHI2,np1_PI2,np1_U,np1_a2,&Nx,np1_JSxx);


            compute_jsthth_(np1_PHI2,np1_PI2,np1_U,np1_a2,np1_b1,&Nx,np1_JSthth);


            compute_ts_(np1_PHI2,np1_PI2,np1_U,np1_a2,&Nx,np1_TS);


            compute_sx_(np1_a1,np1_iPHI,np1_iPI,np1_rPHI,np1_rPI,&Nx,np1_Sx);


            gettimeofday(&wstop_time,NULL);

            timersub(&wstop_time,&wstart_time,&welapsed_time);
            wtime += welapsed_time.tv_sec*1000000+welapsed_time.tv_usec;




            gettimeofday(&rstart_time,NULL);

            tot_res = 0.0;
            if ( flat_spacetime ==0 ) {

                tres = 0.0;
                //resid_dgtheta_(ctfm,ctfmp,n_Athth,n_B,n_alpha,n_beta,n_divbeta,np1_Athth,np1_B,np1_alpha,np1_beta,np1_divbeta,x,&Nx,&ht,&hx,&myzero,&vee,&zepsdis,phys_bdy,&tres);

                //if (norm_B > evol_tol) { tres = tres/norm_B; }
                //tot_res += tres;

                resid_dgx_(ctfmp,n_A,n_Axx,n_alpha,n_beta,n_divbeta,np1_A,np1_Axx,np1_alpha,np1_beta,np1_divbeta,x,&Nx,&ht,&hx,&myzero,&vee,&zepsdis,phys_bdy,&tres);

                if (norm_A > evol_tol) { tres = tres/norm_A; }
                tot_res += tres;

                //resid_daatheta_(n_Athth,n_DDLthetatheta,n_DDLxx,n_JSthth,n_JSxx,n_Pthetatheta,n_RRthetatheta,n_RRxx,n_Uthetatheta,n_alpha,n_em4phi,n_metthetatheta,n_metxx,np1_Athth,np1_DDLthetatheta,np1_DDLxx,np1_JSthth,np1_JSxx,np1_Pthetatheta,np1_RRthetatheta,np1_RRxx,np1_Uthetatheta,np1_alpha,np1_em4phi,np1_metthetatheta,np1_metxx,x,&Nx,&ht,&hx,&myzero,&zepsdis,phys_bdy,&tres);
                //if (norm_Athth > evol_tol) { tres = tres/norm_Athth; }
                //tot_res += tres;

                resid_daax_(n_Axx,n_DDLthetatheta,n_DDLxx,n_JSthth,n_JSxx,n_Pxx,n_RRthetatheta,n_RRxx,n_Uxx,n_alpha,n_em4phi,n_metthetatheta,n_metxx,np1_Axx,np1_DDLthetatheta,np1_DDLxx,np1_JSthth,np1_JSxx,np1_Pxx,np1_RRthetatheta,np1_RRxx,np1_Uxx,np1_alpha,np1_em4phi,np1_metthetatheta,np1_metxx,x,&Nx,&ht,&hx,&myzero,&zepsdis,phys_bdy,&tres);
                if (norm_Axx > evol_tol) { tres = tres/norm_Axx; }
                tot_res += tres;

                resid_dlamx_(ctfm,ctfmp,n_A,n_Athth,n_Axx,n_B,n_K,n_Lamx,n_Sx,n_alpha,n_beta,n_divbeta,n_phi,np1_A,np1_Athth,np1_Axx,np1_B,np1_K,np1_Lamx,np1_Sx,np1_alpha,np1_beta,np1_divbeta,np1_phi,octfmp,x,&Nx,&ht,&hx,&myzero,&vee,&zepsdis,phys_bdy,&tres);


                if (norm_Lamx > evol_tol) { tres = tres/norm_Lamx; }
                tot_res += tres;

                resid_dph_(ctfmp,n_K,n_alpha,n_beta,n_divbeta,n_phi,np1_K,np1_alpha,np1_beta,np1_divbeta,np1_phi,x,&Nx,&ht,&hx,&myzero,&zepsdis,phys_bdy,&tres);

                if (norm_phi > evol_tol) { tres = tres/norm_phi; }
                tot_res += tres;

                resid_dtk_(ctfmp,n_A,n_Athth,n_Axx,n_B,n_DDLthetatheta,n_DDLxx,n_K,n_TS,n_alpha,n_beta,n_metthetatheta,n_metxx,n_rho,np1_A,np1_Athth,np1_Axx,np1_B,np1_DDLthetatheta,np1_DDLxx,np1_K,np1_TS,np1_alpha,np1_beta,np1_metthetatheta,np1_metxx,np1_rho,x,&Nx,&ht,&hx,&myzero,&zepsdis,phys_bdy,&tres);
                if (norm_K > evol_tol) { tres = tres/norm_K; }
                tot_res += tres;

                if (slicing_c == 1) {
                    resid_alphadyopl_(n_K,n_alpha,np1_K,np1_alpha,x,&Nx,&ht,&myzero,phys_bdy,&tres);
                    if (norm_alpha > evol_tol) { tres = tres / norm_alpha; }
                    tot_res += tres;

                }

                if (slicing_c == 2) {
                    resid_alphaopl_(n_K,n_alpha,np1_K,np1_alpha,x,&Nx,&gma,&gmb,&ht,&myzero,phys_bdy,&tres);
                    if (norm_alpha > evol_tol) { tres = tres / norm_alpha; }
                    tot_res += tres;
                }

                if (slicing_c == 3) {
                    resid_alphakd_(n_K,n_alpha,np1_K,np1_alpha,x,&Nx,&ck,&epsal,&ht,&myzero,phys_bdy,&tres);
                    if (norm_alpha > evol_tol) { tres = tres / norm_alpha; }
                    tot_res += tres;
                }
                if (slicing_c == 4) {
                    resid_alphadyopl_(n_K,n_alpha,np1_K,np1_alpha,x,&Nx,&ht,&myzero,phys_bdy,&tres);
                    if (norm_alpha > evol_tol) { tres = tres / norm_alpha; }
                    tot_res += tres;

                    resid_betagd_(n_Lamx,n_beta,np1_Lamx,np1_beta,x,&Nx,&eta,&ht,&mus,&myzero,phys_bdy,&tres);

                    if (norm_beta > evol_tol) { tres = tres/ norm_beta ;}
                    tot_res += tres;

                }
                if (slicing_c == 5) {
                    resid_alphaopl_(n_K,n_alpha,np1_K,np1_alpha,x,&Nx,&gma,&gmb,&ht,&myzero,phys_bdy,&tres);
                    if (norm_alpha > evol_tol) { tres = tres / norm_alpha; }
                    tot_res += tres;


                    resid_betagd_(n_Lamx,n_beta,np1_Lamx,np1_beta,x,&Nx,&eta,&ht,&mus,&myzero,phys_bdy,&tres);
                    if (norm_beta > evol_tol) { tres = tres/ norm_beta ;}
                    tot_res += tres;


                }

                if (slicing_c == 6 ) {
                    resid_alphadyopl_(n_K,n_alpha,np1_K,np1_alpha,x,&Nx,&ht,&myzero,phys_bdy,&tres);
                    if (norm_alpha > evol_tol) { tres = tres / norm_alpha; }
                    tot_res += tres;

                    resid_betakd_(n_K,n_beta,np1_K,np1_beta,x,&Nx,&eta,&ht,&hx,&mus,&myzero,phys_bdy,&tres);

                    if (norm_beta > evol_tol) { tres = tres/ norm_beta ;}
                    tot_res += tres;

                }

                if (slicing_c == 7) {
                    resid_alphaopl_(n_K,n_alpha,np1_K,np1_alpha,x,&Nx,&gma,&gmb,&ht,&myzero,phys_bdy,&tres);
                    if (norm_alpha > evol_tol) { tres = tres / norm_alpha; }
                    tot_res += tres;

                    resid_betagd2_(n_Lamx,n_beta,np1_Lamx,np1_beta,x,&Nx,&ck,&epsal,&ht,&myzero,phys_bdy,&tres);
                    if (norm_beta > evol_tol) { tres = tres/ norm_beta ;}
                    tot_res += tres;


                }
                if (slicing_c == 8) {

                    resid_bb_hgd_adv_(ctfmp,n_BB,n_DLamx,n_beta,np1_BB,np1_DLamx,np1_beta,x,&Nx,&advc,&eta,&ht,&hx,&myzero,phys_bdy,&tres);
                    if (norm_BB > evol_tol) { tres = tres / norm_BB; }
                    tot_res += tres;

                    resid_alpha_hgd_adv_(ctfmp,n_K,n_alpha,n_beta,np1_K,np1_alpha,np1_beta,x,&Nx,&advc,&ht,&hx,&myzero,phys_bdy,&tres);
                    if (norm_alpha > evol_tol) { tres = tres / norm_alpha; }
                    tot_res += tres;

                    resid_beta_hgd_adv_(ctfmp,n_BB,n_beta,np1_BB,np1_beta,x,&Nx,&advc,&ht,&hx,&mus,&myzero,phys_bdy,&tres);
                    if (norm_beta > evol_tol) { tres = tres / norm_beta; }
                    tot_res += tres;
                }



            }


            if (freeze_matter == 0 ) {

                //resid_eqipib2_(ctfm,ctfmp,n_a1,n_alpha,n_b1,n_beta,n_iPIb2,n_ipsi,np1_a1,np1_alpha,np1_b1,np1_beta,np1_iPIb2,np1_ipsi,octfmp,x,&Nx,&ht,&hx,&mass,&myzero,&zepsdis,phys_bdy,&tres);
                //if (norm_iPIb2 > evol_tol) { tres = tres/norm_iPIb2; }
                //tot_res += tres;

                //resid_eqipsi_(n_a1,n_alpha,n_beta,n_iPHI,n_iPI,n_ipsi,np1_a1,np1_alpha,np1_beta,np1_iPHI,np1_iPI,np1_ipsi,x,&Nx,&ht,&myzero,&zepsdis,phys_bdy,&tres);
                //if (norm_ipsi > evol_tol) { tres = tres/norm_ipsi; }
                //tot_res += tres;

                resid_eqrpib2_(ctfm,ctfmp,n_a1,n_alpha,n_b1,n_beta,n_rPIb2,n_rpsi,np1_a1,np1_alpha,np1_b1,np1_beta,np1_rPIb2,np1_rpsi,octfmp,x,&Nx,&ht,&hx,&mass,&myzero,&zepsdis,phys_bdy,&tres);
                if (norm_rPIb2 > evol_tol) { tres = tres/norm_rPIb2; }
                tot_res += tres;

                resid_eqrpsi_(n_a1,n_alpha,n_beta,n_rPHI,n_rPI,n_rpsi,np1_a1,np1_alpha,np1_beta,np1_rPHI,np1_rPI,np1_rpsi,x,&Nx,&ht,&myzero,&zepsdis,phys_bdy,&tres);
                if (norm_rpsi > evol_tol) { tres = tres/norm_rpsi; }
                tot_res += tres;

            }
            gettimeofday(&rstop_time,NULL);

            timersub(&rstop_time,&rstart_time,&relapsed_time);
            rtime += relapsed_time.tv_sec*1000000+relapsed_time.tv_usec;

            if (j > 3000) {
                printf("iteration did not converge\n");
                printf("res = %1.14e\n",tot_res);
                printf("step = %d\n",i/(int)pow(2.0,(double)level));
                printf("norm_A = %1.14e\n",norm_A);
                printf("norm_B = %1.14e\n",norm_B);
                printf("norm_phi = %1.14e\n",norm_phi);
                printf("norm_K = %1.14e\n",norm_K);
                printf("norm_Axx = %1.14e\n",norm_Axx);
                printf("norm_Athth = %1.14e\n",norm_Athth);
                printf("norm_rpsi = %1.14e\n",norm_rpsi);
                printf("norm_ipsi = %1.14e\n",norm_ipsi);
                printf("norm_iPIb2 = %1.14e\n",norm_iPIb2);
                printf("norm_rPIb2 = %1.14e\n",norm_rPIb2);
                printf("norm_alpha = %1.14e\n",norm_alpha);
                gft_close_all();
                exit(1);
            }
        }
        // End of one step evolution 
        if ( compute_ire == 1) {

            ire_hs_(ctfm,ctfmp,np1_A,np1_Athth,np1_Axx,np1_B,np1_K,np1_Lamx,np1_psi,np1_rho,octfmp,&Nx,&hx,n_reshs);
            ire_val_lamx_(ctfm,ctfmp,np1_A,np1_B,np1_Lamx,&Nx,&hx,res_lamx_val);
            ire_mk_(ctfm,ctfmp,np1_A,np1_Athth,np1_Axx,np1_B,np1_K,np1_Sx,np1_psi,x,&Nx,&hx,&myzero,phys_bdy,n_mk);


            if ( (i > 3) ) {
                // Compute IREs 
                ire_resthth_(ctfm,ctfmp,n_U,n_a1,n_alpha,n_b1,n_beta,n_ipsi,n_rpsi,nm1_a1,nm1_alpha,nm1_b1,nm1_beta,nm1_ipsi,nm1_rpsi,np1_a1,np1_alpha,np1_b1,np1_beta,np1_ipsi,np1_rpsi,octfmp,&Nx,&ht,&hx,n_resthth);
                ire_restt_(ctfm,ctfmp,n_U,n_a1,n_alpha,n_b1,n_beta,n_ipsi,n_rpsi,nm1_a1,nm1_b1,nm1_ipsi,nm1_rpsi,np1_a1,np1_b1,np1_ipsi,np1_rpsi,octfmp,&Nx,&ht,&hx,n_restt);
                ire_restx_(ctfm,ctfmp,n_a1,n_alpha,n_b1,n_beta,n_ipsi,n_rpsi,nm1_a1,nm1_b1,nm1_ipsi,nm1_rpsi,np1_a1,np1_b1,np1_ipsi,np1_rpsi,octfmp,&Nx,&ht,&hx,n_restx);
                ire_resxx_(ctfm,ctfmp,n_U,n_a1,n_alpha,n_b1,n_beta,n_ipsi,n_rpsi,nm1_a1,nm1_alpha,nm1_b1,nm1_beta,nm1_ipsi,nm1_rpsi,np1_a1,np1_alpha,np1_b1,np1_beta,np1_ipsi,np1_rpsi,&Nx,&ht,&hx,n_resxx);
                ire_rpsi_direct_(ctfm,ctfmp,n_a1,n_alpha,n_b1,n_beta,n_rpsi,nm1_a1,nm1_alpha,nm1_b1,nm1_beta,nm1_rpsi,np1_a1,np1_alpha,np1_b1,np1_beta,np1_rpsi,octfmp,&Nx,&ht,&hx,&mass,res_rpsi_direct);

                tot_ire = 0.0;
                if (flat_spacetime == 0 ) {
                    //ire_daatheta_(n_DDLthetatheta,n_DDLxx,n_JSthth,n_JSxx,n_Pthetatheta,n_RRthetatheta,n_RRxx,n_Uthetatheta,n_alpha,n_em4phi,n_metthetatheta,n_metxx,nm1_Athth,np1_Athth,&Nx,&ht,&tres);

                    //tot_ire += tres;

                    //ire_daax_(n_DDLthetatheta,n_DDLxx,n_JSthth,n_JSxx,n_Pxx,n_RRthetatheta,n_RRxx,n_Uxx,n_alpha,n_em4phi,n_metthetatheta,n_metxx,nm1_Axx,np1_Axx,&Nx,&ht,&tres);
                    //tot_ire += tres;

                    //ire_dlamx_(ctfm,ctfmp,n_A,n_Athth,n_Axx,n_B,n_K,n_Lamx,n_Sx,n_alpha,n_beta,n_divbeta,n_phi,nm1_Lamx,np1_Lamx,octfmp,&Nx,&ht,&hx,&vee,&tres);
                    //tot_ire += tres;

                    //ire_dph_(ctfmp,n_K,n_alpha,n_beta,n_divbeta,n_phi,nm1_phi,np1_phi,&Nx,&ht,&hx,&tres);
                    //tot_ire += tres;

                    //ire_dtk_(ctfmp,n_A,n_Athth,n_Axx,n_B,n_DDLthetatheta,n_DDLxx,n_K,n_TS,n_alpha,n_beta,n_metthetatheta,n_metxx,n_rho,nm1_K,np1_K,&Nx,&ht,&hx,&tres);
                    //tot_ire += tres;

                    //ire_dgtheta_(ctfm,ctfmp,n_Athth,n_B,n_alpha,n_beta,n_divbeta,nm1_B,np1_B,&Nx,&ht,&hx,&vee,&tres);
                    //tot_ire += tres;

                    //ire_dgx_(ctfmp,n_A,n_Axx,n_alpha,n_beta,n_divbeta,nm1_A,np1_A,&Nx,&ht,&hx,&vee,&tres);
                    //tot_ire += tres;


                }

                //ire_eqipib2_(ctfm,ctfmp,n_a1,n_alpha,n_b1,n_beta,n_iPIb2,n_ipsi,nm1_iPIb2,np1_iPIb2,octfmp,&Nx,&ht,&hx,&mass,&tres);
                //tot_ire += tres;

                //ire_eqipsi_(n_a1,n_alpha,n_beta,n_iPHI,n_iPI,nm1_ipsi,np1_ipsi,&Nx,&ht,&tres);
                //tot_ire += tres;

                //ire_eqrpib2_(ctfm,ctfmp,n_a1,n_alpha,n_b1,n_beta,n_rPIb2,n_rpsi,nm1_rPIb2,np1_rPIb2,octfmp,&Nx,&ht,&hx,&mass,&tres);
                //tot_ire += tres;

                //ire_eqrpsi_(n_a1,n_alpha,n_beta,n_rPHI,n_rPI,nm1_rpsi,np1_rpsi,&Nx,&ht,&tres);
                //tot_ire += tres;


                //irev_daatheta_(n_DDLthetatheta,n_DDLxx,n_JSthth,n_JSxx,n_Pthetatheta,n_RRthetatheta,n_RRxx,n_Uthetatheta,n_alpha,n_em4phi,n_metthetatheta,n_metxx,nm1_Athth,np1_Athth,&Nx,&ht,n_ire_daath);

                //irev_daax_(n_DDLthetatheta,n_DDLxx,n_JSthth,n_JSxx,n_Pxx,n_RRthetatheta,n_RRxx,n_Uxx,n_alpha,n_em4phi,n_metthetatheta,n_metxx,nm1_Axx,np1_Axx,&Nx,&ht,n_ire_daax);

                //irev_dlamx_(ctfm,ctfmp,n_A,n_Athth,n_Axx,n_B,n_K,n_Lamx,n_Sx,n_alpha,n_beta,n_divbeta,n_phi,nm1_Lamx,np1_Lamx,octfmp,&Nx,&ht,&hx,&vee,n_ire_dlamx);

                //irev_dph_(ctfmp,n_K,n_alpha,n_beta,n_divbeta,n_phi,nm1_phi,np1_phi,&Nx,&ht,&hx,n_ire_dph);

                //irev_dtk_(ctfmp,n_A,n_Athth,n_Axx,n_B,n_DDLthetatheta,n_DDLxx,n_K,n_TS,n_alpha,n_beta,n_metthetatheta,n_metxx,n_rho,nm1_K,np1_K,&Nx,&ht,&hx,n_ire_dtk);

                //irev_dgtheta_(ctfm,ctfmp,n_Athth,n_B,n_alpha,n_beta,n_divbeta,nm1_B,np1_B,&Nx,&ht,&hx,&vee,n_ire_dgth);

                //irev_dgx_(ctfmp,n_A,n_Axx,n_alpha,n_beta,n_divbeta,nm1_A,np1_A,&Nx,&ht,&hx,&vee,n_ire_dgx);

            }
        }

        time = time + ht;


        compute_tt_(np1_PHI2,np1_PI2,np1_a2,&Nx,np1_tt);

        if (mabs(np1_tt[0]) > max_tt) {
            max_tt = mabs(np1_tt[0]);
        }

        compute_th_(ctfm,ctfmp,n_a1,n_alpha,n_b1,n_beta,np1_a1,np1_alpha,np1_b1,np1_beta,x,&Nx,&ht,&hx,&myzero,phys_bdy,nph_th);

        find_min(nph_th,Nx,&min_theta);

        compute_mass_(np1_A,np1_B,np1_psi,ctfm,ctfmp,np1_rho,np1_b1,n_a1,n_mass,n_mass2,n_tmr,&maxtmr,&bhmass,&hx,&Nx);

        compute_dmdr_(ctfm,n_rho,&Nx,dmdr);

        tau = tau + ht*(n_alpha[0] + np1_alpha[0])/2.0;


        find_inner_radious(n_mass2,Nx,&in_rd);

        if ( !(wave_reached_org) ) {
            if ( maxtmr > 0.35 ) {
                wave_reached_org = true;
            }
        }

        if ((i + 1) % output_freq_phi  == 0) {
            fprintf(phi_0_h,"%1.14e %1.14e %1.14e\n",tau,(np1_rpsi[0]+n_rpsi[0])/2.0,np1_tt[0]);
        }

        if (check_bh == 1) {
            if (min_theta < trap_threshold ) {
                if ( in_rd > 0.0 ) {
                    printf("False Blackhole detected\n");
                    fclose(mf);
                    fclose(phi_0_h);
                    gft_close_all();
                    exit(1);
                } else {
                    printf("Blackhole detected\n");
                    bh_h=fopen("BH","w");
                    fprintf(bh_h,"1");
                    printf("step: %d time: %f in_rd: %f min_theta: %f tmr: %f iter: %d res: %1.14e bhmass: %1.14e\n",i+1,time,in_rd,min_theta,maxtmr,j,tot_res,bhmass);
                    fprintf(bh_h,"step: %d time: %f in_rd: %f min_theta: %f tmr: %f iter: %d res: %1.14e bhmass: %1.14e\n",i+1,time,in_rd,min_theta,maxtmr,j,tot_res,bhmass);
                    gft_out_bbox("mass2.sdf",time,shape,dim,bbox,n_mass2);
                    gft_out_bbox("rpsi.sdf",time,shape,dim,bbox,np1_rpsi);
                    gft_out_bbox("alpha.sdf",time,shape,dim,bbox,np1_alpha);
                    gft_out_bbox("A.sdf",time,shape,dim,bbox,np1_A);
                    gft_out_bbox("B.sdf",time,shape,dim,bbox,np1_B);
                    gft_out_bbox("phi.sdf",time,shape,dim,bbox,np1_phi);
                    gft_out_bbox("K.sdf",time,shape,dim,bbox,np1_K);
                    gft_out_bbox("Lamx.sdf",time,shape,dim,bbox,np1_Lamx);
                    gft_out_bbox("Axx.sdf",time,shape,dim,bbox,np1_Axx);
                    gft_out_bbox("Athth.sdf",time,shape,dim,bbox,np1_Athth);
                    gft_out_bbox("tmr.sdf",time,shape,dim,bbox,n_tmr);



                    fprintf(bh_h,"%1.14e\n",max_tt);
                    fclose(bh_h);
                    fclose(mf);
                    fclose(phi_0_h);
                    gft_close_all();
                    exit(1);
                }
            }
        }

        if (check_dis == 1 ) {

            if ( in_rd > ( (x_max-i_in_rd)/10.0 + i_in_rd ) ) {
                printf("dispersal detected\n");
                dis_h=fopen("DISPERSAL","w");
                printf("step: %d time: %f in_rd: %f min_theta: %f tmr: %f iter: %d res: %1.14e bhmass: %1.14e\n",i+1,time,in_rd,min_theta,maxtmr,j,tot_res,bhmass);
                fprintf(dis_h,"step: %d time: %f in_rd: %f min_theta: %f tmr: %f iter: %d res: %1.14e bhmass: %1.14e\n",i+1,time,in_rd,min_theta,maxtmr,j,tot_res,bhmass);
                gft_out_bbox("mass2.sdf",time,shape,dim,bbox,n_mass2);
                gft_out_bbox("rpsi.sdf",time,shape,dim,bbox,np1_rpsi);
                gft_out_bbox("alpha.sdf",time,shape,dim,bbox,np1_alpha);
                gft_out_bbox("A.sdf",time,shape,dim,bbox,np1_A);
                gft_out_bbox("B.sdf",time,shape,dim,bbox,np1_B);
                gft_out_bbox("phi.sdf",time,shape,dim,bbox,np1_phi);
                gft_out_bbox("K.sdf",time,shape,dim,bbox,np1_K);
                gft_out_bbox("Lamx.sdf",time,shape,dim,bbox,np1_Lamx);
                gft_out_bbox("Axx.sdf",time,shape,dim,bbox,np1_Axx);
                gft_out_bbox("Athth.sdf",time,shape,dim,bbox,np1_Athth);
                gft_out_bbox("tmr.sdf",time,shape,dim,bbox,n_tmr);


                fprintf(dis_h,"%1.14e\n",max_tt);
                fclose(dis_h);
                fclose(mf);
                fclose(phi_0_h);
                gft_close_all();
                exit(1);
            }

            if ( wave_reached_org ) {
                if ( in_rd > i_in_rd  ) {
                    printf("dispersal detected\n");
                    dis_h=fopen("DISPERSAL","w");
                    printf("step: %d time: %f in_rd: %f min_theta: %f tmr: %f iter: %d res: %1.14e bhmass: %1.14e\n",i+1,time,in_rd,min_theta,maxtmr,j,tot_res,bhmass);
                    fprintf(dis_h,"step: %d time: %f in_rd: %f min_theta: %f tmr: %f iter: %d res: %1.14e bhmass: %1.14e\n",i+1,time,in_rd,min_theta,maxtmr,j,tot_res,bhmass);
                    gft_out_bbox("mass2.sdf",time,shape,dim,bbox,n_mass2);
                    gft_out_bbox("rpsi.sdf",time,shape,dim,bbox,np1_rpsi);
                    gft_out_bbox("alpha.sdf",time,shape,dim,bbox,np1_alpha);
                    gft_out_bbox("A.sdf",time,shape,dim,bbox,np1_A);
                    gft_out_bbox("B.sdf",time,shape,dim,bbox,np1_B);
                    gft_out_bbox("phi.sdf",time,shape,dim,bbox,np1_phi);
                    gft_out_bbox("K.sdf",time,shape,dim,bbox,np1_K);
                    gft_out_bbox("Lamx.sdf",time,shape,dim,bbox,np1_Lamx);
                    gft_out_bbox("Axx.sdf",time,shape,dim,bbox,np1_Axx);
                    gft_out_bbox("Athth.sdf",time,shape,dim,bbox,np1_Athth);
                    gft_out_bbox("tmr.sdf",time,shape,dim,bbox,n_tmr);


                    fprintf(dis_h,"%1.14e\n",max_tt);
                    fclose(dis_h);
                    fclose(mf);
                    fclose(phi_0_h);
                    gft_close_all();
                    exit(1);
                }

                if (maxtmr < dis_threshold ) {
                    printf("dispersal detected\n");
                    dis_h=fopen("DISPERSAL","w");
                    printf("step: %d time: %f in_rd: %f min_theta: %f tmr: %f iter: %d res: %1.14e bhmass: %1.14e\n",i+1,time,in_rd,min_theta,maxtmr,j,tot_res,bhmass);
                    fprintf(dis_h,"step: %d time: %f in_rd: %f min_theta: %f tmr: %f iter: %d res: %1.14e bhmass: %1.14e\n",i+1,time,in_rd,min_theta,maxtmr,j,tot_res,bhmass);
                    gft_out_bbox("mass2.sdf",time,shape,dim,bbox,n_mass2);
                    gft_out_bbox("rpsi.sdf",time,shape,dim,bbox,np1_rpsi);
                    gft_out_bbox("alpha.sdf",time,shape,dim,bbox,np1_alpha);
                    gft_out_bbox("A.sdf",time,shape,dim,bbox,np1_A);
                    gft_out_bbox("B.sdf",time,shape,dim,bbox,np1_B);
                    gft_out_bbox("phi.sdf",time,shape,dim,bbox,np1_phi);
                    gft_out_bbox("K.sdf",time,shape,dim,bbox,np1_K);
                    gft_out_bbox("Lamx.sdf",time,shape,dim,bbox,np1_Lamx);
                    gft_out_bbox("Axx.sdf",time,shape,dim,bbox,np1_Axx);
                    gft_out_bbox("Athth.sdf",time,shape,dim,bbox,np1_Athth);
                    gft_out_bbox("tmr.sdf",time,shape,dim,bbox,n_tmr);


                    fprintf(dis_h,"%1.14e\n",max_tt);
                    fclose(dis_h);
                    fclose(mf);
                    fclose(phi_0_h);
                    gft_close_all();
                    exit(1);

                }

            }
        }

        if ((i + 1) % (output_freq*(int)pow(2.0,(double)level))  == 0) {

            gft_out_bbox("dmdr.sdf",time,shape,dim,bbox,dmdr);
            fprintf(mf,"%1.14e %1.14e\n",time, n_mass[Nx-2]);
            gft_out_bbox("A.sdf",time,shape,dim,bbox,np1_A);
            gft_out_bbox("B.sdf",time,shape,dim,bbox,np1_B);
            gft_out_bbox("phi.sdf",time,shape,dim,bbox,np1_phi);
            gft_out_bbox("K.sdf",time,shape,dim,bbox,np1_K);
            gft_out_bbox("Lamx.sdf",time,shape,dim,bbox,np1_Lamx);
            gft_out_bbox("Axx.sdf",time,shape,dim,bbox,np1_Axx);
            gft_out_bbox("Athth.sdf",time,shape,dim,bbox,np1_Athth);
            gft_out_bbox("tmr.sdf",time,shape,dim,bbox,n_tmr);

            gft_out_bbox("alpha.sdf",time,shape,dim,bbox,np1_alpha);
            gft_out_bbox("beta.sdf",time,shape,dim,bbox,np1_beta);
            gft_out_bbox("ipsi.sdf",time,shape,dim,bbox,np1_ipsi);
            gft_out_bbox("rpsi.sdf",time,shape,dim,bbox,np1_rpsi);
            gft_out_bbox("a1.sdf",time,shape,dim,bbox,np1_a1);
            gft_out_bbox("b1.sdf",time,shape,dim,bbox,np1_b1);
            gft_out_bbox("rho.sdf",time,shape,dim,bbox,np1_rho);
            gft_out_bbox("RRxx.sdf",time,shape,dim,bbox,np1_RRxx);
            gft_out_bbox("mass.sdf",time,shape,dim,bbox,n_mass);
            gft_out_bbox("mass2.sdf",time,shape,dim,bbox,n_mass2);
            gft_out_bbox("Sx.sdf",time,shape,dim,bbox,np1_Sx);
            gft_out_bbox("psi.sdf",time,shape,dim,bbox,np1_psi);
            gft_out_bbox("TS.sdf",time,shape,dim,bbox,np1_TS);
            gft_out_bbox("JSxx.sdf",time,shape,dim,bbox,np1_JSxx);
            gft_out_bbox("JSthth.sdf",time,shape,dim,bbox,np1_JSthth);

            gft_out_bbox("tt.sdf",time,shape,dim,bbox,np1_tt);

            if (compute_ire == 1) {
                gft_out_bbox("restx.sdf",time,shape,dim,bbox,n_restx);
                gft_out_bbox("resxx.sdf",time,shape,dim,bbox,n_resxx);
                gft_out_bbox("restt.sdf",time,shape,dim,bbox,n_restt);
                gft_out_bbox("resthth.sdf",time,shape,dim,bbox,n_resthth);
                gft_out_bbox("reshs.sdf",time,shape,dim,bbox,n_reshs);
                gft_out_bbox("resrpsi.sdf",time,shape,dim,bbox,res_rpsi_direct);
                gft_out_bbox("resmk.sdf",time,shape,dim,bbox,n_mk);
                gft_out_bbox("reslamx.sdf",time,shape,dim,bbox,res_lamx_val);
                gft_out_bbox("iredgx.sdf",time,shape,dim,bbox,n_ire_dgx);
                gft_out_bbox("iredgth.sdf",time,shape,dim,bbox,n_ire_dgth);
                gft_out_bbox("iredph.sdf",time,shape,dim,bbox,n_ire_dph);
                gft_out_bbox("iredtk.sdf",time,shape,dim,bbox,n_ire_dtk);
                gft_out_bbox("iredlamx.sdf",time,shape,dim,bbox,n_ire_dlamx);
                gft_out_bbox("iredaax.sdf",time,shape,dim,bbox,n_ire_daax);
                gft_out_bbox("iredaath.sdf",time,shape,dim,bbox,n_ire_daath);
            }



            /*gft_out_bbox("em4phi.sdf",time,shape,dim,bbox,np1_em4phi);
              gft_out_bbox("metxx.sdf",time,shape,dim,bbox,np1_metxx);
              gft_out_bbox("metthetatheta.sdf",time,shape,dim,bbox,np1_metthetatheta);
              gft_out_bbox("RRthetatheta.sdf",time,shape,dim,bbox,np1_RRthetatheta);
              gft_out_bbox("Pxx.sdf",time,shape,dim,bbox,np1_Pxx);
              gft_out_bbox("Pthetatheta.sdf",time,shape,dim,bbox,np1_Pthetatheta);
              gft_out_bbox("Uxx.sdf",time,shape,dim,bbox,np1_Uxx);
              gft_out_bbox("Uthetatheta.sdf",time,shape,dim,bbox,np1_Uthetatheta);
              gft_out_bbox("DDLxx.sdf",time,shape,dim,bbox,np1_DDLxx);
              gft_out_bbox("DDLthetatheta.sdf",time,shape,dim,bbox,np1_DDLthetatheta);
              gft_out_bbox("U.sdf",time,shape,dim,bbox,np1_U);
              gft_out_bbox("iPHI.sdf",time,shape,dim,bbox,np1_iPHI);
              gft_out_bbox("rPHI.sdf",time,shape,dim,bbox,np1_rPHI);

              gft_out_bbox("iPI.sdf",time,shape,dim,bbox,np1_iPI);
              gft_out_bbox("rPI.sdf",time,shape,dim,bbox,np1_rPI);


              gft_out_bbox("iPIb2.sdf",time,shape,dim,bbox,np1_iPIb2);
              gft_out_bbox("rPIb2.sdf",time,shape,dim,bbox,np1_rPIb2);
              */


            dtphi_center = mabs(np1_rpsi[0] - n_rpsi[0] ) / ht;
            //printf("step: %d time: %f in_rd: %f min_theta: %f tmr: %f iter: %d res: %1.14e ire: %1.14e\n",i+1,time,in_rd,min_theta,maxtmr,j,tot_res,tot_ire);
            printf("step: %d time: %f in_rd: %f min_theta: %f tmr: %f iter: %d res: %1.14e dtphc: %1.14e\n bhmass: %1.14e\n",i+1,time,in_rd,min_theta,maxtmr,j,tot_res,dtphi_center,bhmass);
        }

        if (compute_ire == 1) {
            rdvcpy(nm1_A,n_A,Nx);
            rdvcpy(nm1_B,n_B,Nx);
            rdvcpy(nm1_phi,n_phi,Nx);
            rdvcpy(nm1_K,n_K,Nx);
            rdvcpy(nm1_Axx,n_Axx,Nx);
            rdvcpy(nm1_Athth,n_Athth,Nx);
            rdvcpy(nm1_Lamx,n_Lamx,Nx);
            rdvcpy(nm1_ipsi,n_ipsi,Nx);
            rdvcpy(nm1_rpsi,n_rpsi,Nx);
            rdvcpy(nm1_iPIb2,n_iPIb2,Nx);
            rdvcpy(nm1_rPIb2,n_rPIb2,Nx);
            rdvcpy(nm1_a1,n_a1,Nx);
            rdvcpy(nm1_b1,n_b1,Nx);
            rdvcpy(nm1_alpha,n_alpha,Nx);
            rdvcpy(nm1_beta,n_beta,Nx);
        }

        if (adaptive_time == 1) {
            dtphi_center = mabs( np1_rpsi[0] - n_rpsi[0] ) / ht;
            if (dtphi_center <= 0.1) { 
                ht = ht0; 
                if (first_adt_warn == 0) {
                    printf("Adaptive time mesh disabled\n");
                    first_adt_warn = 1;
                }
            }
            if (dtphi_center > 0.1) { 
                ht = ht0 / (10.0*dtphi_center); 
                if (first_adt_warn == 1) {
                    printf("Adaptive time enabeled: ht0/ht = %f\n",ht0/ht);
                    first_adt_warn = 0;
                }
            } 
        }

        swap_levels(&np1_alpha,&n_alpha);
        swap_levels(&np1_beta,&n_beta);


        swap_levels(&np1_A,&n_A);
        swap_levels(&np1_B,&n_B);
        swap_levels(&np1_phi,&n_phi);
        swap_levels(&np1_K,&n_K);
        swap_levels(&np1_Axx,&n_Axx);
        swap_levels(&np1_Athth,&n_Athth);
        swap_levels(&np1_Lamx,&n_Lamx);

        swap_levels(&np1_divbeta,&n_divbeta);

        swap_levels(&np1_em4phi,&n_em4phi);

        swap_levels(&np1_metxx,&n_metxx);
        swap_levels(&np1_metthetatheta,&n_metthetatheta);

        swap_levels(&np1_Uxx,&n_Uxx);
        swap_levels(&np1_Uthetatheta,&n_Uthetatheta);

        swap_levels(&np1_RRxx,&n_RRxx);
        swap_levels(&np1_RRthetatheta,&n_RRthetatheta);

        swap_levels(&np1_Pxx,&n_Pxx);
        swap_levels(&np1_Pthetatheta,&n_Pthetatheta);

        swap_levels(&np1_DDLxx,&n_DDLxx);
        swap_levels(&np1_DDLthetatheta,&n_DDLthetatheta);


        swap_levels(&np1_a1,&n_a1);
        swap_levels(&np1_a2,&n_a2);

        swap_levels(&np1_b1,&n_b1);
        swap_levels(&np1_psi,&n_psi);

        swap_levels(&np1_ipsi,&n_ipsi);
        swap_levels(&np1_rpsi,&n_rpsi);
        swap_levels(&np1_iPIb2,&n_iPIb2);
        swap_levels(&np1_rPIb2,&n_rPIb2);
        swap_levels(&np1_iPI,&n_iPI);
        swap_levels(&np1_rPI,&n_rPI);
        swap_levels(&np1_iPHI,&n_iPHI);
        swap_levels(&np1_rPHI,&n_rPHI);

        swap_levels(&np1_PHI2,&n_PHI2);
        swap_levels(&np1_PI2,&n_PI2);

        swap_levels(&np1_U,&n_U);

        swap_levels(&np1_rho,&n_rho);
        swap_levels(&np1_TS,&n_TS);
        swap_levels(&np1_Sx,&n_Sx);
        swap_levels(&np1_JSxx,&n_JSxx);
        swap_levels(&np1_JSthth,&n_JSthth);

    }
    gettimeofday(&stop_time,NULL);
    timersub(&stop_time,&start_time,&elapsed_time);
    printf("Total time was %f seconds.\n",elapsed_time.tv_sec+
            elapsed_time.tv_usec/1000000.0);

    printf("Evolve time = %f\n",etime/1000000.0);
    printf("Resid time = %f\n",rtime/1000000.0);
    printf("Work time = %f\n",wtime/1000000.0);
    gft_close_all();

}

