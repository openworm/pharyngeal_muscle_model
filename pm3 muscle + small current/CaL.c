/* Created by Language version: 6.2.0 */
/* VECTORIZED */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib_ansi.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define _threadargscomma_ _p, _ppvar, _thread, _nt,
#define _threadargs_ _p, _ppvar, _thread, _nt
 
#define _threadargsprotocomma_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt,
#define _threadargsproto_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 /* Thread safe. No static _p or _ppvar. */
 
#define t _nt->_t
#define dt _nt->_dt
#define gmax _p[0]
#define oo _p[1]
#define c1 _p[2]
#define c2 _p[3]
#define c3 _p[4]
#define c4 _p[5]
#define i1 _p[6]
#define i2 _p[7]
#define i3 _p[8]
#define eca _p[9]
#define ica _p[10]
#define alpha _p[11]
#define beta _p[12]
#define gamma _p[13]
#define Doo _p[14]
#define Dc1 _p[15]
#define Dc2 _p[16]
#define Dc3 _p[17]
#define Dc4 _p[18]
#define Di1 _p[19]
#define Di2 _p[20]
#define Di3 _p[21]
#define v _p[22]
#define _g _p[23]
#define _ion_eca	*_ppvar[0]._pval
#define _ion_ica	*_ppvar[1]._pval
#define _ion_dicadv	*_ppvar[2]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 static Datum* _extcall_thread;
 static Prop* _extcall_prop;
 /* external NEURON variables */
 /* declaration of user functions */
 static void _hoc_settables(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _extcall_prop = _prop;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_CaL", _hoc_setdata,
 "settables_CaL", _hoc_settables,
 0, 0
};
 
static void _check_settables(double*, Datum*, Datum*, _NrnThread*); 
static void _check_table_thread(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, int _type) {
   _check_settables(_p, _ppvar, _thread, _nt);
 }
 /* declare global and static user variables */
#define usetable usetable_CaL
 double usetable = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "usetable_CaL", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "gmax_CaL", "mho/cm2",
 0,0
};
 static double c40 = 0;
 static double c30 = 0;
 static double c20 = 0;
 static double c10 = 0;
 static double delta_t = 0.01;
 static double i30 = 0;
 static double i20 = 0;
 static double i10 = 0;
 static double oo0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "usetable_CaL", &usetable_CaL,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 static void nrn_cur(_NrnThread*, _Memb_list*, int);
static void  nrn_jacob(_NrnThread*, _Memb_list*, int);
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(_NrnThread*, _Memb_list*, int);
static void _ode_matsol(_NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[3]._i
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "6.2.0",
"CaL",
 "gmax_CaL",
 0,
 0,
 "oo_CaL",
 "c1_CaL",
 "c2_CaL",
 "c3_CaL",
 "c4_CaL",
 "i1_CaL",
 "i2_CaL",
 "i3_CaL",
 0,
 0};
 static Symbol* _ca_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 24, _prop);
 	/*initialize range parameters*/
 	gmax = 1.65e-006;
 	_prop->param = _p;
 	_prop->param_size = 24;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 4, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_ca_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* eca */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ica */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dicadv */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*f)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _CaL_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("ca", -10000.);
 	_ca_sym = hoc_lookup("ca_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 1);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
     _nrn_thread_table_reg(_mechtype, _check_table_thread);
  hoc_register_prop_size(_mechtype, 24, 4);
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 CaL C:/nrn/work_muscle_model_v4/CaL.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double *_t_alpha;
 static double *_t_beta;
 static double *_t_gamma;
static int _reset;
static char *modelname = "Calcium L channel";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int _f_settables(_threadargsprotocomma_ double);
static int settables(_threadargsprotocomma_ double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static void _n_settables(_threadargsprotocomma_ double _lv);
 static int _slist1[8], _dlist1[8];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset = 0; {
   settables ( _threadargscomma_ v ) ;
   Doo = alpha * c4 - 4.0 * beta * oo + 0.0005 * i1 - gamma * oo + 0.001 * ( alpha * i2 - 13.0 * oo ) ;
   Dc1 = - 4.0 * alpha * c1 + beta * c2 ;
   Dc2 = 4.0 * alpha * c1 - beta * c2 + 2.0 * beta * c3 - 3.0 * alpha * c2 ;
   Dc3 = 3.0 * alpha * c2 - 2.0 * beta * c3 + 3.0 * beta * c4 - 2.0 * alpha * c3 ;
   Dc4 = 2.0 * alpha * c3 - 3.0 * beta * c4 + 4.0 * beta * oo - alpha * c4 + 0.01 * ( 4.0 * 0.0005 * beta * i1 - alpha * gamma * c4 ) + 0.002 * ( 4.0 * beta * i2 - 13.0 * c4 ) + 4.0 * beta * 0.0005 * i3 - gamma * 13.0 * c4 ;
   Di1 = gamma * oo - 0.0005 * i1 + 0.001 * ( alpha * i3 - 13.0 * i1 ) + 0.01 * ( alpha * gamma * c4 - 4.0 * beta * 0.0005 * i1 ) ;
   Di2 = 0.001 * ( 13.0 * oo - alpha * i2 ) + 0.0005 * i3 - gamma * i2 + 0.002 * ( 13.0 * c4 - 4.0 * beta * i2 ) ;
   Di3 = 0.001 * ( 13.0 * i1 - alpha * i3 ) + gamma * i2 - 0.0005 * i3 + gamma * 13.0 * c4 - 4.0 * beta * 0.0005 * i3 ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
 settables ( _threadargscomma_ v ) ;
 Doo = Doo  / (1. - dt*( ( - (4.0 * beta)*(1.0) ) - (gamma)*(1.0) + (0.001)*(( ( - (13.0)*(1.0) ) )) )) ;
 Dc1 = Dc1  / (1. - dt*( (- 4.0 * alpha)*(1.0) )) ;
 Dc2 = Dc2  / (1. - dt*( ( - (beta)*(1.0) ) - (3.0 * alpha)*(1.0) )) ;
 Dc3 = Dc3  / (1. - dt*( ( - (2.0 * beta)*(1.0) ) - (2.0 * alpha)*(1.0) )) ;
 Dc4 = Dc4  / (1. - dt*( ( - (3.0 * beta)*(1.0) ) - (alpha)*(1.0) + (0.01)*(( ( - (alpha * gamma)*(1.0) ) )) + (0.002)*(( ( - (13.0)*(1.0) ) )) - (gamma * 13.0)*(1.0) )) ;
 Di1 = Di1  / (1. - dt*( ( - (0.0005)*(1.0) ) + (0.001)*(( ( - (13.0)*(1.0) ) )) + (0.01)*(( ( - (4.0 * beta * 0.0005)*(1.0) ) )) )) ;
 Di2 = Di2  / (1. - dt*( (0.001)*(( ( - (alpha)*(1.0) ) )) - (gamma)*(1.0) + (0.002)*(( ( - (4.0 * beta)*(1.0) ) )) )) ;
 Di3 = Di3  / (1. - dt*( (0.001)*(( ( - (alpha)*(1.0) ) )) - (0.0005)*(1.0) - (4.0 * beta * 0.0005)*(1.0) )) ;
 return 0;
}
 /*END CVODE*/
 static int states (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) { {
   settables ( _threadargscomma_ v ) ;
    oo = oo + (1. - exp(dt*(( - (4.0 * beta)*(1.0) ) - (gamma)*(1.0) + (0.001)*(( ( - (13.0)*(1.0) ) )))))*(- ( (alpha)*(c4) + (0.0005)*(i1) + (0.001)*(( (alpha)*(i2) )) ) / ( ( - ((4.0)*(beta))*(1.0)) - (gamma)*(1.0) + (0.001)*(( ( - (13.0)*(1.0)) )) ) - oo) ;
    c1 = c1 + (1. - exp(dt*((- 4.0 * alpha)*(1.0))))*(- ( (beta)*(c2) ) / ( ((- 4.0)*(alpha))*(1.0) ) - c1) ;
    c2 = c2 + (1. - exp(dt*(( - (beta)*(1.0) ) - (3.0 * alpha)*(1.0))))*(- ( ((4.0)*(alpha))*(c1) + ((2.0)*(beta))*(c3) ) / ( ( - (beta)*(1.0)) - ((3.0)*(alpha))*(1.0) ) - c2) ;
    c3 = c3 + (1. - exp(dt*(( - (2.0 * beta)*(1.0) ) - (2.0 * alpha)*(1.0))))*(- ( ((3.0)*(alpha))*(c2) + ((3.0)*(beta))*(c4) ) / ( ( - ((2.0)*(beta))*(1.0)) - ((2.0)*(alpha))*(1.0) ) - c3) ;
    c4 = c4 + (1. - exp(dt*(( - (3.0 * beta)*(1.0) ) - (alpha)*(1.0) + (0.01)*(( ( - (alpha * gamma)*(1.0) ) )) + (0.002)*(( ( - (13.0)*(1.0) ) )) - (gamma * 13.0)*(1.0))))*(- ( ((2.0)*(alpha))*(c3) + ((4.0)*(beta))*(oo) + (0.01)*(( (((4.0)*(0.0005))*(beta))*(i1) )) + (0.002)*(( ((4.0)*(beta))*(i2) )) + (((4.0)*(beta))*(0.0005))*(i3) ) / ( ( - ((3.0)*(beta))*(1.0)) - (alpha)*(1.0) + (0.01)*(( ( - ((alpha)*(gamma))*(1.0)) )) + (0.002)*(( ( - (13.0)*(1.0)) )) - ((gamma)*(13.0))*(1.0) ) - c4) ;
    i1 = i1 + (1. - exp(dt*(( - (0.0005)*(1.0) ) + (0.001)*(( ( - (13.0)*(1.0) ) )) + (0.01)*(( ( - (4.0 * beta * 0.0005)*(1.0) ) )))))*(- ( (gamma)*(oo) + (0.001)*(( (alpha)*(i3) )) + (0.01)*(( ((alpha)*(gamma))*(c4) )) ) / ( ( - (0.0005)*(1.0)) + (0.001)*(( ( - (13.0)*(1.0)) )) + (0.01)*(( ( - (((4.0)*(beta))*(0.0005))*(1.0)) )) ) - i1) ;
    i2 = i2 + (1. - exp(dt*((0.001)*(( ( - (alpha)*(1.0) ) )) - (gamma)*(1.0) + (0.002)*(( ( - (4.0 * beta)*(1.0) ) )))))*(- ( (0.001)*(( (13.0)*(oo) )) + (0.0005)*(i3) + (0.002)*(( (13.0)*(c4) )) ) / ( (0.001)*(( ( - (alpha)*(1.0)) )) - (gamma)*(1.0) + (0.002)*(( ( - ((4.0)*(beta))*(1.0)) )) ) - i2) ;
    i3 = i3 + (1. - exp(dt*((0.001)*(( ( - (alpha)*(1.0) ) )) - (0.0005)*(1.0) - (4.0 * beta * 0.0005)*(1.0))))*(- ( (0.001)*(( (13.0)*(i1) )) + (gamma)*(i2) + ((gamma)*(13.0))*(c4) ) / ( (0.001)*(( ( - (alpha)*(1.0)) )) - (0.0005)*(1.0) - (((4.0)*(beta))*(0.0005))*(1.0) ) - i3) ;
   }
  return 0;
}
 static double _mfac_settables, _tmin_settables;
  static void _check_settables(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  if (!usetable) {return;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_settables =  - 100.0 ;
   _tmax =  100.0 ;
   _dx = (_tmax - _tmin_settables)/200.; _mfac_settables = 1./_dx;
   for (_i=0, _x=_tmin_settables; _i < 201; _x += _dx, _i++) {
    _f_settables(_p, _ppvar, _thread, _nt, _x);
    _t_alpha[_i] = alpha;
    _t_beta[_i] = beta;
    _t_gamma[_i] = gamma;
   }
  }
 }

 static int settables(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _lv) { 
#if 0
_check_settables(_p, _ppvar, _thread, _nt);
#endif
 _n_settables(_p, _ppvar, _thread, _nt, _lv);
 return 0;
 }

 static void _n_settables(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 _f_settables(_p, _ppvar, _thread, _nt, _lv); return; 
}
 _xi = _mfac_settables * (_lv - _tmin_settables);
 if (isnan(_xi)) {
  alpha = _xi;
  beta = _xi;
  gamma = _xi;
  return;
 }
 if (_xi <= 0.) {
 alpha = _t_alpha[0];
 beta = _t_beta[0];
 gamma = _t_gamma[0];
 return; }
 if (_xi >= 200.) {
 alpha = _t_alpha[200];
 beta = _t_beta[200];
 gamma = _t_gamma[200];
 return; }
 _i = (int) _xi;
 _theta = _xi - (double)_i;
 alpha = _t_alpha[_i] + _theta*(_t_alpha[_i+1] - _t_alpha[_i]);
 beta = _t_beta[_i] + _theta*(_t_beta[_i+1] - _t_beta[_i]);
 gamma = _t_gamma[_i] + _theta*(_t_gamma[_i+1] - _t_gamma[_i]);
 }

 
static int  _f_settables ( _threadargsprotocomma_ double _lv ) {
   alpha = 0.40 * exp ( ( _lv + 20.0 + 12.0 ) / 10.0 ) ;
   beta = 0.05 * exp ( - ( _lv - 20.0 + 12.0 ) / 13.0 ) ;
   gamma = 0.11662 * 0.136058 / ( 10.0 + 0.136058 ) ;
    return 0; }
 
static void _hoc_settables(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 
#if 1
 _check_settables(_p, _ppvar, _thread, _nt);
#endif
 _r = 1.;
 settables ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 8;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  eca = _ion_eca;
     _ode_spec1 (_p, _ppvar, _thread, _nt);
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 8; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  eca = _ion_eca;
 _ode_matsol1 (_p, _ppvar, _thread, _nt);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_ca_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 2, 4);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  int _i; double _save;{
  c4 = c40;
  c3 = c30;
  c2 = c20;
  c1 = c10;
  i3 = i30;
  i2 = i20;
  i1 = i10;
  oo = oo0;
 {
   settables ( _threadargscomma_ v ) ;
   oo = 0.0 ;
   c1 = 1.0 ;
   c2 = 0.0 ;
   c3 = 0.0 ;
   c4 = 0.0 ;
   i1 = 0.0 ;
   i2 = 0.0 ;
   i3 = 0.0 ;
   }
 
}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];

#if 0
 _check_settables(_p, _ppvar, _thread, _nt);
#endif
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
  eca = _ion_eca;
 initmodel(_p, _ppvar, _thread, _nt);
 }}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   ica = gmax * oo * ( v - eca ) ;
   }
 _current += ica;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
  eca = _ion_eca;
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ double _dica;
  _dica = ica;
 _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
  _ion_dicadv += (_dica - ica)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ica += ica ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type) {
 double _break, _save;
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 _break = t + .5*dt; _save = t;
 v=_v;
{
  eca = _ion_eca;
 { {
 for (; t < _break; t += dt) {
   states(_p, _ppvar, _thread, _nt);
  
}}
 t = _save;
 } }}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(oo) - _p;  _dlist1[0] = &(Doo) - _p;
 _slist1[1] = &(c1) - _p;  _dlist1[1] = &(Dc1) - _p;
 _slist1[2] = &(c2) - _p;  _dlist1[2] = &(Dc2) - _p;
 _slist1[3] = &(c3) - _p;  _dlist1[3] = &(Dc3) - _p;
 _slist1[4] = &(c4) - _p;  _dlist1[4] = &(Dc4) - _p;
 _slist1[5] = &(i1) - _p;  _dlist1[5] = &(Di1) - _p;
 _slist1[6] = &(i2) - _p;  _dlist1[6] = &(Di2) - _p;
 _slist1[7] = &(i3) - _p;  _dlist1[7] = &(Di3) - _p;
   _t_alpha = makevector(201*sizeof(double));
   _t_beta = makevector(201*sizeof(double));
   _t_gamma = makevector(201*sizeof(double));
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif
