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
#define ii _p[5]
#define ek _p[6]
#define ik _p[7]
#define alpha_0 _p[8]
#define beta_0 _p[9]
#define k_f _p[10]
#define k_b _p[11]
#define alpha_1 _p[12]
#define beta_1 _p[13]
#define alpha_i _p[14]
#define beta_i _p[15]
#define alpha_i3 _p[16]
#define psi _p[17]
#define Doo _p[18]
#define Dc1 _p[19]
#define Dc2 _p[20]
#define Dc3 _p[21]
#define Dii _p[22]
#define v _p[23]
#define _g _p[24]
#define _ion_ek	*_ppvar[0]._pval
#define _ion_ik	*_ppvar[1]._pval
#define _ion_dikdv	*_ppvar[2]._pval
 
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
 "setdata_KvEXP2", _hoc_setdata,
 "settables_KvEXP2", _hoc_settables,
 0, 0
};
 /* declare global and static user variables */
#define p16 p16_KvEXP2
 double p16 = 1;
#define p15 p15_KvEXP2
 double p15 = 1;
#define p14 p14_KvEXP2
 double p14 = 1;
#define p13 p13_KvEXP2
 double p13 = 1;
#define p12 p12_KvEXP2
 double p12 = 1;
#define p11 p11_KvEXP2
 double p11 = 1;
#define p10 p10_KvEXP2
 double p10 = 1;
#define p9 p9_KvEXP2
 double p9 = 1;
#define p8 p8_KvEXP2
 double p8 = 1;
#define p7 p7_KvEXP2
 double p7 = 1;
#define p6 p6_KvEXP2
 double p6 = 1;
#define p5 p5_KvEXP2
 double p5 = 1;
#define p4 p4_KvEXP2
 double p4 = 1;
#define p3 p3_KvEXP2
 double p3 = 1;
#define p2 p2_KvEXP2
 double p2 = 1;
#define p1 p1_KvEXP2
 double p1 = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "gmax_KvEXP2", "S/cm2",
 0,0
};
 static double c30 = 0;
 static double c20 = 0;
 static double c10 = 0;
 static double delta_t = 0.01;
 static double ii0 = 0;
 static double oo0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "p1_KvEXP2", &p1_KvEXP2,
 "p2_KvEXP2", &p2_KvEXP2,
 "p3_KvEXP2", &p3_KvEXP2,
 "p4_KvEXP2", &p4_KvEXP2,
 "p5_KvEXP2", &p5_KvEXP2,
 "p6_KvEXP2", &p6_KvEXP2,
 "p7_KvEXP2", &p7_KvEXP2,
 "p8_KvEXP2", &p8_KvEXP2,
 "p9_KvEXP2", &p9_KvEXP2,
 "p10_KvEXP2", &p10_KvEXP2,
 "p11_KvEXP2", &p11_KvEXP2,
 "p12_KvEXP2", &p12_KvEXP2,
 "p13_KvEXP2", &p13_KvEXP2,
 "p14_KvEXP2", &p14_KvEXP2,
 "p15_KvEXP2", &p15_KvEXP2,
 "p16_KvEXP2", &p16_KvEXP2,
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
"KvEXP2",
 "gmax_KvEXP2",
 0,
 0,
 "oo_KvEXP2",
 "c1_KvEXP2",
 "c2_KvEXP2",
 "c3_KvEXP2",
 "ii_KvEXP2",
 0,
 0};
 static Symbol* _k_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 25, _prop);
 	/*initialize range parameters*/
 	gmax = 0.0001;
 	_prop->param = _p;
 	_prop->param_size = 25;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 4, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_k_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* ek */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ik */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dikdv */
 
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

 void _kvherg_mgwmn_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("k", -10000.);
 	_k_sym = hoc_lookup("k_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 1);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
  hoc_register_prop_size(_mechtype, 25, 4);
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 KvEXP2 C:/nrn/work_muscle_model_v4/kvherg_mgwmn.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "K+ HERG channel - MGWMN model";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int settables(_threadargsprotocomma_ double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[5], _dlist1[5];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset = 0; {
   settables ( _threadargscomma_ v ) ;
   Dc1 = - alpha_0 * c1 + beta_0 * c2 ;
   Dc2 = alpha_0 * c1 - beta_0 * c2 - k_f * c2 + k_b * c3 ;
   Dc3 = k_f * c2 - k_b * c3 - alpha_1 * c3 + beta_1 * oo + psi * ii - alpha_i3 * c3 ;
   Doo = alpha_1 * c3 - beta_1 * oo - alpha_i * oo + beta_i * ii ;
   Dii = alpha_i * oo - beta_i * ii + alpha_i3 * c3 - psi * ii ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
 settables ( _threadargscomma_ v ) ;
 Dc1 = Dc1  / (1. - dt*( (- alpha_0)*(1.0) )) ;
 Dc2 = Dc2  / (1. - dt*( ( - (beta_0)*(1.0) ) - (k_f)*(1.0) )) ;
 Dc3 = Dc3  / (1. - dt*( ( - (k_b)*(1.0) ) - (alpha_1)*(1.0) - (alpha_i3)*(1.0) )) ;
 Doo = Doo  / (1. - dt*( ( - (beta_1)*(1.0) ) - (alpha_i)*(1.0) )) ;
 Dii = Dii  / (1. - dt*( ( - (beta_i)*(1.0) ) - (psi)*(1.0) )) ;
 return 0;
}
 /*END CVODE*/
 static int states (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) { {
   settables ( _threadargscomma_ v ) ;
    c1 = c1 + (1. - exp(dt*((- alpha_0)*(1.0))))*(- ( (beta_0)*(c2) ) / ( (- alpha_0)*(1.0) ) - c1) ;
    c2 = c2 + (1. - exp(dt*(( - (beta_0)*(1.0) ) - (k_f)*(1.0))))*(- ( (alpha_0)*(c1) + (k_b)*(c3) ) / ( ( - (beta_0)*(1.0)) - (k_f)*(1.0) ) - c2) ;
    c3 = c3 + (1. - exp(dt*(( - (k_b)*(1.0) ) - (alpha_1)*(1.0) - (alpha_i3)*(1.0))))*(- ( (k_f)*(c2) + (beta_1)*(oo) + (psi)*(ii) ) / ( ( - (k_b)*(1.0)) - (alpha_1)*(1.0) - (alpha_i3)*(1.0) ) - c3) ;
    oo = oo + (1. - exp(dt*(( - (beta_1)*(1.0) ) - (alpha_i)*(1.0))))*(- ( (alpha_1)*(c3) + (beta_i)*(ii) ) / ( ( - (beta_1)*(1.0)) - (alpha_i)*(1.0) ) - oo) ;
    ii = ii + (1. - exp(dt*(( - (beta_i)*(1.0) ) - (psi)*(1.0))))*(- ( (alpha_i)*(oo) + (alpha_i3)*(c3) ) / ( ( - (beta_i)*(1.0)) - (psi)*(1.0) ) - ii) ;
   }
  return 0;
}
 
static int  settables ( _threadargsprotocomma_ double _lv ) {
   alpha_0 = 0.0069 * p1 * exp ( 0.0272 * p2 * _lv ) ;
   beta_0 = 0.0227 * p3 * exp ( - 0.0431 * p4 * _lv ) ;
   k_f = 0.0266 * p5 ;
   k_b = 0.1348 * p6 ;
   alpha_1 = 0.0218 * p7 * exp ( 0.0262 * p8 * _lv ) ;
   beta_1 = 0.0009 * p9 * exp ( - 0.0269 * p10 * _lv ) ;
   alpha_i = 0.0622 * p11 * exp ( 0.0120 * p12 * _lv ) ;
   beta_i = 0.0059 * p13 * exp ( - 0.0443 * p14 * _lv ) ;
   alpha_i3 = 1.29e-5 * p15 * exp ( 2.71e-6 * p16 * _lv ) ;
   psi = ( beta_1 * beta_i * alpha_i3 ) / ( alpha_1 * alpha_i ) ;
    return 0; }
 
static void _hoc_settables(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r = 1.;
 settables ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 5;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ek = _ion_ek;
     _ode_spec1 (_p, _ppvar, _thread, _nt);
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 5; ++_i) {
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
  ek = _ion_ek;
 _ode_matsol1 (_p, _ppvar, _thread, _nt);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_k_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_k_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_k_sym, _ppvar, 2, 4);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  int _i; double _save;{
  c3 = c30;
  c2 = c20;
  c1 = c10;
  ii = ii0;
  oo = oo0;
 {
   settables ( _threadargscomma_ v ) ;
   oo = 0.0 ;
   c1 = 1.0 ;
   c2 = 0.0 ;
   c3 = 0.0 ;
   ii = 0.0 ;
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
  ek = _ion_ek;
 initmodel(_p, _ppvar, _thread, _nt);
 }}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   ik = gmax * oo * ( v - ek ) ;
   }
 _current += ik;

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
  ek = _ion_ek;
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ double _dik;
  _dik = ik;
 _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
  _ion_dikdv += (_dik - ik)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ik += ik ;
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
  ek = _ion_ek;
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
 _slist1[0] = &(c1) - _p;  _dlist1[0] = &(Dc1) - _p;
 _slist1[1] = &(c2) - _p;  _dlist1[1] = &(Dc2) - _p;
 _slist1[2] = &(c3) - _p;  _dlist1[2] = &(Dc3) - _p;
 _slist1[3] = &(oo) - _p;  _dlist1[3] = &(Doo) - _p;
 _slist1[4] = &(ii) - _p;  _dlist1[4] = &(Dii) - _p;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif
