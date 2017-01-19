TITLE K+ HERG channel - MGWMN model

UNITS {
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

NEURON {
	SUFFIX KvEXP2
	USEION k READ ek WRITE ik
	RANGE gmax
}

PARAMETER { 
	gmax = 0.0001 (S/cm2) 
	p1 = 1.0
	p2 = 1.0
	p3 = 1.0
	p4 = 1.0
	p5 = 1.0
	p6 = 1.0
	p7 = 1.0
	p8 = 1.0
	p9 = 1.0
	p10 = 1.0
	p11 = 1.0
	p12 = 1.0
	p13 = 1.0
	p14 = 1.0
	p15 = 1.0
	p16 = 1.0	
}

ASSIGNED { 
	v (mV)
	ek (mV)
	ik (mA/cm2)
	alpha_0 (/ms)
	beta_0 (/ms)
	k_f (/ms)
	k_b (/ms)
	alpha_1 (/ms)
	beta_1 (/ms)
	alpha_i (/ms)
	beta_i (/ms)
	alpha_i3 (/ms)
	psi (/ms)
}

STATE {
	oo c1 c2 c3 ii
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ik  = gmax*oo*(v-ek)
}

INITIAL {
	settables(v)
	oo = 0
	c1 = 1
	c2 = 0
	c3 = 0
	ii = 0
}

DERIVATIVE states {  
	settables(v)
	c1' = -alpha_0 * c1 + beta_0 * c2
	c2' =  alpha_0 * c1 - beta_0 * c2 - k_f * c2 + k_b * c3
	c3' =  k_f * c2 - k_b * c3 - alpha_1 * c3 + beta_1 * oo + psi * ii - alpha_i3 * c3
	oo' =  alpha_1 * c3 - beta_1 * oo - alpha_i * oo + beta_i * ii
	ii' =  alpha_i * oo - beta_i * ii + alpha_i3 * c3 - psi * ii
}

UNITSOFF

PROCEDURE settables(v (mV)) {
	alpha_0   = 0.0069   * p1  * exp(  0.0272  * p2  * v )
	beta_0    = 0.0227   * p3  * exp( -0.0431  * p4  * v )
	k_f       = 0.0266   * p5
	k_b       = 0.1348   * p6
	alpha_1   = 0.0218   * p7  * exp(  0.0262  * p8  * v )
	beta_1    = 0.0009   * p9  * exp( -0.0269  * p10 * v )
	alpha_i   = 0.0622   * p11 * exp(  0.0120  * p12 * v ) 
	beta_i    = 0.0059   * p13 * exp( -0.0443  * p14 * v )
	alpha_i3  = 1.29e-5  * p15 * exp(  2.71e-6 * p16 * v ) 
	psi       = ( beta_1 * beta_i * alpha_i3 ) / (alpha_1 * alpha_i)
}

UNITSON