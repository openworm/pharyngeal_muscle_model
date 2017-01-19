TITLE Voltage-gated K channel

UNITS {
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

NEURON {
	SUFFIX Kv
	USEION k READ ek WRITE ik
	RANGE gbar, g, i
}

PARAMETER { 
	gbar = 0.036 (S/cm2) 
}

ASSIGNED {
	v (mV)
	ek (mV)
	ik (mA/cm2)
	i (mA/cm2)
	g (S/cm2)
}

STATE { 
	n 
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gbar * n^4
	i = g * (v - ek)
	ik = i
}

INITIAL {
	n = alpha(v)/(alpha(v) + beta(v))
}

DERIVATIVE states {
	n' = (1-n)*alpha(v) - n*beta(v)
}

FUNCTION alpha(Vm (mV)) (/ms) {
	LOCAL x
	UNITSOFF
	x = (Vm+55)/10
	if (fabs(x) > 1e-6) {
		alpha = 0.1*x/(1 - exp(-x))
	} 
	else {
		alpha = 0.1/(1 - 0.5*x)
	} 
	UNITSON
}

FUNCTION beta(Vm (mV)) (/ms) {
	UNITSOFF
	beta = 0.125*exp(-(Vm+65)/80)
	UNITSON
}