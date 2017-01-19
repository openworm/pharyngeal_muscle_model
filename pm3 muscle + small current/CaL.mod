COMMENT
A model of graded calcium release and L-type 
Ca2+ channel inactivation in cardiac muscle
Am J Physiol Heart Circ Physiol 286: H1154-H1169, 2004;
10.1152/ajpheart.00168.2003.
ENDCOMMENT

TITLE Calcium L channel

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
}

NEURON {
    SUFFIX CaL
    USEION ca READ eca WRITE ica
    RANGE gmax
}

PARAMETER {
    gmax  = 0.00000165 (mho/cm2)
	}

ASSIGNED { 
    v (mV)
    eca (mV)
    ica (mA/cm2)
    alpha (/ms)
    beta (/ms)
    gamma (/ms)
}

STATE {
    oo c1 c2 c3 c4 i1 i2 i3
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    ica  = gmax*oo*(v-eca)
}

INITIAL {
	settables(v)
	oo = 0
	c1 = 1
	c2 = 0
	c3 = 0
	c4 = 0
	i1 = 0
	i2 = 0
	i3 = 0
}

DERIVATIVE states {  
	settables(v)
	oo'  = alpha*c4 - 4*beta*oo + 0.0005*i1 - gamma*oo + 0.001*(alpha*i2 - 13.0*oo)
	c1' = -4*alpha*c1 + beta*c2
	c2' = 4*alpha*c1 - beta*c2 + 2*beta*c3 - 3*alpha*c2
	c3' = 3*alpha*c2 - 2*beta*c3 + 3*beta*c4 - 2*alpha*c3
	c4' = 2*alpha*c3 - 3*beta*c4 + 4*beta*oo - alpha*c4 + 0.01*(4*0.0005*beta*i1 - alpha*gamma*c4) + 0.002*(4*beta*i2 - 13.0*c4) + 4*beta*0.0005*i3 - gamma*13.0*c4
	i1' = gamma*oo - 0.0005*i1 + 0.001*(alpha*i3 - 13.0*i1) + 0.01*(alpha*gamma*c4 - 4*beta*0.0005*i1)
	i2' = 0.001*(13.0*oo - alpha*i2) + 0.0005*i3 - gamma*i2 + 0.002*(13.0*c4 - 4*beta*i2)
	i3' = 0.001*(13.0*i1 - alpha*i3) + gamma*i2 - 0.0005*i3 + gamma*13.0*c4 - 4*beta*0.0005*i3
}

UNITSOFF

PROCEDURE settables(v (mV)) {
    TABLE alpha, beta, gamma
          FROM -100 TO 100 WITH 200

    alpha = 0.40 * exp( (v+20.0+12.0)/10.0)
    beta  = 0.05 * exp(-(v-20.0+12.0)/13.0)
    gamma = 0.11662 * 0.136058 /(10.0+0.136058)
}

UNITSON
