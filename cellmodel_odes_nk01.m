
% This file implements the right hand side of the ODE.

%%
function dydt= cellmodel_odes(t, y, rates, parameters)

	dm= rates(1);       % mrna degradation rate [min-1]
	kb= rates(2);       % mRNA-ribosome binding rate [cell/min molecs]
	ku= rates(3);       % mRNA-ribosome unbinding rate [min-1]

    
	thetar= parameters(1);      % ribosome transcription threshold [molecs/cell]
	s0= parameters(2);          % external nutrients [molesc]
	gmax= parameters(3);        % maximal translational elongation rate
	thetax= parameters(4);      % non-ribosomal transcription threshold [molecs/cell]
	Kt= parameters(5);          % nutrients import threshold [molecs]
	M= parameters(6);           % total cell mass [aa]
	we= parameters(7);          % max enzyme transcription rate [molecs/min cell]
	Km= parameters(8);         % enzymatic threshold [molecs/cell]
	vm= parameters(9);         % max enzymatic rate [min-1]
	nx= parameters(10);         % length of non-ribosomal proteins [aa/molecs]
	Kq= parameters(11);         % q-autoinhibition threshold
	Kp= parameters(12);
	vt= parameters(13);         % max nutrient import rate [min-1]
	wr= parameters(14);         % max ribosome transcription rate [molecs/min cell]
	wq= parameters(15);         % max q-transcription rate [molecs/min cell]
	wp= parameters(16);         % max enzyme transcription rate [molecs/min cell]
	nq= parameters(17);         % q-autoinhibition hill coeff. [none]
	nr= parameters(18);         % ribosome length [aa/molecs]
	ns= parameters(19);         % nutrient efficiency

	rmr= y(1);  % mrna-ribo complex of ribosomal proteins
	em= y(2);   % protein metabolic enzyme
	rmq= y(3);  % mrna-ribo complex of housekeeping proteins
	rmt= y(4);  % mrna-ribo complex of transporter enzyme proteins
	et= y(5);   % protein transporter enzyme
	rmm= y(6);  % mrna-ribo complex of metabolic enzyme proteins
	mt= y(7);  % mrna of transporter enzyme
	mm= y(8);  % mrna of metabolic enzyme
	q= y(9);   % protein housekeeping
	si= y(10);  % internal nutrients
	mq= y(11);  % mrna of house-keeping protein
	mr= y(12);  % mrna of ribosomes
	r= y(13);   % ribosomes
	a= y(14);   % ATP
    
    mp= y(15);  % mrna of gratuitous protein
    rmp= y(16);  % mrna-ribo complex of gratuitous protein
	p= y(17);   % gratuitous protein
    
	Kgamma= gmax/Kp;    % translational elongation threshold [moles/cell]   % essential (part equation 3)
	gamma= gmax*a/(Kgamma + a); % rate of translational elongation          % essential (eq. 3)
	ttrate= (rmq + rmr + rmp + rmt + rmm)*gamma;    % total number of ribosomes * rate of translational elongation      % essential (eq.9b)
	lam= ttrate/M;  % growth rate       % essential (eq.9b)
	fr= nr*(r + rmr + rmp + rmt + rmm + rmq) / ( nr*(r + rmr + rmp + rmt + rmm + rmq) + nx * (p + q + et + em)); % ribosome fraction        % equation 10 (essential???) 
	nucat= em*vm*si/(Km + si);  % rate of metabolism of nutrients       % essential (part of eq.1)

	dydt(size(y,1),1)= 0;
	dydt(1)= +kb*r*mr-ku*rmr-gamma/nr*rmr-f*rmr-lam*rmr;  % d(rmr)/dt     % essential (equation 6A) + non-essential terms
	dydt(2)= +gamma/nx*rmm-lam*em;                              % d(em)/dt      
	dydt(3)= +kb*r*mq-ku*rmq-gamma/nx*rmq-f*rmq-lam*rmq;  % d(rmq)/dt     % essential (equation 6B) + non-essential terms
	dydt(4)= +kb*r*mt-ku*rmt-gamma/nx*rmt-f*rmt-lam*rmt;  % d(rmt)/dt     % essential (equation 6C) + non-essential terms
	dydt(5)= +gamma/nx*rmt-lam*et;                              % d(et)/dt
	dydt(6)= +kb*r*mm-ku*rmm-gamma/nx*rmm-f*rmm-lam*rmm;  % d(rmm)/dt     % essential (equation 6D) + non-essential terms
	dydt(7)= +(we*a/(thetax + a))+ku*rmt+gamma/nx*rmt-kb*r*mt-dm*mt-lam*mt;    % d(mt)/dt      % essential (eq.5D & 8)
	dydt(8)= +(we*a/(thetax + a))+ku*rmm+gamma/nx*rmm-kb*r*mm-dm*mm-lam*mm;    % d(mm)/dt      % essential (eq.5C & 8)
	dydt(9)= +gamma/nx*rmq-lam*q;                              % d(q)/dt
	dydt(10)= +(et*vt*s0/(Kt + s0))-nucat-lam*si;               % d(si)/dt  nutrient import into cell % essential (eq1)
	dydt(11)= +(wq*a/(thetax + a)/(1 + (q/Kq)^nq))+ku*rmq+gamma/nx*rmq-kb*r*mq-dm*mq-lam*mq;     % d(mq)/dt  transcription, translation and degratation of q-proteins   % essential (eq.5A & 8)
	dydt(12)= +(wr*a/(thetar + a))+ku*rmr+gamma/nr*rmr-kb*r*mr-dm*mr-lam*mr;     %  d(mr)/dt  transcription, translation and degratation of ribosomal-proteins      % essential (eq.5B & 8)
	dydt(13)= +ku*rmr+ku*rmt+ku*rmm+ku*rmp+ku*rmq+gamma/nr*rmr+gamma/nr*rmr+gamma/nx*rmt+gamma/nx*rmm+gamma/nx*rmp+gamma/nx*rmq-kb*r*mr-kb*r*mt-kb*r*mm-kb*r*mp-kb*r*mq-lam*r; % d(r)/dt    % essential (eq.7)
	dydt(14)= +ns*nucat-ttrate-lam*a;       % d(a)/dt % essential (eq.2)
    
    
    dydt(15)= +(wp*a/(thetax + a))+ku*rmp+gamma/nx*rmp-kb*r*mp-dm*mp-lam*mp;    % d(mp)/dt   transcription, translation and degratation of the mRNA of gratuitous-proteins
    dydt(16)= +kb*r*mp-ku*rmp-gamma/nx*rmp-f*rmp-lam*rmp;                        % d(rmp)/dt
    dydt(17)= +gamma/nx*rmp-lam*p;                                              % d(p)/dt
    
    