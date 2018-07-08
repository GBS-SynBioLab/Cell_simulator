function dydt= cellmodel_odes(t, y)
    %% whole cell model rates, parameters and species assigments
    global a r dm kb ku thetax gamma lam;
    
    % parameters
    thetar = 426.8693338968694;
    s0= 1.0e4;
    gmax= 1260.0;
    thetax = 4.379733394834643;
    Kt= 1.0e3;
    M= 1.0e8;
    we= 4.139172187824451;
    Km= 1.0e3;
    vm= 5800.0;
    nx= 300.0;
    Kq= 1.522190403737490e+05;
    Kp= 180.1378030928276;
    vt= 726.0;
    wr= 929.9678874564831;
    wq= 948.9349882947897;
    nq= 4;
    nr= 7549.0;
    ns= 0.5;


    % define rate constants
    dm= 0.1;
    kb= 1;
    ku= 1.0;

    
    % whole cell model molecular species
    
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
    
    %% synthetic gene circuit parameters & species assignments
    
    % synthetic genetic circuit  parameters
    wp= 1000;
    
    % synthetic gene circuit species assignments
    mp= y(15);  % mrna of gratuitous protein
    rmp= y(16);  % mrna-ribo complex of gratuitous protein
	p = y(17);   % gratuitous protein
    
    
    %% whole cell model odes and fluxes
    
    % fluxes
	Kgamma= gmax/Kp;    % translational elongation threshold [moles/cell]   % essential (part equation 3)
	gamma= gmax*a/(Kgamma + a); % rate of translational elongation          % essential (eq. 3)
	ttrate= (rmq + rmr + rmp + rmt + rmm)*gamma;    % total number of ribosomes * rate of translational elongation      % essential (eq.9b)
	lam= ttrate/M;  % growth rate       % essential (eq.9b)
	fr= nr*(r + rmr + rmp + rmt + rmm + rmq) / ( nr*(r + rmr + rmp + rmt + rmm + rmq) + nx * (p + q + et + em)); % ribosome fraction        % equation 10 (essential???) 
	nucat= em*vm*si/(Km + si);  % rate of metabolism of nutrients       % essential (part of eq.1)

	dydt(size(y,1),1)= 0;
    % odes
	dydt(1)= +kb*r*mr-ku*rmr-gamma/nr*rmr-lam*rmr;  % d(rmr)/dt     % essential (equation 6A) + non-essential terms
	dydt(2)= +gamma/nx*rmm-lam*em;                              % d(em)/dt      
	dydt(3)= +kb*r*mq-ku*rmq-gamma/nx*rmq-lam*rmq;  % d(rmq)/dt     % essential (equation 6B) + non-essential terms
	dydt(4)= +kb*r*mt-ku*rmt-gamma/nx*rmt-lam*rmt;  % d(rmt)/dt     % essential (equation 6C) + non-essential terms
	dydt(5)= +gamma/nx*rmt-lam*et;                              % d(et)/dt
	dydt(6)= +kb*r*mm-ku*rmm-gamma/nx*rmm-lam*rmm;  % d(rmm)/dt     % essential (equation 6D) + non-essential terms
	dydt(7)= +(we*a/(thetax + a))+ku*rmt+gamma/nx*rmt-kb*r*mt-dm*mt-lam*mt;    % d(mt)/dt      % essential (eq.5D & 8)
	dydt(8)= +(we*a/(thetax + a))+ku*rmm+gamma/nx*rmm-kb*r*mm-dm*mm-lam*mm;    % d(mm)/dt      % essential (eq.5C & 8)
	dydt(9)= +gamma/nx*rmq-lam*q;                              % d(q)/dt
	dydt(10)= +(et*vt*s0/(Kt + s0))-nucat-lam*si;               % d(si)/dt  nutrient import into cell % essential (eq1)
	dydt(11)= +(wq*a/(thetax + a)/(1 + (q/Kq)^nq))+ku*rmq+gamma/nx*rmq-kb*r*mq-dm*mq-lam*mq;     % d(mq)/dt  transcription, translation and degratation of q-proteins   % essential (eq.5A & 8)
	dydt(12)= +(wr*a/(thetar + a))+ku*rmr+gamma/nr*rmr-kb*r*mr-dm*mr-lam*mr;     %  d(mr)/dt  transcription, translation and degratation of ribosomal-proteins      % essential (eq.5B & 8)
	dydt(13)= +ku*rmr+ku*rmt+ku*rmm+ku*rmp+ku*rmq+gamma/nr*rmr+gamma/nr*rmr+gamma/nx*rmt+gamma/nx*rmm+gamma/nx*rmp+gamma/nx*rmq-kb*r*mr-kb*r*mt-kb*r*mm-kb*r*mp-kb*r*mq-lam*r; % d(r)/dt    % essential (eq.7)
	dydt(14)= +ns*nucat-ttrate-lam*a;       % d(a)/dt % essential (eq.2)
    
    %% synthetic gene circuit odes
    
    % odes
    dydt(15)= +transcription(wp) + ribosome_unbind(rmp) + translation(rmp,nx) - ribosome_bind(mp) - degradation(mp) - dilution(mp);    % d(mp)/dt  
    dydt(16)= +ribosome_bind(mp) - ribosome_unbind(rmp) - translation(rmp,nx) - dilution(rmp);                             % d(rmp)/dt
    dydt(17)= +translation(rmp,nx) - dilution(p) ;                                                                 % d(p)/dt
    
end