%% Whole cell model
% parameters
thetar= 426.8693338968694;
s0= 1.0e4;
gmax= 1260.0;
thetax= 4.379733394834643;
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
parameters = [thetar s0 gmax thetax Kt M we Km vm nx Kq Kp vt wr wq nq nr ns];

% define rate constants
dm= 0.1;
kb= 1;
ku= 1.0;
rates= [dm kb ku];

% define initial conditions
rmr_0= 0;
em_0= 0;
rmq_0= 0;
rmt_0= 0;
et_0= 0;
rmm_0= 0;
mt_0= 0;
mm_0= 0;
q_0= 0;
si_0= 0;
mq_0= 0;
mr_0= 0;
r_0= 10.0;
a_0= 1000.0;
init= [rmr_0 em_0 rmq_0 rmt_0 et_0 rmm_0 mt_0 mm_0 q_0 si_0 mq_0 mr_0 r_0 a_0];


%% Synthetic gene circuit

% parameters and rate constants
wp= 0.0;
% add parameters into the parameters array (includes parameters of whole cell model)
parameters = [parameters wp];

% initial amounts of species:
mp_0= 0;
rmp_0= 0;
p_0= 0;
% add initial amounts of species into the init array (includes initial amounts of species of whole cell model)
init = [init mp_0 rmp_0 p_0];


%% call solver routine 
t0= 0;
tf= 1e9;
[t,y]= ode15s(@(t,y) cellmodel_odes_nk01(t, y, rates, parameters), [t0 tf], init);

% whole cell model molecular species
rmr= y(:,1);    % mrna-ribo complex of ribosomal proteins
em= y(:,2);     % protein metabolic enzyme
rmq= y(:,3);    % mrna-ribo complex of housekeeping proteins
rmt= y(:,4);    % mrna-ribo complex of transporter enzyme proteins
et= y(:,5);     % protein transporter enzyme
rmm= y(:,6);    % mrna-ribo complex of metabolic enzyme proteins
mt= y(:,7);     % mrna of transporter enzyme
mm= y(:,8);     % mrna of metabolic enzyme
q= y(:,9);      % protein housekeeping
si= y(:,10);    % internal nutrients
mq= y(:,11);    % mrna of house-keeping protein
mr= y(:,12);    % mrna of ribosomal proteins
r= y(:,13);     % free ribosomes
a= y(:,14);     % energy

% synthetic gene circuits molecular species 
mp= y(:,15);
rmp= y(:,16);
p= y(:,17);

