
% This file initializes parameters and calls the solver routine.

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
wp= 0.0;
nq= 4;
nr= 7549.0;
ns= 0.5;
parameters= [thetar s0 gmax thetax Kt M we Km vm nx Kq Kp vt wr wq wp nq nr ns];

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

mp_0= 0;
rmp_0= 0;
p_0= 0;

init= [rmr_0 em_0 rmq_0 rmt_0 et_0 rmm_0 mt_0 mm_0 q_0 si_0 mq_0 mr_0 r_0 a_0 mp_0 rmp_0 p_0];


%% call solver routine 
t0= 0;
tf= 1e9;
[t,y]= ode15s(@(t,y) cellmodel_odes(t, y, rates, parameters), [t0 tf], init);
rmr= y(:,1);
em= y(:,2);
rmq= y(:,3);
rmt= y(:,4);
et= y(:,5);
rmm= y(:,6);
mt= y(:,7);
mm= y(:,8);
q= y(:,9);
si= y(:,10);
mq= y(:,11);
mr= y(:,12);
r= y(:,13);
a= y(:,14);

mp= y(:,15);
rmp= y(:,16);
p= y(:,17);
