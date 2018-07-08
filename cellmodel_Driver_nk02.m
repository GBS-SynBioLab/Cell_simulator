%% Whole cell model initial values of species

% define initial conditions
rmr_0; em_0; rmq_0; rmt_0; et_0; rmm_0; mt_0; mm_0; q_0; si_0; mq_0; mr_0 = 0;
r_0= 10.0;
a_0= 1000.0;

init= [rmr_0 em_0 rmq_0 rmt_0 et_0 rmm_0 mt_0 mm_0 q_0 si_0 mq_0 mr_0 r_0 a_0];


%% Synthetic gene circuit initial values of species

% initial amounts of species:
mp_0= 0;
rmp_0= 0;
p_0= 0;
% add initial amounts of species into the init array (includes initial amounts of species of whole cell model)
init = [init mp_0 rmp_0 p_0];


%% call solver routine 
t0= 0;
tf= 1e9;
[t,y]= ode15s(@(t,y) cellmodel_odes_nk02(t, y), [t0 tf], init);

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

%% plotting

close all;

f1 = figure;
subplot(3,3,1);
loglog(t,si);
xlabel('time');
ylabel('internal nutrients');
title('Whole cell model species');
ylim([1e-15 1e5]);
xlim([1e-5 1e510]);
hold on

subplot(3,3,2);
loglog(t,a);
xlabel('time');
ylabel('energy');
title('Whole cell model species');
ylim([1e-15 1e5]);
xlim([1e-5 1e510]);
hold on

subplot(3,3,3);
loglog(t,r);
xlabel('time');
ylabel('free ribosomes');
title('Whole cell model species');
ylim([1e-15 1e5]);
xlim([1e-5 1e510]);
hold on

subplot(3,3,4);
loglog(t,et);
xlabel('time');
ylabel('transporter enzymes');
title('Whole cell model species');
ylim([1e-15 1e5]);
xlim([1e-5 1e510]);
hold on

subplot(3,3,5);
loglog(t,em);
xlabel('time');
ylabel('metabolic enzymes');
title('Whole cell model species');
ylim([1e-15 1e5]);
xlim([1e-5 1e510]);
hold on

subplot(3,3,6);
loglog(t,q);
xlabel('time');
ylabel('housekeping proteins');
title('Whole cell model species');
ylim([1e-15 1e10]);
xlim([1e-5 1e510]);
hold on


subplot(3,3,7);
loglog(t,mm);
xlabel('time');
ylabel('mRNAs of m, t & q');
title('Whole cell model species');
ylim([1e-15 1e5]);
xlim([1e-5 1e510]);
hold on
loglog(t,mt);
hold on
loglog(t,mq);
hold on

subplot(3,3,8);
loglog(t,mr);
xlabel('time');
ylabel('mRNA of ribosomal proteins');
title('Whole cell model species');
ylim([1e-15 1e5]);
xlim([1e-5 1e510]);
hold on

act_tl_r = rmr + rmq + rmt + rmm;
subplot(3,3,9);
loglog(t,act_tl_r);
xlabel('time');
ylabel('translating ribosomes');
title('Whole cell model species');
ylim([1e-15 1e5]);
xlim([1e-5 1e510]);
hold off


f2 = figure;

subplot(3,3,1);
loglog(t,mp);
xlabel('time');
ylabel('mRNA');
title('synthetic circuit species');
ylim([1e-15 1e5]);
xlim([1e-5 1e510]);
hold on

subplot(3,3,2);
loglog(t,rmp);
xlabel('time');
ylabel('ribocomplex');
title('synthetic circuit species');
ylim([1e-15 1e5]);
xlim([1e-5 1e510]);
hold on

subplot(3,3,3);
loglog(t,p);
xlabel('time');
ylabel('protein');
title('synthetic circuit species');
ylim([1e-15 1e10]);
xlim([1e-5 1e510]);
hold off
