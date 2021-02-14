%% Prerequisites for the Question 3

% Load the diffusion signal
fid = fopen('isbi2015_data_normalised.txt', 'r', 'b');
fgetl(fid); % Read in the header
D = fscanf(fid, '%f', [6, inf])'; % Read in the data
fclose(fid);

% Select the first of the 6 voxels
meas = D(:,1);
% Load the protocol
fid = fopen('isbi2015_protocol.txt', 'r', 'b');
fgetl(fid);
A = fscanf(fid, '%f', [7, inf]);
fclose(fid);

% Create the protocol
grad_dirs = A(1:3,:);
G = A(4,:);
delta = A(5,:);
smalldel = A(6,:);
TE = A(7,:);


GAMMA = 2.675987E8;
bvals = ((GAMMA*smalldel.*G).^2).*(delta-smalldel/3);
% convert bvals units from s/m^2 to s/mm^2
 
bvals = bvals/10^6;


qhat = grad_dirs;

%%

% A=Gx we know A and we can work out G to find the missing values for x

G = [ones(1,3612); -bvals.*qhat(1,:).^2; -2*bvals.*qhat(1,:).*qhat(2,:); -2*bvals.*qhat(1,:).*qhat(3,:); -bvals.*qhat(2,:).^2; -2*bvals.*qhat(2,:).*qhat(3,:); -bvals.*qhat(3,:).^2]';

x = G\log(meas);


%% Q1.3.1 Ball and stick model 


% Produce starting matrix

startx = [x(1), x(2), x(3), x(4), x(5)];

h=optimset('MaxFunEvals',20000,'Algorithm','quasi-newton','TolX',1e-10,'TolFun',1e-10, 'Display', 'iter');
oldopts = optimset('TolFun',1e-6);
options = optimset(oldopts,'PlotFcns','optimplotfval','TolX',1e-8);

[parameter_hat,RESNORMq131,EXITFLAG,OUTPUT]=fminunc('q3BallStick_noAvox_Positive',startx,h,meas,bvals,qhat);

Resnorm_q1_3_1 = RESNORMq131;

S0 = parameter_hat(1); % initial signal
d  = parameter_hat(2); % who knows what this is
f  = parameter_hat(3); % this guy is a joke
theta = parameter_hat(4); % I like him
phi = parameter_hat(5); % He does nothing for my life

parameter_hat = [ S0^2 d^2 Siggy(f) theta phi];

% pd stands for PreDicted Model hence the pd

pd_q131 = BallStick_noAvox(parameter_hat, bvals, qhat);

%% Plot the data points of the measurements and the model
figure(1);
plot(meas, ' bs', 'MarkerSize', 2, 'LineWidth', 4); % Actual data 
hold on
plot(pd_q131, ' rx', 'MarkerSize', 2, 'LineWidth', 2) % Model predictions
xlabel('qindex')
ylabel('S')
title('Q1.3.1')


%% Q1.3.2 Two compartment model Nothing Added

% Produce starting matrix

startx = [x(1), x(2), x(3), x(4), x(5)];


[parameter_hat,RESNORMq132i,EXITFLAG,OUTPUT]=fminunc('q3twoCOMPSSD_Positive',startx,h,meas,bvals,qhat);

% Resnorm_q1_3_2i = RESNORMq132i;

S0_q132i = parameter_hat(1); % initial signal
d_q132i  = parameter_hat(2); % who knows what this is
f_q132i  = parameter_hat(3); % this guy is a joke
theta_q132i = parameter_hat(4); % I like him
phi_q132i = parameter_hat(5); % He does nothing for my life

parameter_hat_q132i = [ S0_q132i^2 d_q132i^2 Siggy(f_q132i) theta_q132i phi_q132i];

% pd stands for PreDicted Model hence the pd

pd_q132i = q3twoCOMPSSD_noAvox(parameter_hat, bvals, qhat);

disp(startx);


%% Add Different Starting Points

R = RESNORMq132i;

% Housing all the integers
numIter = 100;
Resnorm_numIter=zeros(1,numIter);
x1_rand_index=zeros(1,numIter);
x2_rand_index=zeros(1,numIter);
x3_rand_index=zeros(1,numIter);
x4_rand_index=zeros(1,numIter);
x5_rand_index=zeros(1,numIter);

counter = zeros(1,numIter);

for i = 1:numIter
    
    x1_rand = normrnd(x(1),0.5);
    x1_rand_index(:,i)=x1_rand;
    x2_rand = normrnd(x(2),0.5e-4);
    x2_rand_index(:,i)=x2_rand;
    x3_rand = normrnd(x(3),0.5e-6);
    x3_rand_index(:,i)=x3_rand;
    x4_rand = normrnd(x(4),0.5e-6);
    x4_rand_index(:,i)=x4_rand;
    x5_rand = normrnd(x(5),0.5e-6);
    x5_rand_index(:,i)=x5_rand;
    
    startx2 = [x1_rand_index(:,i),x2_rand_index(:,i),x3_rand_index(:,i),x4_rand_index(:,i),x5_rand_index(:,i)];
    
   [parameter_hat_q132i,RESNORMq132i,EXITFLAG,OUTPUT]=fminunc('q3twoCOMPSSD_Positive',startx2,h,meas,bvals,qhat);
   
   Resnorm_numIter(:,i) = R;
   
   if RESNORMq132i<R
       counter(:,i)= RESNORMq132i;
   else
       continue
   end
       
   
    
   
end
%%

N_counter = nnz(counter);
N_Resnorm_numIter=nnz(Resnorm_numIter);
Perc_q1_1_4 = (N_counter/N_Resnorm_numIter)*100;
display(Perc_q1_1_4);
% disp(R);

%% Plot the data points of the measurements and the model

figure(2);
plot(meas, ' bs', 'MarkerSize', 2, 'LineWidth', 4); % Actual data 
hold on
plot(pd_q132i, ' rx', 'MarkerSize', 2, 'LineWidth', 2) % Model predictions

