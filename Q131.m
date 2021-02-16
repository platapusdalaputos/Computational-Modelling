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
legend('Measured')
hold on
plot(pd_q131, ' rx', 'MarkerSize', 2, 'LineWidth', 2) % Model predictions
xlabel('q-index')
ylabel('S')
title('Q1.3.1')
legend('Measured','Modelled')


%% Q1.3.2i Two compartment model Nothing Added

% Produce starting matrix

startx = [x(1), x(2), x(3), x(4), x(5)];

% parameter_hat =  [-0.944815499367264   0.045855710911745  60.528238481733418  -1.544932176062478  -0.069107594560283];


[parameter_hat,RESNORMq132i,EXITFLAG,OUTPUT]=fminunc('q3twoCOMPSSD_Positive',startx,h,meas,bvals,qhat);

% Resnorm_q1_3_2i = RESNORMq132i;

S0_q132i = parameter_hat(1); % initial signal
d_q132i  = parameter_hat(2); % who knows what this is
f_q132i  = parameter_hat(3); % this guy is a joke
theta_q132i = parameter_hat(4); % I like him
phi_q132i = parameter_hat(5); % He does nothing for my life

parameter_hat = [ S0_q132i^2 d_q132i^2 Siggy(f_q132i) theta_q132i phi_q132i];

% pd stands for PreDicted Model hence the pd

pd_q132i = q3twoCOMPSSD_noAvox(parameter_hat, bvals, qhat);

% disp(startx);


%% Add Different Starting Points

R = RESNORMq132i;

% Housing all the integers
numIter = 10;

Resnorm_numIter=zeros(1,numIter);

x1_rand_index=zeros(1,numIter);

x2_rand_index=zeros(1,numIter);

x3_rand_index=zeros(1,numIter);

x4_rand_index=zeros(1,numIter);

x5_rand_index=zeros(1,numIter);

x6_rand_index=zeros(1,numIter);

x7_rand_index=zeros(1,numIter);

resnorm_q2_keep = zeros(1,numIter);

param_hat_best = zeros(1,5,numIter);

s0_boot = zeros(1, numIter);

d_boot = zeros(1, numIter);

f_boot = zeros(1, numIter);

counter = zeros(1,numIter);

j = 0;
k = 0;

for i = 2:numIter
    
    k=k+1;
   
    x1_rand = mvnrnd(startx(1),20);
    x1_rand_index(:,i)=x1_rand;
    x2_rand = mvnrnd(startx(2),0.5e-1);
    x2_rand_index(:,i)=x2_rand;
    x3_rand = mvnrnd(startx(3),0.5e-1);
    x3_rand_index(:,i)=x3_rand;
    x4_rand = mvnrnd(startx(4),pi);
    x4_rand_index(:,i)=x4_rand;
    x5_rand = mvnrnd(startx(5),pi);
    x5_rand_index(:,i)=x5_rand;
    
    startx2 = [x1_rand_index(:,i),x2_rand_index(:,i),x3_rand_index(:,i),x4_rand_index(:,i),x5_rand_index(:,i)];
    
   [parameter_hat_q132i,RESNORMq132i,EXITFLAG,OUTPUT]=fminunc('q3twoCOMPSSD_Positive',startx2,h,meas,bvals,qhat);
   
%    Resnorm_numIter(:,i) = R;
%    
%    if RESNORMq132i<R
%        counter(:,i)= RESNORMq132i;
%    else
%        continue
%    end
   
   Resnorm_numIter(:,1)= RESNORMq132i;
   
    
   Resnorm_numIter(:,i)= RESNORMq132i;
       
   if Resnorm_numIter(:,i)<= min(Resnorm_numIter(i-k:i))
        
        j=j+1;
        
        resnorm_q2_keep(:,i) = Resnorm_numIter(:,i);
        index = find(resnorm_q2_keep);
        
        param_hat_best(:,:,j) = parameter_hat_q132i;
        s0_boot(:,j) = param_hat_best(1,1,j).^2;
        d_boot(:,j) = param_hat_best(1,2,j).^2;
        f_boot(:,j) = Siggy(param_hat_best(1,3,j));
        
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
%%

P = param_hat_best(:,:,5);

S0_q132ii = P(1); % initial signal
d_q132ii  = P(2); % who knows what this is
f_q132ii  = P(3); % this guy is a joke
theta_q132ii = P(4); % I like him
phi_q132ii = P(5); % He does nothing for my life

Ps = [1.0049^2 0.0378^2 Siggy(0.3020) 1.5447 3.0587];


%%
% pd stands for PreDicted Model hence the pd

pd_q132ii = q3twoCOMPSSD_noAvox(Ps, bvals, qhat);

% disp(startx);


%% Plot the data points of the measurements and the model

figure(2);
plot(meas, ' bs', 'MarkerSize', 2, 'LineWidth', 4); % Actual data 
hold on
plot(pd_q132ii, ' rx', 'MarkerSize', 2, 'LineWidth', 2) % Model predictions
xlabel('q-index')
ylabel('S')
title('Q1.3.2i')
legend('Measured','Modelled')

%% Q1.3.2ii

N = length(meas);

[QX, QY, QZ] = deal(qhat(1,:),qhat(2,:),qhat(3,:));

G = [ones(1,N); -bvals .* QX.^2; -2*bvals.*QX.*QY; -2*bvals.*QX.*QZ; -bvals.*QY.^2; -2*bvals.*QY.*QZ; -bvals.*QZ.^2]';

x = G\log(meas); 

% [logS0, Dxx, Dxy, Dxz, Dyy, Dyz, Dzz] = deal(x(1),x(2),x(3),x(4),x(5),x(6),x(7));

Diffusion_tensor= [x(1),x(2),x(3); x(3),x(5),x(6); x(4),x(6),x(7) ];

predicted_diffusion = exp(G * x);
SSD_diffusion = sum((predicted_diffusion - meas) .^2);

%%
figure(3);
plot(meas, ' bs', 'MarkerSize', 2, 'LineWidth', 4); % Actual data 
legend('Measured')
hold on
plot(predicted_diffusion, ' rx', 'MarkerSize', 2, 'LineWidth', 2) % Model predictions
xlabel('q-index')
ylabel('S')
title('Q1.3.2diffusion')
legend('Measured','Modelled')

%% Add Different Starting Points


% Housing all the integers
numIter = 20;

x1_rand_index=zeros(1,numIter);
x2_rand_index=zeros(1,numIter);
x3_rand_index=zeros(1,numIter);
x4_rand_index=zeros(1,numIter);
x5_rand_index=zeros(1,numIter);
x6_rand_index=zeros(1,numIter);
x7_rand_index=zeros(1,numIter);


% 
param_hat_best = zeros(1,5,numIter);
s0_boot = zeros(1, numIter);
d_boot = zeros(1, numIter);
f_boot = zeros(1, numIter);
counter = zeros(1,numIter);



for i = 1:numIter
    
   
    x1_rand = mvnrnd(x(1),0.1);
    x2_rand = mvnrnd(x(2),0.0001);
    x3_rand = mvnrnd(x(3),0.00001);
    x4_rand = mvnrnd(x(4),0.00001);
    x5_rand = mvnrnd(x(5),0.00001);
    x6_rand = mvnrnd(x(6),0.00001);
    x7_rand = mvnrnd(x(7),0.00001);

    x1_rand_index(:,i)=x1_rand;
    x2_rand_index(:,i)=x2_rand;
    x3_rand_index(:,i)=x3_rand;
    x4_rand_index(:,i)=x4_rand;
    x5_rand_index(:,i)=x5_rand;
    x6_rand_index(:,i)=x6_rand;
    x7_rand_index(:,i)=x7_rand;  
    
    startxdiff = [x1_rand_index(:,i),x2_rand_index(:,i),x3_rand_index(:,i),x4_rand_index(:,i),x5_rand_index(:,i), x6_rand_index(:,i),x7_rand_index(:,i)];
    
%     startxdiff = [x(1),x(2),x(3),x(4),x(5),x(6),x(7)];

%       startxdiff = [x1_rand,x2_rand,x3_rand,x4_rand,x5_rand,x6_rand,x7_rand];

    
   [parameter_hat_diffusion,RESNORM_diffusion,EXITFLAG,OUTPUT]=fminunc('ledZepplinpart1',startxdiff,h,meas,bvals,qhat);
   
   
   
end

%%



Predicted_Zep = ledZepplinpart1_Sig(parameter_hat_diffusion,bvals,qhat);

figure(3);
plot(meas, ' bs', 'MarkerSize', 2, 'LineWidth', 4); % Actual data 
legend('Measured')
hold on
plot(Predicted_Zep, ' rx', 'MarkerSize', 2, 'LineWidth', 2) % Model predictions
xlabel('q-index')
ylabel('S')
title('Q1.3.2.Zepplin')
legend('Measured','Modelled')



%%

numIter = 10;

for i = 1:numIter
    
   
%     x1_rand = mvnrnd(x(1),0.1);
%     x2_rand = mvnrnd(x(2),0.0001);
%     x3_rand = mvnrnd(x(3),0.00001);
%     x4_rand = mvnrnd(x(4),0.00001);
%     x5_rand = mvnrnd(x(5),0.00001);
%     x6_rand = mvnrnd(x(6),0.00001);
%     x7_rand = mvnrnd(x(7),0.00001);
% 
%     x1_rand_index(:,i)=x1_rand;
%     x2_rand_index(:,i)=x2_rand;
%     x3_rand_index(:,i)=x3_rand;
%     x4_rand_index(:,i)=x4_rand;
%     x5_rand_index(:,i)=x5_rand;
%     x6_rand_index(:,i)=x6_rand;
%     x7_rand_index(:,i)=x7_rand;  
%     
%     startxdiff = [x1_rand_index(:,i),x2_rand_index(:,i),x3_rand_index(:,i),x4_rand_index(:,i),x5_rand_index(:,i), x6_rand_index(:,i),x7_rand_index(:,i)];
    
    startxdiff = [x(1),x(2),x(3),x(4),x(5),x(6),x(7)];

%       startxdiff = [x1_rand,x2_rand,x3_rand,x4_rand,x5_rand,x6_rand,x7_rand];
%     bootstrap_sample(i,:) = predicted + normrnd(predicted, sigma);
% %     bootstrap_sample(i,:) = normrnd(Avox, sigma);
% 
%     
%     [parameter_hat,RESNORM_boot,EXITFLAG,OUTPUT]=fminunc('BallStickSSD_Positive_q114',startx,h,bootstrap_sample(i,:),bvals,qhat);
    
    
   [parameter_hat_diffusion,RESNORM_diffusion,EXITFLAG,OUTPUT]=fminunc('ledZepplinpart2',startxdiff,h,meas,bvals,qhat);
   
   
   
end
%%
Predicted_Zep2 = ledZepplinpart2_Sig(parameter_hat_diffusion,bvals,qhat);

figure(3);
plot(meas, ' bs', 'MarkerSize', 2, 'LineWidth', 4); % Actual data 
legend('Measured')
hold on
plot(Predicted_Zep2, ' rx', 'MarkerSize', 2, 'LineWidth', 2) % Model predictions
xlabel('q-index')
ylabel('S')
title('Q1.3.2.ZepplinT')
legend('Measured','Modelled')

