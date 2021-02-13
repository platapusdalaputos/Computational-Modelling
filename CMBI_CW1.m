%% Coursework 1
%% Q1.1.1




load('data.mat');
dwis=double(dwis);
dwis=permute(dwis,[4,1,2,3]);

% Middle slice of the lst image volume, which has b =0
imshow(flipud(squeeze(dwis(1,:,:,72))'),[]);

% Middle slice of the 2nd image volume, which has b=1000
imshow(flipud(squeeze(dwis(2,:,:,72))'), []);

qhat = load('bvecs');
size(dwis)

%%
bvals = 1000*sum(qhat.*qhat);

Avox = dwis(:,92,65,72); % A single voxel element
%% Computing the value of the Matric x on one voxel

G = [ones(1,108); -bvals.*qhat(1,:).^2; -2*bvals.*qhat(1,:).*qhat(2,:); -2*bvals.*qhat(1,:).*qhat(3,:); -bvals.*qhat(2,:).^2; -2*bvals.*qhat(2,:).*qhat(3,:); -bvals.*qhat(3,:).^2]';

x = G\log(Avox);

%% Computing Diffusion Tensor

diffusionTensor = [[x(2) x(3) x(4)]; [x(3) x(5) x(6)]; [x(4) x(6) x(7)]];

Mean_diffusivity = trace(diffusionTensor); % Computing the trace of the diffusion tensor




%% Computing Fractional Anisotrophy

[a, b] = eig(diffusionTensor); %  Computing eigenvalues 

% Two columns, Matrix a = the eigenvectors of diffusion tensor   

% Matrix b = eigenvalues of the diffusion tensor, largest values is best guess of fibre orientation


%% Mapping Diffusion tensor of a whole slice

% there are seven values because of the 7 unknowns we need to calculate 

dt_map72 = zeros(7, 145, 174); 


for i = 1:145
    for j =1:174
        A = dwis(:,i,j,72);
        if min(A)>0
            dt_map72(:,i,j) = G\log(A);
            
        end
    end
end

% imshow(flipud(squeeze(dt_map72(1,:,:))'), []);
imshow(flipud(squeeze(exp(dt_map72(7,:,:)))'), []);
% imshow(flipud(squeeze(dt_map72(1,:,:))'), []);

%% Mean Diffusivity & Fractional Anisotrpy


MD_map72 = zeros(3, 3, 145, 174);
Mean_diffusivity_map72 = zeros(1,145, 174);
FA_map72 = zeros(1,145, 174);
% dwis = diffusion weighted images

for i = 1:145
    for j =1:174        
         MD_map72(:,:,i,j) = [[dt_map72(2,i,j) dt_map72(3,i,j) dt_map72(4,i,j)]; [dt_map72(3,i,j) dt_map72(5,i,j) dt_map72(6,i,j)]; [dt_map72(4,i,j) dt_map72(6,i,j) dt_map72(7,i,j)]];                   
         Mean_diffusivity_map72(:,i,j) = trace(MD_map72(:,:,i,j));
         FA_map72(:,i,j) = sqrt(0.5*(3 - (1/(trace(( MD_map72(:,:,i,j)/trace(MD_map72(:,:,i,j)))^2)))));
    end
end
%% plot Anisotropy

Red = zeros(1,145, 174);
Green = zeros(1,145, 174);
Blue = zeros(1,145, 174);


for i = 1:145
    for j =1:174        
         MD_map72(:,:,i,j) = [[dt_map72(2,i,j) dt_map72(3,i,j) dt_map72(4,i,j)]; [dt_map72(3,i,j) dt_map72(5,i,j) dt_map72(6,i,j)]; [dt_map72(4,i,j) dt_map72(6,i,j) dt_map72(7,i,j)]];                   
         Mean_diffusivity_map72(:,i,j) = trace(MD_map72(:,:,i,j));
         FA_map72(:,i,j) = sqrt(0.5*(3 - (1/(trace(( MD_map72(:,:,i,j)/trace(MD_map72(:,:,i,j)))^2)))));
         Red(:,i,j) = FA_map72(:,i,j).*cos(x(2));
         Green(:,i,j) = FA_map72(:,i,j).*cos(x(5));
         Blue(:,i,j) = FA_map72(:,i,j).*cos(x(7));
    end
end



%% Plot Mean Diffusivity

figure(2);
subplot(1,4,1)
imshow(flipud(squeeze(Mean_diffusivity_map72(1,:,:))'), []);

subplot(1,4,2)
imshow(flipud(squeeze(exp(Mean_diffusivity_map72(1,:,:)))'), []);

subplot(1,4,3)
imshow(flipud(squeeze(dt_map72(1,:,:))'), []);

subplot(1,4,4)
imshow(flipud(squeeze(exp(dt_map72(1,:,:)))'), []);


%% plot fractional Anisotrophy
figure(3);

imshow(flipud(squeeze((FA_map72(1,:,:)))'), []);
% hold on
% imshow(flipud(squeeze((Green(1,:,:)))'), []);


% data = Red';                                 % Original Data
% datascl = bsxfun(@times, data, [100 50 1]);         % ‘Weigh’ Data
% datadisp = reshape(sum(datascl,2), 10, 10);         % Sum & Reshape
% colormap(jet);                                      % Set Appropriate Colormap
% % figure(1)
% image(datadisp)  

%% Q1.1.2

Avox = dwis(:,92,65,72);


% Define a starting point for the non-linear fit
startx1 = 3.5e+00;  
startx2 = 3e-03;   
startx3 = 2.5e-01;
startx4 = 0;
startx5 = 0;

startx =[startx1, startx2, startx3, startx4, startx5];





% algorithm.
h=optimset('MaxFunEvals',20000,'Algorithm','quasi-newton','TolX',1e-10,'TolFun',1e-10, 'Display', 'iter');
oldopts = optimset('TolFun',1e-6);
options = optimset(oldopts,'PlotFcns','optimplotfval','TolX',1e-8);

% Now run the fitting
[parameter_hat1,RESNORM,EXITFLAG,OUTPUT]=fminunc('BallStickSSD',startx,h,Avox,bvals,qhat);


% fun = @BallStickSSD;
% x0 = [startx1, startx2, startx3, startx4, startx5, Avox, bvals, qhat]; 
% [x,fval] = fminsearch(fun,x0,options);

% RESNORM is the value of the function at the solution found (parameter_hat)
pd1 = BallStickSSD(parameter_hat1, Avox, bvals, qhat);
% Predict measurements from model
pd = BallStick_noAvox(parameter_hat1, bvals, qhat);



%% Plot results for Q1.1.2


figure(4);

plot(Avox, ' bs', 'MarkerSize', 16, 'LineWidth', 4); % Actual data 

hold on;

plot(pd, ' rx', 'MarkerSize', 16, 'LineWidth', 4) % Model predictions
set(gca, 'FontSize', 15);
xlabel('q index');
ylabel('S');
legend('Data', 'Model');
title('q1.1.2');

%% Q1.1.3



[parameter_hat,RESNORM,EXITFLAG,OUTPUT]=fminunc('BallStickSSD_Positive',startx,h,Avox,bvals,qhat);

Resnorm_q1_1_3 = RESNORM;
S0 = parameter_hat(1);
d  = parameter_hat(2);
f  = parameter_hat(3);
theta = parameter_hat(4);
phi = parameter_hat(5);

parameter_hat = [ S0^2 d^2 atan(f)/(pi+0.5) theta phi];


pd_positive = BallStick_noAvox(parameter_hat, bvals, qhat);


%% Plot results for Q1.1.3

figure(5);

plot(Avox, ' bs', 'MarkerSize', 16, 'LineWidth', 4); % Actual data 
hold on;
plot(pd_positive, ' rx', 'MarkerSize', 16, 'LineWidth', 4) % Model predictions
set(gca, 'FontSize', 15);
xlabel('q index');
ylabel('S');
legend('Data', 'Model');
title('q1.1.3');


%% Q1.1.4


x1 = 3.5e+00;   
x2 = 3e-03;
x3 = 2.5e-01;
x4 = 0;
x5 = 0;


startx1 = [x1, x2, x3, x4, x5];

numIter = 500;
Resnorm_numIter=zeros(1,numIter);

x1_rand_index=zeros(1,numIter);
x2_rand_index=zeros(1,numIter);
x3_rand_index=zeros(1,numIter);
x4_rand_index=zeros(1,numIter);
x5_rand_index=zeros(1,numIter);


% Computing RESNORM as a reference with the same starting point as 1.1.2

[parameter_hat1_1_4,RESNORM_Q4,EXITFLAG,OUTPUT]=fminunc('BallStickSSD_Positive',startx1,h,Avox,bvals,qhat);

R = RESNORM_Q4;
display(R);

%%


counter = zeros(1,numIter);

for i = 1:numIter
    
    x1_rand = normrnd(x1,0.5);
    x1_rand_index(:,i)=x1_rand;
    x2_rand = normrnd(x2,1e-1);
    x2_rand_index(:,i)=x2_rand;
    x3_rand = normrnd(x3,1);
    x3_rand_index(:,i)=x3_rand;
    x4_rand = normrnd(x4,pi);
    x4_rand_index(:,i)=x4_rand;
    x5_rand = normrnd(x5,pi);
    x5_rand_index(:,i)=x5_rand;
    
    startx2 = [x1_rand_index(:,i),x2_rand_index(:,i),x3_rand_index(:,i),x4_rand_index(:,i),x5_rand_index(:,i)];
    
   [parameter_hat1_1_4,R1,EXITFLAG,OUTPUT]=fminunc('BallStickSSD_Positive',startx2,h,Avox,bvals,qhat);
   
   Resnorm_numIter(:,i) = R;
   
   if R1<=R
       counter(:,i)= R1;
   else
       continue
   end
       
   
    
   
end

%% 

% Calculate the percentage of iterations that are less than or are below
% the minimum value

N_counter = nnz(counter);
N_Resnorm_numIter=nnz(Resnorm_numIter);
Perc_q1_1_4 = (N_counter/N_Resnorm_numIter)*100;
display(Perc_q1_1_4);


%% Q1.1.5

% % Attempt parameter maps S0
% 
% S0_map72 = zeros(1,145, 174);
% f_map72 = zeros(1,145, 174);
% d_map72 = zeros(1,145, 174);
% theta_map72 = zeros(1,145, 174);
% phi_map72 = zeros(1,145, 174);
% RESNORM_map72 = zeros(1,145, 174);
% 
% parameter_hat_data = zeros(1, 5, 145, 174); 
% 
% for i = 1:145
%     for j =1:174  
%         Avox_multi = dwis(:,i,j,72);
%         [parameter_hat,RESNORM,EXITFLAG,OUTPUT]=fminunc('BallStickSSD_Positive',startx,h,Avox_multi,bvals,qhat);
%          parameter_hat_data(:,:,i,j) = parameter_hat;
% 
% 
%     end
%     
% end

% S0_map72 = zeros(1,145, 174);
% f_map72 = zeros(1,145, 174);
% d_map72 = zeros(1,145, 174);
% theta_map72 = zeros(1,145, 174);
% phi_map72 = zeros(1,145, 174);
% RESNORM_map72 = zeros(1,145, 174);

load('parameter_hat_data.mat')
for i = 1:145
    for j =1:174  
        Avox_multi = dwis(:,i,j,72);
         parameter_hat_data(:,:,i,j) = parameter_hat;
    end
    
end




%% Plot image for Q1.1.5



load('parameter_hat_data.mat')

S0_map72 = zeros(1,145, 174);
f_map72 = zeros(1,145, 174);
d_map72 = zeros(1,145, 174);
theta_map72 = zeros(1,145, 174);
phi_map72 = zeros(1,145, 174);
RESNORM_map72 = zeros(1,145, 174);

for i = 1:145
    for j =1:174  
         S0_map72(:,i,j) = parameter_hat_data(1,1,i,j);
         f_map72(:,i,j) = parameter_hat_data(1,2,i,j);
         d_map72(:,i,j) = parameter_hat_data(1,3,i,j);
         theta_map72(:,i,j) = parameter_hat_data(1,4,i,j);
         phi_map72(:,i,j) = parameter_hat_data(1,5,i,j);
    end
end

S0_map72_sq = S0_map72.^2;
f_map72_sq = 1/(1+(-1.*(f_map72)));
d_map72_sq = d_map72.^2;




% figure(5);
% subplot(2,3,1)
% imshow(flipud(squeeze(S0_map72(:,:,:))'), []);
% subplot(2,3,2)
% imshow(flipud(squeeze(f_map72(:,:,:))'), []);
% subplot(2,3,3)
% imshow(flipud(squeeze(d_map72(:,:,:))'), []);
% subplot(2,3,4)
% imshow(flipud(squeeze(theta_map72(:,:,:))'), []);
% subplot(2,3,5)
% imshow(flipud(squeeze(phi_map72(:,:,:))'), []);


%% Q1.2.1 Non - parametric modelling


% Non - parametric booting
% Bootstrap method requires sampling the input randomly 
orig_data = Avox';
num_bs_iterations = 100;
N = 108;

bootstrap_sample = zeros(num_bs_iterations, N);
s0_boot = zeros(1, num_bs_iterations);
d_boot = zeros(1, num_bs_iterations);
f_boot = zeros(1, num_bs_iterations);

cols = size(orig_data,2);
resnorm_q2 = zeros(1,num_bs_iterations);

for i = 1:num_bs_iterations
    
    P = randperm(cols);
    
    bootstrap_sample(i,:) = orig_data(:,P);
    
    [parameter_hat,RESNORM_boot,EXITFLAG,OUTPUT]=fminunc('BallStickSSD_Positive',startx,h,bootstrap_sample(i,:),bvals,qhat);

    resnorm_q2(:,i) = RESNORM_boot;
    s0_boot(:,i) = parameter_hat(1);
    d_boot(:,i) = parameter_hat(2);
    f_boot(:,i) = parameter_hat(3);
    
    
end 

%% Q1.2.1 Parametric Bootstrapping


% Use the value sigma (and try sigma^2 too to find the
[sumResSquared, S] = BallStickSSD_Positive(parameter_hat,Avox, bvals, qhat);
K = 108 ;% Number of measurements
N = 5 ;% Number of Parameters

sigma = sqrt((1/(K-N))*sumResSquared);


%% 2_Sigma Value for Parameters

% [Mean_s0, two_std_s0] = normfit(s0_boot);

% s0
Mean_s0 = mean(s0_boot);
two_sigma_s0 = 2.*(std(s0_boot));
histogram(s0_boot);

% d


% f

%% 95% Values for Parameters


% s0


% d


% f

%% Q1.2 Plot
figure(5);
subplot(2,2,1)
histogram(resnorm_q2,10);
subplot(2,2,2)

% histogram(s0_boot);
% subplot(2,2,3)
% 
% histogram(f_boot);
% subplot(2,2,4)
% 
% histogram(d_boot);

%% Q1.2.2





%% Q1.3.1





%% Q1.3.2

