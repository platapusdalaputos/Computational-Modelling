function S = BallStick_noAvox(x, bvals, qhat)
% Extract the parameters
S0 = x(1); 
diff = x(2);
f = x(3);
theta = x(4);
phi = x(5);
% Synthesize the signals according to the model
fibdir = [cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];

fibdotgrad = sum(qhat.*repmat(fibdir, [length(qhat) 1])');

S = S0*(f*exp(-bvals*diff.*(fibdotgrad.^2)) + (1-f)*exp(-bvals*diff));
