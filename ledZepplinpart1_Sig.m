function S = ledZepplinpart1_Sig(x, bvals, qhat)

S0 = x(1);
diff = x(2);
f = x(3);
theta = x(4);
phi = x(5);
lambda1 = x(6); 
lambda2 = x(7);

fibDir = [cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];

fibDotGrad_s = (sum(qhat .* repmat(fibDir, [length(qhat) 1])')).^2;

Si = exp(-bvals * diff .* fibDotGrad_s); 
Se = exp(-bvals.*(lambda2 + (lambda1 - lambda2)*fibDotGrad_s)); 

S = S0*(f*Si + (1-f)*Se);
