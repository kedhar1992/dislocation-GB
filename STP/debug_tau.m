clc
clear all
close all

sigma_a = [0 1 0; 1 0 0; 0 0 0];

n = [-1,1,0];
d = [1,1,1];
n = n/norm(n);
d = d/norm(d);
x = 1:180;

for i = 1:180
    eulers = [i 54.7403 45];
    g = eulers2g(eulers);
    
    sigma_rot = g * (sigma_a * transpose(g));
    tau(i,:) = abs(dot((sigma_rot * n'),d));    % n - slip plane normal;  d - burger vector
    
end

plot(x,tau)
