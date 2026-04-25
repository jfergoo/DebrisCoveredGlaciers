function [results] = DEBISO(tf,dt,restart,DT,filename,opts)

% DEBISO is a function that solves the coupled ice flow-debris transport problem
% as described in the paper "Modelling steady states and the transient response
% of debris-covered glaciers" by Ferguson and Vieli (in The Cryosphere).
%
% In order to run this code, a configuration file is needed. For an example
% of how to set up a model run, see the example code provided in the file
% DebrisExample.m which came with this code.
%
% This version was written on 15-Feb-2021 by James C. Ferguson and is
% protected by the license found in the accompanying file LICENSE.txt.


config = str2func(filename);
[x,b,cons] = config();

if restart == 1
    H0 = opts.H0;
    D = opts.D;
else
    H0 = 10*ones(size(x));
    D = zeros(size(x));
end

N = DT/dt;
t = [dt:dt:tf];
tt = [1:DT:tf];
NN = length(tt);

i_ct = 0;
j_ct = 1;

K = length(x);

HH = zeros(K,NN); 
DD = zeros(K,NN);
uu = zeros(K,NN);
aa = zeros(K,NN);

for i = 1:length(t) 
    
    i_ct = i_ct + 1;
    
    [H,D,x,u,a] = icedebrissolver(dt,H0,b,D,x,cons);
    
    H0 = H;
    
    if i_ct == N
        
        HH(:,j_ct) = H;
        uu(:,j_ct) = u;
        aa(:,j_ct) = a;
        DD(:,j_ct) = D;
        j_ct = j_ct + 1;
        i_ct = 0;
        
        figure(1)
        plot(x,b+H)
        hold on
        if restart == 1
            plot(x,b+opts.H0,':r')
        end
        plot(x,b,'k')
        hold off
        title("Time = " + floor(i*dt) + " years")
        
    end

end

results.H = HH;
results.D = DD;
results.u = uu;
results.a = aa;
results.t = tt;
results.x = x;
results.b = b;
results.cons = cons;

end

function [H,D,x,u,a] = icedebrissolver(tf,H,b,D,x,cons)

% 1D ice-debris solver

% Inputs: 
%   dt
%   H
%   b
%   D
%   x
%   cons
%
% Outputs:
%   H
%   D
%   x
%   x_t

Gamma  = 2 * cons.A * (cons.rho * cons.g)^3 / 5;

t = 0.0;
CFL = 0.25;

while t < tf

    % find dt from CFL condition for D and H evo
    u = u_solve(H,D,b,x,Gamma);
    udx = (u(2:end)-u(1:end-1))./(x(2:end)-x(1:end-1));
    dt_D = CFL / max(udx);
    [Dstag,dt_H] = Mahaffy(H,b,x,Gamma,tf);
    dt = min(min(dt_D,dt_H),tf-t);
    
    % compute forcing, accounting for BCs
    [a,a_H,f_D] = massbalance(x,b,H,D,dt,cons);
    
    % evolve H and D
    D = advect_upwind(D,u,x,dt,f_D);
    H = diffusion(H,Dstag,b,x,dt,a);
    
    % advance time
    t = t + dt;
    
    % find terminus position
    if max(D) ~= 0
        indy = find(H<=0.1,1);
        [dummy,indymax] = max(H);
        x_t = interp1(H(indymax:indy),x(indymax:indy),cons.H_t);
        indx = find(x>x_t,1);
        D(indx:end) = 0;
    else
        indy = find(H<=0.1,1);
        x_t = x(indy);
    end
        
    
end
end

    function u = u_solve(H,D,b,x,Gamma)
        
        % upwind differentiation
        dx = x(2:end)-x(1:end-1);
        K = length(H);
        n = 3;
        s = H+b;
        sgrad = 0*s;
        sgrad(2:K-1) = (s(3:K) - s(1:K-2))./(2*dx(1));
        sgrad(1) = sgrad(2);
        sgrad(end) = sgrad(end-1);

        u = -1.25*Gamma*H.^(n+1).*sgrad.^2.*sgrad;
                         
        if max(D) ~= 0
            
            HH = max(H);
            ind2 = find(H <= 0.1,1);
            indd = ceil(HH/dx(1));
            umin = mean(u(ind2-10*indd-2:ind2-2));
            u(ind2-2:end) = umin;
             
        else
                u = u;     
        end    
    end

    function [Dstag,dt0] = Mahaffy(H,b,x,Gamma,tf)
        
        % compute Mahaffy diffusivity
        
        Hmid = (H(1:end-1) + H(2:end))/2;
        h = H+b;
        gradh = (h(2:end) - h(1:end-1))./(x(2:end) - x(1:end-1));
        Dstag = Gamma * Hmid.^5 .* gradh.^2;
        
        max_Ddx = max(Dstag(2:end)./(x(3:end) - x(2:end-1)),Dstag(1:end-1)./(x(2:end-1) - x(1:end-2)));  
        max_Ddx = max(2./(x(3:end)-x(1:end-2)).*max_Ddx);
        
        
        if max_Ddx <= 0.0  % for zero thickness ice glaciers
            dt0 = tf;
        else
            dt0 = 0.25 / max_Ddx;
        end
    end

    function T = advect_upwind(T,u,x,dt,F)
        
        J = length(T);
        
        T(2:J) = T(2:J) - dt*(u(2:J).*T(2:J) - u(1:J-1).*T(1:J-1))./(x(2:J)-x(1:J-1));
        
        T = T + dt*F;
        
        T(end) = 0;
        
    end
   
    function H = diffusion(H0,DS,b,x,dt,a)
        
        % Solves dH/dt = d/dx(D * dh/dx) + a
        % with Dirichlet BCs at x(1) and x(end)
        % such that H(1) = H0(1) and H(end) = H0(end)

        J = length(x)-1;
        h = H0+b;
        H = H0;
        j = 2:J; 

        mu = dt*2./(x(j+1)-x(j-1));     

        H(j) = H(j) + ...
            mu.*(DS(j).*(h(j+1)-h(j))./(x(j+1)-x(j)) - ...
            DS(j-1).*(h(j)-h(j-1))./(x(j)-x(j-1)));

        a(1) = 0;
        a(end) = 0;

        H(j) = H(j) + dt*a(j);
        H = max(H,0.0);

    end
    
    function [a,a_H,f_D] = massbalance(x,b,H,D,dt,c)

% this function computes the surface mass balance (SMB) and debris source
% term for the icedebsolver
%
% inputs:
%
% x - the spatial grid
% b - the bed profile
% H - the ice thickness profile
% D - the debris thickness profile
% dt - the time step
% c - any constants needed for computation
%
% outputs:
%
% a - the SMB for the entire glacier
% a_H - the debris-free SMB
% f_D - the debris evolution source term
%
% these computations can be adjusted to suit the needs of the model

% Elevation dependent debris-free SMB


a_H = c.gamma*(H+b - c.ELA);
a_H(a_H > c.amax) = c.amax; %truncated at c.amax

% SMB for the entire glacier
a = a_H*c.D0./(c.D0 + D);

% forcing for the debris evolution equation

f_D = -a*c.c;
f_D(f_D<0) = 0; % negative source term not allowed
f_Dmax = c.c*H/dt;
f_D = min(f_D,f_Dmax); % cannot create more debris than stored in ice

% terminus correction

if max(D) ~= 0
%             
            indy = find(H<=0.1,1);
            [dummy,indymax] = max(H);
            x_tt = interp1(H(indymax:indy),x(indymax:indy),c.H_t);
            indx = find(x >= x_tt,1);
                  
            if c.cryo == 1
                         
                dx = x(2) - x(1);        
                K = length(H);        
                s = H+b;
                sgrad = 0*s;
                sgrad(2:K-1) = (s(3:K) - s(2:K-1))/dx(1);
                sgrad(1) = sgrad(2);
                tau_d = -c.rho*c.g*H.*sgrad;
                tau_d(H<=0.0) = 0.0;
                
                fac = c.lambdam*(c.maxth - tau_d)/(c.maxth-c.minth);
                fac(tau_d > c.maxth) = 0;
                fac(tau_d < c.minth) = c.lambdam;
                fac(H <=c.H_t) = 0;
         
                a = (1-fac).*a + fac.*a_H;
                
            end
            
            raty = (x(indx)-x_tt)/(x(indx)-x(indx-1)); % fraction of clean ice between grid points     
            
            if (raty < 0.5)

                a(indx) = (a(indx)*(0.5-raty)+a_H(indx)*raty) + a_H(indx)/2;
                a(indx+1:end) = a_H(indx+1:end);
                
            elseif (raty > 0.5)
                a(indx-1) = a_H(indx-1)*(raty-0.5) + a(indx-1)*(1-raty) + a(indx-1)/2;
                a(indx:end) = a_H(indx:end);
                
            else a(indx:end) = a_H(indx:end);
                
            end
end

    end
        
        

          
