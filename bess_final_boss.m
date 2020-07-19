N = 51; 
nu = 4;
b = 30*10^(-6);                  % maks raduis [my m]
r = linspace(0, b, N);           % radius change over time
t = linspace(0, 50, N);          % time 
kappa = 1.36 * 10^(-18)                % permabilitet 
T = 50;                          % Total tid [s]
omega = (2*pi)/T;                % [1/s]
    
A = 1 ;                        % p(a) = A pressure at start
B = 0.1 ;                          % p(b) = B pressure at end

[R_n,J,z] = Bessel2(nu,kappa,omega,N,b,r,B,A,1);
[p,J,z,R_n] = pressure(t,omega,R_n,J,z);

txt = {};
for i=0:nu
    txt{end+1} = sprintf('J_%d',i);
end

figure1=figure('Position', [100, 100, 1024, 1200]);
subplot(2,1,1);
plot(z,J);
grid on
legend(txt,'Location','Best')
title(['Bessel Functions of the First Kind for $\nu \in $[0,', num2str(nu),']'],'interpreter','latex')
xlabel('z','interpreter','latex')
ylabel('$J_\nu(z)$','interpreter','latex')
%xlim([0,10]);

subplot(2,1,2);
plot(t,p);
grid on
%legend(txt2)
title('p(t,r)')
xlabel('r','interpreter','latex')
ylabel('p(t)','interpreter','latex')
%xlim([0,10]);
