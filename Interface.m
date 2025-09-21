clear all;

pi = 3.14; 
 G = 0.5; % Greshof number
  B = 0.01; % Heat Source Parameter

h = @(x) 1+0.5 * sin(2*pi*x);
%1 + 0.5 * sin(2 * pi * x);
       
q =-1.025;
q1 = @(mu, q) (-(0.8)^8 * (mu-1) * G*B/120*mu + (0.8)^6 * (mu-1) * G*B/20*mu - (0.8)^5 * (mu-5/4) * G*B/30*mu - (0.8)^4 * (mu-1) * (1+G*B/24*mu) + (0.8)^3 * (mu-3/2) * (q + 1+ G*B/30*mu) + (0.8) * (G*B/(120*mu) + 3*q/2 + 1/2))/((mu-1)*(0.8)^3 + 1);

x_values = linspace(0, 1, 11)


% Values of mu to consider
mu_values = [0.1, 1, 10];
colors = {'-r', '-g', '-b'};

figure(5);
hold on;

% Loop over different values of mu
for j = 1:length(mu_values)

    mu = mu_values(j);
    x_all = [];
    
    h_all = [];
    q1_value = q1(mu, q);
    for i = 1:length(x_values)
        x = x_values(i);
        h_value = h(x);
        q1_value = q1(mu, q);
        
        q2 = q - q1_value;

        A = @(mu) (mu - 1) * G*B/(120*mu);
        B1 = @(mu, h) -(mu - 1) * G*B*h^2/(20*mu);
        C = @(mu, h) (mu - 5/4) * G*B*h^3/(30*mu);
        D = @(mu, h) (mu -1) * (1 + G*B*h^4/(24*mu));
        E = @(mu, h, q1) - ((mu - 3/2) * (q + h + G*B*h^5/(30*mu)) - q1*(mu - 1));
        F = @(h) -((q + h) * 3*h^2/2 + G*B*h^7/120 - h^3);
        G1 = @(q1, h) q1*h^3;
        
        a = A(mu);
        b1 = B1(mu, h_value);
        c = C(mu, h_value);
        d = D(mu, h_value);
        e = E(mu, h_value, q1_value);
        f = F(h_value);
        g1 = G1(q1_value, h_value);

        % Calculate the roots
        roots_array = roots([a 0 b1 c d e 0 f g1]);
        
        % Filter only the real roots
        real_roots = roots_array(imag(roots_array) == 0 & real(roots_array) < h_value & real(roots_array) > 0);
        
disp(real_roots);

        % Store x values and corresponding real roots
        x_all = [x_all, repmat(x, 1, length(real_roots))];
         h_all = [h_all, real(real_roots)'];
    end
    
    % Plot the results
    plot(x_all, h_all, colors{j},'Linewidth',1.6, 'DisplayName', sprintf('mu = %1.3f', mu));
end

% Plot h(x) versus x
h_values = arrayfun(h, x_values);
plot(x_values, h_values,'-k', 'Linewidth',1.6, 'DisplayName', 'h(x)');

xlabel('z');
ylabel('r');
title('Plot of x versus Different \mu Values and h(x)');
legend show;
grid off;
hold off;


