load("24jun19_highQ.mat");

% Load BTCs and time vector from time hour zero to 30 hr and convert to
% concentration
time_hr = zero_time_R4_US;
n = length(time_hr);
n_wanted=1024;

cout_1A = (R4_US_cond_no_bg(1:n).*0.3);
%cout_1A(cout_1A<0) = 0; % Replace negative values with zero

cout_1B = (R4_DS_cond_no_bg(1:n).*0.3); 
%cout_1B(cout_1B<0) = 0;

cout_2A = (R3_US_cond_no_bg(1:n).*0.3);
%cout_2A(cout_2A<0) = 0;

cout_2B = (R3_DS_cond_no_bg(1:n).*0.3); 
%cout_2B(cout_2B<0) = 0;

in=cout_1A(1:n_wanted); out=cout_1B(1:n_wanted); time=time_hr(1:n_wanted);
save ('highQ_R1.mat', 'time', 'in', 'out')

figure(44)
plot(time,in,'o-')
hold on 
plot(time,out)
axis([-1 5 -1 20])
legend
hold off

in=cout_2A(1:n_wanted); out=cout_2B(1:n_wanted); 
save ('highQ_R2.mat', 'time', 'in', 'out')

figure(45)
plot(time,in,'o-')
hold on 
plot(time,out)
axis([-1 5 -1 20])
legend
hold off
