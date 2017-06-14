clear all;
close all;
measurements = [-0.577883;0.117166;-0.460717];
xlm = [1;1];
t = [0:1:100];
for i = 1:100%length(B) 
      [state,xlm] = LevMarFunc(measurements,xlm);
%       x(1, i) = state(1, 1);
%       y(1, i) = state(2, 1); 
      mag(1, i) = norm(xlm(:,i));
end
[state, cov] = TDOA([1;1], eye(2), measurements);
for i = 1:100
     [state, cov] = TDOA(state,cov,measurements);
     x2(1, i) = state(1, 1);
     y2(1, i) = state(2, 1); 
     mag(2, i) = norm(state);
end
x2 = horzcat(0.995,x2);
y2 = horzcat(0.995,y2);

[state, cov] = TDOA_EKF([1;1], eye(2), measurements);
for i = 1:100
     [state, cov] = TDOA_EKF(state,cov,measurements);
     x3(1, i) = state(1, 1);
     y3(1, i) = state(2, 1); 
     mag(3, i) = norm(state);
end
x3 = horzcat(1.005,x3);
y3 = horzcat(1.005,y3);

mag = horzcat([1;1;1],mag);
scatter(xlm(1,:), xlm(2,:), 'b');
hold on
scatter(x2, y2, 'g');
hold on
scatter(x3, y3, 'r');
xlim([0 6])
ylim([0 6])
title('Levenberg-Marquardt vs. IEKF vs. EKF Scatter Plot (10 iterations)')
xlabel('X Position (cm)') % x-axis label
ylabel('Y Position (cm)') % y-axis label

figure, plot(t, mag);
xlim([0 100])
ylim([0 7.5])
title('Levenberg-Marquardt vs. IEKF vs. EKF Magnitude Convergence (10 iterations)')
xlabel('Iterations (n)') % x-axis label
ylabel('Magnitude (cm)') % y-axis label
legend('LevMar','IEKF', 'EKF')