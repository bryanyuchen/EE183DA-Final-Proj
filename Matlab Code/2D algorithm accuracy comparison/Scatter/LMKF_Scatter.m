clear all;
close all;
%measurements = [-0.577;0.117;-0.46];
 B = [[-35; 0; -38], [-35; -4; -38], [-36; 0; -39], [-36; 0; -38], [-33; 0; -37], [-35; 0; -37], [-39; 0; -37], [-34; -4; -40], [-35; -5; -40], [-35; -4; -39], [-33; -5; -38], [-36; -6; -40], [-34; -5; -39], [-34; 0; -37], [-33; -5; -39], [-32; -5; -40], [-34; -5; -39], [-33; -6; -40], [-34; -4; -39], [-33; 0; -37], [-35; 0; -40], [-34; 0; -40], [-32; -5; -39], [-34; 0; -37], [-34; -4; -39], [-34; 0; -38], [-35; -5; -38], [-33; -4; -38], [-35; 0; -41], [-33; 0; -37], [-34; -5; -39], [-35; 0; -39], [-35; 0; -37], [-36; 0; -38], [-35; 0; -38], [-34; 0; -37], [-35; 0; -36], [-34; 0; -36], [-35; 0; -37], [-35; 0; -36], [-35; 0; -39], [-35; 0; -37], [-34; 0; -38], [-34; 0; -35], [-35; 0; -36], [-35; 0; -37], [-35; 0; -37], [-35; 0; -38], [-36; 0; -38], [-35; 0; -38], [-35; -4; -38], [-35; 0; -37], [-35; 0; -38], [-33; 0; -36], [-36; 0; -39], [-34; -4; -40], [-36; -4; -36], [-35; 0; -38], [-34; 0; -37], [-34; 0; -38], [-36; 0; -39], [-34; 0; -38], [-34; 0; -36], [-35; 0; -39], [-35; 0; -39], [-35; 3; -40], [-35; 0; -39], [-35; 0; -38], [-37; -5; -40], [-36; 0; -41], [-35; 0; -42], [-36; -4; -40], [-34; 0; -40], [-35; -5; -40], [-35; -4; -39], [-34; -4; -39], [-34; 0; -38], [-34; -4; -40], [-35; -4; -40], [-34; -4; -39], [-35; 0; -40], [-35; -4; -39], [-36; -5; -40], [-36; -5; -41], [-36; -3; -39], [-36; -5; -40], [-35; -5; -40], [-36; -5; -39], [-35; -4; -39], [-37; 0; -41], [-35; 0; -38], [-34; 0; -40], [-35; -6; -40], [-33; -5; -40], [-31; -7; -35], [-34; -5; -40], [-35; 0; -40], [-34; -7; -40], [-32; -8; -39], [-33; -6; -40]];
 lengthB = length(B);
 random = randn(3,lengthB);
 B = B + random;
  
   %LPF
   for i = 1:lengthB-8
      temp = B(:,i:i+8); 
      B(:,i+4) = mean(temp,2);
   end
 B(1,:) = B(1,:)-4;
 B(2,:) = B(2,:)+2;
 B = 0.0143 * B;
 cov = eye(2);
 state = [5;5];
 
for i = 1:length(B) 
      state_est = LevMarFunc(B(:,i),state);
      if abs(state_est) > 3*abs(state) % catch off cases
          continue;
      end
      x(1, i) = state_est(1, 1);
      y(1, i) = state_est(2, 1); 
      t(1, i) = i;

      [state, cov] = Kalmf(state,cov,[x(i);y(i)]);
      x2(1, i) = state(1, 1);
      y2(1, i) = state(2, 1); 
 end
% [state, cov] = TDOA_EKF([1;1], eye(2), measurements);
% for i = 1:10
%      [state, cov] = TDOA_EKF(state,cov,measurements);
%      x3(1, i) = state(1, 1);
%      y3(1, i) = state(2, 1); 
%      mag(3, i) = norm(state);
% end

% for i=1:length(B) 
    %plot(x2(i),y2(i),'or','MarkerSize',5,'MarkerFaceColor','r')
    
    scatter(x, y,'b');
    hold on
    scatter(x2,y2,'m');
    hold on
    scatter(x2(100),y2(100),'g');
    hold on
    scatter(0,15,'r');
    %set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    axis([-20 20 0 40])
    xlabel('x position (cm)');
    ylabel('y position (cm)');
    title('LM Real Data Scatter, Constant State Space KF (mnoise=3cm), 100 samples, LPF');
    legend('LM Estimation', 'LM after Kalman Filtering', 'Final', 'Expected ~ (0,15)');
    pause(.1)
    hold off
% end

