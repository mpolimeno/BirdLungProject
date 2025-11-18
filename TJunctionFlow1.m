function [p_A2A1,p_A3A1,phi_A2A1,phi_A3A1] = TJunctionFlow1(u_minus,u_plus,zeta_vort,cf,RuttaCondition)

% TJunctionFlow1 - Computes the potential flow solution for flow inside a T-junction with a single vortex
%
%   This function models the potential flow around a T-junction with a vortex positioned at a given location 
%   in the complex plane. 
%
%   INPUTS:
%     u_minus      - Inflow velocity on the left horizontal branch (scalar).
%     u_plus       - Outflow velocity on the right horizontal branch (scalar).
%                    Possible values for (u_minus,u_plus):
%                    Case 1: (1,0.7);
%                    Case 2: (1,0.75);
%                    Case 3: (-1,-0.75);
%     zeta_vort    - Complex number representing the position of the vortex in the complex plane.
%                    Possible values for zeta_vort:
%                    Case 1: zeta_vort = 0 + 1i*pi/2   -> Position for fixed Gamma=0.5
%                    Case 2: zeta_vort = -0.25 + 1i*1.25 -> Vortex in the vertical branch
%                    Case 3: zeta_vort = -0.75 + 1i*0.75 -> Vortex in the left horizontal branch
%     cf           - Integer representing the figure object number for plotting.
%     RuttaCondition - Boolean flag. If True, the Rutta condition is applied at the left corner, 
%                      and gamma is computed. If False, gamma is set to 0.5

% NOTES
% Case 1 applies for RuttaCondition="False", while the other two cases are for RuttaCondition="True"

k = sqrt(5);
options_fsolve = optimset('display','off');

% T-junction coordinates: 
% junction is between -2<x<2 and 0<y<3, diameter 1

% "1" matrices are 99x401; "2" matrices are 200x99
%load('../AnandFiles/TJunctionPotential_3-1a.mat'); % coarse solution for Z,Zeta
%load('TJunctionPotential_3-2.mat'); % fine solution for Z,Zeta

% Matteo: trying to see if I can get the same results by modifying
% TJunctionPotential_3.m so that I will not have to load the data already
[X1,Y1,Z1,X2,Y2,Z2,Zeta1,Zeta2] = TJunctionPotential_3;
ReZeta1 = real(Zeta1)'; % transpose to allow for concatenation
ReZeta2 = real(Zeta2);
ImZeta1 = imag(Zeta1)'; % transpose to allow for concatenation
ImZeta2 = imag(Zeta2);

% assign correct matrix (1 or 2) in which to find vortex based on its location
if (imag(zeta_vort) > 1) 
    ReZeta_vort = ReZeta2;
    ImZeta_vort = ImZeta2;
    Z_vort = Z2;
else
    ReZeta_vort = ReZeta1;
    ImZeta_vort = ImZeta1;
    Z_vort = Z1;
end;

% determine range of zeta values you have precomputed
ReZeta_min = ReZeta_vort(1,1); 
ReZeta_max = ReZeta_vort(1,end);
ImZeta_min = ImZeta_vort(1,1);
ImZeta_max = ImZeta_vort(end,1);

% compute vortex position in z-plane by interpolation or extrapolation
if (and(and(and(real(zeta_vort) > ReZeta_min, real(zeta_vort) < ReZeta_max), imag(zeta_vort) > ImZeta_min), imag(zeta_vort) < ImZeta_max))
    x_vort = interp2(ReZeta_vort,ImZeta_vort,real(Z_vort),real(zeta_vort),imag(zeta_vort));
    y_vort = interp2(ReZeta_vort,ImZeta_vort,imag(Z_vort),real(zeta_vort),imag(zeta_vort));
else % use extrapolation by spline
    x_vort = interp2(ReZeta_vort,ImZeta_vort,real(Z_vort),real(zeta_vort),imag(zeta_vort),'spline');
    y_vort = interp2(ReZeta_vort,ImZeta_vort,imag(Z_vort),real(zeta_vort),imag(zeta_vort),'spline');
end;
z_vort = x_vort + 1i*y_vort;
disp(['vortex z-position is ' num2str(z_vort)])

% Matteo: I added this to either fix gamma or compute it
if RuttaCondition=="True"
    % vortex circulation computed using Kutta condition at left corner
    gamma = abs(1 + z_vort)^2/imag(z_vort)*(u_minus/(1 - 1/k) - u_plus/(1 + 1/k)); 
elseif RuttaCondition=="False"
    gamma = 0.5;
else
    error("Must Specify RuttaCondition as either True or False")
end
disp(['gamma = ' num2str(gamma)])

% disp(['velocity at left corner is ' num2str(ComplexVelocity(-1,k,u_minus,u_plus,gamma,z_vort))])
% plot velocity along line (in z-plane) approaching left corner
% figure(cf)
% yv = 0:.00001:.001;
% plot(yv,abs(ComplexVelocity(-1+1i*yv,k,u_minus,u_plus,gamma,z_vort)),'.-k','MarkerSize',16)
% pause

%%% compute and plot stagnation points %%%
z_stag_boundary = 1/k*(u_minus + u_plus)/(u_minus - u_plus);
zeta_stag_boundary = TJunctionPotential_Zeta_Func(z_stag_boundary,k);
% disp(['stagnation point on boundary (without vortex) is at ' num2str(zeta_stag_boundary)])

% formula for complex velocity without prefactor as function of z
stag_func = @(z) -u_minus./(z + 1/k) + u_plus./(z - 1/k) - gamma*imag(z_vort)./abs(z - z_vort).^2;
% % formula with gamma-value substituted in
% % stag_func = @(z) -u_minus*(1./(z + 1/k) + 1/(1 - 1/k)*abs(1 + z_vort)^2./abs(z - z_vort).^2) + u_plus*(1./(z - 1/k) + 1/(1 + 1/k)*abs(1 + z_vort)^2./abs(z - z_vort).^2);

% METHOD 1: find stagnation point using fsolve
[z_stag_boundary_2,~,eflag] = fsolve(stag_func,1 + .1*1i,options_fsolve);
zeta_stag_boundary_2 = TJunctionPotential_Zeta_Func(z_stag_boundary_2,k);
disp(['eflag = ' num2str(eflag) '; stagnation point on boundary (with vortex) using fsolve is at zeta = ' num2str(zeta_stag_boundary_2) '; z = ' num2str(z_stag_boundary_2)])

% METHOD 2: find stagnation point by assuming its on right L-boundary
z_stag_vec = [1/k+.01:.01:300]; % 1:.01:300;
stag_func_vec = stag_func(z_stag_vec);

figure(cf+4)
title(['stagnation point identification'])
subplot(1,2,1); 
grid on;
plot(z_stag_vec,stag_func_vec,'.-k','MarkerSize',16)
xlabel('Re($z$)','interpreter','latex')
ylabel('velocity (Im($z$)=0)','interpreter','latex')

%%METHOD 3: find stagnation point as intersection of contours
[X_stag,Y_stag] = meshgrid(-2:pi/300:2);
ZMat_stag = X_stag + 1i*Y_stag;
subplot(1,2,2)
hold on
[M_stagR,c_LstagR] = contour(X_stag,Y_stag,real(stag_func(ZMat_stag)),[0;0]);
c_LstagR.LineWidth = 3;
c_LstagR.Color = 'r';
[M_stagI,c_LstagI] = contour(X_stag,Y_stag,imag(ZMat_stag),[0;0]);
c_LstagI.LineWidth = 3;
c_LstagI.Color = 'b';
hold off
xlabel('Re($z$)','interpreter','latex')
ylabel('Im($z$)','interpreter','latex')

% % use METHOD 2 to identify stagnation point
stag_loc = find(stag_func_vec(2:end).*stag_func_vec(1:end-1) < 0);
if (~isempty(stag_loc))
    stag_loc = stag_loc(1);
    z_stag_boundary_3 = z_stag_vec(stag_loc);
    zeta_stag_boundary_3 = TJunctionPotential_Zeta_Func(z_stag_boundary_3,k);
else
    zeta_stag_boundary_3 = NaN;
end;

% assemble zeta and z-coordinates by combining two branches of T-junction
size1 = size(Z1);
size2 = size(Z2);
num_col = (size1(2) - size2(2))/2;
ReZeta = [ReZeta1; NaN*ones(size2(1),num_col) ReZeta2 NaN*ones(size2(1),num_col)];
ImZeta = [ImZeta1; NaN*ones(size2(1),num_col) ImZeta2 NaN*ones(size2(1),num_col)];
Zeta = ReZeta + 1i*ImZeta;
Zeta1 = ReZeta1 + 1i*ImZeta1;
Zeta2 = ReZeta2 + 1i*ImZeta2;
Z = [Z1; NaN*ones(size2(1),num_col) Z2 NaN*ones(size2(1),num_col)];

% assemble x- and y-vectors from which Zeta matrices are derived
ImZeta1_vec = ImZeta1(:,1);
ReZeta1_vec = ReZeta1(1,:);
ImZeta2_vec = ImZeta2(:,1);
ReZeta2_vec = ReZeta2(1,:);

% complex velocities
CVel = ComplexVelocity(Z,k,u_minus,u_plus,gamma,z_vort);
U = real(CVel);
V = -imag(CVel);

CVel1 = ComplexVelocity(Z1,k,u_minus,u_plus,gamma,z_vort);
U1 = real(CVel1);
V1 = -imag(CVel1);

CVel2 = ComplexVelocity(Z2,k,u_minus,u_plus,gamma,z_vort);
U2 = real(CVel2);
V2 = -imag(CVel2);

% complex velocities in absence of vortex
CVel_free = ComplexVelocity(Z,k,u_minus,u_plus,0,z_vort);
U_free = real(CVel_free);
V_free = -imag(CVel_free);

CVel1_free = ComplexVelocity(Z1,k,u_minus,u_plus,0,z_vort);
U1_free = real(CVel1_free);
V1_free = -imag(CVel1_free);

CVel2_free = ComplexVelocity(Z2,k,u_minus,u_plus,0,z_vort);
U2_free = real(CVel2_free);
V2_free = -imag(CVel2_free);

% velocity potential
wpot1 = VelocityPotential(Z1,k,u_minus,u_plus,gamma,z_vort);
wpot2 = VelocityPotential(Z2,k,u_minus,u_plus,gamma,z_vort);
wpot = VelocityPotential(Z,k,u_minus,u_plus,gamma,z_vort);

% pressure coefficients according to Bernoulli equation
% assume pressure is zero at negative infinity on horizontal branch
P = 1 - abs(CVel).^2/u_minus^2;
P1 = 1 - abs(CVel1).^2/u_minus^2;
P2 = 1 - abs(CVel2).^2/u_minus^2;

% pressure coefficients in absence of vortex
P_free = 1 - abs(CVel_free).^2/u_minus^2;
P1_free = 1 - abs(CVel1_free).^2/u_minus^2;
P2_free = 1 - abs(CVel2_free).^2/u_minus^2;

% pressure differences assuming Bernoulli (so u^2-diff), averaged over pipe diameter
A1_loc = find(ReZeta1(end,:) == (-0.5-1.5));
A2_loc = find(ReZeta1(end,:) == (0.5+1.5));
A3_loc = find(ImZeta2(:,end) == (1+1.5));
SpeedSq1 = mean(abs(CVel1).^2,1);
SpeedSq2 = mean(abs(CVel2).^2,2);
p_A2A1 = 1/2*(SpeedSq1(A1_loc) - SpeedSq1(A2_loc));
p_A3A1 = 1/2*(SpeedSq1(A1_loc) - SpeedSq2(A3_loc));

% velocity potential differences
PhiAvg1 = real(mean(wpot1,1));
PhiAvg2 = real(mean(wpot2,2));
phi_A2A1 = PhiAvg1(A2_loc) - PhiAvg1(A1_loc);
phi_A3A1 = PhiAvg2(A3_loc) - PhiAvg1(A1_loc);

%%%%%%%% CONSTRUCT STREAMLINES %%%%%%%%
% if (u_minus > 0) % if flow is coming in from left horizontal branch
%     % start with streamlines originating from left horizontal branch
%     y_space = .025;
%     y_stream = y_space:y_space:(1-y_space);
%     figure(cf)
%     strm_line = streamline(ReZeta1,ImZeta1,U1,V1,-1.9*ones(1,length(y_stream)),y_stream);
% 
%     %%%% arrange streamlines into matrix format %%%
% 
%     % determine max number of points in streamline
%     num_pts = 0; % number of points in each streamline
%     for ind = 1:length(y_stream)
%         lth = length(strm_line(ind).XData);
%         if (lth > num_pts)
%             num_pts = lth;
%         end;
%     end;
% 
%     loc = []; % indices where streamline exits into vertical branch (rather than horizontal) 
%     x_stream = []; % x-coordinates where these streamlines exit
%     %strm_Mat_X = NaN*ones(length(y_stream),num_pts);
%     %strm_Mat_Y = NaN*ones(length(y_stream),num_pts);
%     for ind = 1:length(y_stream)
%         xs = strm_line(ind).XData;
%         ys = strm_line(ind).YData;
%         %strm_Mat_X(ind,1:length(xs)) = xs;
%         %strm_Mat_Y(ind,1:length(xs)) = ys;
% 
%         if (~(abs(xs - 2) < 0.1)) % detect if streamline does not exit right horizontal branch
%             loc = [loc ind];
%             x_stream = [x_stream xs(end)];
%         end;
%     end;
% 
%     % continue streamlines into top vertical branch
%     strm_line2 = streamline(ReZeta2,ImZeta2,U2,V2,x_stream,1.01*ones(1,length(x_stream)));
%     strm_line3 = streamline(ReZeta2,ImZeta2,U2,V2,-.49:.05:.49,1.1*ones(1,length(-.49:.05:.49)));
%     strm_line4 = streamline(ReZeta2,ImZeta2,U2,V2,-.49:.05:.49,1.5*ones(1,length(-.49:.05:.49)));
% else % if flow is going out from left horizontal branch
%     % start with streamlines originating from top vertical branch
%     y_space = .05;
%     y_stream = y_space:y_space:(1-y_space);
%     figure(cf)
%     strm_line = streamline(ReZeta1,ImZeta1,U1,V1,1.9*ones(1,length(y_stream)),y_stream);
%     
%     % then do streamlines originating from right horizontal branch
%     x_space = 0.05;
%     x_stream = (-.5+x_space):y_space:(.5-x_space);
%     strm_line2 = streamline(ReZeta2,ImZeta2,U2,V2,x_stream,2.9*ones(1,length(x_stream)));
%     
%     % then do streamlines near left corner
%     y_space = .025;
%     y_stream = .5:y_space:(1-y_space);
%     strm_line3 = streamline(ReZeta1,ImZeta1,U1,V1,-0.5*ones(1,length(y_stream)),y_stream);
% end;
% 
% hold on
% % construct lines showing T-junction boundary
% line([-0.5 -0.5],[1 3],'LineWidth',3,'Color','k')
% line([0.5 0.5],[1 3],'LineWidth',3,'Color','k')
% line([-2 -.5],[1 1],'LineWidth',3,'Color','k')
% line([.5 2],[1 1],'LineWidth',3,'Color','k')
% line([-2 2],[0 0],'LineWidth',3,'Color','k')
% plot(real(zeta_vort),imag(zeta_vort),'.k','MarkerSize',36)
% hold off
% xlabel('$x$','interpreter','latex')
% ylabel('$y$','interpreter','latex')
% set(gca,'FontName','Times','FontSize',24)
% title(['$u_- = $ ' num2str(u_minus) '; $u_+ = $ ' num2str(u_plus) '; $\gamma =$ ' num2str(gamma)],'interpreter','latex')

% figure(1)
% plot(strm_Mat_X',strm_Mat_Y','.-','MarkerSize',16)

%%%%%%%% STREAMLINE CONSTRUCTION AND PLOTTING OVER %%%%%%%%

%%%% CODE SNIPPET SOLVES STREAMLINE EQUATIONS IN ODE45: DOESN'T WORK WELL
% figure(cf+1); subplot(1,3,1); hold on; subplot(1,3,2); hold on; subplot(1,3,3); hold on;
% for ind = 1:length(y_stream)
%     [T,Xi] = ode45(@(t,xi) StreamRHS(t,xi,k,u_minus,u_plus,Z1,Z2,Zeta1,Zeta2),[0 50],[-1.5; y_stream(ind)]);
% 
%     figure(cf+1)
%     subplot(1,3,1); plot(T,Xi(:,1),'.-k','MarkerSize',16)
%     subplot(1,3,2); plot(T,Xi(:,2),'.-k','MarkerSize',16)
%     subplot(1,3,3); plot(Xi(:,1),Xi(:,2),'.-k','MarkerSize',16)
% end;
% subplot(1,3,1); hold off; subplot(1,3,2); hold off; subplot(1,3,3); hold off;
% subplot(1,3,3); xlim([-2 2]); ylim([0 3]);
% 
% plot streamlines using Im(w) = 0
figure(cf+1)
[M1,c1] = contour(ReZeta1,ImZeta1,imag(wpot1),50);
c1.LineWidth = 3;
hold on
[M2,c2] = contour(ReZeta2,ImZeta2,imag(wpot2),50); % use 25 if no vortex; 50 if vortex
%size(M2)
c2.LineWidth = 3;
% construct lines showing T-junction boundary
line([-0.5 -0.5],[1 3],'LineWidth',3,'Color','k')
line([0.5 0.5],[1 3],'LineWidth',3,'Color','k')
line([-2 -.5],[1 1],'LineWidth',3,'Color','k')
line([.5 2],[1 1],'LineWidth',3,'Color','k')
line([-2 2],[0 0],'LineWidth',3,'Color','k')
% plot vortex
plot(real(zeta_vort),imag(zeta_vort),'.k','MarkerSize',24)
% plot stagnation point
plot(real(zeta_stag_boundary_3),imag(zeta_stag_boundary_3),'.m','MarkerSize',24)
xlabel('$x$','interpreter','latex'); ylabel('$y$','interpreter','latex')
set(gca,'FontName','Times','FontSize',24)
title(['$u_- = $ ' num2str(u_minus) '; $u_+ = $ ' num2str(u_plus) '; $\gamma =$ ' num2str(gamma)],'interpreter','latex')
% 
% % %%%%% compute and plot stagnation streamline on streamline plot %%%%%
% % % THIS SNIPPET IS IRRELEVANT IF STAGNATION POINT IS FIXED TO BE AT CORNER
% % % determine stagnation point, if it exists
[Re_stag,Im_stag] = StagPointFind1(ReZeta2_vec,ImZeta2_vec,U2,V2);
zeta_stag = Re_stag + 1i*Im_stag;
zeta_stag = zeta_stag(abs(zeta_stag - zeta_vort) > 0.02); % remove artificial stagnation point detected at vortex position

if (~isempty(zeta_stag))

    % plot stagnation point on streamline plot
    figure(cf+1)
    hold on
    plot(real(zeta_vort),imag(zeta_vort),'.k','MarkerSize',36)
    plot(real(zeta_stag),imag(zeta_stag),'.r','MarkerSize',36)
    hold off

    % compute stagnation point in z-plane by interpolation
    x_stag = interp2(ReZeta2,ImZeta2,real(Z2),real(zeta_stag),imag(zeta_stag));
    y_stag = interp2(ReZeta2,ImZeta2,imag(Z2),real(zeta_stag),imag(zeta_stag));
    z_stag = x_stag + 1i*y_stag;

    % compute velocity potential at stagnation point
    w_stag = VelocityPotential(z_stag,k,u_minus,u_plus,gamma,z_vort);

    figure(cf+1)
    hold on
    [M_stag,c_stag] = contour(ReZeta2,ImZeta2,imag(wpot2),imag(w_stag)*[1;1]);
    c_stag.LineWidth = 3;
    c_stag.Color = 'r';
    hold off
end;
% 
% compute and plot streamline emerging from left corner
w_LTcorner = VelocityPotential(-1,k,u_minus,u_plus,gamma,z_vort);
% disp(['velocity potential at left corner is ' num2str(w_LTcorner)])
figure(cf+1)
hold on
[M_LTcorner1,c_LTcorner1] = contour(ReZeta1,ImZeta1,imag(wpot1),imag(w_LTcorner)*[1;1]);
c_LTcorner1.LineWidth = 3;
c_LTcorner1.Color = 'r';
[M_LTcorner2,c_LTcorner2] = contour(ReZeta2,ImZeta2,imag(wpot2),imag(w_LTcorner)*[1;1]);
c_LTcorner2.LineWidth = 3;
c_LTcorner2.Color = 'r';
hold off

% compute and plot streamline emerging from right corner
w_RTcorner = VelocityPotential(1,k,u_minus,u_plus,gamma,z_vort);
% disp(['velocity potential at right corner is ' num2str(w_RTcorner)])
figure(cf+1)
hold on
[M_RTcorner1,c_RTcorner1] = contour(ReZeta1,ImZeta1,imag(wpot1),imag(w_RTcorner)*[1;1]);
c_RTcorner1.LineWidth = 3;
c_RTcorner1.Color = 'm';
[M_RTcorner2,c_RTcorner2] = contour(ReZeta2,ImZeta2,imag(wpot2),imag(w_RTcorner)*[1;1]);
c_RTcorner2.LineWidth = 3;
c_RTcorner2.Color = 'm';
hold off

if (u_minus > 0) % hack to differentiate forward and reverse cycles
    disp(['vertical position of pink streamline is ' num2str(M_RTcorner1(2,2))]);
else
    disp(['vertical position of pink streamline is ' num2str(M_RTcorner1(2,end))]);
end;

% %%%%% plot zero-contours of velocity %%%%%
% figure(cf+2)
% [MU2,cu2] = contour(ReZeta2,ImZeta2,U2,[0;0]);
% cu2.LineWidth = 3; cu2.Color = 'b';
% hold on
% [MV2,cv2] = contour(ReZeta2,ImZeta2,V2,[0;0]);
% cv2.LineWidth = 3; cv2.Color = 'r';
% % construct lines showing T-junction boundary
% line([-0.5 -0.5],[1 3],'LineWidth',3,'Color','k')
% line([0.5 0.5],[1 3],'LineWidth',3,'Color','k')
% line([-2 -.5],[1 1],'LineWidth',3,'Color','k')
% line([.5 2],[1 1],'LineWidth',3,'Color','k')
% line([-2 2],[0 0],'LineWidth',3,'Color','k')
% plot(real(zeta_vort),imag(zeta_vort),'.k','MarkerSize',36)
% hold off
% xlabel('$x$','interpreter','latex')
% ylabel('$y$','interpreter','latex')
% set(gca,'FontName','Times','FontSize',24)
% title(['$u_- = $ ' num2str(u_minus) '; $u_+ = $ ' num2str(u_plus) '; $\gamma =$ ' num2str(gamma)],'interpreter','latex')
% 
% %%%%%%%%%%

%%%%% plot velocity components and magnitude in surface plot %%%%%
% figure(cf+3)
% 
% subplot(1,3,1)
% surf(ReZeta,ImZeta,U,'edgecolor','none'); view([0 90]); colorbar;
% xlabel('$x$','interpreter','latex'); ylabel('$y$','interpreter','latex')
% h = colorbar;
% ylabel(h,'$u$','interpreter','latex');
% set(gca,'FontName','Times','FontSize',24)
% 
% subplot(1,3,2)
% surf(ReZeta,ImZeta,V,'edgecolor','none'); view([0 90]); colorbar;
% xlabel('$x$','interpreter','latex'); ylabel('$y$','interpreter','latex')
% h = colorbar;
% ylabel(h,'$v$','interpreter','latex');
% set(gca,'FontName','Times','FontSize',24)
% 
% subplot(1,3,3)
% surf(ReZeta,ImZeta,sqrt(U.^2 + V.^2),'edgecolor','none'); view([0 90]); colorbar;
% xlabel('$x$','interpreter','latex'); ylabel('$y$','interpreter','latex')
% h = colorbar;
% ylabel(h,'$\sqrt{u^2+v^2}$','interpreter','latex');
% set(gca,'FontName','Times','FontSize',24)

%%% plot pressure coefficients averaged along pipe diamter %%%
% figure(cf+2)
% subplot(1,2,1)
% plot(ReZeta1_vec,mean(P1,1),'.-r','MarkerSize',16)
% hold on
% plot(ReZeta1_vec,mean(P1_free,1),'.-k','MarkerSize',16)
% hold off
% xlabel('Re($\zeta$)','interpreter','latex')
% ylabel('$\langle C_p\rangle_y$','interpreter','latex')
% set(gca,'FontName','Times','FontSize',24)
% set(gca,'XTick',-2:.5:2)
% grid on
% 
% ind_neg = find(ReZeta1(end,:) == -.5);
% ind_pos = find(ReZeta1(end,:) == 0.5);
% subplot(1,2,2)
% plot([ImZeta1_vec; ImZeta2_vec],[mean(P1(:,ind_neg+1:ind_pos),2); mean(P2,2)],'.-r','MarkerSize',16)
% hold on
% plot([ImZeta1_vec; ImZeta2_vec],[mean(P1_free(:,ind_neg+1:ind_pos),2); mean(P2_free,2)],'.-k','MarkerSize',16)
% xlabel('Im($\zeta$)','interpreter','latex')
% ylabel('$\langle C_p\rangle_x$','interpreter','latex')
% set(gca,'FontName','Times','FontSize',24)
% xlim([0 3])
% set(gca,'XTick',0:.5:3)
% grid on

%%%%% plot streamfunction = Im(w) and Re(w) as colored surfaces %%%%% 
figure(cf+5)
subplot(1,2,1); surf(ReZeta,ImZeta,imag(wpot),'edgecolor','none'); view([0 90]); colorbar;
xlabel('$x$','interpreter','latex'); ylabel('$y$','interpreter','latex')
h = colorbar;
ylabel(h,'Im($w$)','interpreter','latex');
set(gca,'FontName','Times','FontSize',24)

subplot(1,2,2); surf(ReZeta,ImZeta,real(wpot),'edgecolor','none'); view([0 90]); colorbar;
xlabel('$x$','interpreter','latex'); ylabel('$y$','interpreter','latex')
h = colorbar;
ylabel(h,'Re($w$)','interpreter','latex');
set(gca,'FontName','Times','FontSize',24)

%%%%% plot complex velocity on boundary to identify stagnation points %%%%%
% z_vec1 = (-1/k+.01):.01:(1/k-.01); % lower horizontal boundary of junction
% z_vec2 = (-1+.01):.01:(-1/k-.01); % left horizontal boundary
% z_vec3 = (1/k+.01):.01:(1-.01); % right horizontal boundary
% z_vec4 = -5:.01:(-1-.01); % left vertical boundary
% z_vec5 = (1+.01):.01:5; % right vertical boundary
% 
% % calculate corresponding zeta-values
% zeta_vec1 = Zeta_Vec_Calc(z_vec1,k);
% zeta_vec2 = Zeta_Vec_Calc(z_vec2,k);
% zeta_vec3 = Zeta_Vec_Calc(z_vec3,k);
% zeta_vec4 = Zeta_Vec_Calc(z_vec4,k);
% zeta_vec5 = Zeta_Vec_Calc(z_vec5,k);
% 
% % calculate corresponding complex velocity
% vel1 = ComplexVelocity(z_vec1,k,u_minus,u_plus,gamma,z_vort);
% vel2 = ComplexVelocity(z_vec2,k,u_minus,u_plus,gamma,z_vort);
% vel3 = ComplexVelocity(z_vec3,k,u_minus,u_plus,gamma,z_vort);
% vel4 = ComplexVelocity(z_vec4,k,u_minus,u_plus,gamma,z_vort);
% vel5 = ComplexVelocity(z_vec5,k,u_minus,u_plus,gamma,z_vort);

% figure(cf+3)
% subplot(1,5,1)
% plot(real(zeta_vec1),real(vel1),'.-k','MarkerSize',16)
% xlabel('Re($\zeta$)','interpreter','latex')
% ylabel('$u$','interpreter','latex')
% set(gca,'FontName','Times','FontSize',24)
% title(['bottom horiz. boundary'])
% 
% subplot(1,5,2)
% plot(real(zeta_vec2),real(vel2),'.-k','MarkerSize',16)
% xlabel('Re($\zeta$)','interpreter','latex')
% ylabel('$u$','interpreter','latex')
% set(gca,'FontName','Times','FontSize',24)
% title(['left horiz. boundary'])
% 
% subplot(1,5,3)
% plot(real(zeta_vec3),real(vel3),'.-k','MarkerSize',16)
% xlabel('Re($\zeta$)','interpreter','latex')
% ylabel('$u$','interpreter','latex')
% set(gca,'FontName','Times','FontSize',24)
% title(['right horiz. boundary'])
% 
% subplot(1,5,4)
% plot(imag(zeta_vec4),-imag(vel4),'.-k','MarkerSize',16)
% xlabel('Im($\zeta$)','interpreter','latex')
% ylabel('$v$','interpreter','latex')
% set(gca,'FontName','Times','FontSize',24)
% title(['left vert. boundary'])
% 
% subplot(1,5,5)
% plot(imag(zeta_vec5),-imag(vel5),'.-k','MarkerSize',16)
% xlabel('Im($\zeta$)','interpreter','latex')
% ylabel('$v$','interpreter','latex')
% set(gca,'FontName','Times','FontSize',24)
% title(['right vert. boundary'])

function out = StreamRHS(t,xi,k,u_minus,u_plus,Z1,Z2,Zeta1,Zeta2)

% find correct z-value by interpolation
if (xi(2) < 1) % if coordinate is in horizontal branch of T-junction
    Zeta = Zeta1;
    Z = Z1;
else % if coordinate is in vertical branch of T-junction
    Zeta = Zeta2;
    Z = Z2;
end;
x = interp2(real(Zeta),imag(Zeta),real(Z),xi(1),xi(2));
y = interp2(real(Zeta),imag(Zeta),imag(Z),xi(1),xi(2));
z = x + 1i*y;

CVel = ComplexVelocity(z,k,u_minus,u_plus,0,NaN);
u = real(CVel);
v = -imag(CVel);
out = [u; v];

function out = VelocityPotential(Z,k,u_minus,u_plus,gamma,z_vort)

out = 1/pi*(u_minus*log((Z + 1/k)) - u_plus*log(Z - 1/k)) - 1i*gamma/(2*pi)*log((Z - z_vort)./(Z - conj(z_vort)));
%out = 1/pi*(u_minus*log((Z + 1/k)) - u_plus*log(Z - 1/k)) - 1i*gamma/(2*pi)*(log(Z - z_vort) - log(Z - conj(z_vort)));

function out = ComplexVelocity(Z,k,u_minus,u_plus,gamma,z_vort) % velocity u-i*v = dw/dzeta

% term1 is "velocity" term (no vortex)
term1 = (-u_minus*(Z - 1/k) + u_plus*(Z + 1/k))./sqrt(1 - Z.^2);

% term2 incluces the vortex
term2 = -1i*gamma/(2*k^2)*(1 - (k*Z).^2)./sqrt(1 - Z.^2).*(1./(Z - z_vort) - 1./(Z - conj(z_vort)));

out = term1 + term2;

function out = Zeta_Vec_Calc(zvec,k)

out = NaN*ones(1,length(zvec));
for ind = 1:length(zvec)
    z = zvec(ind);
    out(ind) = TJunctionPotential_Zeta_Func(z,k);
end;