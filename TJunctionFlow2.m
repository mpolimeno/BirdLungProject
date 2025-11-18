function TJunctionFlow2(cf)

% potential flow solution for flow inside T-junction with two vortices

k = sqrt(5);

% vortex (complex) position; used 0 + 1i*pi/2 before
% configuration for "forward" phase of cycle
%zeta_vort_minus = -.25 + 1i*1.25; % vortex in vertical branch
%zeta_vort_plus = 0.75 + 0.75*1i; % vortex in right horizontal branch

% configuration for "reverse" phase of cycle
%Matteo % Fig 12b
%zeta_vort_minus = -0.75 + 1i*0.75; % vortex in left horizontal branch
%zeta_vort_plus = 0.25 + 1i*0.75; % vortex right-of-center in horizontal branch

% Matteo % Fig 12a
zeta_vort_minus = -0.25 + 1i*1.25; % vortex in left horizontal branch
zeta_vort_plus = 0.75 + 1i*0.75; % vortex right-of-center in horizontal branch

u_minus = 1; % inflow velocity on left horizontal branch
u_plus = 0.75; % outflow velocity on right horizontal branch

% T-junction coordinates: junction is between -2<x<2 and 0<y<3, diameter 1
% "1" matrices are 99x401; "2" matrices are 200x99
load('../AnandFiles/TJunctionPotential_3-2.mat')

% assign correct matrix (1 or 2) in which to find vortex based on its location
if (imag(zeta_vort_minus) > 1) 
    ReZeta_vort_minus = ReZeta2;
    ImZeta_vort_minus = ImZeta2;
    Z_vort_minus = Z2;
else
    ReZeta_vort_minus = ReZeta1;
    ImZeta_vort_minus = ImZeta1;
    Z_vort_minus = Z1;
end;

% do the same for zeta_vort_plus
if (imag(zeta_vort_plus) > 1) 
    ReZeta_vort_plus = ReZeta2;
    ImZeta_vort_plus = ImZeta2;
    Z_vort_plus = Z2;
else
    ReZeta_vort_plus = ReZeta1;
    ImZeta_vort_plus = ImZeta1;
    Z_vort_plus = Z1;
end;

% compute vortex positions in z-plane by interpolation
x_vort_minus = interp2(ReZeta_vort_minus,ImZeta_vort_minus,real(Z_vort_minus),real(zeta_vort_minus),imag(zeta_vort_minus));
y_vort_minus = interp2(ReZeta_vort_minus,ImZeta_vort_minus,imag(Z_vort_minus),real(zeta_vort_minus),imag(zeta_vort_minus));
z_vort_minus = x_vort_minus + 1i*y_vort_minus;
disp(['vortex z-position is ' num2str(z_vort_minus)])
x_vort_plus = interp2(ReZeta_vort_plus,ImZeta_vort_plus,real(Z_vort_plus),real(zeta_vort_plus),imag(zeta_vort_plus));
y_vort_plus = interp2(ReZeta_vort_plus,ImZeta_vort_plus,imag(Z_vort_plus),real(zeta_vort_plus),imag(zeta_vort_plus));
z_vort_plus = x_vort_plus + 1i*y_vort_plus;
disp(['vortex z-position is ' num2str(z_vort_plus)])

% vortex circulation computed using Kutta condition at both sharp corners
Mat(1,1) = imag(z_vort_minus)/abs(-1-z_vort_minus)^2;
Mat(1,2) = imag(z_vort_plus)/abs(-1-z_vort_plus)^2;
Mat(2,1) = imag(z_vort_minus)/abs(1-z_vort_minus)^2;
Mat(2,2) = imag(z_vort_plus)/abs(1-z_vort_plus)^2;
vel_delta = u_minus - u_plus;
vel_sigma = u_minus + u_plus;
rhs_vec = 1/(1 - 1/k^2)*[vel_delta + vel_sigma/k; -vel_delta + vel_sigma/k];
gamma_vec = inv(Mat)*rhs_vec;
gamma_minus = gamma_vec(1);
gamma_plus = gamma_vec(2);

z_stag_boundary = 1/k*(u_minus + u_plus)/(u_minus - u_plus);
zeta_stag_boundary = TJunctionPotential_Zeta_Func(z_stag_boundary,k);
disp(['stagnation point on boundary (without vortex) is at ' num2str(zeta_stag_boundary)])

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
CVel = ComplexVelocity(Z,k,u_minus,u_plus,gamma_minus,z_vort_minus,gamma_plus,z_vort_plus);
U = real(CVel);
V = -imag(CVel);

CVel1 = ComplexVelocity(Z1,k,u_minus,u_plus,gamma_minus,z_vort_minus,gamma_plus,z_vort_plus);
U1 = real(CVel1);
V1 = -imag(CVel1);

CVel2 = ComplexVelocity(Z2,k,u_minus,u_plus,gamma_minus,z_vort_minus,gamma_plus,z_vort_plus);
U2 = real(CVel2);
V2 = -imag(CVel2);

% complex velocities in absence of vortex
CVel_free = ComplexVelocity(Z,k,u_minus,u_plus,0,z_vort_minus,0,z_vort_plus);
U_free = real(CVel_free);
V_free = -imag(CVel_free);

CVel1_free = ComplexVelocity(Z1,k,u_minus,u_plus,0,z_vort_minus,0,z_vort_plus);
U1_free = real(CVel1_free);
V1_free = -imag(CVel1_free);

CVel2_free = ComplexVelocity(Z2,k,u_minus,u_plus,0,z_vort_minus,0,z_vort_plus);
U2_free = real(CVel2_free);
V2_free = -imag(CVel2_free);

% velocity potential
wpot1 = VelocityPotential(Z1,k,u_minus,u_plus,gamma_minus,z_vort_minus,gamma_plus,z_vort_plus);
wpot2 = VelocityPotential(Z2,k,u_minus,u_plus,gamma_minus,z_vort_minus,gamma_plus,z_vort_plus);
wpot = VelocityPotential(Z,k,u_minus,u_plus,gamma_minus,z_vort_minus,gamma_plus,z_vort_plus);

% pressure coefficients according to Bernoulli equation
% assume pressure is zero at negative infinity on horizontal branch
P = 1 - abs(CVel).^2/u_minus^2;
P1 = 1 - abs(CVel1).^2/u_minus^2;
P2 = 1 - abs(CVel2).^2/u_minus^2;

% pressure coefficients in absence of vortex
P_free = 1 - abs(CVel_free).^2/u_minus^2;
P1_free = 1 - abs(CVel1_free).^2/u_minus^2;
P2_free = 1 - abs(CVel2_free).^2/u_minus^2;

%%%%%%%% CONSTRUCT STREAMLINES %%%%%%%%
if (u_minus > 0) % if flow is coming in from left horizontal branch
    % start with streamlines originating from left horizontal branch
    y_space = .025;
    y_stream = y_space:y_space:(1-y_space);
    figure(cf)
    strm_line = streamline(ReZeta1,ImZeta1,U1,V1,-1.9*ones(1,length(y_stream)),y_stream);

    %%%% arrange streamlines into matrix format %%%

    % determine max number of points in streamline
    num_pts = 0; % number of points in each streamline
    for ind = 1:length(y_stream)
        lth = length(strm_line(ind).XData);
        if (lth > num_pts)
            num_pts = lth;
        end;
    end;

    loc = []; % indices where streamline exits into vertical branch (rather than horizontal) 
    x_stream = []; % x-coordinates where these streamlines exit
    %strm_Mat_X = NaN*ones(length(y_stream),num_pts);
    %strm_Mat_Y = NaN*ones(length(y_stream),num_pts);
    for ind = 1:length(y_stream)
        xs = strm_line(ind).XData;
        ys = strm_line(ind).YData;
        %strm_Mat_X(ind,1:length(xs)) = xs;
        %strm_Mat_Y(ind,1:length(xs)) = ys;

        if (~(abs(xs - 2) < 0.1)) % detect if streamline does not exit right horizontal branch
            loc = [loc ind];
            x_stream = [x_stream xs(end)];
        end;
    end;

    % continue streamlines into top vertical branch
    strm_line2 = streamline(ReZeta2,ImZeta2,U2,V2,x_stream,1.01*ones(1,length(x_stream)));
    strm_line3 = streamline(ReZeta2,ImZeta2,U2,V2,-.49:.05:.49,1.1*ones(1,length(-.49:.05:.49)));
    strm_line4 = streamline(ReZeta2,ImZeta2,U2,V2,-.49:.05:.49,1.5*ones(1,length(-.49:.05:.49)));
else % if flow is going out from left horizontal branch
    % start with streamlines originating from top vertical branch
    y_space = .05;
    y_stream = y_space:y_space:(1-y_space);
    figure(cf)
    strm_line = streamline(ReZeta1,ImZeta1,U1,V1,1.9*ones(1,length(y_stream)),y_stream);
    
    % then do streamlines originating from right horizontal branch
    x_space = 0.05;
    x_stream = (-.5+x_space):y_space:(.5-x_space);
    strm_line2 = streamline(ReZeta2,ImZeta2,U2,V2,x_stream,2.9*ones(1,length(x_stream)));
    
    % then do streamlines near left corner
    y_space = .025;
    y_stream = .5:y_space:(1-y_space);
    strm_line3 = streamline(ReZeta1,ImZeta1,U1,V1,-0.5*ones(1,length(y_stream)),y_stream);
end;

hold on
% construct lines showing T-junction boundary
line([-0.5 -0.5],[1 3],'LineWidth',3,'Color','k')
line([0.5 0.5],[1 3],'LineWidth',3,'Color','k')
line([-2 -.5],[1 1],'LineWidth',3,'Color','k')
line([.5 2],[1 1],'LineWidth',3,'Color','k')
line([-2 2],[0 0],'LineWidth',3,'Color','k')
plot(real(zeta_vort_minus),imag(zeta_vort_minus),'.k','MarkerSize',36)
plot(real(zeta_vort_plus),imag(zeta_vort_plus),'.k','MarkerSize',36)
hold off
xlabel('$x$','interpreter','latex')
ylabel('$y$','interpreter','latex')
set(gca,'FontName','Times','FontSize',24)
title(['$u_- = $ ' num2str(u_minus) '; $u_+ = $ ' num2str(u_plus) '; $\gamma_- =$ ' num2str(gamma_minus) '; $\gamma_+ =$ ' num2str(gamma_plus)],'interpreter','latex')

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
xlabel('$x$','interpreter','latex'); ylabel('$y$','interpreter','latex')
set(gca,'FontName','Times','FontSize',24)
title(['$u_- = $ ' num2str(u_minus) '; $u_+ = $ ' num2str(u_plus) '; $\gamma_- =$ ' num2str(gamma_minus) '; $\gamma_+ =$ ' num2str(gamma_plus)],'interpreter','latex')

% %%%%% compute and plot stagnation streamline on streamline plot %%%%%
% % THIS SNIPPET IS IRRELEVANT IF STAGNATION POINT IS FIXED TO BE AT CORNER
% % determine stagnation point, if it exists
% [Re_stag,Im_stag] = StagPointFind1(ReZeta2_vec,ImZeta2_vec,U2,V2);
% zeta_stag = Re_stag + 1i*Im_stag;
% zeta_stag = zeta_stag(abs(zeta_stag - zeta_vort) > 0.02); % remove artificial stagnation point detected at vortex position
% 
% if (~isempty(zeta_stag))
%     
%     % plot stagnation point on streamline plot
%     figure(cf+1)
%     hold on
%     plot(real(zeta_vort),imag(zeta_vort),'.k','MarkerSize',36)
%     plot(real(zeta_stag),imag(zeta_stag),'.r','MarkerSize',36)
%     hold off
%     
%     % compute stagnation point in z-plane by interpolation
%     x_stag = interp2(ReZeta2,ImZeta2,real(Z2),real(zeta_stag),imag(zeta_stag));
%     y_stag = interp2(ReZeta2,ImZeta2,imag(Z2),real(zeta_stag),imag(zeta_stag));
%     z_stag = x_stag + 1i*y_stag;
%     
%     % compute velocity potential at stagnation point
%     w_stag = VelocityPotential(z_stag,k,u_minus,u_plus,gamma,z_vort);
% 
%     figure(cf+1)
%     hold on
%     [M_stag,c_stag] = contour(ReZeta2,ImZeta2,imag(wpot2),imag(w_stag)*[1;1]);
%     c_stag.LineWidth = 3;
%     c_stag.Color = 'r';
%     hold off
% end;

% compute and plot streamline emerging from left corner
w_LTcorner = VelocityPotential(-1,k,u_minus,u_plus,gamma_minus,z_vort_minus,gamma_plus,z_vort_plus);
disp(['velocity potential at left corner is ' num2str(w_LTcorner)])
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
w_RTcorner = VelocityPotential(1,k,u_minus,u_plus,gamma_minus,z_vort_minus,gamma_plus,z_vort_plus);
disp(['velocity potential at right corner is ' num2str(w_RTcorner)])
figure(cf+1)
hold on
[M_RTcorner1,c_RTcorner1] = contour(ReZeta1,ImZeta1,imag(wpot1),imag(w_RTcorner)*[1;1]);
c_RTcorner1.LineWidth = 3;
c_RTcorner1.Color = 'm';
[M_RTcorner2,c_RTcorner2] = contour(ReZeta2,ImZeta2,imag(wpot2),imag(w_RTcorner)*[1;1]);
c_RTcorner2.LineWidth = 3;
c_RTcorner2.Color = 'm';

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
figure(cf+4)

subplot(2,4,1)
surf(ReZeta,ImZeta,U_free,'edgecolor','none'); view([0 90]); colorbar;
xlabel('$x$','interpreter','latex'); ylabel('$y$','interpreter','latex')
h = colorbar;
ylabel(h,'$u$','interpreter','latex');
set(gca,'FontName','Times','FontSize',24)

subplot(2,4,2)
surf(ReZeta,ImZeta,V_free,'edgecolor','none'); view([0 90]); colorbar;
xlabel('$x$','interpreter','latex'); ylabel('$y$','interpreter','latex')
h = colorbar;
ylabel(h,'$v$','interpreter','latex');
set(gca,'FontName','Times','FontSize',24)

subplot(2,4,3)
surf(ReZeta,ImZeta,sqrt(U_free.^2 + V_free.^2),'edgecolor','none'); view([0 90]); colorbar;
xlabel('$x$','interpreter','latex'); ylabel('$y$','interpreter','latex')
h = colorbar;
ylabel(h,'$\sqrt{u^2+v^2}$','interpreter','latex');
set(gca,'FontName','Times','FontSize',24)

subplot(2,4,4)
surf(ReZeta,ImZeta,P_free,'edgecolor','none'); view([0 90]); colorbar;
xlabel('$x$','interpreter','latex'); ylabel('$y$','interpreter','latex')
h = colorbar;
ylabel(h,'$C_p$','interpreter','latex');
set(gca,'FontName','Times','FontSize',24)

subplot(2,4,5)
surf(ReZeta,ImZeta,U,'edgecolor','none'); view([0 90]); colorbar;
xlabel('$x$','interpreter','latex'); ylabel('$y$','interpreter','latex')
h = colorbar;
ylabel(h,'$u$','interpreter','latex');
set(gca,'FontName','Times','FontSize',24)

subplot(2,4,6)
surf(ReZeta,ImZeta,V,'edgecolor','none'); view([0 90]); colorbar;
xlabel('$x$','interpreter','latex'); ylabel('$y$','interpreter','latex')
h = colorbar;
ylabel(h,'$v$','interpreter','latex');
set(gca,'FontName','Times','FontSize',24)

subplot(2,4,7)
surf(ReZeta,ImZeta,sqrt(U.^2 + V.^2),'edgecolor','none'); view([0 90]); colorbar;
xlabel('$x$','interpreter','latex'); ylabel('$y$','interpreter','latex')
h = colorbar;
ylabel(h,'$\sqrt{u^2+v^2}$','interpreter','latex');
set(gca,'FontName','Times','FontSize',24)

subplot(2,4,8)
surf(ReZeta,ImZeta,P,'edgecolor','none'); view([0 90]); colorbar;
xlabel('$x$','interpreter','latex'); ylabel('$y$','interpreter','latex')
h = colorbar;
ylabel(h,'$C_p$','interpreter','latex');
set(gca,'FontName','Times','FontSize',24)

% plot line-averaged values of pressure coefficients
figure(cf+2)
subplot(2,2,1)
plot(ReZeta1_vec,mean(P1,1),'.-r','MarkerSize',16)
hold on
plot(ReZeta1_vec,mean(P1_free,1),'.-k','MarkerSize',16)
hold off
xlabel('Re($\zeta$)','interpreter','latex')
ylabel('$\langle C_p\rangle_y$','interpreter','latex')
set(gca,'FontName','Times','FontSize',24)
set(gca,'XTick',-2:.5:2)
grid on

ind_neg = find(ReZeta1(end,:) == -.5);
ind_pos = find(ReZeta1(end,:) == 0.5);
subplot(2,2,2)
plot([ImZeta1_vec; ImZeta2_vec],[mean(P1(:,ind_neg+1:ind_pos),2); mean(P2,2)],'.-r','MarkerSize',16)
hold on
plot([ImZeta1_vec; ImZeta2_vec],[mean(P1_free(:,ind_neg+1:ind_pos),2); mean(P2_free,2)],'.-k','MarkerSize',16)
xlabel('Im($\zeta$)','interpreter','latex')
ylabel('$\langle C_p\rangle_x$','interpreter','latex')
set(gca,'FontName','Times','FontSize',24)
xlim([0 3])
set(gca,'XTick',0:.5:3)
grid on

% subplot(2,2,3)
% plot(ImZeta1_vec,P1_free(:,1:20:end),'.-','MarkerSize',16)
% xlabel('Im($\zeta$)','interpreter','latex')
% ylabel('$C_p$','interpreter','latex')
% set(gca,'FontName','Times','FontSize',24)
% xlim([0 1])
% set(gca,'XTick',0:.25:1)
% grid on
% 
% subplot(2,2,4)
% plot([ReZeta1_vec(152:250)],[P1_free(1:10:end,152:250); P2_free(1:10:end,:)],'.-','MarkerSize',16)
% xlabel('Re($\zeta$)','interpreter','latex')
% ylabel('$C_p$','interpreter','latex')
% set(gca,'FontName','Times','FontSize',24)
% xlim([-.5 .5])
% set(gca,'XTick',-.5:.25:.5)
% grid on

%%%%% plot streamfunction = Im(w) as a colored surface %%%%% 
% figure(cf+5)
% surf(ReZeta,ImZeta,imag(wpot),'edgecolor','none'); view([0 90]); colorbar;
% xlabel('$x$','interpreter','latex'); ylabel('$y$','interpreter','latex')
% h = colorbar;
% ylabel(h,'Im($w$)','interpreter','latex');
% set(gca,'FontName','Times','FontSize',24)

% plot complex velocity on boundary to identify stagnation points
z_vec1 = (-1/k+.001):.001:(1/k-.001); % lower horizontal boundary of junction
z_vec2 = (-1+.001):.001:(-1/k-.001); % left horizontal boundary
z_vec3 = (1/k+.001):.001:(1-.001); % right horizontal boundary
z_vec4 = -5:.001:(-1-.001); % left vertical boundary
z_vec5 = (1+.001):.001:5; % right vertical boundary

% calculate corresponding zeta-values
zeta_vec1 = Zeta_Vec_Calc(z_vec1,k);
zeta_vec2 = Zeta_Vec_Calc(z_vec2,k);
zeta_vec3 = Zeta_Vec_Calc(z_vec3,k);
zeta_vec4 = Zeta_Vec_Calc(z_vec4,k);
zeta_vec5 = Zeta_Vec_Calc(z_vec5,k);

% calculate corresponding complex velocity
vel1 = ComplexVelocity(z_vec1,k,u_minus,u_plus,gamma_minus,z_vort_minus,gamma_plus,z_vort_plus);
vel2 = ComplexVelocity(z_vec2,k,u_minus,u_plus,gamma_minus,z_vort_minus,gamma_plus,z_vort_plus);
vel3 = ComplexVelocity(z_vec3,k,u_minus,u_plus,gamma_minus,z_vort_minus,gamma_plus,z_vort_plus);
vel4 = ComplexVelocity(z_vec4,k,u_minus,u_plus,gamma_minus,z_vort_minus,gamma_plus,z_vort_plus);
vel5 = ComplexVelocity(z_vec5,k,u_minus,u_plus,gamma_minus,z_vort_minus,gamma_plus,z_vort_plus);

figure(cf+3)
subplot(1,5,1)
plot(real(zeta_vec1),real(vel1),'.-k','MarkerSize',16)
title(['bottom horizontal'])
subplot(1,5,2)
plot(real(zeta_vec2),real(vel2),'.-k','MarkerSize',16)
title(['left horizontal'])
subplot(1,5,3)
plot(real(zeta_vec3),real(vel3),'.-k','MarkerSize',16)
title(['right horizontal'])
subplot(1,5,4)
plot(imag(zeta_vec4),-imag(vel4),'.-k','MarkerSize',16)
title(['left vertical'])
subplot(1,5,5)
plot(imag(zeta_vec5),-imag(vel5),'.-k','MarkerSize',16)
title(['right vertical'])

function out = VelocityPotential(Z,k,u_minus,u_plus,gamma_minus,z_vort_minus,gamma_plus,z_vort_plus)

out = 1/pi*(u_minus*log((Z + 1/k)) - u_plus*log(Z - 1/k)) - 1i*gamma_minus/(2*pi)*log((Z - z_vort_minus)./(Z - conj(z_vort_minus))) - 1i*gamma_plus/(2*pi)*log((Z - z_vort_plus)./(Z - conj(z_vort_plus)));

function out = ComplexVelocity(Z,k,u_minus,u_plus,gamma_minus,z_vort_minus,gamma_plus,z_vort_plus) % velocity u-i*v = dw/dzeta

% term1 is "velocity" term (no vortex)
term1 = (-u_minus*(Z - 1/k) + u_plus*(Z + 1/k))./sqrt(1 - Z.^2);

% term2 incluces the vortex
term2_minus = -1i*gamma_minus/(2*k^2)*(1 - (k*Z).^2)./sqrt(1 - Z.^2).*(1./(Z - z_vort_minus) - 1./(Z - conj(z_vort_minus)));
term2_plus = -1i*gamma_plus/(2*k^2)*(1 - (k*Z).^2)./sqrt(1 - Z.^2).*(1./(Z - z_vort_plus) - 1./(Z - conj(z_vort_plus)));

out = term1 + term2_minus + term2_plus;

function out = Zeta_Vec_Calc(zvec,k)

out = NaN*ones(1,length(zvec));
for ind = 1:length(zvec)
    z = zvec(ind);
    out(ind) = TJunctionPotential_Zeta_Func(z,k);
end;