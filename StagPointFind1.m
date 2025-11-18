function [Re_stag,Im_stag] = StagPointFind1(ReZeta,ImZeta,U,V)
% computes stagnation points of 2D flow defined by (U,V) by computing
% ... intersections of zero-level curves of U and V

% INPUTS: ReZeta and ImZeta are vectors of the x- and y-coordinates
% U and V are the two velocity components
% OUTPUTS: Re_stag and Im_stag are the real- and imaginary-parts (or x-
% ... and y-coordinates) of the stagnation points

addpath('/Users/Anand/Dropbox/Nematic'); % folder containing intersections_dnld.m
Re_stag = [];
Im_stag = [];

MU = contourc(ReZeta,ImZeta,U,[0;0]);
MV = contourc(ReZeta,ImZeta,V,[0;0]);

%size(MU)
%size(MV)
%MU(:,1:5)
%MV(:,1:5)
%MV(:,30:34)

contour_lengths_U = contour_length(MU);
contour_lengths_V = contour_length(MV);

%contour_lengths_U
%contour_lengths_V

start_U = 2;
for ind_U = 1:length(contour_lengths_U)
    start_V = 2;
    end_U = start_U + contour_lengths_U(ind_U) - 1;
    for ind_V = 1:length(contour_lengths_V)
        end_V = start_V + contour_lengths_V(ind_V) - 1;
        U_curve = MU(:,start_U:end_U);
        V_curve = MV(:,start_V:end_V);
        [re_stag,im_stag] = intersections_dnld(U_curve(1,:),U_curve(2,:),V_curve(1,:),V_curve(2,:));
        start_V = end_V + 2;
        Re_stag = [Re_stag re_stag];
        Im_stag = [Im_stag im_stag];
    end;
    start_U = end_U + 2;
end;

function out = contour_length(M)
% computes lengths of each disconnected curve comprising contour
% INPUT is matrix M from contourc function; see Matlab help to see
% how data in M is arranged.
% OUTPUT is a vector of these lengths

out = []; % vector showing lengths of each disconnected curve comprising contour
ind = 1;
while (ind <= length(M(1,:)))
    out = [out M(2,ind)];
    ind = ind + M(2,ind) + 1;
end;