%
% Data in the Trajs variable
%
%%
clear 

setcctpath() 
nK = 10;
nBoot = 15;
nLhood = zeros(nK,nBoot);
nspread = zeros(nK,nBoot);

for K = 1:nK
for iBoot = 1:nBoot
disp([K,iBoot])
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% (1) load the .xls data into cell forWPmat
f = readtable('track_PostFiltering.ERAI.IC.xlsx');
track = f{:,:};

track_lon = track(:,2);
track_lat = track(:,3);

% id 
id = track(:,1);
% Number of trajectories [n_uni_track]
uni_track = unique(track(:,1));
n_uni_track = size(uni_track);
n_uni_track = n_uni_track(1,1);

uni_track_random = randsample(uni_track,n_uni_track);
uni_track_random = uni_track_random(:,1);

Trajs.Y = {};
for i = 1:n_uni_track
  itrack =  find( id==uni_track_random(i) );
  Trajs.Y(i,:) = { [track_lon(itrack), track_lat(itrack)] } ;
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% lrm_model
% (2) Set the model settings
ops.zero    = 'nozero';    % see CELL2ARRAY for possible values
ops.method  = 'lrm';       % see above for possible values
ops.K       = K;           % number of clusters
ops.order   = 2;           % 1-linear, 2-quadratic, 3-cubic, etc.
    
%
% (3) Set the general EM settings
ops.NumEMStarts    = 1;  % number of times to run EM (best model chosen)
ops.IterLimit      = 20; % maximum allowable iterations for any one EM start
ops.ValidateLhood  = 0;  % run a large-sample Lhood calculation at end
                         % of each EM start; used for time-based 
                         % (approximate) alignment. This will give a more
                         % accurate approximation of the log-likelihood.
                         % You only need this if you are running comparison
                         % experiments between different models based on
                         % log likelihood scores.
ops.MsgHnd         = -1;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Cluster the data according to the above settings
% set the warning state
% warning off MATLAB:singularMatrix;
warning off all;
warning on MATLAB:divideByZero;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% perform the clustering
model(K,iBoot) = curve_clust(Trajs,ops);
nLhood(K,iBoot) = model(K,iBoot).TrainLhood_ppt;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% calculate the mean regression trajectories 
t = 0:1:200;
[P, K, D] = size(model(K,iBoot).Mu);
XX = regmat(t,P-1);
YY = regmat(t,P-1);

for ik = 1:K %ik: ith cluster
    % mean regression trajectories (xb,yb)
    xb{K,iBoot}(:,ik) = XX*model(K,iBoot).Mu(:,ik,1);
    yb{K,iBoot}(:,ik) = YY*model(K,iBoot).Mu(:,ik,2);
    
    % move the mean regression trajectories to (0,0) (rxb,ryb)
    %rxb{K,iBoot}(:,ik) = xb{K,iBoot}(:,ik) - xb{K,iBoot}(1,ik);
    %ryb{K,iBoot}(:,ik) = yb{K,iBoot}(:,ik) - yb{K,iBoot}(1,ik);
end

for i = 1:n_uni_track
    length_itrack = length( Trajs.Y{i,:}(:,1) );
    clval = model(K,iBoot).C(i);
    
    difflon = Trajs.Y{i,:}(:,1) - xb{K,iBoot}(1:length_itrack,clval);
    difflat = Trajs.Y{i,:}(:,2) - yb{K,iBoot}(1:length_itrack,clval);
    difftot = sum(difflon.^2+difflat.^2);
    nspread(K,iBoot) = nspread(K,iBoot)+difftot;    
end
nspread(K,iBoot) = nspread(K,iBoot)/n_uni_track;

end 
end

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
p1 = plot(linspace(1,nK,nK),max(nLhood,[],2),'k-o')
axis([0 11 -3.8 -2.75])
xlabel('Number of clusters')
ylabel('Log likelihood')
set(gca,'FontSize',10)

set(gcf,'position',[3,3,300,320])
saveas(gcf,'log.pdf')


%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot(linspace(1,nK,nK),min(nspread,[],2),'k-s')
axis([0 11 500 6500])
xlabel('Number of clusters')
ylabel('Within cluster error')
set(gca,'FontSize',10)
set(gcf,'position',[3,3,300,320])
saveas(gcf,'log.pdf')

%%
csvwrite('ClusterSpread.ERAI.txt', [max(nLhood,[],2),min(nspread,[],2)]);

