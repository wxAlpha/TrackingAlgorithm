%
% Data in the Trajs variable
%
%%
clear

setcctpath() 
K = 3;
nBoot = 500;
nLhood = zeros(nBoot,1);

for iBoot = 1:nBoot
disp([K,iBoot])
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% (1) load the .xls data into cell format
f = readtable('track_PostFiltering.ERAI.IC.xlsx');
track = f{:,:};

track_lon = track(:,2);
track_lat = track(:,3);

% id 
id = track(:,1);
% Number of trajectories [n_uni_track]
uni_track = unique(id);
n_uni_track = size(uni_track);
n_uni_track = n_uni_track(1,1);

uni_track_random = randsample(uni_track,n_uni_track);
uni_track_random = uni_track_random(:,1);
uni_track_random_record(iBoot,:) = uni_track_random;

Trajs.Y = {};
for i = 1:n_uni_track
  itrack =  find( id==uni_track_random(i) );
  Trajs.Y(i,:) = { [track_lon(itrack), track_lat(itrack)] } ;
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% lrm_model
%
% (2) Set the model settings
%
ops.zero             = 'nozero'; % see CELL2ARRAY for possible values
ops.method           = 'lrm';    % see above for possible values
ops.K                = K;        % number of clusters
ops.order            = 2;        % 1-linear, 2-quadratic, 3-cubic, etc.
    
%
% (3) Set the general EM settings
%

ops.NumEMStarts       = 1;  % number of times to run EM (best model chosen)
ops.IterLimit         = 30; % maximum allowable iterations for any one EM start
ops.ValidateLhood     = 0;  % run a large-sample Lhood calculation at end
                            % of each EM start; used for time-based 
                            % (approximate) alignment. This will give a more
                            % accurate approximation of the log-likelihood.
                            % You only need this if you are running comparison
                            % experiments between different models based on
                            % log likelihood scores.
ops.MsgHnd            = -1; 

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Cluster the data according to the above settings
% set the warning state
% warning off MATLAB:singularMatrix;
warning off all;
warning on MATLAB:divideByZero;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% perform the clustering
model(iBoot) = curve_clust(Trajs,ops);
nLhood(iBoot) = model(iBoot).TrainLhood_ppt;

end 
%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% convert back to .xls format
maxrun = find(nLhood==max(nLhood))
max_uni_track = uni_track_random_record(maxrun,:);

cluster = id;
cluster = nan;
for i = 1:n_uni_track
  itrack =  find(id==max_uni_track(i));
  cluster(itrack,1) = model(maxrun).C(i);
end
track =[track cluster];

%xlswrite('wCluster_track_PostFiltering_IC_C3.xlsx',track)

%% 
t = 0:1:78;
[P,K,D] = size(model(maxrun).Mu);
XX = regmat(t,P-1);
YY = regmat(t,P-1);

for ik=1:K
    % mean regression trajectories (xb,yb)
    xb(:,ik) = XX*model(maxrun).Mu(:,ik,1);
    yb(:,ik) = YY*model(maxrun).Mu(:,ik,2);
    
    % move the mean regression trajectories to (0,0) (rxb,ryb)
    rxb(:,ik) = xb(:,ik) - xb(1,ik);
    ryb(:,ik) = yb(:,ik) - yb(1,ik);
end

%% Plot the mean regression trajectories 
plot(rxb(1:16,1),ryb(1:16,1),'-ko','LineWidth',1.2)
hold on
plot(rxb(1:16,2),ryb(1:16,2),'-ks','LineWidth',1.2)
plot(rxb(1:16,3),ryb(1:16,3),'-k^','LineWidth',1.2)
hold off
xlabel('Rel. Longitude')
ylabel('Rel. Latitude')
lgd = legend('A','B','C','Location','northwest','FontSize',20)
set(lgd,'position',[0.15 0.7 0.2 0.2])
set(gca,'FontSize',14)
set(gcf,'position',[5,5,500,400])

%% 
addpath /home/yujia/software/matlab/m_map
%m_proj('lambert','lon',[100 128],'lat',[18 42])
m_proj('miller','lon',[103 128],'lat',[17 43])
m_coast('color','k');
hold on
m_line(xb(1:16,1),yb(1:16,1),'linewi',1,'color','k','LineStyle','-','Marker','o');
m_line(xb(1:16,2),yb(1:16,2),'linewi',1,'color','k','LineStyle','-','Marker','s');
m_line(xb(1:16,3),yb(1:16,3),'linewi',1,'color','k','LineStyle','-','Marker','^');
m_grid('linestyle','none','tickdir','out');
hold off
xlabel('Longitude')
ylabel('Latitude')
set(gca,'FontSize',10)
set(gcf,'position',[5,5,290,400])

%%
csvwrite('Reg_trajectory.C1.txt', [yb(1:16,1),xb(1:16,1)]);
csvwrite('Reg_trajectory.C2.txt', [yb(1:16,2),xb(1:16,2)]);
csvwrite('Reg_trajectory.C3.txt', [yb(1:16,3),xb(1:16,3)]);






