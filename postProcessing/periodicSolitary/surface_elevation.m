datas="surfaceElevation.dat";

raw_datas=table2array(readtable(datas));

indices=find(raw_datas(:,1)>0);
y=raw_datas(indices,:);
timesteps=y(:,1);
free_surface_elevation=y(:,2:end);
plot(timesteps,free_surface_elevation(:,1))

xlabel('timestep');
ylabel('waveheight');


