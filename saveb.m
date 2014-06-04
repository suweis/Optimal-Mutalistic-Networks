%% Function to save outputs of the main simulation.

function saveb(fname,bold,b,A,dim,s,beta)
name=['SET YOUR DIRECTORY/resultHTI','_S=',num2str(dim),'_omega=',num2str(beta),'_sigma=',num2str(s),fname];
save(name, 'bold','bN','b','xeq0','xeqN','xeq','maxEVold','maxEVN','maxEV0')
%save([name,'_bold','.dat'],'bold','-ascii') %If you want also .dat ouptut
%save([name,'_b','.dat'],'b','-ascii')
%save([name,'_A','.dat'],'A','-ascii')
end
