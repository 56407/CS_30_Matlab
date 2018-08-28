clc;clear; close all
 
eta0=1;
g=10;
h=10;
f0=1e-4;
l=sqrt(g*h)/f0;
dt=10;
ef=0.05;
pint=50;

c=sqrt(g*h);

x=-400000:100:400000;
dx=x(2)-x(1);

eta=zeros(3,length(x));

eta(1,:)=-eta0*sign(x);
eta(2,:)=eta(1,:);
 
xip=eta0/h*sign(x); % potential vorticity

fid=fopen('output_geoadj.dat','w');
fprintf(fid,'%+f ',eta(1,:));
fprintf(fid,'\n');
% fprintf(fid,'%f ',eta(2,:));

wb=waitbar(0,'Simulation is running...');
for t=1:10000
    for i=2:size(eta,2)-1
        eta(3,i)=2*eta(2,i)-eta(1,i)+(dt)^2*(g*h*(eta(2,i+1)-2*eta(2,i)+eta(2,i-1))/(dx)^2-h*(f0^2)*(xip(i)+eta(2,i)/h));
    end
    
    eta(3,1)=1;
    eta(3,end)=-1;
    
    eta(3,1)=eta(2,1)+c*dt/dx*(eta(2,2)-eta(2,1));
    eta(3,end)=eta(2,end)-c*dt/dx*(eta(2,end)-eta(2,end-1));
    
    etaf=eta(3,:);
    for i=2:size(eta,2)-1
        eta(3,i)=(1-ef)*etaf(i)+0.5*ef*(etaf(i-1)+etaf(i+1));
    end
    
    if rem(t,pint)==0
    fprintf(fid,'%+f ',eta(3,:));
    fprintf(fid,'\n');
    end

    
    eta(1,:)=eta(2,:);
    eta(2,:)=eta(3,:);
    eta(3,:)=0;
    
    waitbar(t/10000,wb)
end
close(wb);
fclose(fid);

fid=fopen('output_geoadj.dat');
for t=1:inf
    dat=fgetl(fid);
    if length(dat)~=1
        et=str2num(dat);
        plot(x(2:end-1),et(2:end-1));
        title(['t=' num2str((t-1)*dt*pint,'%4.0f')]);
        ylim([-1.5 1.5])
        set(gcf,'color','w')
        frame = getframe(1);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if t == 1;
            imwrite(imind,cm,'geoadj.gif','gif','DelayTime',0.02, 'Loopcount',inf);
        else
            imwrite(imind,cm,'geoadj.gif','gif','DelayTime',0.02,'WriteMode','append');
        end
        
    else
        break
    end
end
fclose(fid);