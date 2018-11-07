clc
clear all
close all
% Cargando matrices
load('A.mat')
load('b.mat')
load('XDOM.mat')
load('YDOM.mat')
load('ZDOM.mat')
load('u.mat')
load('tiempo.mat')
load('Su.mat')
load('XDO.mat')
load('YDO.mat')
load('ZSOL.mat')
arica=[];
iquique=[];
mejillones=[];
patache=[];
dart=[];
boyas= [0] %[0 1 2 3 4]% input('Ingrese con que boyas quiere invertir como vector (Arica 0, Iquique 1, Patache 2, Mejillones 3, DART 4): ');
% for tinv=[30] %tiempos en que se genera una nueva inversi�n
tinv = 70
%tinv=input('Ingrese los minutos que desea invertir: '); %por si se
%necesita una inversi�n exacta
tini=[100 80 100 100 50];% desde que tiempo comienza a tomar datos, excluyendo as� los del terremoto
sigm=[];
AA=[];
BB=[];
load('A.mat')
load('b.mat')
tlong=[];
%genero ahora matrices nuevas de A y b para excluir los datos en que aun no
%hay variaci�n debido al tsunami
for i=boyas
    AA=[AA;A(i*961+1+tini(i+1):i*961+tinv*4,:)];
    BB=[BB;b(i*961+1+tini(i+1):i*961+tinv*4)];
    tlong(i+1)=length(b(i*961+1+tini(i+1):i*961+tinv*4));    
end
%%
tau=2;
time=[];
tlim=240; % Tiempo en min en el que se actualiza el valor de varianza de la diagonal (240 no actualiza).
for i=boyas
    COP=[];
    if tinv>tlim %actualizando covarianza con los resultados de la inversi�n anterior
        resid=abs(AA*mm-BB);
        COP=cov(resid(sum(tlong(1:i))+1:sum(tlong(1:i+1)))); 
    else %utilizando s�lo los datos previos al tsunami
        COP=cov(b(i*961+1:i*961+tini(i+1)));
    end
    %Correlaci�n
    for j=1:1:tinv*4-tini(i+1)
        for ii=1:1:tinv*4-tini(i+1)
          time(ii+sum(tlong(1:i)),j+sum(tlong(1:i)))=COP*exp(-abs(tiempo(ii)-tiempo(j))/tau);
        end 
    end 
end
%Dimensiones matriciales
M=length(A(1,:));
N=length(BB);
b=BB;
A=AA;
save('time.mat','time')
L=chol(time,'lower');%cholesky
Linv=L^-1;

%%
% Golden Section Search (Iterando)
a=0;% start of interval
en=500;
epsilon=0.001;               % accuracy value
iter=1000;                       % maximum number of iterations
tau=double((sqrt(5)-1)/2);      % golden proportion coefficient, around 0.618
k=0;                            % number of iterations
x1=a+(1-tau)*(en-a);             % computing x values
x2=a+tau*(en-a);
aux=1;
for uu=[x1 x2]
    
bb=[Linv*b; linspace(0,0,length(Su(:,1)))'];
G=[Linv*A; sqrt(uu)*Su];
i=i+1;
Q=[];
R=[];
for j=1:M
    v=G(:,j);
    for i=1:j-1
        R(i,j)=Q(:,i)'*G(:,j);
        v=v-R(i,j)*Q(:,i);
    end
    R(j,j)=norm(v);
    Q(:,j)=v/R(j,j);
end

mm=mldivide(R,Q'*bb); %Calcula m
Umin(aux)=norm(bb-G*mm)^2;
sigabic2(aux)=Umin(aux)/N;
tt(aux)=0;
for j=1:1:M
    tt(aux)=tt(aux)+2*log(abs(R(j,j)));
end
aux=aux+1;

end

%%
f_x1=N*log(2*pi*Umin(1)/N)-M*log(x1)+tt(1)+2*sum(log(diag(Linv)))+N+4;                     % computing values in x points

f_x2=N*log(2*pi*Umin(2)/N)-M*log(x2)+tt(2)+2*sum(log(diag(Linv)))+N+4;
plot(x1,f_x1,'rx')              % plotting x
plot(x2,f_x2,'rx')

while ((abs(en-a)>epsilon) && (k<iter))
    k=k+1;
    if(f_x1<f_x2)
        en=x2;
        x2=x1;
        x1=a+(1-tau)*(en-a);
        aux=1
         for uu=[x1 x2] 
            bb=[Linv*b; linspace(0,0,length(Su(:,1)))'];
            G=[Linv*A; sqrt(uu)*Su];
            i=i+1;
            Q=[];
            R=[];
            for j=1:M
            v=G(:,j);
                    for i=1:j-1
                    R(i,j)=Q(:,i)'*G(:,j);
                     v=v-R(i,j)*Q(:,i);
                    end
            R(j,j)=norm(v);
            Q(:,j)=v/R(j,j);
            end
mm=mldivide(R,Q'*bb);
Umin(aux)=(norm(bb-G*mm))^2;
sigabic2(aux)=Umin(aux)/(length(b)-2);
tt(aux)=0;
for j=1:1:length(mm)
    tt(aux)=tt(aux)+2*log(abs(R(j,j)));
end
aux=aux+1;


end
        
     f_x1=N*log(2*pi*Umin(1)/N)-M*log(x1)+tt(1)+2*sum(log(diag(Linv)))+N+4;                     % computing values in x points
    f_x2=N*log(2*pi*Umin(2)/N)-M*log(x2)+tt(2)+2*sum(log(diag(Linv)))+N+4;
        plot(x1,f_x1,'rx');
    else
        a=x1;
        x1=x2;
        aux=1
        x2=a+tau*(en-a);
             for uu=[x1 x2] 
bb=[Linv*b; linspace(0,0,length(Su(:,1)))'];
G=[Linv*A; sqrt(uu)*Su];
i=i+1;
Q=[];
R=[];
for j=1:M
    v=G(:,j);
    for i=1:j-1
        R(i,j)=Q(:,i)'*G(:,j);
        v=v-R(i,j)*Q(:,i);
    end
    R(j,j)=norm(v);
    Q(:,j)=v/R(j,j);
end
mm=mldivide(R,Q'*bb);
Umin(aux)=norm(bb-G*mm)^2;
sigabic2(aux)=Umin(aux)/N;
tt(aux)=0;
for j=1:1:M
    tt(aux)=tt(aux)+2*log(abs(R(j,j)));
end
aux=aux+1;


end
f_x1=N*log(2*pi*Umin(1)/N)-M*log(x1)+tt(1)+2*sum(log(diag(Linv)))+N+4;                     % computing values in x points

f_x2=N*log(2*pi*Umin(2)/N)-M*log(x2)+tt(2)+2*sum(log(diag(Linv)))+N+4;       
        plot(x2,f_x2,'rx')
    end
    aux=1;
    k=k+1;
end


% chooses minimum point
if(f_x1<f_x2)
    ABIC=f_x1;
    uu=x1;%Hyperparameter mu
    Umin=Umin(1);
    sprintf('x_min=%f', x1)
    sprintf('f(x_min)=%f ', f_x1)
    Cm=(diag(sigabic2(1)*R^(-1)*(R')^-1)).^0.5 
    plot(x1,f_x1,'ro')
    sigabic2=sigabic2(1) %Hyperparameter sigma
else
    ABIC=f_x2;
    uu=x2
    Umin=Umin(2)
    sprintf('x_min=%f', x2)
    sprintf('f(x_min)=%f ', f_x2)
    plot(x2,f_x2,'ro')
    sigabic2=sigabic2(2)
    Cm=(diag(sigabic2*R^(-1)*(R')^-1)).^0.5
end
Varpost=sigabic2*(G'*G)^-1;
plot(x1,f_x1,'rx')              % plotting x
%%
plot(x2,f_x2,'rx')
%Hasta aqu� termina la inversi�n
%% Reconstruyendo la superficie inicial
close all
load('XDOM.mat')
load('XDO.mat')
load('YDO.mat')
load('YDOM.mat')
load('ZSOL.mat')
load('u.mat')
load('ZDOM.mat')

 s=0;
 can=[];
 eta=zeros*XDOM;

 for a=1:1:length(A(1,:));
  ETA=zeros*XDOM;
  lam=0.1; 
  i= u(a,1); % punto medio de la onda gaussiana
  j= u(a,2);
  
  ETA=exp(-((XDOM-i)/lam).^2-((YDOM-j)/lam).^2);
  idx=find(ETA<0.001);
  ETA(idx)=0;
  eta=ETA*mm(a)+eta;
 
 end
    %% Reconstruyendo la superficie inicial de una posible soluci�n cualquiera
close all
load('XDOM.mat')
load('XDO.mat')
load('YDO.mat')
load('YDOM.mat')
load('ZSOL.mat')
load('u.mat')
load('ZDOM.mat')
h=mvnrnd(mm,Varpost)'; %elecci�n de m por muestreo random en posterior
 s=0;
 can=[];
 etaCh=zeros*XDOM;

 for a=1:1:length(A(1,:));
  ETA=zeros*XDOM;
  lam=0.1; 
  i= u(a,1); % punto medio de la onda gaussiana
  j= u(a,2);
  
  ETA=exp(-((XDOM-i)/lam).^2-((YDOM-j)/lam).^2);
  idx=find(ETA<0.001);
  ETA(idx)=0;
  etaCh=ETA*h(a)+etaCh;
 
    end
%% Gr�fico Comparaci�n alturas en Iquique con diferentes soluciones
load('A.mat')
load('b.mat')
%try
close all
clf
clc
figure (2)
aux3=1;
title('Alturas en Iquique')
subplot 311
hold on
title('Inversión ABIC, Datos reales')
xlabel('tiempo [hrs]')
ylabel('altura [m]')
plot(tiempo(tini(aux3+1)+1:240*4+1)/60,A(961+1+tini(aux3+1):2*961,:)*mm,'m')
plot(tiempo(tini(aux3+1)+1:240*4+1)./60,b(961+1+tini(aux3+1):2*961),'k','LineWidth',1)
%%
subplot 312
%load('Mareogramas.mat')
title('Resultados de Hayes, Datos reales')
hold on
xlabel('tiempo [hrs]')
ylabel('altura [m]')
plot(tiempo(tini(aux3+1)+1:240*4+1)/60,b(962+tini(aux3+1):2*961),'k','LineWidth',1)
plot(TG(2).time(tini(aux3+1)/4*60+1:240*60)./3600,TG(2).eta(tini(aux3+1)/4*60+1:240*60),'Color',[0 0.5 0.5])
%%
subplot 313
%load('AnMareo.mat')
hold on
title('Resultados de An, Datos reales')
xlabel('tiempo [hrs]')
ylabel('altura [m]')
plot(tiempo(tini(aux3+1)+1:240*4+1)/60,b(962+tini(aux3+1):2*961),'k','LineWidth',1)
plot(TG(2).time(tini(aux3+1)/4*60+1:240*60)./3600,TG(2).eta(tini(aux3+1)/4*60+1:240*60),'Color',[0 0.5 0.5])
tinvh=num2str(tinv/60);
%end
% tinvm=num2str(tinv);
% cd ..
% cd ..
% 
% cd ImagTesis
% cd IdentMod
% cd Caso3
% if tinv<60
%     print(['Comparacion_Iquique_en_' tinvm 'min'],'-dpng')
% else
%     print(['Comparacion_Iquique_en_' tinvh 'hr'],'-dpng')
% end
% cd ..
% cd ..
% cd ..
% cd ABIC
% cd ABIC
% end
%% Gr�fico soluciones de deformaci�n inicial
hold off
load('XDP.mat')
load('YDP.mat')
load('CHAO.mat')
figure (1)
clear h aux3
aux4=0;
for i=0:1:4
    if find(boyas==i)>0
        aux4=aux4+1;
    else
       aux3=i;
    end
end
if aux4<5
if aux3==0;
    a='Arica';
elseif aux3==1
    a='Iquique';
elseif aux3==2
    a='Patache';
elseif aux3==3
    a='Mejillones'; 
elseif aux3==4
    a='DART';
end
end
subplot 411
hold on
title('a)              Promedio Modelo, usando ABIC')
s=surf(XDOM,YDOM,eta),shading interp;
contour(XDOM,YDOM,ZDOM, [0 0],'w');
plot(-70.817,-19.642,'*');
ylim([-21 -18])
xlim([-72 -70])
colorbar
hold off
subplot 412
contour(XDOM,YDOM,ZDOM, [0 0],'w');
hold on
title('b)              Soluci�n de Hayes')
h=surf(XDO,YDO,ZSOL),shading interp;
plot(-70.817,-19.642,'*');
ylim([-21 -18])
xlim([-72 -70])
colorbar
hold off
subplot 413
contour(XDOM,YDOM,ZDOM, [0 0],'w');
hold on
title('c)              Soluci�n de An')
h=surf(XD,YD,CHAO),shading interp;
plot(-70.817,-19.642,'*');
ylim([-21 -18])
xlim([-72 -70])
colorbar
hold off
subplot 414
contour(XDOM,YDOM,ZDOM, [0 0],'w');
hold on
title('d)              Desviaci�n est�ndar Inversi�n')
h=surf(XDOM,YDOM,etaCh),shading interp;
plot(-70.817,-19.642,'*');
ylim([-21 -18])
xlim([-72 -70])
colorbar
hold off


tinvh=num2str(tinv/60);
tinvm=num2str(tinv);
% cd ..
% cd ..
% cd ImagTesis
% cd IdentMod
% cd Caso3
% if tinv<60
%     print(['Soluci�n_sin_' a '_en_' tinvm 'min'],'-dpng')
% else
%     print(['Soluci�n_sin_' a '_en_' tinvh 'hr'],'-dpng')
% end
% cd ..
% cd ..
% cd ..
% cd ABIC
% cd ABIC
%%
% clear Y
% X=min(mm)-min(Cm):0.01:max(mm)+max(Cm);
% figure (3)
% YY=[];
% for i=1:1:17
%     kk=find(u(:,4)<i+1&u(:,4)>i-1);
%     ax=0;
%    for j=u(kk(1),3)+4:1:u(kk(length(kk)),3)+4
%        ax=ax+1;
%     Y(i,:)=normpdf(X,mm(kk(ax)),Cm(kk(ax)));
% subplot(17,12,(17-i)*12+j)
% plot(X,Y(i,:))
%    end
%    YY=[YY;Y];
% end



%% Gr�fico del valor medio de los coeficientes y de sus respectivas desviaciones est�ndar
% figure (4)
% subplot 211
% hold on
% title('Promedio Coeficientes')
% imagesc(XDOM(:,1),YDOM(1,:),flipud(imrotate(ZDOM,90)))
% set(gca,'Ydir','Normal');
% contour(XDOM,YDOM,ZDOM, [0 0],'b');
% h(1)=scatter(u(:,1),u(:,2),[],mm,'filled');
% caxis([-1 1])
% colorbar;
% ylim([-21 -18])
% xlim([-73 -70])
% hold off
% % subplot 312
% % hold on
% % title('Desvest/coef.')
% % imagesc(XDOM(:,1),YDOM(1,:),flipud(imrotate(ZDOM,90)))
% % set(gca,'Ydir','Normal');
% % contour(XDOM,YDOM,ZDOM, [0 0],'b');
% % hh=scatter(u(:,1),u(:,2),[],abs(Cm)./abs(mm),'filled');
% % caxis([min(abs(Cm)./abs(mm)) max(abs(Cm)./abs(mm))])
% % colorbar;
% % ylim([-21 -18])
% % xlim([-73 -70])
% % hold off
% subplot 212
% hold on
% title('Desviaci�n Est�ndar')
% imagesc(XDOM(:,1),YDOM(1,:),flipud(imrotate(ZDOM,90)))
% set(gca,'Ydir','Normal');
% contour(XDOM,YDOM,ZDOM, [0 0],'b');
% hh=scatter(u(:,1),u(:,2),[],Cm,'filled');
% caxis([min(Cm) max(Cm)])
% colorbar;
% ylim([-21 -18])
% xlim([-73 -70])
% hold off
% cd ..
% cd ..
% cd ImagTesis
% if tinv<60
%     print(['Scatter_' a '_en_' tinvm 'min'],'-dpng')
% else
%     print(['Scatter_' a '_en_' tinvh 'hr'],'-dpng')
% end
% cd ..
% cd ABIC
% cd ABIC
% 
% 

%% Varianza prior y posterior
% figure (5)
% Varprior=Su'*Su;
% Cprior=diag(Varprior).^0.5;
% subplot 211
% title('Desv. Est Posterior')
% scatter(u(:,1),u(:,2),[],Cm,'filled');
% subplot 212
% title('Desv. Est Prior')
% scatter(u(:,1),u(:,2),[],Cprior,'filled');
% figure (6)
% clear Y
% X=min(mm)-min(min(Cm,Cprior)):0.01:max(mm)+max(max(Cm,Cprior));
% YY=[];
% YYprior=[];
% Yprior=[];
% for i=1:1:17 %Para todas las fuentes unitarias
%      kk=find(u(:,4)<i+1&u(:,4)>i-1);
%     ax=0;
%    for j=u(kk(1),3)+4:1:u(kk(length(kk)),3)+4
%        ax=ax+1;
%         Y(i,:)=normpdf(X,mm(kk(ax)),Cm(kk(ax)));
%         Yprior(i,:)=normpdf(X,x(kk(ax)),Cprior(kk(ax)));
%         subplot(17,12,(17-i)*12+j)
%         plot(X,Y(i,:),X,Yprior(i,:))
%    end
%    YYprior=[YYprior;Yprior];
%    YY=[YY;Y];
%  end
%% Calcular para mare�grafo faltante
% clear h aux3
% aux4=0;
% for i=0:1:4
%     if find(boyas==i)>0
%         aux4=aux4+1;
%     else
%        aux3=i;
%     end
% end
% if aux4<5
% if aux3==0;
%     a='Arica';
% elseif aux3==1
%     a='Iquique';
% elseif aux3==2
%     a='Patache';
% elseif aux3==3
%     a='Mejillones'; 
% elseif aux3==4
%     a='DART';
% end
% end
% if aux4==4;
% A1=A;
% load('A.mat')
% load('b.mat')
% tinvh=num2str(tinv/60);
% tinvm=num2str(tinv);
% tpron=240;
% %tpron=input('tiempo de pronostico [min]: ')
% figure (7)
% if tinv<60
%     title(['Predicci�n de ' a ' con inversi�n de ' tinvm ' minutos de datos'])
% elseif tinv==60
%     title(['Predicci�n de ' a ' con inversi�n de ' tinvh ' hr de datos'])
% elseif tinv>60
%     title(['Predicci�n de ' a ' con inversi�n de ' tinvh ' hrs de datos'])
% end
% 
% for j=1:1:1000
% h=mvnrnd(mm,sigabic2*(G'*G)^-1);
% A2=A(961*aux3+tini(aux3+1)+1:aux3*961+tpron*4,:);
% PronBoy=A2*h';
% hold on
% p1=plot(tiempo(1+tini(aux3+1):tpron*4)/4,PronBoy,'Color',[0.8 0.55 1])
% end
% 
% p2=plot(tiempo(1+tini(aux3+1):tpron*4)/4,b(aux3*961+1+tini(aux3+1):aux3*961+tpron*4),'b','LineWidth',2)
% p3=plot(tiempo(1+tini(aux3+1):tpron*4)/4,A2*mm,'m','LineWidth',1.5)
% legend([p1 p2 p3],'Incertidumbre','Datos en tiempo real','Modelo mas probable')
% xlabel('Tiempo [min]')
% ylabel('Altura [m]')
% 
% 
% cd ..
% cd ..
% cd ImagTesis
% cd IdentMod
% if tinv<60
% print(['Prediccion_' a '_inversi�n_de_' tinvm 'min'],'-dpng')
% elseif tinv==60
%     print(['Prediccion_' a '_inversi�n_de_' tinvh 'hr'],'-dpng')
% elseif tinv>60
%     print(['Prediccion_' a '_inversi�n_de_' tinvh 'hrs'],'-dpng')
% end
% cd ..
% cd ..
% cd ABIC
% cd ABIC
% end
% 
% 
% %% Resultados mareogramas y deformaci�n inicial 
% if aux4==5
% au=0;
% figure(9)
% load('A.mat')
% load('b.mat')
% 
% for aux3=[0 2 3 4]
%     if aux3==3
%         au=au+2
%     else
%     au=1+au;
%     end
% subplot(3,3,au+1)
% if aux3==0;
%     a='Mare�grafo de Arica';
% elseif aux3==1
%     a='Mare�grafo de Iquique';
% elseif aux3==2
%     a='Mare�grafo de Patache';
% elseif aux3==3
%     a='Mare�grafo de Mejillones'; 
% elseif aux3==4
%     a='DART';
% end
% title([a])
% tpron=240;
% maPro=0;
% miPro=0;
% for j=1:1:1000
% h=mvnrnd(mm,sigabic2*(G'*G)^-1);
% A2=A(961*aux3+tini(aux3+1)+1:aux3*961+tpron*4,:);
% PronBoy=A2*h';
% if miPro>min(PronBoy)
%     miPro=min(PronBoy);
% end
% if maPro<max(PronBoy)
%     maPro=max(PronBoy);
% end
% 
% hold on
% p1=plot(tiempo(1+tini(aux3+1):tpron*4),PronBoy,'Color',[0.8 0.55 1]);
% end
% p2=plot(tiempo(1+tini(aux3+1):tpron*4),b(aux3*961+1+tini(aux3+1):aux3*961+tpron*4),'b','LineWidth',1);
% p3=plot(tiempo(1+tini(aux3+1):tpron*4),A2*mm,'m','LineWidth',1);
% xlim([0 tpron])
% ylim([miPro maPro])
% end  
% 
% 
% subplot(3,3,[1,4,7])
% hold on
% s=surf(XDOM,YDOM,eta),shading interp;
% contour(XDOM,YDOM,ZDOM, [0 0],'w');
% 
% ylim([-21 -18])
% xlim([-72 -70])
% colorbar
% hold off
% 
% subplot(3,3,[8,9])
% hold on
% title('Mare�grafo de Iquique')
% tpron=240;
% aux3=1;
% for j=1:1:1000
% h=mvnrnd(mm,sigabic2*(G'*G)^-1);
% A2=A(961*aux3+tini(aux3+1)+1:aux3*961+tpron*4,:);
% PronBoy=A2*h';
% hold on
% p1=plot(tiempo(1+tini(aux3+1):tpron*4),PronBoy,'Color',[0.8 0.55 1])
% if miPro>min(PronBoy)
%     miPro=min(PronBoy);
% end
% if maPro<max(PronBoy)
%     maPro=max(PronBoy);
% end
% end
% p2=plot(tiempo(1+tini(aux3+1):tpron*4),b(aux3*961+1+tini(aux3+1):aux3*961+tpron*4),'b','LineWidth',1)
% p3=plot(tiempo(1+tini(aux3+1):tpron*4),A2*mm,'m','LineWidth',1)
% xlim([0 tpron])
% ylim([miPro maPro])
% hold off
% cd .. 
% cd ..
% cd ImagTesis
% cd IdentMod
% if tinv<60
% print(['Prediccion_inversi�n_de_' tinvm 'min'],'-dpng')
% elseif tinv==60
%     print(['Prediccion_inversi�n_de_' tinvh 'hr'],'-dpng')
% elseif tinv>60
%     print(['Prediccion_inversi�n_de_' tinvh 'hrs'],'-dpng')
% end
% cd ..
% cd ..
% cd ABIC
% cd ABIC
% end
% cami=x1
% save('mm.mat','mm')
% %%
% figure (10)
% clear bb G Umin sig2 uu q 
% i=0;
% rr=0;
% ee=1;
% q=[];
% for sig2=0:0.5:10
%  for uu=0:0.5:10
%  %uu=cami;
% bb=[Linv*BB; linspace(0,0,length(Su(:,1)))'];
% G=[Linv*AA; sqrt(uu)*Su];
% Q=[];
% R=[];
% for j=1:M
%     v=G(:,j);
%     for i=1:j-1
%         R(i,j)=Q(:,i)'*G(:,j);
%         v=v-R(i,j)*Q(:,i);
%     end
%     R(j,j)=norm(v);
%     Q(:,j)=v/R(j,j);
% end
% 
% % calculando coeficientes
% mm=mldivide(R,Q'*bb);
% Umin=norm(bb-G*mm)^2;
% tt=0;
% for j=1:1:length(mm)
%     tt=tt+2*log(abs(R(j,j)));
% end
% rr=rr+1;
%  q(rr,ee)=N*log(2*sig2*pi())-M*log(uu)+tt+2*sum(log(diag(Linv)))+Umin/sig2+4;  
%  end
%  ee=ee+1;
%  rr=0;
% end
%  [uu,sig2]=meshgrid(0:0.5:10);
%  surf(uu,sig2,q), shading interp;
% %end
% sig21=0:0.5:10;
% ux=find(abs(uu(1,:)-cami)==min(abs(uu(1,:)-cami)));
% siy=find(abs(sig2(:,1)-sigabic2)== min(abs(sig2(:,1)-sigabic2)));
% plot(sig2(:,1),q(ux,:))
% figure(17)
% plot(uu(1,:),q(:,siy))

%% skill
% sk=[];
% load('A.mat')
% load('b.mat')
% for aux3=[0 1 2 3 4]
%     for j=1:1:1000
%      hh=mvnrnd(mm,sigabic2*(G'*G)^-1)';
%    for i=1:1:179
%     if tini(aux3+1)<61+i
%         
%         if tini(aux3+1)<i
%             O=b(aux3*961+tini(aux3+1):aux3*961+i+60);
%             P=A(aux3*961+tini(aux3+1):aux3*961+i+60,:)*hh;
%         else
%             O=b(aux3*961+i:aux3*961+i+60);
%             P=A(aux3*961+i:aux3*961+i+60,:)*hh;
%         end
%     sk(aux3*1000+j,i)=1-(sum((P-O).^2)/sum(O.^2))^0.5;
%     end
%     end
% end
% end
% %  skup=floor(sk)+ceil((sk-floor(sk))/0.1)*0.1;
% % arica=[skup(1:1000,tini(1)-60:length(sk(1,:)))'];
% % iquique=[skup(1001:2000,:)'];
% % patache=[skup(2001:3000,tini(3)-60:length(sk(3,:)))'];
% % mejillones=[skup(3001:4000,tini(4)-60:length(sk(4,:)))'];
% % dart=[skup(4001:5000,tini(5)-60:length(sk(5,:)))'];
% 
% arica2=[sk(1:1000,tini(1)-60:length(sk(1,:)))'];
% iquique2=[sk(1001:2000,:)'];
% patache2=[sk(2001:3000,tini(3)-60:length(sk(3,:)))'];
% mejillones2=[sk(3001:4000,tini(4)-60:length(sk(4,:)))'];
% dart2=[sk(4001:5000,tini(5)-60:length(sk(5,:)))'];
% cd ..
% cd ..
% cd ImagTesis
% cd ImagCorrec
% tinvm=num2str(tinv);
% % figure(21)
% % histogram(arica,[-4:0.1:1])
% % ylim([0 5*10^4]);
% % print(['Histograma_Arica_' tinvm 'min_3'],'-dpng')
% figure (26)
% boxplot(arica2','Colors','k','Notch','off','BoxStyle','filled','Symbol','w')
% print(['BP_Arica_' tinvm 'min_3'],'-dpng')
% ylim([-12 1])
% % figure(22)
% % histogram(iquique,[-0.5:0.1:1])
% % ylim([0 10*10^4]);
% % print(['Histograma_Iquique_' tinvm 'min_3'],'-dpng')
% figure (27)
% boxplot(iquique2','Colors','k','Notch','off','BoxStyle','filled','Symbol','w')
% print(['BP_Iquique_' tinvm 'min_3'],'-dpng')
% % figure(23)
% % % histogram(patache,[-2:0.1:1])
% % ylim([0 12*10^4]);
% % print(['Histograma_Patache_' tinvm 'min_3'],'-dpng')
% figure (28)
% boxplot(patache2','Colors','k','Notch','off','BoxStyle','filled','Symbol','w')
% print(['BP_Patache_' tinvm 'min_3'],'-dpng')
% % figure(24)
% % histogram(mejillones,[-4:0.1:1])
% % ylim([0 5*10^4]);
% % print(['Histograma_Mejillones_' tinvm 'min_3'],'-dpng')
% figure (29)
% boxplot(mejillones2','Colors','k','Notch','off','BoxStyle','filled','Symbol','w')
% print(['BP_Mejillones_' tinvm 'min_3'],'-dpng')
% % figure(25)
% % histogram(dart,[-0.5:0.1:1])
% % ylim([0 8*10^4]);
% % print(['Histograma_Dart_' tinvm 'min_3'],'-dpng')
% figure (30)
% boxplot(dart2','Colors','k','Notch','off','BoxStyle','filled','Symbol','w')
% print(['BP_Dart_' tinvm 'min_3'],'-dpng')
% cd .. 
% cd ..
% cd ABIC
% cd ABIC
% % skill
% % sk=[];
% % load('A.mat')
% % load('b.mat')
% % if tinv>30
% %     load('Mapacalor.mat')
% % end
% % for aux3=[0 1 2 3 4]
% % for i=1:1:179
% %     if tini(aux3+1)<61+i
% %         if tini(aux3+1)<i
% %             O=b(aux3*961+tini(aux3+1):aux3*961+i+60);
% %             P=A(aux3*961+tini(aux3+1):aux3*961+i+60,:)*mm;
% %         else
% %             O=b(aux3*961+i:aux3*961+i+60);
% %             P=A(aux3*961+i:aux3*961+i+60,:)*mm;
% %         end
% %     sk(aux3+1,i)=1-(sum((P-O).^2)/sum(O.^2))^0.5;
% %     end
% % end
% % end
% % arica=[arica sk(1,tini(1)-60:length(sk(1,:)))'];
% % iquique=[iquique sk(2,:)'];
% % patache=[patache sk(3,tini(3)-60:length(sk(3,:)))'];
% % mejillones=[mejillones sk(4,tini(4)-60:length(sk(4,:)))'];
% % dart=[dart sk(5,tini(5)-60:length(sk(5,:)))'];
% % save('Mapacalor.mat','arica','iquique','patache','mejillones','dart')
% % if tinv==240
% %     tinv2=[30 45 60 90 120 150 180 210 240]
% %     HeatMap(arica,'ColumnLabels',tinv2,'colormap',hot,'DisplayRange',2)
% % HeatMap(iquique,'ColumnLabels',tinv2,'colormap',hot,'DisplayRange',0.8)
% % HeatMap(patache,'ColumnLabels',tinv2,'colormap',hot,'DisplayRange',1)
% % HeatMap(mejillones,'ColumnLabels',tinv2,'colormap',hot,'DisplayRange',1)
% % HeatMap(dart,'ColumnLabels',tinv2,'colormap',hot,'DisplayRange',0.5)
% % 
% % end

% end