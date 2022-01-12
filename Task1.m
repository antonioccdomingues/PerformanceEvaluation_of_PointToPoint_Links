%% a)
lambda = [400 800 1200 1600 2000];
C=10;
f=1000000;
P=10000;

N=50; %Run Simulator1 50 times
per1 = zeros(5,N); % 5 arrays de 50 posições
per2 = zeros(5,N);%average packet delay storage
per3 = zeros(5,N);
per4 = zeros(5,N);


for i = 1 : 5
    for it = 1 : N
        [per1(i,it), per2(i,it), per3(i,it), per4(i,it)] = Simulator1(lambda(i),C,f,P);
       
        
    end
end
alfa = 0.1;


%avgPacketLoss = mean(per1);
%term = norminv(1-alfa/2)*sqrt(var(per1)/N);
%fprintf('PacketLoss             = %.2e +- %.2e\n',avgPacketLoss,term)
term1=zeros(1,5);
avgpacketDelay= zeros(1,5);
errHigh= zeros(1,5);
errLow = zeros(1,5);
for i = 1 : 5
    avgpacketDelay(i) = mean(per2(i, 1:N));
    term1(i) = norminv(1-alfa/2)*sqrt(var(per2(i, 1:N))/N);
    fprintf('Av. Packet Delay (ms) with λ= %d  = %.2e +- %.2e\n',lambda(i), avgpacketDelay(i),term1(i))
    errHigh(i) = term1(i);
    errLow(i) = -(term1(i));
end

bar(lambda,avgpacketDelay)
title('Average packet delays')
ylabel('Average packet delay (ms)')
xlabel('Packet Rate (packets/sec)')

hold on

er = errorbar(lambda,avgpacketDelay,errLow,errHigh);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

hold off

% avgpacketDelay = mean(per2);
% term1 = norminv(1-alfa/2)*sqrt(var(per2)/N);
% fprintf('Av. Packet Delay (ms)  = %.2e +- %.2e\n',avgpacketDelay,term1)
% 
% 
% maxPacketDelay = mean(per3);
% term2 = norminv(1-alfa/2)*sqrt(var(per3)/N);
% fprintf('Max. Packet Delay (ms) = %.2e +- %.2e\n',maxPacketDelay,term2)
% 
% throughput = mean(per4);
% term3 = norminv(1-alfa/2)*sqrt(var(per4)/N);
% fprintf('Throughput (Mbps)      = %.2e +- %.2e\n',throughput,term3)

%% b)
C = 10; %MBbps
f = [100000 20000 10000 2000]; %tamanho da fila
P = 10000; %pacotes (stop criterium)
lambda = 1800; %pps
N = 50; %numero de vezes que vamos chamar o simulator1
size = length(f); %numero de valores de f
alfa = 0.1;

%matrizes para guardar os valores do simulator1
%como chamamos N vezes, entao serao guardados N valores
%apesar de haver interesse apenas no PL e no PD e preciso
%criar matrizes para os 4 parametros pois o simulator
%retorna sempre os 4 parametros de performance
PLmatriz = zeros(size,N); 
APDmatriz = zeros(size,N);
MPDmatriz = zeros(size,N);
TTmatriz = zeros(size,N);

PLmean = zeros(1,size);
APDmean = zeros(1,size);
PLterm = zeros(1,size);
APDterm = zeros(1,size);

for j = 1:size %para cada valor de f
   for i = 1:N %chamamos o simulator1 50vezes para cada valor de f
       [PLmatriz(j,i), APDmatriz(j,i), MPDmatriz(j,i), TTmatriz(j,i)] = Simulator1(lambda,C,f(j),P); %f(0)=100000, f(1)=20000, ...
   end
   
   PLmean(j) = mean(PLmatriz(j,:)); %media dos valores de cada valor de f (f(0) fica no PLmean(0), f(1) fica no PLmean(1)...
   APDmean(j) = mean(APDmatriz(j,:));
   PLterm(j) = norminv(1-alfa/2)*sqrt(var(PLmatriz(j,:))/N);
   APDterm(j) = norminv(1-alfa/2)*sqrt(var(APDmatriz(j,:))/N);
   
end

%grafico do packet loss
figure(1)
bar(f, PLmean)
ylabel('Average packet loss (%)')
xlabel('queue size(bytes)')
hold on
er = errorbar(f, PLmean, PLterm);
er.Color = [0 1 0];
er.LineStyle = 'none';
title('Packet Loss (%)');
hold off

%grafico do packet delay
figure(2)
bar(f, APDmean)
ylabel('Average packet delay (ms)')
xlabel('queue size(bytes)')
hold on
er = errorbar(f, APDmean, APDterm);
er.Color = [0 1 0];
er.LineStyle = 'none';
title('Avg Packet Delay (ms)');
hold off


%% c)

lambda = 1800;
C=[10 20 30 40];
f=1000000;
P=10000;

N=50; %Run Simulator1 50 times
per1 = zeros(4,N); % 5 arrays de 50 posições
per2 = zeros(4,N);%average packet delay storage
per3 = zeros(4,N);
per4 = zeros(4,N);


for i = 1 : 4
    for it = 1 : N
        [per1(i,it), per2(i,it), per3(i,it), per4(i,it)] = Simulator1(lambda,C(i),f,P);
       
        
    end
end
alfa = 0.1;


%avgPacketLoss = mean(per1);
%term = norminv(1-alfa/2)*sqrt(var(per1)/N);
%fprintf('PacketLoss             = %.2e +- %.2e\n',avgPacketLoss,term)
term1=zeros(1,4);
avgpacketDelay= zeros(1,4);
errHigh= zeros(1,4);
errLow = zeros(1,4);
for i = 1 : 4
    avgpacketDelay(i) = mean(per2(i, 1:N));
    term1(i) = norminv(1-alfa/2)*sqrt(var(per2(i, 1:N))/N);
    fprintf('Av. Packet Delay (ms) with C= %d  = %.2e +- %.2e\n',C(i), avgpacketDelay(i),term1(i))
    errHigh(i) = term1(i);
    errLow(i) = -(term1(i));
end

bar(C,avgpacketDelay)
title('Average packet delays')
ylabel('Average packet delay (ms)')
xlabel('link bandwidth (Mbps)')

hold on

er = errorbar(C,avgpacketDelay,errLow,errHigh);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

hold off

%% d)
lambda = 1800; %pps
packetSize = [65:1:109 111:1:1517]; %tamanho que os pacotes podem ter, para alem de 64, 110 e 1518
C = [40000000 30000000 20000000 10000000]; %valores de C
size = length(C); %numero de valores de C
mg1 = zeros(1, size); %matriz para guardar os valores teoricos de acordo com o valor de C

for j = 1:size %percorre cada posicao de C
    
    Sresto2 = 0; %quadrado do remaining (para usar no calculo do ES2)

    S64 = (64*8)/C(j);
    S110 = (110*8)/C(j);
    S1518 = (1518*8)/C(j);
    Sresto = (mean(packetSize)*8)/C(j);
    for i = [65:109 111:1517]
        Sresto2 = Sresto2 + ((i*8)/C(j))^2; %valor total do S do restantes packetSize
    end

    Sresto2 = Sresto2/length(packetSize); %valor medio do S dos restantes packetSize
    
    ES = 0.19 * S64 + 0.23 * S110 + 0.17 * S1518 + (1 - 0.19 - 0.23 - 0.17) * Sresto;
    ES2 = 0.19 * S64^2 + 0.23 * S110^2 + 0.17 * S1518^2 + (1 - 0.19 - 0.23 - 0.17) * Sresto2;
    mg1(j) = (((lambda*ES2) / (2*(1-(lambda*ES)))) + ES)*1000; %APD em ms
   
end
   
disp(mg1);
%grafico do packet loss
figure(1)
bar(C, mg1)
ylabel('Average packet delay (ms)')
xlabel('link bandwidth (Mbps)')
hold on
title('Avg Packet Delay (ms) - theoretical values');
hold off

%% e)

lambda = 1800;
C=[10 20 30 40];
f=1000000;
P=10000;

N=50; %Run Simulator1 50 times
per1 = zeros(4,N); % 4 arrays de N posições
per2 = zeros(4,N);%average packet delay storage
per3 = zeros(4,N);
per4 = zeros(4,N);
per5 = zeros(4,N);  %average packet delay for packetsize = 64
per6 = zeros(4,N);  %average packet delay for packetsize = 110
per7 = zeros(4,N);  %average packet delay for packetsize = 1518



for i = 1 : 4
    for it = 1 : N
        [per1(i,it), per2(i,it), per3(i,it), per4(i,it), per5(i,it), per6(i,it), per7(i,it)] = Simulator1e(lambda,C(i),f,P);
       
        
    end
end
alfa = 0.1;

term1=zeros(1,4);
term2=zeros(1,4);
term3=zeros(1,4);
term4=zeros(1,4);
avgpacketDelay= zeros(1,4);
avgpacketDelay64= zeros(1,4);
avgpacketDelay110= zeros(1,4);
avgpacketDelay1518= zeros(1,4);

error= [];
for i = 1 : 4
%     avgpacketDelay(i) = mean(per2(i, 1:N));
%     term1(i) = norminv(1-alfa/2)*sqrt(var(per2(i, 1:N))/N);
%     fprintf('Av. Packet Delay (ms)                      with C= %d      = %.2e +- %.2e\n',C(i), avgpacketDelay(i),term1(i))

    avgpacketDelay64(i) = mean(per5(i, 1:N));
    term2(i) = norminv(1-alfa/2)*sqrt(var(per5(i, 1:N))/N);
    fprintf('Av. Packet Delay (ms) of packetSize = 64   with C= %d      = %.2e +- %.2e\n',C(i), avgpacketDelay64(i),term2(i))

    avgpacketDelay110(i) = mean(per6(i, 1:N));
    term3(i) = norminv(1-alfa/2)*sqrt(var(per6(i, 1:N))/N);
    fprintf('Av. Packet Delay (ms) of packetSize = 110  with C= %d      = %.2e +- %.2e\n',C(i), avgpacketDelay110(i),term3(i))

    avgpacketDelay1518(i) = mean(per7(i, 1:N));
    term4(i) = norminv(1-alfa/2)*sqrt(var(per7(i, 1:N))/N);
    fprintf('Av. Packet Delay (ms) of packetSize = 1518 with C= %d      = %.2e +- %.2e\n',C(i), avgpacketDelay1518(i),term4(i))
    fprintf('---------------------------------------------------------------------------------\n')
    error =[error; term2(i) term3(i) term4(i)];

end
y = [];

%determining y values (average packet delay per packet size)
for i =1:4
    y = [y; avgpacketDelay64(i) avgpacketDelay110(i) avgpacketDelay1518(i)];
end
%creating bar graph
b = bar(C,y, 'grouped');
hold on
%determining error
[ngroups,nbars] = size(y);
x = nan(nbars, ngroups);

for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end

errorbar(x',y,error,'k','linestyle','none');
title('Average packet delays')
ylabel('Average packet delay (ms)')
xlabel('link bandwidth (Mbps)')
legend('64 Bytes', '110 Bytes', '1518 Bytes', Location='northeast')
hold off

