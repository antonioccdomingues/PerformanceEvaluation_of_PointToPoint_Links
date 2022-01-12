function [PL , PLVoIP , APD , APDVoIP , MPD , MPDVoIP ,  TT] = Simulator4f(lambda,C,f,P,n)
% INPUT PARAMETERS:
%  lambda - packet rate (packets/sec)
%  C      - link bandwidth (Mbps)
%  f      - queue size (Bytes)
%  P      - number of packets (stopping criterium)
%  n      - number of VoIP packets
% OUTPUT PARAMETERS:
%  PLdata  - Packet Loss of data packets (%)
%  PLVoIP  - Packet Loss of VoIP packets (%)
%  APDdata - Average Delay of data packets (milliseconds)
%  APDVoIP - Average Delay of VoIP packets (milliseconds)
%  MPDdata - Maximum Delay of data packets (milliseconds)
%  MPDVoIP - Maximum Delay of VoIP packets (milliseconds)
%  TT - Transmitted Throughput (data + VoIP) (Mbps)


%Events:
ARRIVAL= 0;       % Arrival of a packet            
DEPARTURE= 1;     % Departure of a packet

DADOS=0;          % Packet ser de dados
VOIP=1;           % Packet ser de VoIP

%State variables:
STATE = 0;          % 0 - connection free; 1 - connection bysy
QUEUEOCCUPATION= 0; % Occupation of the queue (in Bytes)
QUEUE= [];          % Size and arriving time instant of each packet in the queue

%Statistical Counters:
TOTALPACKETS= 0;            % No. of packets arrived to the system
TOTALVOIPPACKETS= 0;        % No. of VoIP packets arrived to the system
LOSTPACKETS= 0;         % No. of packets dropped due to buffer overflow
LOSTVOIPPACKETS= 0;             % No. of VoIP packets dropped due to buffer overflow
TRANSMITTEDPACKETS= 0;      % No. of transmitted packets
TRANSMITTEDVOIPPACKETS= 0;  % No. of VoIP transmitted packets
TRANSMITTEDBYTES= 0;        % Sum of the Bytes of transmitted packets
DELAYS= 0;                  % Sum of the delays of transmitted packets
DELAYSVOIP= 0;              % Sum of the delays of transmitted VoIP packets
MAXDELAY= 0;                % Maximum delay among all transmitted packets
MAXDELAYVOIP= 0;            % Maximum delay among all transmitted VoIP packets

% Initializing the simulation clock:
Clock= 0;

% Initializing the List of Events with the first ARRIVAL:
tmp= Clock + exprnd(1/lambda);
tmpVoIP = Clock + (randi([16 24])) * 0.001;  % clock para os packets VoIP

EventList = [ARRIVAL, tmp, GeneratePacketSize(), tmp , DADOS];
for i=1:n
    EventList = [EventList; ARRIVAL , tmpVoIP, GenerateVoIPPacketSize(), tmpVoIP , VOIP ];
end

%Similation loop:
while TRANSMITTEDPACKETS + TRANSMITTEDVOIPPACKETS<P               % Stopping criterium
    EventList= sortrows(EventList,2);    % Order EventList by time
    Event= EventList(1,1);               % Get first event and 
    Clock= EventList(1,2);               %   and
    PacketSize= EventList(1,3);          %   associated
    ArrivalInstant= EventList(1,4);      %   parameters
    PacketType= EventList(1,5);          % Tipo de packet, VoIP ou Dados
    EventList(1,:)= [];                  % Eliminate first event
    switch Event
        case ARRIVAL                     % If first event is an ARRIVAL
            if PacketType==0             % Caso o Packet seja de Dados
                TOTALPACKETS= TOTALPACKETS+1;
                tmp= Clock + exprnd(1/lambda);
                EventList = [EventList; ARRIVAL, tmp, GeneratePacketSize(), tmp, DADOS];
                if STATE==0
                    STATE= 1;
                    EventList = [EventList; DEPARTURE, Clock + 8*PacketSize/(C*10^6), PacketSize, Clock, DADOS];
                else
                    if QUEUEOCCUPATION + PacketSize <= 0.9*f %se a queue mais o tamanho do pacote for menor que 90% da capacidade o pacote é aceite
                        QUEUE= [QUEUE;PacketSize , Clock, DADOS];
                        QUEUEOCCUPATION= QUEUEOCCUPATION + PacketSize;
                    else
                        LOSTPACKETS= LOSTPACKETS + 1; %caso contrario o pacote e descartado
                    end
                end
            end
            if PacketType==1            % Caso o Packet seja de VoIP
                TOTALVOIPPACKETS= TOTALVOIPPACKETS+1;
                tmp= Clock + (randi([16 24])) * 0.001;
                EventList = [EventList; ARRIVAL, tmp, GenerateVoIPPacketSize(), tmp , VOIP];
                if STATE==0
                    STATE= 1;
                    EventList = [EventList; DEPARTURE, Clock + 8*PacketSize/(C*10^6), PacketSize, Clock, VOIP];
                else
                    if QUEUEOCCUPATION + PacketSize <= f
                        QUEUE= [QUEUE;PacketSize , Clock, VOIP];
                        QUEUEOCCUPATION= QUEUEOCCUPATION + PacketSize;
                    else
                        LOSTVOIPPACKETS= LOSTVOIPPACKETS + 1;
                    end
                end
            end
        case DEPARTURE                     % If first event is a DEPARTURE
            if PacketType==0               % Caso o Packet seja de Dados
                TRANSMITTEDBYTES= TRANSMITTEDBYTES + PacketSize;
                DELAYS= DELAYS + (Clock - ArrivalInstant);
                if Clock - ArrivalInstant > MAXDELAY
                    MAXDELAY= Clock - ArrivalInstant;
                end
                TRANSMITTEDPACKETS= TRANSMITTEDPACKETS + 1;
                if QUEUEOCCUPATION > 0
                    QUEUE = sortrows(QUEUE, -3);
                    EventList = [EventList; DEPARTURE, Clock + 8*QUEUE(1,1)/(C*10^6), QUEUE(1,1), QUEUE(1,2), QUEUE(1,3)];
                    QUEUEOCCUPATION= QUEUEOCCUPATION - QUEUE(1,1);
                    QUEUE(1,:)= [];
                else
                    STATE= 0;
                end
            end

            if PacketType==1            % Caso o Packet seja de VoIP
                TRANSMITTEDBYTES= TRANSMITTEDBYTES + PacketSize;
                DELAYSVOIP= DELAYSVOIP + (Clock - ArrivalInstant);
                if Clock - ArrivalInstant > MAXDELAYVOIP
                    MAXDELAYVOIP= Clock - ArrivalInstant;
                end
                TRANSMITTEDVOIPPACKETS= TRANSMITTEDVOIPPACKETS + 1;
                if QUEUEOCCUPATION > 0
                    QUEUE = sortrows(QUEUE, -3);
                    EventList = [EventList; DEPARTURE, Clock + 8*QUEUE(1,1)/(C*10^6), QUEUE(1,1), QUEUE(1,2), QUEUE(1,3)];
                    QUEUEOCCUPATION= QUEUEOCCUPATION - QUEUE(1,1);
                    QUEUE(1,:)= [];
                else
                    STATE= 0;
                end
            end
    end
end

%Performance parameters determination:
PL= 100*LOSTPACKETS/TOTALPACKETS;               % in %
PLVoIP= 100*LOSTVOIPPACKETS/TOTALVOIPPACKETS;   % in %
APD= 1000*DELAYS/TRANSMITTEDPACKETS;            % in milliseconds
APDVoIP=1000*DELAYSVOIP/TRANSMITTEDVOIPPACKETS; % in milliseconds
MPD= 1000*MAXDELAY;                             % in milliseconds
MPDVoIP=1000*MAXDELAYVOIP;                      % in milliseconds
TT= 10^(-6)*TRANSMITTEDBYTES*8/Clock;           % in Mbps

end

function out= GeneratePacketSize()
    aux= rand();
    aux2= [65:109 111:1517];
    if aux <= 0.19
        out= 64;
    elseif aux <= 0.19 + 0.23
        out= 110;
    elseif aux <= 0.19 + 0.23 + 0.17
        out= 1518;
    else
        out = aux2(randi(length(aux2)));
    end
end

function out= GenerateVoIPPacketSize()
    out = randi([110, 130]);
end