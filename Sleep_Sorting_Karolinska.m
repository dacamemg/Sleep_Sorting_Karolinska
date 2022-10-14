%***AUTOMATIC POLYSONOGRAPHY SCRIPT HOUR PER HOUR ANALYSIS***
% This script performs the states sorting of 24h recording, along with
% macro (awake/sleep features) and micro (EEG oscillations) analysis.

%These variables must be addressed
% - Original sample OS - 2000 or 6250
% - Downsampling factor DSF - 2 or 5 (according with the OS: 2000<-2 and
    % 6250<-5
% - Channels 1:2 or 1:3
% - Presence of Cortex channel
% - Variables channels HIP 1 - EMG 2 or 3


%INPUTS
% Sequentials one hour length .m files


%OUTPUTS
%States_H: Percentage of each state per hour (AWA, SWS, REM)

%States_Trans: Number of transitions (bouts) per state (AWA, SWS, REM) with minimal window of 10s.

%States_Total_Length: Total bouts length per state (AWA, SWS, REM) and
%mininal 10s window.

%States_Bouts_Length: Average bouts length calculated by the total legth
%divided by number of bouts (AWA, SWS, REM - minimal 10s window).

%States_Cont_Nor_H: Percentage from total state of continuous state (AWA, SWS, REM - minimal window of 20s)

%States_Period_Mean_H: Mean duration of continuous state (AWA, SWS, REM - minimal window of 20s)

%States_Period_Max_H: Max duration of continous state (AWA, SWS, REM - minimal window of 20s)

%States_Cont_Number: Number of continous episodes (AWA, SWS, REM - minimal window of 20s)

%States_Cont_Prev_Post_H: State previous and posterior continous states
%(AWA, SWS, REM - minimal window of 20s)

%States_FFT_PLFP_norm: 24h FFT of normalized signal separated by states
%(AWA, SWS and REM: 24 hours/lines per state)

%States_FFT_PLFP: 24h FFT signal separated by states (AWA, SWS and REM: 24 hours/lines per state)

%States_FFT_PLFP_All_norm: FFT of normalized signal separated mean per hour
%(disregarding the states).

%States_FFT_PLFP_All: FFT separated mean per hour(disregarding the states).

%States_PB_H_norm: normalized LFP power organized as 24 hours (lines) and each 3 colums (AWA,
    ...SWS, REM) is one band, Delta, Theta, Beta, Gamma

%States_PB_H: LFP power organized as 24 hours (lines) and each 4 colums (AWA,
    ...SWS, REM) is one band, Delta, Theta, Beta, Gamma

%States_MIraw_AV_H: Raw MI value per hour (AWA, SWS, REM)

%States_Spindle: Spindle analysis during SWS period divided by night and day period (Total number, Duration[s], Episodes/min and frequency) 



clear 
close all



% #####Basal parameters#####
%Block time in seconds
Time=10;
%Downsampling factor
DSF=5;
%Original Sample
OS=6250;
%Frequency sample
FS=OS/DSF;
%Window analysis
FST=FS*Time;
%Hour by hour analysis (60min * 60sec)
Hour=(60*60)/Time;


%Designing a bandstop filter for 50Hz noise IIR
BSfilter_IIR=designfilt('bandstopiir','FilterOrder',20,'HalfPowerFrequency1',49,'HalfPowerFrequency2',51,'SampleRate',FS);

LPfilter_IIR = designfilt('lowpassiir','FilterOrder',8,'PassbandFrequency',FS/2,'PassbandRipple',0.2,'SampleRate',OS);



%     ***Choose the animal file***
% Loading the the choosen animal files

% Animal 1
% channels=(1:3);

% Animal 2
channels=(4:6);

%Loading and concatenating the .m files (each file is one hour duration) 
matfiles = dir('*.mat') ;
N = length(matfiles) ;
filenames = cell(N,1) ;
for i = 1:N
   filenames{i} = matfiles(i).name ;
   temp = load(matfiles(i).name) ;
   
   if i==2
       data=temp.data(:,channels);
   end
   
   if i>2
       data=[data; temp.data(:,channels)];
   end
       
end
data=data';



%   ***VARIABLES*** 
%   ***New filter before downsampling - 22/02/15***
%Naming Variable Hippocampus 
HIP=data(1,:); %not filtered signal
HIP=filtfilt(LPfilter_IIR,HIP); %Digital Low pass filter
HIP=downsample(HIP,DSF); %downsampling
%Naming Variable Cortex 
CTX=data(2,:); %not filtered signal
CTX=filtfilt(LPfilter_IIR,CTX);%Digital Low pass filter
CTX=downsample(CTX,DSF); %downsampling
%Naming Variable EMG
EMG=data(3,:); %not filtered signal
%EMG=filtfilt(LPfilter_IIR,EMG);%Digital Low pass filter
EMG=downsample(EMG,DSF); %downsampling

%Hippocampus analysis
% LFP=HIP;

%Cortex analysis
LFP=CTX;


clear data temp N i


%EMG 
%Taking only the envelop of the signal
[yupper,ylower] = envelope(EMG,FST,'rms');
clear yupper
%Getting only the derivate of the signal
EMG_D = diff(ylower);
clear ylower

%Pwelch - Getting the EMG power
m=1;
for i=1:FST:(length(EMG_D)-FST)
    pxx= pwelch(EMG_D(i:(i+FST)),FS,0,FS);
    EMG_P(:,m)=pxx';
    m=m+1;
end
clear m

%EMG
%Selecting only the power from the spectrogram
PEMG=EMG_P(1:500,:);
%Normalizing the EMG power
PEMG_mean=mean(PEMG');
PEMG_mean=PEMG_mean';
for i=1:(length(PEMG(:,1))) %i = Frequencia
    for j=1:(length(PEMG(1,:))) %j = Tempo
        PEMG_norm(i,j)=PEMG(i,j)./PEMG_mean(i,1);
    end
    i=0;
end
% Selecting the EMG power frequencies
RMS=mean(PEMG_norm(12:400,:));
%Calculating the EMG zscore
ZRMS=zscore(RMS);

%Adjusting the threshold as 25% from min zcored PEMG_norm power (EMG power).
threshold=0.75*min(ZRMS);



%LFP 
% LFP Filtering noise
LFP_F=filtfilt(BSfilter_IIR,LFP);


%Pwelch - Getting the LFP power
m=1;
for i=1:FST:(length(LFP_F)-FST)
    pxx= pwelch(LFP_F(i:(i+FST)),FS,0,FS);
    PLFP(:,m)=pxx';
    m=m+1;
end
% Normalizing the Power by the total power of the bin
for i=1:(length(PLFP(1,:))) %i = Time
    for j=1:(length(PLFP(:,1))) %j = Frequency
    PLFP_norm(j,i)=PLFP(j,i)./sum(PLFP(:,i));
    end
    j=0;
end

%Selecting the LFP power band frequencies
Ptheta=mean(PLFP_norm(7:9,:));
Pdelta=mean(PLFP_norm(2:4,:));
%Calculating the LFP theta/delta power
HRMS=Ptheta./Pdelta;
ZHRMS=zscore(HRMS);%Calculating the LFP zscore



%FIGURE
%Checking the Theta/Delta Zscore
figure;
plot(ZHRMS);
line([0 (length(ZHRMS))],[0 0],'color','r','LineWidth',1);
xsteps=(720:720:8640);
xticks(xsteps);
xticklabels({'2','4','6','8','10','12','14','16','18','20','22','24'}) 
ylabel('Theta/Delta Power Zscore')
set(gca,'FontSize',16)
xlim([0 (length(ZHRMS))])
%FIGURE
%Checking the threshold
figure;
plot(ZRMS);
line([0 (length(ZRMS))],[threshold threshold],'color','r','LineWidth',1);
xsteps=(720:720:8640);
xticks(xsteps);
xticklabels({'2','4','6','8','10','12','14','16','18','20','22','24'})  
ylabel('EMG Power Zscore')
set(gca,'FontSize',16)
xlim([0 (length(ZRMS))])
clear xsteps
%FIGURE
%Checking the four states separation
figure; plot(ZRMS,ZHRMS,'b.');
line([threshold threshold],[min(ZHRMS) max(ZHRMS)],'color','r','LineWidth',1);
line([min(ZRMS) max(ZRMS)],[0 0],'color','r','LineWidth',1);
xlabel('EMG Power Zscore') 
ylabel('Theta/Delta Power Zscore')
set(gca,'FontSize',16)









%   ***CATEGORIZATION***

% Categoring awake and sleep according with zscore
clear i j
for i=1:length(ZRMS)
if ZRMS(1,i)>threshold
M_I(1,i)=3; % Moviment
else 
    M_I(1,i)=0; % Immobility
end
end
clear i j

%Sleep and delta
for i=1:length(M_I)
    if M_I(1,i)==3
        M_I_delta(1,i)=0;
    else
        M_I_delta(1,i)=mean(PLFP_norm(2:4,i));
    end
end

j=1;
for i=1:length(M_I_delta)
    if M_I_delta(1,i)>0
        M_I_delta_Z_temp(1,j)=M_I_delta(1,i);
        j=j+1;
    end
end


%Sleep and theta
for i=1:length(M_I)
    if M_I(1,i)==3
        M_I_theta(1,i)=0;
    else
        M_I_theta(1,i)=mean(PLFP_norm(7:9,i));
    end
end

j=1;
for i=1:length(M_I_theta)
    if M_I_theta(1,i)>0
        M_I_theta_Z_temp(1,j)=M_I_theta(1,i);
        j=j+1;
    end
end


ZHRMS_temp = M_I_theta_Z_temp./M_I_delta_Z_temp; %theta/delta only in immobil momments 

ZHRMS_temp = zscore(ZHRMS_temp); 



%Categoring the REM (1) and SWS (2)
for i=1:length(ZHRMS_temp)
if ZHRMS_temp(1,i)>0
ZHRMS_temp_final(1,i)=1; 
else
ZHRMS_temp_final(1,i)=2; 
end
end
%Unifying with M_I
j=1;
for i=1:length(M_I)
if M_I(1,i)==3
result(1,i)=3;
else
result(1,i)= ZHRMS_temp_final(1,j);
j=j+1;
end
end
clear ZHRMS_temp ZHRMS_temp_final 


%Categorizando
REM=find(result==1);
SWS=find(result==2);
AWA=find(result==3);



%###OUTPUT### - Sleep time

%Calculating the amount time in each state in seconds in % of total
%recording
rem=size(REM);
rem=rem(1,2);
rem_time=((rem*100)/(length(result)));

sws=size(SWS);
sws=sws(1,2);
sws_time=((sws*100)/(length(result)));

awa=size(AWA);
awa=awa(1,2);
awa_time=((awa*100)/(length(result)));

Sleep(1,1)=awa_time;
Sleep(2,1)=sws_time;
Sleep(3,1)=rem_time;


  




%   ***Separation of states episodes per hour***

%REM episodes separed hour by hour
clear j m i REM_Temp REM_H
j=Hour;
m=1;
for i=1:Hour:(length(result-Hour))%Steps hour-hour
    REM_Temp(m,:)=(REM>=i & REM<=j); %Finding episodes within one hour
    REM_H(m,:)=REM_Temp(m,:).*REM; %Getting the time stamps from REM categorization
    m=m+1;%Changing lines (hours)
    j=j+360; %Steps
end
clear j m i REM_Temp
% Counting the number of episodes per hour
clear j m i REM_Temp REM_H_Count
for m=1:1:(length(REM_H(:,1))) %Steps hour-hour (24 total)
REM_Temp=find(REM_H(m,:)>0); %Finding valid values
REM_H_Count_Temp=size(REM_Temp); %Couting valid values
REM_H_Count(m,1)=REM_H_Count_Temp(1,2); %Getting the valid values 
clear REM_H_Count_Temp
end
REM_H_Count_Nor=(REM_H_Count./Hour)*100; %Normalizing the valid values
clear REM_H_Count_Temp REM_Temp


%SWS episodes separed hour by hour
clear j m i SWS_Temp SWS_H
j=Hour;
m=1;
for i=1:Hour:(length(result-Hour))
    SWS_Temp(m,:)=(SWS>=i & SWS<=j);
    SWS_H(m,:)=SWS_Temp(m,:).*SWS;
    m=m+1;
    j=j+360;
end
clear j m i SWS_Temp
% Counting the number of episodes per hour
clear j m i SWS_Temp SWS_H_Count
for m=1:1:(length(SWS_H(:,1)))
SWS_Temp=find(SWS_H(m,:)>0);
SWS_H_Count_Temp=size(SWS_Temp);
SWS_H_Count(m,1)=SWS_H_Count_Temp(1,2);
clear SWS_H_Count_Temp
end
SWS_H_Count_Nor=(SWS_H_Count./Hour)*100;
clear SWS_H_Count_Temp SWS_Temp


%AWA episodes separed hour by hour
clear j m i AWA_Temp AWA_H
j=Hour;
m=1;
for i=1:Hour:(length(result-Hour))
    AWA_Temp(m,:)=(AWA>=i & AWA<=j);
    AWA_H(m,:)=AWA_Temp(m,:).*AWA;
    m=m+1;
    j=j+360;
end
clear j m i AWA_Temp
% Counting the number of episodes per hour
clear j m i AWA_Temp AWA_H_Count
for m=1:1:(length(AWA_H(:,1)))
AWA_Temp=find(AWA_H(m,:)>0);
AWA_H_Count_Temp=size(AWA_Temp);
AWA_H_Count(m,1)=AWA_H_Count_Temp(1,2);
clear AWA_H_Count_Temp
end
AWA_H_Count_Nor=(AWA_H_Count./Hour)*100;
clear AWA_H_Count_Temp AWA_Temp



%###OUTPUT###

States_H(:,1)=AWA_H_Count_Nor;
States_H(:,2)=SWS_H_Count_Nor;
States_H(:,3)=REM_H_Count_Nor;



%FIGURE (Smoothed)
figure;
plot(smooth(States_H(:,1),5),'LineWidth',2);xlim([1 24]);ylim([0 100]);
hold
plot(smooth(States_H(:,2),5),'LineWidth',2)
plot(smooth(States_H(:,3),5),'LineWidth',2)
legend('AWA','SWS','REM')
xsteps=(2:2:24);
xticks(xsteps)
xlabel('Hour') 
ylabel('Time Percentage')
set(gca,'FontSize',16)


figure;
plot(smooth(States_H(:,1)),'LineWidth',2);xlim([1 24]);ylim([0 100]);
hold
plot(smooth((States_H(:,2))+(States_H(:,3)),5),'LineWidth',2);
legend('Awake','Sleep')
xsteps=(2:2:24);
xticks(xsteps)
xlabel('Hour') 
ylabel('Time Percentage')
set(gca,'FontSize',16)








%***Transitions (number of bouts) and total bouts length - minimal window of 10s*** 

%Total Transitions
transition = zeros(1,(length(result)+1)); %Pre-locate +1 (sum)

for i=1:(length(result)-1)
    if result(i+1)-result(i)==0 %if it is the same state add zero
        transition(i)=0; %no transition
    else
        transition(i)=1; %transition: add 1
    end
end
clear i

Transition_Total = NaN(24,1); 

j=360; %steps 24 hours (hour = 3600 seconds/10)
m=1; %steps for the hour
for i=1:360:((length(transition)-360)) %steps for the sum
    Transition_Total(m,1)=sum(transition(1,i:(i+j))); %sum transition per hour
    m=m+1;
end
clear i j m 



%Transition AWA
trans_awa=zeros(1,(length(result)+2));
trans_awa(1:length(result))=result;
trans_awa(trans_awa==2)=1;

for i=1:(length(trans_awa)-1)
    if trans_awa(i+1)-trans_awa(i)==-2 %if it is any state (1) followed by awa (3)
        transition_awa(i)=1; %transition
    else
        transition_awa(i)=0; %no transition: any other possibility
    end
end
clear i

Transition_AWA = NaN(24,1); 

j=360; %steps 24 hours (hour = 3600 seconds/10)
m=1; %steps for the hour
for i=1:360:((length(transition_awa)-360)) %steps for the sum
    Transition_AWA(m,1)=sum(transition_awa(1,i:(i+j))); %sum transition per hour
    m=m+1;
end
clear i j m



%Transitions SWS
trans_sws=zeros(1,(length(result)+2));
trans_sws(1:length(result))=result;
trans_sws(trans_sws==3)=1;

for i=1:(length(trans_sws)-1)
    if trans_sws(i+1)-trans_sws(i)==-1 %if it is any state (1) followed by awa (3)
        transition_sws(i)=1; %transition
    else
        transition_sws(i)=0; %no transition: any other possibility
    end
end
clear i

Transition_SWS = NaN(24,1); 

j=360; %steps 24 hours (hour = 3600 seconds/10)
m=1; %steps for the hour
for i=1:360:((length(transition_sws)-360)) %steps for the sum
    Transition_SWS(m,1)=sum(transition_sws(1,i:(i+j))); %sum transition per hour
    m=m+1;
end
clear i j m



%Transitions REM
trans_rem=zeros(1,(length(result)+2));
trans_rem(1:length(result))=result;
trans_rem(trans_rem==2)=3;

for i=1:(length(trans_rem)-1)
    if trans_rem(i+1)-trans_rem(i)==2 %if it is any state (1) followed by awa (3)
        transition_rem(i)=1; %transition
    else
        transition_rem(i)=0; %no transition: any other possibility
    end
end
clear i

Transition_REM = NaN(24,1); 

j=360; %steps 24 hours (hour = 3600 seconds/10)
m=1; %steps for the hour
for i=1:360:((length(transition_rem)-360)) %steps for the sum
    Transition_REM(m,1)=sum(transition_rem(1,i:(i+j))); %sum transition per hour
    m=m+1;
end
clear i j m


%OUTPUT
States_Trans(:,1)=Transition_AWA;
States_Trans(:,2)=Transition_SWS;
States_Trans(:,3)=Transition_REM; 








%Bouts length with minal window of 10s
%AWA
length_awa_temp=zeros(1,(length(result)+2));
length_awa_temp(1:length(result))=result;
length_awa_temp(length_awa_temp==2)=1;

j=1;
for i=1:(length(length_awa_temp)-1)
    if length_awa_temp(i+1)+length_awa_temp(i)>5 %if it is any state (1) followed by awa (3)
        length_awa_1(i)=0; 
        length_awa_1(i+1)=j+1;
        j=j+1;
    else
        length_awa_1(i+1)=0; %no transition: any other possibility
        j=1;
    end
end
clear i

 for i=1:(length(length_awa_temp)-2)
    length_awa_2(i)=(length_awa_temp(i+1)-length_awa_temp(i))-(length_awa_temp(i+2));
 end
 
 
 for i=1:(length(length_awa_2))
      if length_awa_2(i)==1
          length_awa_2(i)=1;
      else
          length_awa_2(i)=0;
      end
 end
length_awa_2(8641)=0;
 
length_awa_final=length_awa_1(1,1:8641)+length_awa_2;
 
length_awa_24 = NaN(24,1); 

j=360; %steps 24 hours (hour = 3600 seconds/10)
m=1; %steps for the hour
for i=1:360:((length(length_awa_final)-360)) %steps for the sum
    length_awa_24(m,1)=sum(length_awa_final(1,i:(i+j))); %sum transition per hour
    m=m+1;
end
clear i j m



%SWS
length_sws_temp=zeros(1,(length(result)+2));
length_sws_temp(1:length(result))=result;
length_sws_temp(length_sws_temp==3)=1;
length_sws_temp(length_sws_temp==2)=3;

j=1;
for i=1:(length(length_sws_temp)-1)
    if length_sws_temp(i+1)+length_sws_temp(i)>5 %if it is any state (1) followed by awa (3)
        length_sws_1(i)=0; 
        length_sws_1(i+1)=j+1;
        j=j+1;
    else
        length_sws_1(i+1)=0; %no transition: any other possibility
        j=1;
    end
end
clear i

 for i=1:(length(length_sws_temp)-2)
    length_sws_2(i)=(length_sws_temp(i+1)-length_sws_temp(i))-(length_sws_temp(i+2));
 end
 
 
 for i=1:(length(length_sws_2))
      if length_sws_2(i)==1
          length_sws_2(i)=1;
      else
          length_sws_2(i)=0;
      end
 end
length_sws_2(8641)=0;
 
length_sws_final=length_sws_1(1,1:8641)+length_sws_2;
 
length_sws_24 = NaN(24,1); 

j=360; %steps 24 hours (hour = 3600 seconds/10)
m=1; %steps for the hour
for i=1:360:((length(length_sws_final)-360)) %steps for the sum
    length_sws_24(m,1)=sum(length_sws_final(1,i:(i+j))); %sum transition per hour
    m=m+1;
end
clear i j m



%REM
length_rem_temp=zeros(1,(length(result)+2));
length_rem_temp(1:length(result))=result;
length_rem_temp(length_rem_temp==1)=4;
length_rem_temp(length_rem_temp==3)=1;
length_rem_temp(length_rem_temp==2)=1;
length_rem_temp(length_rem_temp==4)=3;

j=1;
for i=1:(length(length_rem_temp)-1)
    if length_rem_temp(i+1)+length_rem_temp(i)>5 %if it is any state (1) followed by awa (3)
        length_rem_1(i)=0; 
        length_rem_1(i+1)=j+1;
        j=j+1;
    else
        length_rem_1(i+1)=0; %no transition: any other possibility
        j=1;
    end
end
clear i

 for i=1:(length(length_rem_temp)-2)
    length_rem_2(i)=(length_rem_temp(i+1)-length_rem_temp(i))-(length_rem_temp(i+2));
 end
 
 
 for i=1:(length(length_rem_2))
      if length_rem_2(i)==1
          length_rem_2(i)=1;
      else
          length_rem_2(i)=0;
      end
 end
length_rem_2(8641)=0;
 
length_rem_final=length_rem_1(1,1:8641)+length_rem_2;
 
length_rem_24 = NaN(24,1); 

j=360; %steps 24 hours (hour = 3600 seconds/10)
m=1; %steps for the hour
for i=1:360:((length(length_rem_final)-360)) %steps for the sum
    length_rem_24(m,1)=sum(length_rem_final(1,i:(i+j))); %sum transition per hour
    m=m+1;
end
clear i j m



%OUPUT
States_Total_Length(:,1)=length_awa_24;
States_Total_Length(:,2)=length_sws_24;
States_Total_Length(:,3)=length_rem_24;

States_Bouts_Length=States_Total_Length./States_Trans; %average bouts length calculated by the total legth divided by number of bouts.








%   ***STATES CONTINUITY PER HOUR - minimal window of 20s***

%Counting the continuity of each state. 
    ...The higher the value, more stable in time is the state.

%REM
clear i j Int m n 
for m=1:(length(REM_H(:,1))) %Steps hour-hour (24 total)
    j=2;
    for i=1:(length(REM_H)-1) %Throughout 
        if (REM_H(m,j))-(REM_H(m,i))==1 %Checking subsequent values
            Cont_REM_H(m,i)=1; %If is subsequent, value 1
        else
            Cont_REM_H(m,i)=0; %if not subsequent, value 0 
        end
        j=j+1;
    end
end

clear i j Int m n
for m=1:(length(REM_H(:,1)))
    Cont_REM_Total_H(1,m)=sum(Cont_REM_H(m,:)); %Total continuous periods per hour
    Cont_REM_Nor_Total_H(1,m)=(Cont_REM_Total_H(1,m)./(rem))*100; %Normalization by the total state time
    Cont_REM_Nor_H(1,m)=(Cont_REM_Total_H(1,m)./REM_H_Count(m,1))*100; %Normalization by the total hour event state time
end

%Calculating the time period of each block
clear i j Temp m n 
for m=1:(length(REM_H(:,1))) %Steps hour-hour (24 total)
    j=1;
    for i=1:((length(Cont_REM_H))-1) %Covering up the whole length
        if (Cont_REM_H(m,i)+Cont_REM_H(m,(i+1)))>1 %Suming up the period duration if it is sequential
            Period_REM_H(m,i)=0; %The previous value is 0
            Period_REM_H(m,(i+1))=j+(Cont_REM_H(m,(i+1))); %The subsequent value is the sum
            j=j+1;
        else
            Period_REM_H(m,(i+1))=0; %The previous value is 0
            j=1; %If it is not sequential, then mantain the value
        end
    end
end

%Finding the mean REM period per hour
clear i j Temp m n Period_REM_Mean_H Period_REM_MEAN_H
for m=1:(length(REM_H(:,1))) %Steps hour-hour (24 total)
    Temp = find(Period_REM_H(m,:)>0); %Finding sequential periods (over zero)
    if isempty(Temp)==0 %If the period is not zero
        for i=1:(length(Temp))
            Period_REM_MEAN_H(m,i)=Period_REM_H(m,(Temp(1,i))); %Getting the values
        end
        Period_REM_Mean_H(m,1)=(sum(Period_REM_MEAN_H(m,:)))/(length(Temp)); %Getting the mean 
    else
        Period_REM_Mean_H(m,1)=0;
    end
    clear Temp
end

%Finding the max period of REM sleep
for m=1:(length(REM_H(:,1))) %Steps hour-hour (24 total)
    Period_REM_Max_H(m,1)=max(Period_REM_H(m,:)); %Just getting the mean value
end



%Finding the number of REM countinuous episodes
for i=1:length (Period_REM_H(:,1)) %24h steps
    for j=1:length (Period_REM_H(1,:)) %time stamps steps
        if Period_REM_H(i,j)==0
            temp(i,j)=0;
        else
            temp(i,j)=1;
        end
    end
    Cont_Number_REM(i,1)=sum(temp(i,:));
end
clear temp



%Finding the REM countinuous PREVIOUS and POSTERIOR state
for i=1:length(Period_REM_H(:,1)) % steps de 24 Horas
temp=find((Period_REM_H(i,:))>0); %steps timestamps
for j=1:length(temp)
if (REM_H(i,(temp(1,j))))-(Period_REM_H(i,(temp(1,j))))==0 %Evitar inicio menor ou igual a zero
Cont_REM_H_Previous(i,j)=NaN;
else
Cont_REM_H_Previous(i,j)=result((REM_H(i,(temp(1,j))))-(Period_REM_H(i,(temp(1,j))))); %state before continous period
end
if ((REM_H(i,(temp(1,j))))+2)>length (result) %evitar final maior que length (result)
Cont_REM_H_Posterior(i,j)=NaN;
else
Cont_REM_H_Posterior(i,j)=result(1,((REM_H(i,(temp(1,j))))+2)); %state after continous period
end
end
clear temp
end


Cont_REM_H_Previous(Cont_REM_H_Previous==0) = NaN;
Cont_REM_H_Posterior(Cont_REM_H_Posterior==0) = NaN;


%Finding PREVIOUS normalized
for i=1:length (Cont_REM_H_Previous(:,1)) %hours
for j=1:length(Cont_REM_H_Previous(1,:))%Steps
if Cont_REM_H_Previous(i,j)==3
Cont_REM_H_Previous_AWA(i,j)=1;
else
Cont_REM_H_Previous_AWA(i,j)=0;
end
if Cont_REM_H_Previous(i,j)==2
Cont_REM_H_Previous_SWS(i,j)=1;
else
Cont_REM_H_Previous_SWS(i,j)=0;
end
end
end


for i=1:length (Cont_REM_H_Previous(:,1)) %hours
for j=1:length(Cont_REM_H_Previous(1,:))%Steps
if Cont_REM_H_Previous(i,j)==1
    Cont_REM_H_Previous(i,j)=NaN;
end
end
end


temp=(Cont_REM_H_Previous); %Finding the total event at that hour
temp=(isnan(temp));
temp=length(temp(:,1))-sum(temp);
temp=sum(temp); %total valid events
Cont_REM_H_Previous_AWA_Nor=sum(Cont_REM_H_Previous_AWA);
Cont_REM_H_Previous_AWA_Nor=sum(Cont_REM_H_Previous_AWA_Nor);
Cont_REM_H_Previous_AWA_Nor=Cont_REM_H_Previous_AWA_Nor/temp;
Cont_REM_H_Previous_SWS_Nor=sum(Cont_REM_H_Previous_SWS);
Cont_REM_H_Previous_SWS_Nor=sum(Cont_REM_H_Previous_SWS_Nor);
Cont_REM_H_Previous_SWS_Nor=Cont_REM_H_Previous_SWS_Nor/temp;
clear temp

%Finding POSTERIOR normalized
for i=1:length (Cont_REM_H_Posterior(:,1)) %hours
for j=1:length(Cont_REM_H_Posterior(1,:))%Steps
if Cont_REM_H_Posterior(i,j)==3
Cont_REM_H_Posterior_AWA(i,j)=1;
else
Cont_REM_H_Posterior_AWA(i,j)=0;
end
if Cont_REM_H_Posterior(i,j)==2
Cont_REM_H_Posterior_SWS(i,j)=1;
else
Cont_REM_H_Posterior_SWS(i,j)=0;
end
end
end


for i=1:length (Cont_REM_H_Posterior(:,1)) %hours
for j=1:length(Cont_REM_H_Posterior(1,:))%Steps
if Cont_REM_H_Posterior(i,j)==1
    Cont_REM_H_Posterior(i,j)=NaN;
end
end
end


temp=(Cont_REM_H_Posterior); %Finding the total event at that hour
temp=(isnan(temp));
temp=length(temp(:,1))-sum(temp);
temp=sum(temp); %total valid events
Cont_REM_H_Posterior_AWA_Nor=sum(Cont_REM_H_Posterior_AWA);
Cont_REM_H_Posterior_AWA_Nor=sum(Cont_REM_H_Posterior_AWA_Nor);
Cont_REM_H_Posterior_AWA_Nor=Cont_REM_H_Posterior_AWA_Nor/temp;
Cont_REM_H_Posterior_SWS_Nor=sum(Cont_REM_H_Posterior_SWS);
Cont_REM_H_Posterior_SWS_Nor=sum(Cont_REM_H_Posterior_SWS_Nor);
Cont_REM_H_Posterior_SWS_Nor=Cont_REM_H_Posterior_SWS_Nor/temp;
clear temp





%SWS
clear i j Int m n 
for m=1:(length(SWS_H(:,1))) %Steps hour-hour (24 total)
    j=2;
    for i=1:(length(SWS_H)-1) %Throughout 
        if (SWS_H(m,j))-(SWS_H(m,i))==1 %Checking subsequent values
            Cont_SWS_H(m,i)=1; %If is subsequent, value 1
        else
            Cont_SWS_H(m,i)=0; %if not subsequent, value 0 
        end
        j=j+1;
    end
end

clear i j Int m n
for m=1:(length(SWS_H(:,1)))
    Cont_SWS_Total_H(1,m)=sum(Cont_SWS_H(m,:)); %Total continuous periods per hour
    Cont_SWS_Nor_Total_H(1,m)=(Cont_SWS_Total_H(1,m)./(sws))*100; %Normalization by the total state time
    Cont_SWS_Nor_H(1,m)=(Cont_SWS_Total_H(1,m)./SWS_H_Count(m,1))*100; %Normalization by the total hour event state time
end

%Calculating the time period of each block
clear i j Temp m n 
for m=1:(length(SWS_H(:,1))) %Steps hour-hour (24 total)
    j=1;
    for i=1:((length(Cont_SWS_H))-1) %Covering up the whole length
        if (Cont_SWS_H(m,i)+Cont_SWS_H(m,(i+1)))>1 %Suming up the period duration if it is sequential
            Period_SWS_H(m,i)=0; %The previous value is 0
            Period_SWS_H(m,(i+1))=j+(Cont_SWS_H(m,(i+1))); %The subsequent value is the sum
            j=j+1;
        else
            Period_SWS_H(m,(i+1))=0; %The previous value is 0
            j=1; %If it is not sequential, then mantain the value
        end
    end
end

%Finding the mean SWS period per hour
clear i j Temp m n
for m=1:(length(SWS_H(:,1))) %Steps hour-hour (24 total)
    Temp = find(Period_SWS_H(m,:)>0); %Finding sequential periods (over zero)
    if isempty(Temp)==0 %If the period is not zero
        for i=1:(length(Temp))
            Period_SWS_MEAN_H(m,i)=Period_SWS_H(m,(Temp(1,i))); %Getting the values
        end
        Period_SWS_Mean_H(m,1)=(sum(Period_SWS_MEAN_H(m,:)))/(length(Temp)); %Getting the mean 
    else
        Period_SWS_Mean_H(m,1)=0;
    end
    clear Temp
end

%Finding the max period of SWS sleep
for m=1:(length(SWS_H(:,1))) %Steps hour-hour (24 total)
    Period_SWS_Max_H(m,1)=max(Period_SWS_H(m,:)); %Just getting the mean value
end



%Finding the number of SWS countinuous episodes
for i=1:length (Period_SWS_H(:,1)) %24h steps
    for j=1:length (Period_SWS_H(1,:)) %time stamps steps
        if Period_SWS_H(i,j)==0
            temp(i,j)=0;
        else
            temp(i,j)=1;
        end
    end
    Cont_Number_SWS(i,1)=sum(temp(i,:));
end
clear temp





%Finding the SWS countinuous PREVIOUS and POSTERIOR state
for i=1:length(Period_SWS_H(:,1)) % steps de 24 Horas
temp=find((Period_SWS_H(i,:))>0); %steps timestamps
for j=1:length(temp)
if (SWS_H(i,(temp(1,j))))-(Period_SWS_H(i,(temp(1,j))))==0 %Evitar inicio menor ou igual a zero
Cont_SWS_H_Previous(i,j)=NaN;
else
Cont_SWS_H_Previous(i,j)=result((SWS_H(i,(temp(1,j))))-(Period_SWS_H(i,(temp(1,j))))); %state before continous period
end
if ((SWS_H(i,(temp(1,j))))+2)>length (result) %evitar final maior que length (result)
Cont_SWS_H_Posterior(i,j)=NaN;
else
Cont_SWS_H_Posterior(i,j)=result(1,((SWS_H(i,(temp(1,j))))+2)); %state after continous period
end
end
clear temp
end


Cont_SWS_H_Previous(Cont_SWS_H_Previous==0) = NaN;
Cont_SWS_H_Posterior(Cont_SWS_H_Posterior==0) = NaN;


%Finding PREVIOUS normalized
for i=1:length (Cont_SWS_H_Previous(:,1)) %hours
for j=1:length(Cont_SWS_H_Previous(1,:))%Steps
if Cont_SWS_H_Previous(i,j)==3
Cont_SWS_H_Previous_AWA(i,j)=1;
else
Cont_SWS_H_Previous_AWA(i,j)=0;
end
if Cont_SWS_H_Previous(i,j)==1
Cont_SWS_H_Previous_REM(i,j)=1;
else
Cont_SWS_H_Previous_REM(i,j)=0;
end
end
end


for i=1:length (Cont_SWS_H_Previous(:,1)) %hours
for j=1:length(Cont_SWS_H_Previous(1,:))%Steps
if Cont_SWS_H_Previous(i,j)==2
    Cont_SWS_H_Previous(i,j)=NaN;
end
end
end


temp=(Cont_SWS_H_Previous); %Finding the total event at that hour
temp=(isnan(temp));
temp=length(temp(:,1))-sum(temp);
temp=sum(temp); %total valid events
Cont_SWS_H_Previous_AWA_Nor=sum(Cont_SWS_H_Previous_AWA);
Cont_SWS_H_Previous_AWA_Nor=sum(Cont_SWS_H_Previous_AWA_Nor);
Cont_SWS_H_Previous_AWA_Nor=Cont_SWS_H_Previous_AWA_Nor/temp;
Cont_SWS_H_Previous_REM_Nor=sum(Cont_SWS_H_Previous_REM);
Cont_SWS_H_Previous_REM_Nor=sum(Cont_SWS_H_Previous_REM_Nor);
Cont_SWS_H_Previous_REM_Nor=Cont_SWS_H_Previous_REM_Nor/temp;
clear temp

%Finding POSTERIOR normalized
for i=1:length (Cont_SWS_H_Posterior(:,1)) %hours
for j=1:length(Cont_SWS_H_Posterior(1,:))%Steps
if Cont_SWS_H_Posterior(i,j)==3
Cont_SWS_H_Posterior_AWA(i,j)=1;
else
Cont_SWS_H_Posterior_AWA(i,j)=0;
end
if Cont_SWS_H_Posterior(i,j)==1
Cont_SWS_H_Posterior_REM(i,j)=1;
else
Cont_SWS_H_Posterior_REM(i,j)=0;
end
end
end


for i=1:length (Cont_SWS_H_Posterior(:,1)) %hours
for j=1:length(Cont_SWS_H_Posterior(1,:))%Steps
if Cont_SWS_H_Posterior(i,j)==2
    Cont_SWS_H_Posterior(i,j)=NaN;
end
end
end


temp=(Cont_SWS_H_Posterior); %Finding the total event at that hour
temp=(isnan(temp));
temp=length(temp(:,1))-sum(temp);
temp=sum(temp); %total valid events
Cont_SWS_H_Posterior_AWA_Nor=sum(Cont_SWS_H_Posterior_AWA);
Cont_SWS_H_Posterior_AWA_Nor=sum(Cont_SWS_H_Posterior_AWA_Nor);
Cont_SWS_H_Posterior_AWA_Nor=Cont_SWS_H_Posterior_AWA_Nor/temp;
Cont_SWS_H_Posterior_REM_Nor=sum(Cont_SWS_H_Posterior_REM);
Cont_SWS_H_Posterior_REM_Nor=sum(Cont_SWS_H_Posterior_REM_Nor);
Cont_SWS_H_Posterior_REM_Nor=Cont_SWS_H_Posterior_REM_Nor/temp;
clear temp




%AWA
clear i j Int m n 
for m=1:(length(AWA_H(:,1))) %Steps hour-hour (24 total)
    j=2;
    for i=1:(length(AWA_H)-1) %Throughout 
        if (AWA_H(m,j))-(AWA_H(m,i))==1 %Checking subsequent values
            Cont_AWA_H(m,i)=1; %If is subsequent, value 1
        else
            Cont_AWA_H(m,i)=0; %if not subsequent, value 0 
        end
        j=j+1;
    end
end

clear i j Int m n
for m=1:(length(AWA_H(:,1)))
    Cont_AWA_Total_H(1,m)=sum(Cont_AWA_H(m,:)); %Total continuous periods per hour
    Cont_AWA_Nor_Total_H(1,m)=(Cont_AWA_Total_H(1,m)./(awa))*100; %Normalization by the total state time
    Cont_AWA_Nor_H(1,m)=(Cont_AWA_Total_H(1,m)./AWA_H_Count(m,1))*100; %Normalization by the total hour event state time
end

%Calculating the time period of each block
clear i j Temp m n 
for m=1:(length(AWA_H(:,1))) %Steps hour-hour (24 total)
    j=1;
    for i=1:((length(Cont_AWA_H))-1) %Covering up the whole length
        if (Cont_AWA_H(m,i)+Cont_AWA_H(m,(i+1)))>1 %Suming up the period duration if it is sequential
            Period_AWA_H(m,i)=0; %The previous value is 0
            Period_AWA_H(m,(i+1))=j+(Cont_AWA_H(m,(i+1))); %The subsequent value is the sum
            j=j+1;
        else
            Period_AWA_H(m,(i+1))=0; %The previous value is 0
            j=1; %If it is not sequential, then mantain the value
        end
    end
end

%Finding the mean AWA period per hour
clear i j Temp m n
for m=1:(length(AWA_H(:,1))) %Steps hour-hour (24 total)
    Temp = find(Period_AWA_H(m,:)>0); %Finding sequential periods (over zero)
    if isempty(Temp)==0 %If the period is not zero
        for i=1:(length(Temp))
            Period_AWA_MEAN_H(m,i)=Period_AWA_H(m,(Temp(1,i))); %Getting the values
        end
        Period_AWA_Mean_H(m,1)=(sum(Period_AWA_MEAN_H(m,:)))/(length(Temp)); %Getting the mean 
    else
        Period_AWA_Mean_H(m,1)=0;
    end
    clear Temp
end

%Finding the max period of AWA sleep
for m=1:(length(AWA_H(:,1))) %Steps hour-hour (24 total)
    Period_AWA_Max_H(m,1)=max(Period_AWA_H(m,:)); %Just getting the mean value
end
            


%Finding the number of AWA countinuous episodes
for i=1:length (Period_AWA_H(:,1)) %24h steps
    for j=1:length (Period_AWA_H(1,:)) %time stamps steps
        if Period_AWA_H(i,j)==0
            temp(i,j)=0;
        else
            temp(i,j)=1;
        end
    end
    Cont_Number_AWA(i,1)=sum(temp(i,:));
end
clear temp





%Finding the AWA countinuous PREVIOUS and POSTERIOR state
for i=1:length(Period_AWA_H(:,1)) % steps de 24 Horas
temp=find((Period_AWA_H(i,:))>0); %steps timestamps
for j=1:length(temp)
if (AWA_H(i,(temp(1,j))))-(Period_AWA_H(i,(temp(1,j))))==0 %Evitar inicio menor ou igual a zero
Cont_AWA_H_Previous(i,j)=NaN;
else
Cont_AWA_H_Previous(i,j)=result((AWA_H(i,(temp(1,j))))-(Period_AWA_H(i,(temp(1,j))))); %state before continous period
end
if ((AWA_H(i,(temp(1,j))))+2)>length (result) %evitar final maior que length (result)
Cont_AWA_H_Posterior(i,j)=NaN;
else
Cont_AWA_H_Posterior(i,j)=result(1,((AWA_H(i,(temp(1,j))))+2)); %state after continous period
end
end
clear temp
end


Cont_AWA_H_Previous(Cont_AWA_H_Previous==0) = NaN;
Cont_AWA_H_Posterior(Cont_AWA_H_Posterior==0) = NaN;


%Finding PREVIOUS normalized
for i=1:length (Cont_AWA_H_Previous(:,1)) %hours
for j=1:length(Cont_AWA_H_Previous(1,:))%Steps
if Cont_AWA_H_Previous(i,j)==2
Cont_AWA_H_Previous_SWS(i,j)=1;
else
Cont_AWA_H_Previous_SWS(i,j)=0;
end
if Cont_AWA_H_Previous(i,j)==1
Cont_AWA_H_Previous_REM(i,j)=1;
else
Cont_AWA_H_Previous_REM(i,j)=0;
end
end
end


for i=1:length (Cont_AWA_H_Previous(:,1)) %hours
for j=1:length(Cont_AWA_H_Previous(1,:))%Steps
if Cont_AWA_H_Previous(i,j)==3
    Cont_AWA_H_Previous(i,j)=NaN;
end
end
end


temp=(Cont_AWA_H_Previous); %Finding the total event at that hour
temp=(isnan(temp));
temp=length(temp(:,1))-sum(temp);
temp=sum(temp); %total valid events
Cont_AWA_H_Previous_SWS_Nor=sum(Cont_AWA_H_Previous_SWS);
Cont_AWA_H_Previous_SWS_Nor=sum(Cont_AWA_H_Previous_SWS_Nor);
Cont_AWA_H_Previous_SWS_Nor=Cont_AWA_H_Previous_SWS_Nor/temp;
Cont_AWA_H_Previous_REM_Nor=sum(Cont_AWA_H_Previous_REM);
Cont_AWA_H_Previous_REM_Nor=sum(Cont_AWA_H_Previous_REM_Nor);
Cont_AWA_H_Previous_REM_Nor=Cont_AWA_H_Previous_REM_Nor/temp;
clear temp

%Finding POSTERIOR normalized
for i=1:length (Cont_AWA_H_Posterior(:,1)) %hours
for j=1:length(Cont_AWA_H_Posterior(1,:))%Steps
if Cont_AWA_H_Posterior(i,j)==2
Cont_AWA_H_Posterior_SWS(i,j)=1;
else
Cont_AWA_H_Posterior_SWS(i,j)=0;
end
if Cont_AWA_H_Posterior(i,j)==1
Cont_AWA_H_Posterior_REM(i,j)=1;
else
Cont_AWA_H_Posterior_REM(i,j)=0;
end
end
end


for i=1:length (Cont_AWA_H_Posterior(:,1)) %hours
for j=1:length(Cont_AWA_H_Posterior(1,:))%Steps
if Cont_AWA_H_Posterior(i,j)==3
    Cont_AWA_H_Posterior(i,j)=NaN;
end
end
end


temp=(Cont_AWA_H_Posterior); %Finding the total event at that hour
temp=(isnan(temp));
temp=length(temp(:,1))-sum(temp);
temp=sum(temp); %total valid events
Cont_AWA_H_Posterior_SWS_Nor=sum(Cont_AWA_H_Posterior_SWS);
Cont_AWA_H_Posterior_SWS_Nor=sum(Cont_AWA_H_Posterior_SWS_Nor);
Cont_AWA_H_Posterior_SWS_Nor=Cont_AWA_H_Posterior_SWS_Nor/temp;
Cont_AWA_H_Posterior_REM_Nor=sum(Cont_AWA_H_Posterior_REM);
Cont_AWA_H_Posterior_REM_Nor=sum(Cont_AWA_H_Posterior_REM_Nor);
Cont_AWA_H_Posterior_REM_Nor=Cont_AWA_H_Posterior_REM_Nor/temp;
clear temp









%Histogram of number of episodes for each state (minimal 20s window)


inicial=2; %minimum value 
step=2; %range of the histogram bouts

%REM histogram episodes 

%REM Night
hist_REM_temp=Period_REM_H(1:12,:); 
hist_REM_temp = hist_REM_temp(hist_REM_temp>0);

if isempty (hist_REM_temp)==0
    M_Hist = max( hist_REM_temp ); 
    edges=[inicial:step:M_Hist+1]; 
    [hist_REM_Night,edges] = histcounts(hist_REM_temp,edges);
    clear M_Hist edges hist_REM_temp
else
    hist_REM_Night=0;
end

%REM Day
hist_REM_temp = Period_REM_H(13:24,:);
hist_REM_temp = hist_REM_temp(hist_REM_temp>0);

if isempty (hist_REM_temp)==0
    M_Hist = max( hist_REM_temp );
    edges=[inicial:step:M_Hist+1];
    [hist_REM_Day,edges] = histcounts(hist_REM_temp,edges);
    clear M_Hist edges hist_REM_temp
else
    hist_REM_Day=0;
end
   

%SWS histogram episodes 

%SWS Night
hist_SWS_temp=Period_SWS_H(1:12,:);
hist_SWS_temp = hist_SWS_temp(hist_SWS_temp>0);

if isempty (hist_SWS_temp)==0
    M_Hist = max( hist_SWS_temp );
    edges=[inicial:step:M_Hist+1];
    [hist_SWS_Night,edges] = histcounts(hist_SWS_temp,edges);
    clear M_Hist edges hist_SWS_temp
else
    hist_SWS_Night=0;
end

%SWS Day
hist_SWS_temp = Period_SWS_H(13:24,:);
hist_SWS_temp = hist_SWS_temp(hist_SWS_temp>0);

if isempty (hist_SWS_temp)==0
    M_Hist = max( hist_SWS_temp );
    edges=[inicial:step:M_Hist+1];
    [hist_SWS_Day,edges] = histcounts(hist_SWS_temp,edges);
    clear M_Hist edges hist_SWS_temp
else
    hist_SWS_Day=0;
end


%AWA histogram episodes 

%AWA Night
hist_AWA_temp=Period_AWA_H(1:12,:);
hist_AWA_temp = hist_AWA_temp(hist_AWA_temp>0);

if isempty (hist_AWA_temp)==0
    M_Hist = max( hist_AWA_temp );
    edges=[inicial:step:M_Hist+1];
    [hist_AWA_Night,edges] = histcounts(hist_AWA_temp,edges);
    clear M_Hist edges hist_SWS_temp
else
    hist_AWA_Night=0;
end

%AWA Day
hist_AWA_temp = Period_AWA_H(13:24,:);
hist_AWA_temp = hist_AWA_temp(hist_AWA_temp>0);

if isempty (hist_AWA_temp)==0
    M_Hist = max( hist_AWA_temp );
    edges=[inicial:step:M_Hist+1];
    [hist_AWA_Day,edges] = histcounts(hist_AWA_temp,edges);
    clear M_Hist edges hist_AWA_temp
else
    hist_AWA_Day=0;
end





%###OUTPUT###
%Percentage of time in continous pattern
States_Cont_Nor_H(:,1)=Cont_AWA_Nor_H';
States_Cont_Nor_H(:,2)=Cont_SWS_Nor_H';
States_Cont_Nor_H(:,3)=Cont_REM_Nor_H';

States_Cont_Nor_H = changem(States_Cont_Nor_H,NaN,0);%Changing the value zero for NaN

%Mean continous period
States_Period_Mean_H(:,1)=Period_AWA_Mean_H;
States_Period_Mean_H(:,2)=Period_SWS_Mean_H;
States_Period_Mean_H(:,3)=Period_REM_Mean_H;

States_Period_Mean_H = changem(States_Period_Mean_H,NaN,0);%Changing the value zero for NaN

%Maximum continous period
States_Period_Max_H(:,1)=Period_AWA_Max_H;
States_Period_Max_H(:,2)=Period_SWS_Max_H;
States_Period_Max_H(:,3)=Period_REM_Max_H;

States_Period_Max_H = changem(States_Period_Max_H,NaN,0);%Changing the value zero for NaN


%Previous and Posterior state referenced at a continous state.
%AWA Prev and post SWS REM; SWS prev and post AWA REM; REM prev and post
%AWA SWS
States_Cont_Prev_Post_H=[Cont_AWA_H_Previous_SWS_Nor Cont_AWA_H_Previous_REM_Nor ...
    Cont_AWA_H_Posterior_SWS_Nor Cont_AWA_H_Posterior_REM_Nor Cont_SWS_H_Previous_AWA_Nor ...
    Cont_SWS_H_Previous_REM_Nor Cont_SWS_H_Posterior_AWA_Nor Cont_SWS_H_Posterior_REM_Nor ...
    Cont_REM_H_Previous_AWA_Nor Cont_REM_H_Previous_SWS_Nor Cont_REM_H_Posterior_AWA_Nor ...
    Cont_REM_H_Posterior_SWS_Nor];
States_Cont_Prev_Post_H=States_Cont_Prev_Post_H';


% Number of continous episodes
States_Cont_Number = [Cont_Number_AWA Cont_Number_SWS Cont_Number_REM];


%Histogram of Number of epidoes
States_Hist_AWA_N = hist_AWA_Night';
States_Hist_AWA_D = hist_AWA_Day';

States_Hist_SWS_N = hist_SWS_Night';
States_Hist_SWS_D = hist_SWS_Day';

States_Hist_REM_N = hist_REM_Night';
States_Hist_REM_D = hist_REM_Day'; 








% ***POWER and POWER_norm (FFT) per hour***

%FFT of the spectrogram
steps=(length(PLFP))/24; %Aquiring the periods per hour
j=steps;
m=1;
for i=1:steps:(length(PLFP))
    States_FFT_PLFP(m,:)=mean(PLFP(:,(i:j))'); %FFT by the mean of normalized spectrogram
    j=j+steps;
    m=m+1;
end
clear m j i steps

%FFT of the normalized spectrogram
steps=(length(PLFP_norm))/24; %Aquiring the periods per hour
j=steps;
m=1;
for i=1:steps:(length(PLFP_norm))
    States_FFT_PLFP_norm(m,:)=mean(PLFP_norm(:,(i:j))'); %FFT by the mean of normalized spectrogram
    j=j+steps;
    m=m+1;
end
clear m j i steps

%###OUTPUT###
States_FFT_PLFP_All = States_FFT_PLFP(:,1:240);
States_FFT_PLFP_norm_All = States_FFT_PLFP_norm(:,1:240);







% ***POWER_norm in BANDS per hour in each STATE***

%Range Scpetrogram
R1=1;
R2=240;
%Delta Frequencies
DF1=2;
DF2=4;
%Theta Frequencies
TF1=7; %Theta Frequency one
TF2=9; %Theta Frequency two
%Beta Frequency
BF1=15;
BF2=30;
%Gamma Frequencies
GF1=30;
GF2=80;
%High Gamma Frequencies
HGF1=80;
HGF2=140;
%Ripple frequencies
RF1=140;
RF2=220;


%Finding each spectrogram for each REM hour
clear Temp Temp_1 P_REM_H
Pnorm_REM_H=cell((length(REM_H(:,1))),1); %prellocate
for m=1:(length(REM_H(:,1))) %Steps hour-hour (24 total)
    Temp = find(REM_H(m,:)>0);%Getting the matrix positions of values over 0
    for i=1:(length(Temp))
        Temp_1(1,i)=REM_H(m,(Temp(1,i))); %Getting the time stamps
        Pnorm_REM_H{m,i}=PLFP_norm(:,(Temp_1(1,i))); %Getting the normalized spectral power
    end
    clear Temp_1 Temp
end


%Finding avalibles cells per hour
clear Temp Temp_1 temp j
j=1;
Pnorm_REM_H_mean=NaN(24,240); %prellocate
for m=1:(length(REM_H(:,1))) %Steps hour-hour (24 total)
    for i=1:(length(Pnorm_REM_H(m,:))) 
        if isempty(Pnorm_REM_H{m,i})==0 
            temp(:,j)=(Pnorm_REM_H{m,i}(R1:R2));
        else 
            temp(:,j)=0;
        end
        j=j+1; %passo
    end
    
    temp=temp';
    
    n_media=sum(temp)./(nnz((temp(:,1)))); 
    
    Pnorm_REM_H_mean(m,:)=n_media; 
    
    clear temp j n_media
    j=1;
    
end





%Finding each spectrogram for each SWS hour
clear Temp Temp_1 P_SWS_H
Pnorm_SWS_H=cell((length(SWS_H(:,1))),1); %prellocate
for m=1:(length(SWS_H(:,1))) %Steps hour-hour (24 total)
    Temp = find(SWS_H(m,:)>0);%Getting the matrix positions of values over 0
    for i=1:(length(Temp))
        Temp_1(1,i)=SWS_H(m,(Temp(1,i))); %Getting the time stamps
        Pnorm_SWS_H{m,i}=PLFP_norm(:,(Temp_1(1,i))); %Getting the normalized spectral power
    end
    clear Temp_1 Temp
end


%Finding avalibles cells per hour
clear Temp Temp_1 temp j
j=1;
Pnorm_SWS_H_mean=NaN(24,240); %prellocate
for m=1:(length(SWS_H(:,1))) %Steps hour-hour (24 total)
    for i=1:(length(Pnorm_SWS_H(m,:))) %Percorrendo os time-stamps da hora
        if isempty(Pnorm_SWS_H{m,i})==0 %se houver matriz, copie em temp
            temp(:,j)=(Pnorm_SWS_H{m,i}(R1:R2));
        else 
            temp(:,j)=0;
        end
        j=j+1;
    end
    
    temp=temp';
    
    n_media=sum(temp)./(nnz((temp(:,1)))); 
    
    Pnorm_SWS_H_mean(m,:)=n_media; 
    
    clear temp j n_media
    j=1;
    
end





%Finding each spectrogram for each AWA hour
clear Temp Temp_1 P_AWA_H
Pnorm_AWA_H=cell((length(AWA_H(:,1))),1); %prellocate
for m=1:(length(AWA_H(:,1))) %Steps hour-hour (24 total)
    Temp = find(AWA_H(m,:)>0);%Getting the matrix positions of values over 0
    for i=1:(length(Temp))
        Temp_1(1,i)=AWA_H(m,(Temp(1,i))); %Getting the time stamps
        Pnorm_AWA_H{m,i}=PLFP_norm(:,(Temp_1(1,i))); %Getting the normalized spectral power
    end
    clear Temp_1 Temp
end


%Finding avalibles cells per hour
clear Temp Temp_1 temp j
j=1;
Pnorm_AWA_H_mean=NaN(24,240); %prellocate
for m=1:(length(AWA_H(:,1))) %Steps hour-hour (24 total)
    for i=1:(length(Pnorm_AWA_H(m,:))) 
        if isempty(Pnorm_AWA_H{m,i})==0 
            temp(:,j)=(Pnorm_AWA_H{m,i}(R1:R2));
        else 
            temp(:,j)=0;
        end
        j=j+1; 
    end
    
    temp=temp';
    
    n_media=sum(temp)./(nnz((temp(:,1)))); 
    
    Pnorm_AWA_H_mean(m,:)=n_media; 
    
    clear temp j n_media
    j=1;
    
end




%###OUTPUT###
%Separated by bands

Delta_H_norm(:,1)=mean(Pnorm_AWA_H_mean(:,DF1:DF2)');
Delta_H_norm(:,2)=mean(Pnorm_SWS_H_mean(:,DF1:DF2)');
Delta_H_norm(:,3)=mean(Pnorm_REM_H_mean(:,DF1:DF2)');


Theta_H_norm(:,1)=mean(Pnorm_AWA_H_mean(:,TF1:TF2)');
Theta_H_norm(:,2)=mean(Pnorm_SWS_H_mean(:,TF1:TF2)');
Theta_H_norm(:,3)=mean(Pnorm_REM_H_mean(:,TF1:TF2)');


Beta_H_norm(:,1)=mean(Pnorm_AWA_H_mean(:,BF1:BF2)');
Beta_H_norm(:,2)=mean(Pnorm_SWS_H_mean(:,BF1:BF2)');
Beta_H_norm(:,3)=mean(Pnorm_REM_H_mean(:,BF1:BF2)');


Gamma_H_norm(:,1)=mean(Pnorm_AWA_H_mean(:,GF1:GF2)');
Gamma_H_norm(:,2)=mean(Pnorm_SWS_H_mean(:,GF1:GF2)');
Gamma_H_norm(:,3)=mean(Pnorm_REM_H_mean(:,GF1:GF2)');


High_Gamma_H_norm(:,1)=mean(Pnorm_AWA_H_mean(:,HGF1:HGF2)');
High_Gamma_H_norm(:,2)=mean(Pnorm_SWS_H_mean(:,HGF1:HGF2)');
High_Gamma_H_norm(:,3)=mean(Pnorm_REM_H_mean(:,HGF1:HGF2)');


Ripple_H_norm(:,1)=mean(Pnorm_AWA_H_mean(:,RF1:RF2)');
Ripple_H_norm(:,2)=mean(Pnorm_SWS_H_mean(:,RF1:RF2)');
Ripple_H_norm(:,3)=mean(Pnorm_REM_H_mean(:,RF1:RF2)');

States_PB_H_norm=[Delta_H_norm Theta_H_norm Beta_H_norm Gamma_H_norm High_Gamma_H_norm Ripple_H_norm];


%The FFT separated by states
States_FFT_PLFP_norm  = [Pnorm_AWA_H_mean; Pnorm_SWS_H_mean; Pnorm_REM_H_mean];






% ***POWER in BANDS per hour in each STATE***


%Finding each spectrogram for each REM hour
clear Temp Temp_1 P_REM_H
P_REM_H=cell((length(REM_H(:,1))),1); %prellocate
for m=1:(length(REM_H(:,1))) %Steps hour-hour (24 total)
    Temp = find(REM_H(m,:)>0);%Getting the matrix positions of values over 0
    for i=1:(length(Temp))
        Temp_1(1,i)=REM_H(m,(Temp(1,i))); %Getting the time stamps
        P_REM_H{m,i}=PLFP(:,(Temp_1(1,i))); %Getting the normalized spectral power
    end
    clear Temp_1 Temp
end


%Finding avalibles cells per hour
clear Temp Temp_1 temp j
j=1;
P_REM_H_mean=NaN(24,240);%Prellocate
for m=1:(length(REM_H(:,1))) %Steps hour-hour (24 total)
    for i=1:(length(P_REM_H(m,:))) 
        if isempty(P_REM_H{m,i})==0 
            temp(:,j)=(P_REM_H{m,i}(R1:R2));
        else 
            temp(:,j)=0;
        end
        j=j+1; 
    end
    
    temp=temp';
    
    n_media=sum(temp)./(nnz((temp(:,1)))); 
    
    P_REM_H_mean(m,:)=n_media; 
    
    clear temp j n_media
    j=1;
    
end




%Finding each spectrogram for each SWS hour
clear Temp Temp_1 P_SWS_H
P_SWS_H=cell((length(SWS_H(:,1))),1); %prellocate
for m=1:(length(SWS_H(:,1))) %Steps hour-hour (24 total)
    Temp = find(SWS_H(m,:)>0);%Getting the matrix positions of values over 0
    for i=1:(length(Temp))
        Temp_1(1,i)=SWS_H(m,(Temp(1,i))); %Getting the time stamps
        P_SWS_H{m,i}=PLFP(:,(Temp_1(1,i))); %Getting the normalized spectral power
    end
    clear Temp_1 Temp
end


%Finding avalibles cells per hour
clear Temp Temp_1 temp j
j=1;
P_SWS_H_mean=NaN(24,240);%Prellocate
for m=1:(length(SWS_H(:,1))) %Steps hour-hour (24 total)
    for i=1:(length(P_SWS_H(m,:))) 
        if isempty(P_SWS_H{m,i})==0 
            temp(:,j)=(P_SWS_H{m,i}(R1:R2));
        else 
            temp(:,j)=0;
        end
        j=j+1; 
    end
    
    temp=temp';
    
    n_media=sum(temp)./(nnz((temp(:,1)))); 
    
    P_SWS_H_mean(m,:)=n_media; 
    
    clear temp j n_media
    j=1;
    
end




%Finding each spectrogram for each AWA hour
clear Temp Temp_1 P_AWA_H
P_AWA_H=cell((length(AWA_H(:,1))),1); %prellocate
for m=1:(length(AWA_H(:,1))) %Steps hour-hour (24 total)
    Temp = find(AWA_H(m,:)>0);%Getting the matrix positions of values over 0
    for i=1:(length(Temp))
        Temp_1(1,i)=AWA_H(m,(Temp(1,i))); %Getting the time stamps
        P_AWA_H{m,i}=PLFP(:,(Temp_1(1,i))); %Getting the normalized spectral power
    end
    clear Temp_1 Temp
end


%Finding avalibles cells per hour
clear Temp Temp_1 temp j
j=1;
P_AWA_H_mean=NaN(24,240);%Prellocate
for m=1:(length(AWA_H(:,1))) %Steps hour-hour (24 total)
    for i=1:(length(P_AWA_H(m,:))) 
        if isempty(P_AWA_H{m,i})==0 
            temp(:,j)=(P_AWA_H{m,i}(R1:R2));
        else 
            temp(:,j)=0;
        end
        j=j+1; 
    end
    
    temp=temp';
    
    n_media=sum(temp)./(nnz((temp(:,1)))); 
    
    P_AWA_H_mean(m,:)=n_media; 
    
    clear temp j n_media
    j=1;
    
end




%###OUTPUT###
%Separated by bands

Delta_H(:,1)=mean(P_AWA_H_mean(:,DF1:DF2)');
Delta_H(:,2)=mean(P_SWS_H_mean(:,DF1:DF2)');
Delta_H(:,3)=mean(P_REM_H_mean(:,DF1:DF2)');


Theta_H(:,1)=mean(P_AWA_H_mean(:,TF1:TF2)');
Theta_H(:,2)=mean(P_SWS_H_mean(:,TF1:TF2)');
Theta_H(:,3)=mean(P_REM_H_mean(:,TF1:TF2)');


Beta_H(:,1)=mean(P_AWA_H_mean(:,BF1:BF2)');
Beta_H(:,2)=mean(P_SWS_H_mean(:,BF1:BF2)');
Beta_H(:,3)=mean(P_REM_H_mean(:,BF1:BF2)');


Gamma_H(:,1)=mean(P_AWA_H_mean(:,GF1:GF2)');
Gamma_H(:,2)=mean(P_SWS_H_mean(:,GF1:GF2)');
Gamma_H(:,3)=mean(P_REM_H_mean(:,GF1:GF2)');


High_Gamma_H(:,1)=mean(P_AWA_H_mean(:,HGF1:HGF2)');
High_Gamma_H(:,2)=mean(P_SWS_H_mean(:,HGF1:HGF2)');
High_Gamma_H(:,3)=mean(P_REM_H_mean(:,HGF1:HGF2)');


Ripple_H(:,1)=mean(P_AWA_H_mean(:,RF1:RF2)');
Ripple_H(:,2)=mean(P_SWS_H_mean(:,RF1:RF2)');
Ripple_H(:,3)=mean(P_REM_H_mean(:,RF1:RF2)');

States_PB_H=[Delta_H Theta_H Beta_H Gamma_H High_Gamma_H Ripple_H];


%The FFT separated by states
States_FFT_PLFP = [P_AWA_H_mean; P_SWS_H_mean; P_REM_H_mean];








%   ***COMUDOLOGRAM TORT highest theta REM value***

%Comodulogram of the highest theta power bin during REM
% Define phase bins
nbin = 18; % number of phase bins
position=zeros(1,nbin); % this variable will get the beginning (not the center) of each phase bin (in rads)
winsize = 2*pi/nbin;
for j=1:nbin 
    position(j) = -pi+(j-1)*winsize; 
end

%Band analysis definition
%Theta
Pf1 = 7;
Pf2 = 9;
%Gamma
Af1 = 30;
Af2 = 80;

%Getting the maximum theta during immobility (REM)
ZHRMS=zscore(HRMS);%Calculating the LFP zscore
rep=ZRMS.*ZHRMS; 
[Y,I] = min(rep); %representative bim

lfp=LFP(((I-1)*FST):((I)*FST)); %Selectin one bin (10s) before anf after
data_length = length(lfp); %Adapting to TORT script
srate = FS; %Adapting to TORT script
dt = 1/srate; %Adapting to TORT script
t = (1:data_length)*dt; %Adapting to TORT script

%Define the amplitude- and phase-frequencies
PhaseFreqVector=2:2:50;
AmpFreqVector=10:5:200;
PhaseFreq_BandWidth=4;
AmpFreq_BandWidth=20;

% Filtering and Hilbert transform
Comodulogram=single(zeros(length(PhaseFreqVector),length(AmpFreqVector)));
AmpFreqTransformed = zeros(length(AmpFreqVector), data_length);
PhaseFreqTransformed = zeros(length(PhaseFreqVector), data_length);

for ii=1:length(AmpFreqVector)
    Af1 = AmpFreqVector(ii);
    Af2=Af1+AmpFreq_BandWidth;
    AmpFreq=eegfilt(lfp,srate,Af1,Af2); % filtering
    AmpFreqTransformed(ii, :) = abs(hilbert(AmpFreq)); % getting the amplitude envelope
end

for jj=1:length(PhaseFreqVector)
    Pf1 = PhaseFreqVector(jj);
    Pf2 = Pf1 + PhaseFreq_BandWidth;
    PhaseFreq=eegfilt(lfp,srate,Pf1,Pf2); % filtering 
    PhaseFreqTransformed(jj, :) = angle(hilbert(PhaseFreq)); % getting the phase time series
end

%Compute MI and comodulogram
counter1=0;
for ii=1:length(PhaseFreqVector)
counter1=counter1+1;

    Pf1 = PhaseFreqVector(ii);
    Pf2 = Pf1+PhaseFreq_BandWidth;
    
    counter2=0;
    for jj=1:length(AmpFreqVector)
    counter2=counter2+1;
    
        Af1 = AmpFreqVector(jj);
        Af2 = Af1+AmpFreq_BandWidth;
        [MI,MeanAmp]=ModIndex_v2(PhaseFreqTransformed(ii, :), AmpFreqTransformed(jj, :), position);
        Comodulogram(counter1,counter2)=MI;
    end
end


%FIGURE
%Plot comodulogram
figure;
clf
contourf(PhaseFreqVector+PhaseFreq_BandWidth/2,AmpFreqVector+AmpFreq_BandWidth/2,Comodulogram',30,'lines','none')
set(gca,'fontsize',14)
ylabel('Amplitude Frequency (Hz)');ylim([20 100]);
xlabel('Phase Frequency (Hz)');xlim([4 30]);
colorbar








%   ***CFC per hour*** 

... The CFC is calculated on the non-filtered signal.
% Define phase bins
nbin = 18; % number of phase bins
position=zeros(1,nbin); % this variable will get the beginning (not the center) of each phase bin (in rads)
winsize = 2*pi/nbin;
for j=1:nbin 
    position(j) = -pi+(j-1)*winsize; 
end

%Band analysis definition for THETA-gamma coupling
%Theta
Pf1 = 7;
Pf2 = 9;
%Gamma
Af1 = 60;
Af2 = 80;


%REM
MIraw_REM_H=NaN((length(REM_H(:,1))),1);
for m=1:1:(length(REM_H(:,1)))
    for i=1:1:(length(REM_H))
        if (REM_H(m,i))>1
            period=LFP(((REM_H(m,i))*FST)-(FST):(((REM_H(m,i))*FST)));%Selecting the period of 10s
            [MI,MeanAmp] = ModIndex_v1(period,FS,Pf1,Pf2,Af1,Af2,position);
            MIraw_REM_H(m,i)=MI(1,1);
        end
    end
end
%Mean CFC per hour
for m=1:1:(length(REM_H(:,1)))
    Temp=find(MIraw_REM_H(m,:)>0); %Getting the values over 0
    Temp(1,:)=(MIraw_REM_H(m,Temp(1,:))); 
    MIraw_REM_AV_H(m,1)=mean(Temp);
end
clear temp

%SWS
MIraw_SWS_H=NaN((length(SWS_H(:,1))),1);
for m=1:1:(length(SWS_H(:,1)))
    for i=1:1:(length(SWS_H))
        if (SWS_H(m,i))>1
            period=LFP(((SWS_H(m,i))*FST)-(FST):(((SWS_H(m,i))*FST)));%Selecting the period of 10s
            [MI,MeanAmp] = ModIndex_v1(period,FS,Pf1,Pf2,Af1,Af2,position);
            MIraw_SWS_H(m,i)=MI(1,1);
        end
    end
end
%Mean CFC per hour
for m=1:1:(length(SWS_H(:,1)))
    Temp=find(MIraw_SWS_H(m,:)>0); %Getting the values over 0
    Temp(1,:)=(MIraw_SWS_H(m,Temp(1,:))); 
    MIraw_SWS_AV_H(m,1)=mean(Temp);
end
clear temp


%AWA
MIraw_AWA_H=NaN((length(AWA_H(:,1))),1);
for m=1:1:(length(AWA_H(:,1)))
    for i=1:1:(length(AWA_H))
        if (AWA_H(m,i))>1
            period=LFP(((AWA_H(m,i))*FST)-(FST):(((AWA_H(m,i))*FST)));%Selecting the period of 10s
            [MI,MeanAmp] = ModIndex_v1(period,FS,Pf1,Pf2,Af1,Af2,position);
            MIraw_AWA_H(m,i)=MI(1,1);
        end
    end
end
%Mean CFC per hour
for m=1:1:(length(AWA_H(:,1)))
    Temp=find(MIraw_AWA_H(m,:)>0); %Getting the values over 0
    Temp(1,:)=(MIraw_AWA_H(m,Temp(1,:))); 
    MIraw_AWA_AV_H(m,1)=mean(Temp);
end
clear temp



%###OUTPU###

States_MIraw_AV_H_TG(:,1)=MIraw_AWA_AV_H;
States_MIraw_AV_H_TG(:,2)=MIraw_SWS_AV_H;
States_MIraw_AV_H_TG(:,3)=MIraw_REM_AV_H;








%***Spindle analysis*** 211120

%Building a continuous AWA, SWS and REM LFP_F variable.

AWA_LFP=NaN(1,(length(AWA)*FST)); %Prellocating AWA LFP
m=1;
for i=1:length(result)
    if result(i)==3 %checking which time stamp is AWA
        AWA_LFP(((m*FST)-(FST-1)):(m*FST))=LFP (((i*FST)-(FST-1)):(i*FST)); %collecting the LFP
        m=m+1;
    end
end
clear i m 


SWS_LFP=NaN(1,(length(SWS)*FST)); %Prellocating SWS LFP
m=1;
for i=1:length(result)
    if result(i)==2 %checking which time stamp is SWS
        SWS_LFP(((m*FST)-(FST-1)):(m*FST))=LFP_F (((i*FST)-(FST-1)):(i*FST)); %collecting the LFP
        m=m+1;
    end
end
clear i m 


REM_LFP=NaN(1,(length(REM)*FST)); %Prellocating REM LFP
m=1;
for i=1:length(result)
    if result(i)==1 %checking which time stamp is REM
        REM_LFP(((m*FST)-(FST-1)):(m*FST))=LFP (((i*FST)-(FST-1)):(i*FST)); %collecting the LFP
        m=m+1;
    end
end
clear i m 








%SPINDLES 
%based on the work "Validation of an automated sleep spindle detection method for
%mouse electroencephalography" article

%Building a continuous SWS (whole state) from LFP_F variable.
%OBS: This part uses the whole SWS state to create the threshold that is
%going to be used for the night and day periods.

SWS_LFP=NaN(1,(length(SWS)*FST)); %Prellocating SWS LFP
m=1;
for i=1:length(result)
    if result(i)==2 %checking which time stamp is SWS
        SWS_LFP(((m*FST)-(FST-1)):(m*FST))=LFP_F (((i*FST)-(FST-1)):(i*FST)); %collecting the LFP
        m=m+1;
    end
end
clear i m 


Filter_Spindle=designfilt('bandpassiir','FilterOrder',20,'HalfPowerFrequency1',10,'HalfPowerFrequency2',15,'SampleRate',FS);
SWS_LFP_F=filtfilt(Filter_Spindle,SWS_LFP);

min=FS*0.5; %minimal spindle length in seconds
max=FS*5; %max spindle lenght in seconds 

[yupper,ylower] = envelope(SWS_LFP_F,FS*0.7,'rms'); 

env=yupper.^3; 

limiar=2.0*(mean(env)); 

clear yupper env










%NIGHT
result_Night=result(1:length(result)/2); %getting SWS only at Night
m=1;
for i=1:length(result_Night)
    if result_Night(i)==2 %checking which time stamp is SWS
        SWS_LFP_Night(((m*FST)-(FST-1)):(m*FST))=LFP_F (((i*FST)-(FST-1)):(i*FST)); %collecting the LFP
        m=m+1;
    end
end
clear i m 

%SPINDLES Night
SWS_LFP_Night_F=filtfilt(Filter_Spindle,SWS_LFP_Night);

[yupper,ylower] = envelope(SWS_LFP_Night_F,FS*0.7,'rms'); 

env=yupper.^3; 

%The same limiar is going to be used.

%Checking the point over the limiar
for i=1:length(env)
    if env(i)>limiar
        temp_limiar(i)=1;
    else
        temp_limiar(i)=0;
    end
end

%Getting the spinle length by summing the points over the limiar
j=1;
for i=1:((length(temp_limiar))-1) %Covering up the whole length
        if temp_limiar(i)+temp_limiar(i+1)>1 %Suming up the period duration if it is sequential
            
            temp_continuo(i)=0; %The previous value is 0
            temp_continuo(i+1)=j+(temp_limiar(i+1)); %The subsequent value is the sum
            j=j+1;
            
        else
            temp_continuo(i+1)=0; %The previous value is 0
            j=1; %If it is not sequential, then mantain the value
        end
end 

%Selecting the spindles larger than the minimal lenght and calculating the duration and number of occurance
for i=1:length(temp_continuo)
    if temp_continuo(i)>min %limiar de 0.5 segundos
        temp_dur(i)=temp_continuo(i); %Duration is the sum, calculated i the temp_continuo variable
        temp_sum(i)=1; %Number of episodes is the temp_continuo changed to value 1.
    else
        temp_dur(i)=0;
        temp_sum(i)=0;
    end
end

for i=1:length(temp_continuo)
    if temp_dur(i)<max %limiar de 0.5 segundos
        temp_dur(i)=temp_dur(i);
        temp_sum(i)=temp_sum(i);
    else
        temp_dur(i)=0;
        temp_sum(i)=0;
    end
end

for i=1:length(temp_dur)
    if temp_dur(i)>0 %everything over zero
        temp_SWS_F(1,:)=SWS_LFP_Night_F(1,((i+1)-(temp_dur(i))):i); %Selecting the spindle oscillation from the filtered SWS
        temp_peak=findpeaks(temp_SWS_F); %Finding the peaks of the oscillation
        temp_spindle_freq(1,i)=length(temp_peak)/(temp_dur(i)/FS); %Getting the number of the peaks
        clear temp_SWS_F temp_peak
    else
        temp_spindle_freq(1,i)=0;
    end
end
    



%Getting the data
total_sec_period=length(SWS_LFP_Night)/FS; %tempo total do periodo em segundos

spindle_sum=sum(temp_sum); %total number of oscillations episodes
spindle_dur=(mean(temp_dur(temp_dur>1)))/FS; %mean oscillations episode duration in seconds

spindle_den_sec=spindle_sum/total_sec_period; %normalizacao do numero de episodios por tempo total do periodo em segudnos
spindle_den_min=spindle_den_sec*60; %normalizacao do numero de episodios por tempo total do periodo em um minuto

spindle_peaks=(mean(temp_spindle_freq(temp_spindle_freq>1)));    
  

%OUTPUT NIGHT
State_Spindle(1,1)=spindle_sum; %Numero total de spindles
State_Spindle(2,1)=spindle_dur; %Duracao media dos spindles
State_Spindle(3,1)=spindle_den_min; %Spindles por minutos  
State_Spindle(4,1)=spindle_peaks; %Spindles por minutos 

clear temp_spindle_freq temp_peak temp_dur temp_sum temp_continuo j i temp_limiar env yupper temp_SWS_F
clear spindle_sum spindle_dur spindle_den_min spindle_peaks








%DAY
result_Day=result(length(result)/2:end); %Getting the SWS in Day State
result_Day_Temp=zeros(1,(length(result)/2)-1); %Summing up the Night period as zeros
result_Day=[result_Day_Temp result_Day]; %Summing up to have the right potision in the results variable

m=1;
for i=1:length(result_Day)
    if result_Day(i)==2 %checking which time stamp is SWS
        SWS_LFP_Day(((m*FST)-(FST-1)):(m*FST))=LFP_F (((i*FST)-(FST-1)):(i*FST)); %collecting the LFP
        m=m+1;
    end
end
clear i m 


%SPINDLES Day
SWS_LFP_Day_F=filtfilt(Filter_Spindle,SWS_LFP_Day);

[yupper,ylower] = envelope(SWS_LFP_Day_F,FS*0.7,'rms'); %artigo original usou 0.75segundos

env=yupper.^3; %elevacao cubica reduz basal e ressalta os valores positivos

%Definindo pontos acima do limiar
for i=1:length(env)
    if env(i)>limiar
        temp_limiar(i)=1;
    else
        temp_limiar(i)=0;
    end
end

%Definindo o tempo dos spindles (continuous) pela soma dos pontos continuos
%acima do limiar.
j=1;
for i=1:((length(temp_limiar))-1) %Covering up the whole length
        if temp_limiar(i)+temp_limiar(i+1)>1 %Suming up the period duration if it is sequential
            
            temp_continuo(i)=0; %The previous value is 0
            temp_continuo(i+1)=j+(temp_limiar(i+1)); %The subsequent value is the sum
            j=j+1;
            
        else
            temp_continuo(i+1)=0; %The previous value is 0
            j=1; %If it is not sequential, then mantain the value
        end
end 

%Selecionando os spindles que sao maiores que o tempo minimo e
%quantificando a duracao e o numero de ocorrencia.
for i=1:length(temp_continuo)
    if temp_continuo(i)>min %limiar de 0.5 segundos
        temp_dur(i)=temp_continuo(i); %Duration is the sum, calculated i the temp_continuo variable
        temp_sum(i)=1; %Number of episodes is the temp_continuo changed to value 1.
    else
        temp_dur(i)=0;
        temp_sum(i)=0;
    end
end

for i=1:length(temp_continuo)
    if temp_dur(i)<max %limiar de 0.5 segundos
        temp_dur(i)=temp_dur(i);
        temp_sum(i)=temp_sum(i);
    else
        temp_dur(i)=0;
        temp_sum(i)=0;
    end
end

for i=1:length(temp_dur)
    if temp_dur(i)>0 %everything over zero
        temp_SWS_F(1,:)=SWS_LFP_Day_F(1,((i+1)-(temp_dur(i))):i); %Selecting the spindle oscillation from the filtered SWS
        temp_peak=findpeaks(temp_SWS_F); %Finding the peaks of the oscillation
        temp_spindle_freq(1,i)=length(temp_peak)/(temp_dur(i)/FS); %Getting the number of the peaks
        clear temp_SWS_F temp_peak
    else
        temp_spindle_freq(1,i)=0;
    end
end
    



%Getting the data
total_sec_period=length(SWS_LFP_Day)/FS; %tempo total do periodo em segundos

spindle_sum=sum(temp_sum); %total number of oscillations episodes
spindle_dur=(mean(temp_dur(temp_dur>1)))/FS; %mean oscillations episode duration in seconds

spindle_den_sec=spindle_sum/total_sec_period; %normalizacao do numero de episodios por tempo total do periodo em segudnos
spindle_den_min=spindle_den_sec*60; %normalizacao do numero de episodios por tempo total do periodo em um minuto

spindle_peaks=(mean(temp_spindle_freq(temp_spindle_freq>1)));    
  


%OUTPUT Day
State_Spindle(5,1)=spindle_sum; %Numero total de spindles
State_Spindle(6,1)=spindle_dur; %Duracao media dos spindles
State_Spindle(7,1)=spindle_den_min; %Spindles por minutos  
State_Spindle(8,1)=spindle_peaks; %Spindles por minutos









