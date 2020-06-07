
%Align coordinate systems (simscape wrt v3d: x = y, y = -x, z = z)
x = [0 1 0;... % Rotation Matrix
    -1 0 0;...
    0 0 1];

%Format marker variables
LASIavg = x*(mean(LASI{1}(:,1:3))');
RASIavg = x*(mean(RASI{1}(:,1:3))');
LGTRavg = x*(mean(LGTR{1}(:,1:3))');
RGTRavg = x*(mean(RGTR{1}(:,1:3))');
LHCEavg = x*(mean(LHCE{1}(:,1:3))');
RHCEavg = x*(mean(RHCE{1}(:,1:3))');
LHLAavg = x*(mean(LHLA{1}(:,1:3))');
RHLAavg = x*(mean(RHLA{1}(:,1:3))');
LHMEavg = x*(mean(LHME{1}(:,1:3))');
RHMEavg = x*(mean(RHME{1}(:,1:3))');
LLEPavg = x*(mean(LLEP{1}(:,1:3))');
RLEPavg = x*(mean(RLEP{1}(:,1:3))');
LMEPavg = x*(mean(LMEP{1}(:,1:3))');
RMEPavg = x*(mean(RMEP{1}(:,1:3))');
LPATavg = x*(mean(LPAT{1}(:,1:3))');
RPATavg = x*(mean(RPAT{1}(:,1:3))');
LPSIavg = x*(mean(LPSI{1}(:,1:3))');
RPSIavg = x*(mean(RPSI{1}(:,1:3))');
LTCEavg = x*(mean(LTCE{1}(:,1:3))');
RTCEavg = x*(mean(RTCE{1}(:,1:3))');
LTLAavg = x*(mean(LTLA{1}(:,1:3))');
RTLAavg = x*(mean(RTLA{1}(:,1:3))');
LTMEavg = x*(mean(LTME{1}(:,1:3))');
RTMEavg = x*(mean(RTME{1}(:,1:3))');

%Create virtual markers
AJCL = [.5*(LHLAavg(1)+LHMEavg(1)),.5*(LHLAavg(2)+LHMEavg(2)),.5*(LHLAavg(3)+LHMEavg(3))];
KJCL = [.5*(LLEPavg(1)+LMEPavg(1)),.5*(LLEPavg(2)+LMEPavg(2)),.5*(LLEPavg(3)+LMEPavg(3))];
KJCR = [.5*(RLEPavg(1)+RMEPavg(1)),.5*(RLEPavg(2)+RMEPavg(2)),.5*(RLEPavg(3)+RMEPavg(3))];
ASISdist = (abs(RASIavg(1) - LASIavg(1))) + (abs(RASIavg(2) - LASIavg(2))) + (abs(RASIavg(3) - LASIavg(3)));
PSISdist = (abs(RPSIavg(1) - LPSIavg(1))) + (abs(RPSIavg(2) - LPSIavg(2))) + (abs(RPSIavg(3) - LPSIavg(3)));
ASISmidpt = .5*(RASIavg + LASIavg);
PSISmidpt = .5*(RPSIavg + LPSIavg);
HJCL = [LASIavg - (.3*ASISdist), LASIavg + (.14*ASISdist), LASIavg - (.22*ASISdist)]; %Bell and Brand regression model
HJCR = [RASIavg - (.3*ASISdist), RASIavg - (.14*ASISdist), RASIavg - (.22*ASISdist)]; %Bell and Brand regression model
ground = [0,0,0];

%Estimate segment lengths
footXdist = abs(LTCEavg(1) - LHCEavg(1))*100; %Convert from m to cm (matches simscape variable input units)
footYdist = abs(LTLAavg(2) - LTMEavg(2))*100;
footZdist = abs(AJCL(3) - ground(3))*100;

legLengthL = abs(HJCL(3) - KJCL(3))*100;
legLengthR = abs(HJCR(3) - KJCR(3))*100;

shankLengthL = abs(KJCL(3) - AJCL(3))*100;

pelvisXdist = abs(ASISmidpt(1) - PSISmidpt(1))*100;
pelvisYdist = abs(LGTRavg(2)  - RGTRavg(2))*100;
pelvisZdist = (abs(PSISmidpt(3) + - HJCL(3)))*100; 
ASISYdist = (abs(RASIavg(2) - LASIavg(2)))*100;
PSISYdist = (abs(RPSIavg(2) - LPSIavg(2)))*100;

%Save semgent lengths to data structure
subj1.seg.footX = footXdist;
subj1.seg.footY = footYdist;
subj1.seg.footZ = footZdist;
subj1.seg.legLengthL = legLengthL;
subj1.seg.legLengthR = legLengthR;
subj1.seg.shankLengthL = shankLengthL;
subj1.seg.pelvisX = pelvisXdist;
subj1.seg.pelvisY = pelvisYdist;
subj1.seg.pelvisZ = pelvisZdist;
subj1.seg.ASISYdist = ASISYdist;
subj1.seg.PSISYdist = PSISYdist;

clearvars -except subj1