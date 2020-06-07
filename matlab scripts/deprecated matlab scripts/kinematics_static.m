%% Markers
%pelvis
% ASISdist = norm(static.markers.RASI - static.markers.LASI);
ASISmidpt = (static.markers.RASI+static.markers.LASI)/2;
PSISmidpt = (static.markers.RPSI+static.markers.LPSI)/2;
pelvisOrigin = ASISmidpt;
HJCR = [static.markers.RASI(1) - (.36*ASISdist), static.markers.RASI(2) - (.14*ASISdist), static.markers.RASI(3) - (.22*ASISdist)]';
HJCL = [static.markers.LASI(1) - (.36*ASISdist), static.markers.LASI(2) + (.14*ASISdist), static.markers.LASI(3) - (.22*ASISdist)]';
pelvis_j = (pelvisOrigin - static.markers.RASI)/norm(pelvisOrigin - static.markers.RASI);
pelvis_k = cross((ASISmidpt-PSISmidpt)/(norm(ASISmidpt-PSISmidpt)),pelvis_j);
pelvis_i = cross(pelvis_j,pelvis_k);
pelvisR = [pelvis_i, pelvis_j, pelvis_k]';

if det(pelvisR) < .9
warning('determinant of pelvisR is below threshold. Determinant = ');
display(det(pelvisR));
else
end

%Leg
KJCR = (static.markers.RLEP + static.markers.RMEP)/2;
legR_k = (HJCR - KJCR)/(norm(HJCR - KJCR));
legR_i = cross(((static.markers.RMEP - static.markers.RLEP)/norm(static.markers.RMEP - static.markers.RLEP)),legR_k);
legR_j = cross(legR_k,legR_i);
legRR = [legR_i, legR_j, legR_k]';

if det(legRR) < .9
warning('determinant of legRR is below threshold. Determinant = ');
display(det(legRR));
else
end

KJCL = (static.markers.LLEP + static.markers.LMEP)/2;
legL_k = (HJCL - KJCL)/(norm(HJCL - KJCL));
legL_i = cross(((static.markers.LLEP - static.markers.LMEP)/norm(static.markers.LLEP - static.markers.LMEP)),legL_k); 
legL_j = cross(legL_k,legL_i);
legLR = [legL_i, legL_j, legL_k]';
legLRdet = det(legLR);

if det(legLR) < .85
warning('determinant of legLR is below threshold. Determinant = ');
display(det(legLR));
else
end

%Shank
AJCR = static.markers.FAJCR; 
shankR_k = (KJCR - AJCR)/(norm(KJCR - AJCR));
shankR_i = cross(((static.markers.RMEP - static.markers.RLEP)/norm(static.markers.RMEP - static.markers.RLEP)),shankR_k);
% shankR_i = cross(((static.markers.RHME - static.markers.RHLA)/norm(static.markers.RHME - static.markers.RHLA)),shankR_k);
shankR_j = cross(shankR_k,shankR_i);
shankRR = [shankR_i, shankR_j, shankR_k];

if det(shankRR) < .9
warning('determinant of shankRR is below threshold. Determinant = ');
display(det(shankRR));
else
end

AJCL = static.markers.FAJCL; %These vector ss's are low
shankL_k = (KJCL - AJCL)/(norm(KJCL - AJCL));
% static.markers.LHME)/norm(static.markers.LHLA - static.markers.LHME)),shankL_k);
shankL_i = cross(((static.markers.LLEP - static.markers.LMEP)/norm(static.markers.LLEP - static.markers.LMEP)),shankL_k);
shankL_j = cross(shankL_k,shankL_i);
shankLR = [shankL_i, shankL_j, shankL_k];

if det(shankLR) < .9
warning('determinant of shankLR is below threshold. Determinant = ');
display(det(shankLR));
else
end

%Foot
footR_i = (static.markers.RTCE - static.markers.RHCE)/norm(static.markers.RTCE - static.markers.RHCE);
footR_k = cross(footR_i,((static.markers.RTME - static.markers.RTLA)/norm(static.markers.RTME - static.markers.RTLA)));
footR_j = cross(footR_k,footR_i);
footRR = [footR_i, footR_j, footR_k];

if det(footRR) < .9
warning('determinant of footRR is below threshold. Determinant = ');
display(det(footRR));
else
end

footL_i = (static.markers.LTCE - static.markers.LHCE)/norm(static.markers.LTCE - static.markers.LHCE);
footL_k = cross(footL_i,((static.markers.LTLA - static.markers.LTME)/norm(static.markers.LTLA - static.markers.LTME)));
footL_j = cross(footL_k,footL_i);
footLR = [footL_i, footL_j, footL_k];

if det(footLR) < .9
warning('determinant of footLR is below threshold...Determinant =  ');
display(det(footLR));
else
end
%% Calculate segment angles wrt GCS
pelvisAng_GCS = SpinCalc('DCMtoEA123',pelvisR,.2,1);
legRAng_GCS = SpinCalc('DCMtoEA123',legRR,.2,1);
legLAng_GCS = SpinCalc('DCMtoEA123',legLR,.2,1);
shankRAng_GCS = SpinCalc('DCMtoEA123',shankRR,.3,1); %precision tolerance violation 
shankLAng_GCS = SpinCalc('DCMtoEA123',shankLR,.2,1);
footRAng_GCS = SpinCalc('DCMtoEA123',footRR,.2,1);
footLAng_GCS = SpinCalc('DCMtoEA123',footLR,.2,1);

legRAng_pelvis = SpinCalc('DCMtoEA123',((inv(pelvisR))*legRR),.2,1);
legLAng_pelvis = SpinCalc('DCMtoEA123',((inv(pelvisR))*legLR),.2,1);
shankRAng_legR = SpinCalc('DCMtoEA123',((inv(legRR))*shankRR),.3,1); %precision tolerance violation 
shankLAng_legL = SpinCalc('DCMtoEA123',((inv(legLR))*shankLR),.2,1);
footRAng_shankR = SpinCalc('DCMtoEA123',((inv(shankRR))*footRR),.35,1); %precision tolerance violation 
footLAng_shankL = SpinCalc('DCMtoEA123',((inv(shankLR))*footLR),.2,1);

% legRAng_pelvis = SpinCalc('DCMtoEA123',(pelvisR\legRR),.2,1);
% legLAng_pelvis = SpinCalc('DCMtoEA123',(pelvisR\legLR),.2,1);
% shankRAng_legR = SpinCalc('DCMtoEA123',(legRR\shankRR),.3,1); %precision tolerance violation 
% shankLAng_legL = SpinCalc('DCMtoEA123',(legLR\shankLR),.2,1);
% footRAng_shankR = SpinCalc('DCMtoEA123',(shankRR\footRR),.35,1); %precision tolerance violation 
% footLAng_shankL = SpinCalc('DCMtoEA123',(shankLR\footLR),.2,1);

static.GCSOrientations.Pelvis_GCS = pelvisAng_GCS;
static.GCSOrientations.LegLAng_GCS = legLAng_GCS;
static.GCSOrientations.LegRAng_GCS = legRAng_GCS;
static.GCSOrientations.ShankLAng_GCS = shankLAng_GCS;
static.GCSOrientations.ShankRAng_GCS = shankRAng_GCS;
static.GCSOrientations.footRAng_GCS = footRAng_GCS;
static.GCSOrientations.footLAng_GCS = footLAng_GCS;

static.km.Rankle = footRAng_shankR;
static.km.Lankle = footLAng_shankL;
static.km.Rhip = legRAng_pelvis;
static.km.Lhip = legLAng_pelvis;
static.km.Rknee = shankRAng_legR;
static.km.Lknee = shankLAng_legL;
static.km.PelvisAngle = pelvisAng_GCS;
static.km.pelvisCOMPos = pelvisOrigin;




