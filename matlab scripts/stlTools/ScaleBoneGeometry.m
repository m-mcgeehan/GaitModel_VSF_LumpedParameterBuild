addpath(genpath('matlab scripts'),...
        genpath('stlTools'),...
        genpath('Subject Bone Geometry')); %Add dependencies
oldpath = cd;
cd('Subject Bone Geometry');
%% Read in generic geometries
[v, f, n, name] = stlRead('FemurR.stl');
FemurR_sf = subj1.static.segLength.legLengthR/41.58; %Create uniform scaling factor between the generic bone length and the bone length from MoCap
fv.vertices = (v*FemurR_sf)/100; %Convert units
fv.faces = f;
subj1.boneGeom.femurR = fv; 
stlWrite('FemurR_scaled.stl',fv); 
clear v f n name fv

[v, f, n, name] = stlRead('FemurL.stl');
FemurL_sf = subj1.static.segLength.legLengthL/41.58;
fv.vertices = (v*FemurL_sf)/100;
fv.faces = f;
subj1.boneGeom.femurL = fv;
stlWrite('FemurL_scaled.stl',fv); 
clear v f n name fv

[v, f, n, name] = stlRead('ShankL.stl');
ShankL_sf = subj1.static.segLength.shankLengthL/37.25;
fv.vertices = (v*ShankL_sf)/100;
fv.faces = f;
subj1.boneGeom.shankL = fv;
stlWrite('ShankL_scaled.stl',fv); 
clear v f n name fv

[v, f, n, name] = stlRead('ShankR_cut.stl');
ShankR_cut_sf = subj1.static.segLength.shankLengthR/19.0;
fv.vertices = [v(:,1)*ShankL_sf v(:,2)*ShankL_sf v(:,3)*ShankR_cut_sf]/100;
fv.faces = f;
subj1.boneGeom.shankR_cut = fv;
stlWrite('ShankR_cut_scaled.stl',fv); 
clear v f n name fv

[v, f, n, name] = stlRead('Pelvis.stl'); 
% pelvis_sf = (subj1.static.segLength.ASISYdist/25.7+ subj1.static.segLength.PSISYdist/8.63)/2;
pelvis_sfX = 1;
pelvis_sfY = subj1.static.segLength.ASISYdist/26;
pelvis_sfZ = (pelvis_sfX+pelvis_sfY)/2;
fv.vertices = [v(:,1)*pelvis_sfX v(:,2)*pelvis_sfY v(:,3)*pelvis_sfZ]/100;
fv.faces = f;
subj1.boneGeom.pelvis = fv;
stlWrite('Pelvis_scaled.stl',fv); 
clear v f n name fv

[v, f, n, name] = stlRead('HAT.stl');
HAT_sf = (pelvis_sfX + pelvis_sfY + pelvis_sfZ)/3;
fv.vertices = (v*HAT_sf)/100;
fv.faces = f;
subj1.boneGeom.HAT = fv;
stlWrite('HAT_scaled.stl',fv); 
clear v f n name fv


[v, f, n, name] = stlRead('FootL.stl');
FootL_sf = subj1.static.segLength.footLX/21.79;
fv.vertices = (v*FootL_sf)/100;
fv.faces = f;
stlWrite('FootL_scaled.stl',fv); 
clear v f n name fv

cd(oldpath);