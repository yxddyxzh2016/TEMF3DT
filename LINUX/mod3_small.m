clear all;
close all;

% intel MPI library variables
if isunix
    prog_mpi = "mpirun";
elseif ispc
    mpi_root_direc = getenv('I_MPI_ONEAPI_ROOT');
    if (mpi_root_direc == "")
        disp('No Intel MPI library found on Windows system, stop!');
        return
    end
    prog_mpi = append('"',mpi_root_direc,filesep,'bin',filesep,'mpiexec"');
else
    disp('Platform other than Windows and Linux is not supported, stop!')
    return
end
mpi_thread = 1;
omp_thread = 16;

% receivers
coor_all = zeros(76*76,3);
k = 0;
for i=1:76
    for j=1:76
        k = k+1;
        coor_all(k,1) = -750+20*(i-1);
        coor_all(k,2) = -750+20*(j-1);
        if j==1, coor_all(k,3)=95; end
        if j==2, coor_all(k,3)=96; end
        if j==3, coor_all(k,3)=93; end
        if j==4, coor_all(k,3)=86; end
        if j==5, coor_all(k,3)=79; end
        if j==6, coor_all(k,3)=71; end
        if j==7, coor_all(k,3)=67; end
        if j==8, coor_all(k,3)=63; end
        if j==9, coor_all(k,3)=58; end
        if j==10, coor_all(k,3)=53; end
        if j==11, coor_all(k,3)=49; end
        if j==12, coor_all(k,3)=47; end
        if j==13, coor_all(k,3)=33; end
        if j==14, coor_all(k,3)=31; end
        if j==15, coor_all(k,3)=28; end
        if j==16, coor_all(k,3)=26; end
        if j==17, coor_all(k,3)=30; end
        if j==18, coor_all(k,3)=38; end
        if j==19, coor_all(k,3)=44; end
        if j==20, coor_all(k,3)=48; end
        if j==21, coor_all(k,3)=52; end
        if j==22, coor_all(k,3)=50; end
        if j==23, coor_all(k,3)=46; end
        if j==24, coor_all(k,3)=41; end
        if j==25, coor_all(k,3)=34; end
        if j==26, coor_all(k,3)=28; end
        if j==27, coor_all(k,3)=22; end
        if j==28, coor_all(k,3)=17; end
        if j==29, coor_all(k,3)=12; end
        if j==30, coor_all(k,3)=10; end
        if j==31, coor_all(k,3)=13; end
        if j==32, coor_all(k,3)=23; end
        if j==33, coor_all(k,3)=39; end
        if j==34, coor_all(k,3)=50; end
        if j==35, coor_all(k,3)=57; end
        if j==36, coor_all(k,3)=57; end
        if j==37, coor_all(k,3)=56; end
        if j==38, coor_all(k,3)=58; end
        if j==39, coor_all(k,3)=60; end
        if j==40, coor_all(k,3)=61; end
        if j==41, coor_all(k,3)=62; end
        if j==42, coor_all(k,3)=62; end
        if j==43, coor_all(k,3)=64; end
        if j==44, coor_all(k,3)=68; end
        if j==45, coor_all(k,3)=71; end
        if j==46, coor_all(k,3)=72; end
        if j==47, coor_all(k,3)=73; end
        if j==48, coor_all(k,3)=74; end
        if j==49, coor_all(k,3)=73; end
        if j==50, coor_all(k,3)=72; end
        if j==51, coor_all(k,3)=62; end
        if j==52, coor_all(k,3)=56; end
        if j==53, coor_all(k,3)=56; end
        if j==54, coor_all(k,3)=57; end
        if j==55, coor_all(k,3)=58; end
        if j==56, coor_all(k,3)=60; end
        if j==57, coor_all(k,3)=67; end
        if j==58, coor_all(k,3)=72; end
        if j==59, coor_all(k,3)=83; end
        if j==60, coor_all(k,3)=88; end
        if j==61, coor_all(k,3)=92; end
        if j==62, coor_all(k,3)=97; end
        if j==63, coor_all(k,3)=97; end
        if j==64, coor_all(k,3)=101; end
        if j==65, coor_all(k,3)=110; end
        if j==66, coor_all(k,3)=121; end
        if j==67, coor_all(k,3)=126; end
        if j==68, coor_all(k,3)=133; end
        if j==69, coor_all(k,3)=141; end
        if j==70, coor_all(k,3)=147; end
        if j==71, coor_all(k,3)=151; end
        if j==72, coor_all(k,3)=153; end
        if j==73, coor_all(k,3)=154; end
        if j==74, coor_all(k,3)=153; end
        if j==75, coor_all(k,3)=152; end
        if j==76, coor_all(k,3)=149; end
    end
end

coor = zeros(16*16,3);
k = 0;
for i=1:16
    for j=1:16
        k = k+1;
        coor(k,1) = -750+100*(i-1);
        coor(k,2) = -750+100*(j-1);
        if j==1, coor(k,3)=95; end
        if j==2, coor(k,3)=71; end
        if j==3, coor(k,3)=49; end
        if j==4, coor(k,3)=26; end
        if j==5, coor(k,3)=52; end
        if j==6, coor(k,3)=28; end
        if j==7, coor(k,3)=13; end
        if j==8, coor(k,3)=57; end
        if j==9, coor(k,3)=62; end
        if j==10, coor(k,3)=72; end
        if j==11, coor(k,3)=62; end
        if j==12, coor(k,3)=60; end
        if j==13, coor(k,3)=92; end
        if j==14, coor(k,3)=121; end
        if j==15, coor(k,3)=151; end
        if j==16, coor(k,3)=149; end
    end
end

% loop
field = zeros(6,16*16);
for i=162:16*16
    tic
    disp('----------------------------------');
    disp(string(i)+' of '+string(16*16));
    % cntl file
    fid = fopen('mod3_small.cntl','wt');
    fprintf(fid,'receivers: 1\n');
    fprintf(fid,string(coor(i,1))+' '+string(coor(i,2))+' '+string(coor(i,3))+'\n');
    fprintf(fid,'isotropic: 6\n');
    fprintf(fid,'11 1.00d-8 1.0 1.0\n');
    fprintf(fid,'12 2.27d-4 1.0 1.0\n');
    fprintf(fid,'13 5.00d-2 1.0 1.0\n');
    fprintf(fid,'14 2.50d-1 1.0 1.0\n');
    fprintf(fid,'15 8.51d-4 1.0 1.0\n');
    fprintf(fid,'16 2.60d-3 1.0 1.0\n');
    fprintf(fid,'boundary condition: 0\n');
    fprintf(fid,'source: 1\n');
    fprintf(fid,'1    1   1.0e-6    1.0\n');
    fprintf(fid,'sp_mthd: 2\n');
    fprintf(fid,'sp_division: 50\n');
    fprintf(fid,'time_points: 6\n');
    fprintf(fid,'5d-6\n');
    fprintf(fid,'1d-5\n');
    fprintf(fid,'5d-5\n');
    fprintf(fid,'1d-4\n');
    fprintf(fid,'5d-4\n');
    fprintf(fid,'1d-3\n');
    fprintf(fid,'time_maximum: 2e-3\n');
    fprintf(fid,'sp_double: 60\n');
    fprintf(fid,'sp_dblesize: 2\n');
    fprintf(fid,'output:\n');
    fprintf(fid,'0 0 0 0 0 1 0 0 0\n');
    fprintf(fid,'vtkout: 0');
    fclose(fid);
    % geo file
    fid = fopen('mod3_small.geo','wt');
    fprintf(fid,'SetFactory("OpenCASCADE");\n');
%     fprintf(fid,'Mesh.Algorithm3D = 10;\n');
%     fprintf(fid,'General.NumThreads = 16;\n');
%     fprintf(fid,'Mesh.MaxNumThreads1D = 16;\n');
%     fprintf(fid,'Mesh.MaxNumThreads2D = 16;\n');
%     fprintf(fid,'Mesh.MaxNumThreads3D = 16;\n');
    fprintf(fid,'Mesh.SaveElementTagType = 2;\n');
    fprintf(fid,'Mesh.CharacteristicLengthMin = 0.1;\n');
    fprintf(fid,'Mesh.CharacteristicLengthMax = 4000;\n');
    fprintf(fid,'Mesh.CharacteristicLengthFromCurvature = 1;\n');
    fprintf(fid,'Mesh.MinimumElementsPerTwoPi = 6;\n');
    fprintf(fid,'// points\n');
    fprintf(fid,'Point(1) = {-750, -750, 0,300};\n');
    fprintf(fid,'Point(2) = {-750, 750, 0,300};\n');
    fprintf(fid,'Point(3) = {-750, -750, 100,300};\n');
    fprintf(fid,'Point(4) = {-750, -674, 86};\n');
    fprintf(fid,'Point(5) = {-750, -640, 74};\n');
    fprintf(fid,'Point(6) = {-750, -615, 69};\n');
    fprintf(fid,'Point(7) = {-750, -575, 59};\n');
    fprintf(fid,'Point(8) = {-750, -530, 53};\n');
    fprintf(fid,'Point(9) = {-750, -521, 41};\n');
    fprintf(fid,'Point(24) = {-750, -416, 41};\n');
    fprintf(fid,'Point(25) = {-750, -436, 33};\n');
    fprintf(fid,'Point(26) = {-750, -450, 31};\n');
    fprintf(fid,'Point(27) = {-750, -494, 36};\n');
    fprintf(fid,'Point(28) = {-750, -369, 54};\n');
    fprintf(fid,'Point(29) = {-750, -349, 57};\n');
    fprintf(fid,'Point(30) = {-750, -337, 56};\n');
    fprintf(fid,'Point(31) = {-750, -260, 36};\n');
    fprintf(fid,'Point(32) = {-750, -195, 18};\n');
    fprintf(fid,'Point(33) = {-750, -172, 15};\n');
    fprintf(fid,'Point(34) = {-750, -152, 18};\n');
    fprintf(fid,'Point(35) = {-750, -137, 23};\n');
    fprintf(fid,'Point(36) = {-750, -117, 40};\n');
    fprintf(fid,'Point(37) = {-750, -100, 49};\n');
    fprintf(fid,'Point(38) = {-750, -86, 57};\n');
    fprintf(fid,'Point(39) = {-750, -91, 60};\n');
    fprintf(fid,'Point(43) = {-750, -289, 190};\n');
    fprintf(fid,'Point(44) = {-750, -310, 213};\n');
    fprintf(fid,'Point(69) = {-750, -321, 222};\n');
    fprintf(fid,'Point(70) = {-750, -334, 236};\n');
    fprintf(fid,'Point(83) = {-750, -334, 253};\n');
    fprintf(fid,'Point(84) = {-750, -330, 259};\n');
    fprintf(fid,'Point(85) = {-750, -347, 287};\n');
    fprintf(fid,'Point(86) = {-750, -373, 320};\n');
    fprintf(fid,'Point(87) = {-750, -400, 356};\n');
    fprintf(fid,'Point(88) = {-750, -416, 378};\n');
    fprintf(fid,'Point(89) = {-750, -430, 400,200};\n');
    fprintf(fid,'Point(90) = {-750, -391, 400};\n');
    fprintf(fid,'Point(91) = {-750, -380, 378};\n');
    fprintf(fid,'Point(92) = {-750, -365, 358};\n');
    fprintf(fid,'Point(93) = {-750, -365, 400};\n');
    fprintf(fid,'Point(94) = {-750, -359, 345};\n');
    fprintf(fid,'Point(95) = {-750, -348, 333};\n');
    fprintf(fid,'Point(96) = {-750, -337, 308};\n');
    fprintf(fid,'Point(97) = {-750, -329, 293};\n');
    fprintf(fid,'Point(98) = {-750, -264, 226};\n');
    fprintf(fid,'Point(99) = {-750, -265, 194};\n');
    fprintf(fid,'Point(100) = {-750, -277, 205};\n');
    fprintf(fid,'Point(102) = {-750, -300, 241};\n');
    fprintf(fid,'Point(103) = {-750, -310, 255};\n');
    fprintf(fid,'Point(104) = {-750, -265, 180};\n');
    fprintf(fid,'Point(105) = {-750, -225, 131};\n');
    fprintf(fid,'Point(106) = {-750, -177, 98};\n');
    fprintf(fid,'Point(107) = {-750, -126, 76};\n');
    fprintf(fid,'Point(108) = {-750, -64, 62};\n');
    fprintf(fid,'Point(109) = {-750, -34, 61};\n');
    fprintf(fid,'Point(110) = {-750, 2, 64};\n');
    fprintf(fid,'Point(111) = {-750, 36, 66};\n');
    fprintf(fid,'Point(112) = {-750, 69, 67};\n');
    fprintf(fid,'Point(113) = {-750, 101, 71};\n');
    fprintf(fid,'Point(114) = {-750, 111, 74};\n');
    fprintf(fid,'Point(118) = {-750, 20, 101};\n');
    fprintf(fid,'Point(119) = {-750, -1, 127};\n');
    fprintf(fid,'Point(120) = {-750, -26, 152};\n');
    fprintf(fid,'Point(121) = {-750, -45, 178};\n');
    fprintf(fid,'Point(122) = {-750, -68, 205};\n');
    fprintf(fid,'Point(123) = {-750, -108, 257};\n');
    fprintf(fid,'Point(126) = {-750, -153, 300};\n');
    fprintf(fid,'Point(127) = {-750, -122, 288};\n');
    fprintf(fid,'Point(128) = {-750, -93, 258};\n');
    fprintf(fid,'Point(129) = {-750, -54, 204};\n');
    fprintf(fid,'Point(130) = {-750, -20, 156};\n');
    fprintf(fid,'Point(131) = {-750, 7, 139};\n');
    fprintf(fid,'Point(132) = {-750, 30, 132};\n');
    fprintf(fid,'Point(133) = {-750, 44, 138};\n');
    fprintf(fid,'Point(134) = {-750, 49, 147};\n');
    fprintf(fid,'Point(137) = {-750, 30, 152};\n');
    fprintf(fid,'Point(139) = {-750, 28, 234};\n');
    fprintf(fid,'Point(140) = {-750, -6, 272};\n');
    fprintf(fid,'Point(141) = {-750, -36, 302};\n');
    fprintf(fid,'Point(142) = {-750, -50, 313};\n');
    fprintf(fid,'Point(143) = {-750, -51, 326};\n');
    fprintf(fid,'Point(144) = {-750, -96, 340};\n');
    fprintf(fid,'Point(145) = {-750, -130, 354};\n');
    fprintf(fid,'Point(146) = {-750, -173, 390};\n');
    fprintf(fid,'Point(152) = {-750, -177, 400};\n');
    fprintf(fid,'Point(153) = {-750, 88, 91};\n');
    fprintf(fid,'Point(154) = {-750, 71, 111};\n');
    fprintf(fid,'Point(155) = {-750, 64, 135};\n');
    fprintf(fid,'Point(156) = {-750, 62, 165};\n');
    fprintf(fid,'Point(157) = {-750, 57, 178};\n');
    fprintf(fid,'Point(158) = {-750, 49, 176};\n');
    fprintf(fid,'Point(159) = {-750, 43, 167};\n');
    fprintf(fid,'Point(160) = {-750, 40, 152};\n');
    fprintf(fid,'Point(163) = {-750, 44, 202};\n');
    fprintf(fid,'Point(164) = {-750, 46, 197};\n');
    fprintf(fid,'Point(165) = {-750, 48, 192};\n');
    fprintf(fid,'Point(167) = {-750, 107, 80};\n');
    fprintf(fid,'Point(168) = {-750, 118, 95};\n');
    fprintf(fid,'Point(169) = {-750, 132, 123};\n');
    fprintf(fid,'Point(170) = {-750, 132, 198};\n');
    fprintf(fid,'Point(173) = {-750, 161, 179};\n');
    fprintf(fid,'Point(174) = {-750, 156, 185};\n');
    fprintf(fid,'Point(175) = {-750, 147, 187};\n');
    fprintf(fid,'Point(176) = {-750, 140, 193};\n');
    fprintf(fid,'Point(177) = {-750, -39, 325};\n');
    fprintf(fid,'Point(178) = {-750, -36, 319};\n');
    fprintf(fid,'Point(179) = {-750, -7, 312};\n');
    fprintf(fid,'Point(180) = {-750, 24, 293};\n');
    fprintf(fid,'Point(182) = {-750, 71, 263};\n');
    fprintf(fid,'Point(183) = {-750, 52, 273};\n');
    fprintf(fid,'Point(185) = {-750, 240, 77};\n');
    fprintf(fid,'Point(186) = {-750, 160, 158};\n');
    fprintf(fid,'Point(189) = {-750, 156, 81};\n');
    fprintf(fid,'Point(202) = {-750, -185, 254};\n');
    fprintf(fid,'Point(203) = {-750, -180, 255};\n');
    fprintf(fid,'Point(204) = {-750, -177, 257};\n');
    fprintf(fid,'Point(205) = {-750, -190, 284};\n');
    fprintf(fid,'Point(206) = {-750, -210, 313};\n');
    fprintf(fid,'Point(207) = {-750, -216, 313};\n');
    fprintf(fid,'Point(208) = {-750, -218, 307};\n');
    fprintf(fid,'Point(209) = {-750, -200, 283};\n');
    fprintf(fid,'Point(210) = {-750, -299, 334};\n');
    fprintf(fid,'Point(211) = {-750, -301, 313};\n');
    fprintf(fid,'Point(212) = {-750, -303, 310};\n');
    fprintf(fid,'Point(213) = {-750, -291, 290};\n');
    fprintf(fid,'Point(214) = {-750, -262, 264};\n');
    fprintf(fid,'Point(215) = {-750, -227, 215};\n');
    fprintf(fid,'Point(216) = {-750, -202, 192};\n');
    fprintf(fid,'Point(217) = {-750, -178, 164};\n');
    fprintf(fid,'Point(218) = {-750, -160, 148};\n');
    fprintf(fid,'Point(219) = {-750, -144, 142};\n');
    fprintf(fid,'Point(220) = {-750, -129, 141};\n');
    fprintf(fid,'Point(221) = {-750, -128, 144};\n');
    fprintf(fid,'Point(222) = {-750, -133, 151};\n');
    fprintf(fid,'Point(223) = {-750, -134, 161};\n');
    fprintf(fid,'Point(224) = {-750, -149, 178};\n');
    fprintf(fid,'Point(225) = {-750, -159, 193};\n');
    fprintf(fid,'Point(226) = {-750, -168, 209};\n');
    fprintf(fid,'Point(227) = {-750, -179, 221};\n');
    fprintf(fid,'Point(228) = {-750, -192, 229};\n');
    fprintf(fid,'Point(229) = {-750, -193, 226};\n');
    fprintf(fid,'Point(230) = {-750, -191, 223};\n');
    fprintf(fid,'Point(231) = {-750, -181, 216};\n');
    fprintf(fid,'Point(232) = {-750, -167, 194};\n');
    fprintf(fid,'Point(233) = {-750, -157, 179};\n');
    fprintf(fid,'Point(234) = {-750, -145, 167};\n');
    fprintf(fid,'Point(235) = {-750, -140, 160};\n');
    fprintf(fid,'Point(236) = {-750, -142, 151};\n');
    fprintf(fid,'Point(237) = {-750, -161, 159};\n');
    fprintf(fid,'Point(238) = {-750, -183, 185};\n');
    fprintf(fid,'Point(240) = {-750, -202, 211};\n');
    fprintf(fid,'Point(241) = {-750, -219, 238};\n');
    fprintf(fid,'Point(242) = {-750, -248, 275};\n');
    fprintf(fid,'Point(243) = {-750, -274, 307};\n');
    fprintf(fid,'Point(244) = {-750, -135, 122};\n');
    fprintf(fid,'Point(245) = {-750, -187, 121};\n');
    fprintf(fid,'Point(246) = {-750, -199, 126};\n');
    fprintf(fid,'Point(247) = {-750, -217, 139};\n');
    fprintf(fid,'Point(248) = {-750, -213, 145};\n');
    fprintf(fid,'Point(249) = {-750, -197, 136};\n');
    fprintf(fid,'Point(250) = {-750, -188, 135};\n');
    fprintf(fid,'Point(251) = {-750, -192, 139};\n');
    fprintf(fid,'Point(252) = {-750, -201, 147};\n');
    fprintf(fid,'Point(253) = {-750, -204, 153};\n');
    fprintf(fid,'Point(254) = {-750, -202, 156};\n');
    fprintf(fid,'Point(255) = {-750, -173, 135};\n');
    fprintf(fid,'Point(256) = {-750, -165, 131};\n');
    fprintf(fid,'Point(257) = {-750, -162, 133};\n');
    fprintf(fid,'Point(258) = {-750, -183, 143};\n');
    fprintf(fid,'Point(259) = {-750, -194, 151};\n');
    fprintf(fid,'Point(260) = {-750, -199, 153};\n');
    fprintf(fid,'Point(261) = {-750, -175, 143};\n');
    fprintf(fid,'Point(262) = {-750, -188, 156};\n');
    fprintf(fid,'Point(263) = {-750, -198, 171};\n');
    fprintf(fid,'Point(264) = {-750, -206, 181};\n');
    fprintf(fid,'Point(265) = {-750, -212, 184};\n');
    fprintf(fid,'Point(266) = {-750, -214, 188};\n');
    fprintf(fid,'Point(267) = {-750, -209, 192};\n');
    fprintf(fid,'Point(268) = {-750, -130, 125};\n');
    fprintf(fid,'Point(269) = {-750, -135, 132};\n');
    fprintf(fid,'Point(270) = {-750, -143, 131};\n');
    fprintf(fid,'Point(271) = {-750, -161, 144};\n');
    fprintf(fid,'Point(272) = {-750, -175, 157};\n');
    fprintf(fid,'Point(273) = {-750, -185, 170};\n');
    fprintf(fid,'Point(274) = {-750, -193, 181};\n');
    fprintf(fid,'Point(275) = {-750, -201, 188};\n');
    fprintf(fid,'Point(276) = {-750, -62, 76};\n');
    fprintf(fid,'Point(277) = {-750, -103, 87};\n');
    fprintf(fid,'Point(278) = {-750, -154, 94};\n');
    fprintf(fid,'Point(279) = {-750, -158, 97};\n');
    fprintf(fid,'Point(280) = {-750, -154, 102};\n');
    fprintf(fid,'Point(281) = {-750, -58, 79};\n');
    fprintf(fid,'Point(282) = {-750, -62, 84};\n');
    fprintf(fid,'Point(283) = {-750, -89, 93};\n');
    fprintf(fid,'Point(284) = {-750, -117, 100};\n');
    fprintf(fid,'Point(286) = {-750, -10, 94};\n');
    fprintf(fid,'Point(287) = {-750, -11, 101};\n');
    fprintf(fid,'Point(288) = {-750, -7, 91};\n');
    fprintf(fid,'Point(289) = {-750, -2, 92};\n');
    fprintf(fid,'Point(290) = {-750, 3, 97};\n');
    fprintf(fid,'Point(291) = {-750, 1, 102};\n');
    fprintf(fid,'Point(292) = {-750, 0, 106};\n');
    fprintf(fid,'Point(293) = {-750, -8, 113};\n');
    fprintf(fid,'Point(294) = {-750, -19, 125};\n');
    fprintf(fid,'Point(295) = {-750, -26, 133};\n');
    fprintf(fid,'Point(296) = {-750, -32, 140};\n');
    fprintf(fid,'Point(297) = {-750, -38, 140};\n');
    fprintf(fid,'Point(298) = {-750, -41, 136};\n');
    fprintf(fid,'Point(299) = {-750, -36, 129};\n');
    fprintf(fid,'Point(300) = {-750, -25, 114};\n');
    fprintf(fid,'Point(301) = {-750, -125, 400,200};\n');
    fprintf(fid,'Point(302) = {-750, -70, 371};\n');
    fprintf(fid,'Point(303) = {-750, 22, 327};\n');
    fprintf(fid,'Point(304) = {-750, 97, 269};\n');
    fprintf(fid,'Point(305) = {-750, 145, 231};\n');
    fprintf(fid,'Point(306) = {-750, 175, 209};\n');
    fprintf(fid,'Point(308) = {-750, 232, 309};\n');
    fprintf(fid,'Point(314) = {-750, 237, 304};\n');
    fprintf(fid,'Point(315) = {-750, 247, 294};\n');
    fprintf(fid,'Point(316) = {-750, 265, 285};\n');
    fprintf(fid,'Point(317) = {-750, 274, 273};\n');
    fprintf(fid,'Point(318) = {-750, 285, 269};\n');
    fprintf(fid,'Point(319) = {-750, 349, 271};\n');
    fprintf(fid,'Point(320) = {-750, 410, 262};\n');
    fprintf(fid,'Point(321) = {-750, 454, 224};\n');
    fprintf(fid,'Point(322) = {-750, 486, 180};\n');
    fprintf(fid,'Point(323) = {-750, 517, 152};\n');
    fprintf(fid,'Point(324) = {-750, 551, 142};\n');
    fprintf(fid,'Point(336) = {-750, 685, 326};\n');
    fprintf(fid,'Point(337) = {-750, 717, 328};\n');
    fprintf(fid,'Point(338) = {-750, 750, 339};\n');
    fprintf(fid,'Point(339) = {-750, 750, 294};\n');
    fprintf(fid,'Point(340) = {-750, 721, 286};\n');
    fprintf(fid,'Point(341) = {-750, 736, 289};\n');
    fprintf(fid,'Point(345) = {-750, 750, 154,300};\n');
    fprintf(fid,'Point(346) = {-750, 700, 159};\n');
    fprintf(fid,'Point(347) = {-750, 622, 150};\n');
    fprintf(fid,'Point(348) = {-750, 563, 129};\n');
    fprintf(fid,'Point(349) = {-750, 538, 120};\n');
    fprintf(fid,'Point(350) = {-750, 549, 130};\n');
    fprintf(fid,'Point(352) = {-750, 543, 124};\n');
    fprintf(fid,'Point(353) = {-750, 177, 183};\n');
    fprintf(fid,'Point(354) = {-750, 176, 189};\n');
    fprintf(fid,'Point(355) = {-750, 213, 266};\n');
    fprintf(fid,'Point(356) = {-750, 227, 263};\n');
    fprintf(fid,'Point(357) = {-750, 256, 250};\n');
    fprintf(fid,'Point(358) = {-750, 275, 241};\n');
    fprintf(fid,'Point(360) = {-750, 298, 228};\n');
    fprintf(fid,'Point(362) = {-750, 322, 203};\n');
    fprintf(fid,'Point(364) = {-750, 341, 195};\n');
    fprintf(fid,'Point(365) = {-750, 363, 183};\n');
    fprintf(fid,'Point(366) = {-750, 384, 176};\n');
    fprintf(fid,'Point(367) = {-750, 449, 175};\n');
    fprintf(fid,'Point(368) = {-750, 470, 143};\n');
    fprintf(fid,'Point(369) = {-750, 474, 112};\n');
    fprintf(fid,'Point(370) = {-750, 468, 102};\n');
    fprintf(fid,'Point(371) = {-750, 449, 97};\n');
    fprintf(fid,'Point(372) = {-750, 412, 89};\n');
    fprintf(fid,'Point(373) = {-750, 375, 87};\n');
    fprintf(fid,'Point(374) = {-750, 361, 98};\n');
    fprintf(fid,'Point(375) = {-750, 322, 103};\n');
    fprintf(fid,'Point(376) = {-750, 270, 118};\n');
    fprintf(fid,'Point(377) = {-750, 222, 151};\n');
    fprintf(fid,'Point(378) = {-750, 195, 175};\n');
    fprintf(fid,'Point(379) = {-750, -750, 400,300};\n');
    fprintf(fid,'Point(380) = {-750, 750, 400,300};\n');
    fprintf(fid,'Point(387) = {-750, -358, 365};\n');
    fprintf(fid,'Point(388) = {-750, -338, 316};\n');
    fprintf(fid,'Point(389) = {-750, -262, 184};\n');
    fprintf(fid,'Point(390) = {-750, -197, 200};\n');
    fprintf(fid,'Point(391) = {-750, -213, 214};\n');
    fprintf(fid,'Point(392) = {-750, -295, 314};\n');
    fprintf(fid,'Point(393) = {-750, -241, 257};\n');
    fprintf(fid,'Point(394) = {-750, -277, 291};\n');
    fprintf(fid,'Point(395) = {-750, -287, 287};\n');
    fprintf(fid,'Point(396) = {-750, -327, 344};\n');
    fprintf(fid,'Point(397) = {-750, -237, 176};\n');
    fprintf(fid,'Point(398) = {-750, -261, 176};\n');
    fprintf(fid,'Point(399) = {-750, -260, 220};\n');
    fprintf(fid,'Point(400) = {-750, -298, 258};\n');
    fprintf(fid,'Point(402) = {-750, 25, 80};\n');
    fprintf(fid,'Point(403) = {-750, 11, 117};\n');
    fprintf(fid,'Point(404) = {-750, 25, 90};\n');
    fprintf(fid,'Point(405) = {-750, -215, 400};\n');
    fprintf(fid,'Point(406) = {-750, -196, 352};\n');
    fprintf(fid,'Point(407) = {-750, -175, 313};\n');
    fprintf(fid,'Point(408) = {-750, -166, 314};\n');
    fprintf(fid,'Point(409) = {-750, -115, 294};\n');
    fprintf(fid,'Point(410) = {-750, -76, 254};\n');
    fprintf(fid,'Point(411) = {-750, -44, 201};\n');
    fprintf(fid,'Point(412) = {-750, 4, 150};\n');
    fprintf(fid,'Point(413) = {-750, 37, 155};\n');
    fprintf(fid,'Point(414) = {-750, 126, 151};\n');
    fprintf(fid,'Point(415) = {-750, 17, 180};\n');
    fprintf(fid,'Point(416) = {-750, -7, 189};\n');
    fprintf(fid,'Point(417) = {-750, -32, 226};\n');
    fprintf(fid,'Point(418) = {-750, -51, 252};\n');
    fprintf(fid,'Point(419) = {-750, -49, 261};\n');
    fprintf(fid,'Point(420) = {-750, -30, 251};\n');
    fprintf(fid,'Point(421) = {-750, -6, 221};\n');
    fprintf(fid,'Point(422) = {-750, 10, 198};\n');
    fprintf(fid,'Point(423) = {-750, 19, 189};\n');
    fprintf(fid,'Point(424) = {-750, -12, 305};\n');
    fprintf(fid,'Point(425) = {-750, 89, 210};\n');
    fprintf(fid,'Point(426) = {-750, 52, 256};\n');
    fprintf(fid,'Point(427) = {-750, -14, 324};\n');
    fprintf(fid,'Point(428) = {-750, -23, 317};\n');
    fprintf(fid,'Point(429) = {-750, 12, 303};\n');
    fprintf(fid,'Point(431) = {-750, 28, 293};\n');
    fprintf(fid,'Point(432) = {-750, 30, 301};\n');
    fprintf(fid,'Point(433) = {-750, 15, 312};\n');
    fprintf(fid,'Point(434) = {-750, -32, 329};\n');
    fprintf(fid,'Point(435) = {-750, 157, 115};\n');
    fprintf(fid,'Point(437) = {-750, 151, 77};\n');
    fprintf(fid,'Point(438) = {-750, 211, 78};\n');
    fprintf(fid,'Point(439) = {-750, 256, 64};\n');
    fprintf(fid,'Point(440) = {-750, 310, 62};\n');
    fprintf(fid,'Point(441) = {-750, 371, 73};\n');
    fprintf(fid,'Point(442) = {-750, 377, 77};\n');
    fprintf(fid,'Point(443) = {-750, 397, 77};\n');
    fprintf(fid,'Point(444) = {-750, 401, 90};\n');
    fprintf(fid,'Point(446) = {-750, 405, 87};\n');
    fprintf(fid,'Point(447) = {-750, 390, 87};\n');
    fprintf(fid,'Point(448) = {-750, 472, 102};\n');
    fprintf(fid,'Point(449) = {-750, 328, 198};\n');
    fprintf(fid,'Point(450) = {-750, 372, 191};\n');
    fprintf(fid,'Point(451) = {-750, 329, 207};\n');
    fprintf(fid,'Point(452) = {-750, 350, 202};\n');
    fprintf(fid,'Point(453) = {-750, 390, 186};\n');
    fprintf(fid,'Point(454) = {-750, 391, 180};\n');
    fprintf(fid,'Point(455) = {-750, 357, 228};\n');
    fprintf(fid,'Point(456) = {-750, 299, 235};\n');
    fprintf(fid,'Point(457) = {-750, 388, 217};\n');
    fprintf(fid,'Point(458) = {-750, 379, 223};\n');
    fprintf(fid,'Point(459) = {-750, 393, 221};\n');
    fprintf(fid,'Point(460) = {-750, 378, 229};\n');
    fprintf(fid,'Point(461) = {-750, 349, 235};\n');
    fprintf(fid,'Point(462) = {-750, 157, 92};\n');
    fprintf(fid,'Point(463) = {-750, 152, 117};\n');
    fprintf(fid,'Point(464) = {-750, 145, 109};\n');
    fprintf(fid,'Point(465) = {-750, 146, 92};\n');
    fprintf(fid,'Point(466) = {-750, 151, 85};\n');
    fprintf(fid,'Point(467) = {-750, 162, 89};\n');
    fprintf(fid,'Point(468) = {-750, -144, 285};\n');
    fprintf(fid,'Point(71) = {-750, -362, 237};\n');
    fprintf(fid,'Point(72) = {-750, -376, 244};\n');
    fprintf(fid,'Point(73) = {-750, -382, 249};\n');
    fprintf(fid,'Point(74) = {-750, -394, 254};\n');
    fprintf(fid,'Point(381) = {-750, -405, 257};\n');
    fprintf(fid,'Point(382) = {-750, -409, 264};\n');
    fprintf(fid,'Point(383) = {-750, -396, 272};\n');
    fprintf(fid,'Point(384) = {-750, -370, 264};\n');
    fprintf(fid,'Point(385) = {-750, -347, 243};\n');
    fprintf(fid,'Point(386) = {-750, -350, 237};\n');
    fprintf(fid,'Point(192) = {-750, -287, 260};\n');
    fprintf(fid,'Point(193) = {-750, -282, 260};\n');
    fprintf(fid,'Point(194) = {-750, -280, 264};\n');
    fprintf(fid,'Point(195) = {-750, -288, 273};\n');
    fprintf(fid,'Point(196) = {-750, -295, 281};\n');
    fprintf(fid,'Point(198) = {-750, -299, 289};\n');
    fprintf(fid,'Point(199) = {-750, -304, 290};\n');
    fprintf(fid,'Point(200) = {-750, -306, 282};\n');
    fprintf(fid,'Point(201) = {-750, -297, 269};\n');
    fprintf(fid,'// lines\n');
    fprintf(fid,'Line(1) = {1, 2};\n');
    fprintf(fid,'Line(2) = {1, 3};\n');
    fprintf(fid,'Spline(3) = {3, 4, 5, 6, 7, 8};\n');
    fprintf(fid,'Spline(4) = {8, 9, 27, 26, 25, 24, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38};\n');
    fprintf(fid,'Line(6) = {3, 379};\n');
    fprintf(fid,'Line(7) = {379, 89};\n');
    fprintf(fid,'Spline(8) = {89, 88, 87, 86, 85, 84};\n');
    fprintf(fid,'Spline(9) = {84, 83, 70, 69, 44, 43, 104};\n');
    fprintf(fid,'Spline(10) = {382, 383, 384, 385};\n');
    fprintf(fid,'Line(11) = {385, 386};\n');
    fprintf(fid,'Line(12) = {386, 71};\n');
    fprintf(fid,'Line(13) = {381, 382};\n');
    fprintf(fid,'Spline(14) = {71, 72, 73, 74, 381};\n');
    fprintf(fid,'Line(15) = {89, 90};\n');
    fprintf(fid,'Spline(16) = {90, 91, 92};\n');
    fprintf(fid,'Spline(17) = {92, 94, 95, 388};\n');
    fprintf(fid,'Spline(18) = {389, 99, 100, 102, 103, 97, 96, 388};\n');
    fprintf(fid,'Spline(19) = {104, 389};\n');
    fprintf(fid,'Spline(20) = {199, 198, 196, 195, 194};\n');
    fprintf(fid,'Spline(21) = {199, 200, 201, 192, 193, 194};\n');
    fprintf(fid,'Spline(22) = {247, 246, 245};\n');
    fprintf(fid,'Spline(23) = {245, 244};\n');
    fprintf(fid,'Spline(24) = {244, 268};\n');
    fprintf(fid,'Spline(25) = {268, 269, 270};\n');
    fprintf(fid,'Spline(27) = {267, 266, 265, 264};\n');
    fprintf(fid,'Spline(28) = {264, 263, 262, 261, 257};\n');
    fprintf(fid,'Spline(29) = {247, 248, 249, 250};\n');
    fprintf(fid,'Spline(30) = {250, 251, 252, 253, 254};\n');
    fprintf(fid,'Spline(31) = {254, 260, 259, 258, 255, 256};\n');
    fprintf(fid,'Line(32) = {256, 257};\n');
    fprintf(fid,'Spline(33) = {390, 240, 241, 242, 243, 210, 93};\n');
    fprintf(fid,'Spline(34) = {93, 90};\n');
    fprintf(fid,'Spline(38) = {211, 396, 387, 92};\n');
    fprintf(fid,'Line(39) = {400, 388};\n');
    fprintf(fid,'Spline(40) = {400, 98, 399, 397};\n');
    fprintf(fid,'Line(41) = {397, 398};\n');
    fprintf(fid,'Line(42) = {398, 104};\n');
    fprintf(fid,'Spline(44) = {278, 279, 280};\n');
    fprintf(fid,'Spline(45) = {280, 284, 283, 282};\n');
    fprintf(fid,'Spline(46) = {282, 281, 276};\n');
    fprintf(fid,'Spline(47) = {276, 277, 278};\n');
    fprintf(fid,'Spline(49) = {293, 294, 295, 296, 297};\n');
    fprintf(fid,'Spline(50) = {297, 298, 299, 300, 287};\n');
    fprintf(fid,'Spline(51) = {209, 202};\n');
    fprintf(fid,'Spline(52) = {209, 208};\n');
    fprintf(fid,'Spline(53) = {208, 207, 206};\n');
    fprintf(fid,'Spline(54) = {206, 205, 204};\n');
    fprintf(fid,'Spline(55) = {202, 203, 204};\n');
    fprintf(fid,'Line(56) = {236, 222};\n');
    fprintf(fid,'Spline(57) = {236, 235, 234, 233, 232, 231, 230};\n');
    fprintf(fid,'Spline(58) = {222, 223, 224, 225, 226, 227, 228};\n');
    fprintf(fid,'Spline(59) = {230, 229, 228};\n');
    fprintf(fid,'Spline(60) = {390, 238, 237, 236};\n');
    fprintf(fid,'Spline(61) = {222, 221, 220};\n');
    fprintf(fid,'Line(62) = {267, 275};\n');
    fprintf(fid,'Line(63) = {275, 216};\n');
    fprintf(fid,'Spline(64) = {275, 274, 273, 217};\n');
    fprintf(fid,'Spline(65) = {217, 272, 271, 270};\n');
    fprintf(fid,'Spline(66) = {217, 218};\n');
    fprintf(fid,'Spline(67) = {218, 219, 220};\n');
    fprintf(fid,'Spline(69) = {38, 108, 109, 110, 111};\n');
    fprintf(fid,'Spline(70) = {111, 402};\n');
    fprintf(fid,'Line(71) = {402, 404};\n');
    fprintf(fid,'Spline(73) = {407, 406, 405};\n');
    fprintf(fid,'Line(74) = {93, 405};\n');
    fprintf(fid,'Spline(75) = {407, 408, 409, 410, 411, 412, 413};\n');
    fprintf(fid,'Spline(78) = {133, 134};\n');
    fprintf(fid,'Spline(79) = {134, 160, 413};\n');
    fprintf(fid,'Spline(82) = {126, 127, 128, 129, 130, 131, 132, 133};\n');
    fprintf(fid,'Spline(89) = {415, 423, 422, 421, 420, 419};\n');
    fprintf(fid,'Spline(90) = {415, 416, 417, 418, 419};\n');
    fprintf(fid,'Spline(92) = {167, 153, 154, 155, 156, 157};\n');
    fprintf(fid,'Line(93) = {405, 152};\n');
    fprintf(fid,'Spline(94) = {152, 146, 145, 144, 143, 177};\n');
    fprintf(fid,'Spline(95) = {158, 165, 164, 163, 139, 140, 141, 142, 178};\n');
    fprintf(fid,'Spline(96) = {178, 424, 426, 425, 414, 169, 168, 167};\n');
    fprintf(fid,'Spline(97) = {177, 428, 179, 429, 180, 431};\n');
    fprintf(fid,'Spline(98) = {177, 434, 427, 433, 432, 431};\n');
    fprintf(fid,'Spline(99) = {431, 183, 182, 170, 176, 175, 174, 173};\n');
    fprintf(fid,'Spline(100) = {173, 306};\n');
    fprintf(fid,'Spline(101) = {306, 305, 304, 303, 302, 301};\n');
    fprintf(fid,'Line(102) = {301, 152};\n');
    fprintf(fid,'Spline(103) = {111, 112, 113, 114};\n');
    fprintf(fid,'Spline(104) = {114, 167};\n');
    fprintf(fid,'Spline(105) = {114, 437, 438, 185};\n');
    fprintf(fid,'Line(106) = {186, 185};\n');
    fprintf(fid,'Spline(107) = {185, 439,440,441,442,443,446};\n');
    fprintf(fid,'Line(108) = {186, 354};\n');
    fprintf(fid,'Line(109) = {354, 353};\n');
    fprintf(fid,'Spline(110) = {353, 378, 377, 376, 375, 374, 373};\n');
    fprintf(fid,'Line(115) = {373, 447};\n');
    fprintf(fid,'Spline(116) = {447, 444, 446};\n');
    fprintf(fid,'Spline(118) = {446, 372, 371, 370, 448};\n');
    fprintf(fid,'Line(119) = {301, 380};\n');
    fprintf(fid,'Line(120) = {380, 338};\n');
    fprintf(fid,'Spline(121) = {339, 341, 340};\n');
    fprintf(fid,'Spline(122) = {338, 337, 336};\n');
    fprintf(fid,'Line(123) = {340, 336};\n');
    fprintf(fid,'Line(124) = {339, 338};\n');
    fprintf(fid,'Spline(125) = {352, 348, 347, 346, 345};\n');
    fprintf(fid,'Line(126) = {345, 339};\n');
    fprintf(fid,'Line(127) = {345, 2};\n');
    fprintf(fid,'Spline(128) = {448, 349,352};\n');
    fprintf(fid,'Spline(130) = {366, 454, 453};\n');
    fprintf(fid,'Spline(131) = {449, 362, 451};\n');
    fprintf(fid,'Spline(132) = {451, 452, 450, 453};\n');
    fprintf(fid,'Spline(133) = {366, 365, 364, 449};\n');
    fprintf(fid,'Line(134) = {360, 455};\n');
    fprintf(fid,'Line(135) = {360, 456};\n');
    fprintf(fid,'Spline(136) = {455, 458, 457};\n');
    fprintf(fid,'Line(137) = {456, 461};\n');
    fprintf(fid,'Spline(138) = {461, 460, 459};\n');
    fprintf(fid,'Line(139) = {457, 459};\n');
    fprintf(fid,'Spline(140) = {456, 358, 357, 356, 355};\n');
    fprintf(fid,'Line(141) = {355, 308};\n');
    fprintf(fid,'Spline(142) = {308, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324};\n');
    fprintf(fid,'Spline(143) = {459, 367, 368, 369, 448};\n');
    fprintf(fid,'Line(144) = {324, 350};\n');
    fprintf(fid,'Line(145) = {350, 352};\n');
    fprintf(fid,'Spline(146) = {435, 463, 464, 465, 466};\n');
    fprintf(fid,'Line(148) = {462, 435};\n');
    fprintf(fid,'Spline(149) = {216, 215, 214, 395, 213, 212, 211};\n');
    fprintf(fid,'Spline(150) = {211, 392, 394, 393, 391, 390};\n');
    fprintf(fid,'Line(151) = {178, 177};\n');
    fprintf(fid,'Spline(152) = {413, 159, 158};\n');
    fprintf(fid,'Spline(153) = {158, 157};\n');
    fprintf(fid,'Spline(155) = {287, 286, 288, 289, 290};\n');
    fprintf(fid,'Spline(156) = {290, 291, 292, 293};\n');
    fprintf(fid,'Spline(157) = {38, 39, 107, 106, 105};\n');
    fprintf(fid,'Spline(158) = {105, 398};\n');
    fprintf(fid,'Spline(160) = {126, 468, 123};\n');
    fprintf(fid,'Spline(161) = {123, 122, 121, 120};\n');
    fprintf(fid,'Spline(162) = {120, 119, 403, 118, 404};\n');
    fprintf(fid,'Line(163) = {126, 407};\n');
    fprintf(fid,'Line(164) = {466, 189};\n');
    fprintf(fid,'Line(165) = {189, 467};\n');
    fprintf(fid,'Line(166) = {467, 462};\n');
    fprintf(fid,'// Surfaces\n');
    fprintf(fid,'Curve Loop(1) = {14, 13, 10, 11, 12};\n');
    fprintf(fid,'Plane Surface(1) = {1};\n');
    fprintf(fid,'Curve Loop(37) = {3, 4, 157, 158, 42, -9, -8, -7, -6};\n');
    fprintf(fid,'Plane Surface(2) = {37};\n');
    fprintf(fid,'Curve Loop(4) = {1, -127, -125,-128, -118, -107, -105, -103, -69,  -4, -3, -2};\n');
    fprintf(fid,'Plane Surface(3) = {4};\n');
    fprintf(fid,'Curve Loop(5) = {18, -17, -16, -15, 8, 9, 19};\n');
    fprintf(fid,'Plane Surface(4) = {5};\n');
    fprintf(fid,'Curve Loop(6) = {45, 46, 47, 44};\n');
    fprintf(fid,'Plane Surface(5) = {6};\n');
    fprintf(fid,'Curve Loop(7) = {49, 50, 155, 156};\n');
    fprintf(fid,'Plane Surface(6) = {7};\n');
    fprintf(fid,'Curve Loop(8) = {23, 24, 25, -65, -64, -62, 27, 28, -32, -31, -30, -29, 22};\n');
    fprintf(fid,'Plane Surface(7) = {8};\n');
    fprintf(fid,'Curve Loop(9) = {20, -21};\n');
    fprintf(fid,'Plane Surface(8) = {9};\n');
    fprintf(fid,'Curve Loop(10) = {40, 41, 42, 19, 18, -39};\n');
    fprintf(fid,'Plane Surface(9) = {10};\n');
    fprintf(fid,'Curve Loop(11) = {150, 60, 56, 61, -67, -66, -64, 63, 149};\n');
    fprintf(fid,'Plane Surface(10) = {11};\n');
    fprintf(fid,'Curve Loop(12) = {58, -59, -57, 56};\n');
    fprintf(fid,'Plane Surface(11) = {12};\n');
    fprintf(fid,'Curve Loop(13) = {54, -55, -51, 52, 53};\n');
    fprintf(fid,'Plane Surface(12) = {13};\n');
    fprintf(fid,'Curve Loop(38) = {69, 70, 71, -162, -161, -160, 163, 73, -74, -33, 60, 57, 59, -58, 61, -67, -66, 65, -25, -24, -23, -22, 29, 30, 31, 32, -28, -27, 62, 63, 149, 38, 17, -39, 40, 41, -158, -157};\n');
    fprintf(fid,'Plane Surface(13) = {38};\n');
    fprintf(fid,'Curve Loop(19) = {82, 78, 79, 152, 153, -92, -104, -103, 70, 71, -162, -161, -160};\n');
    fprintf(fid,'Plane Surface(14) = {19};\n');
    fprintf(fid,'Curve Loop(20) = {89, -90};\n');
    fprintf(fid,'Plane Surface(15) = {20};\n');
    fprintf(fid,'Curve Loop(21) = {95, 151, -94, -93, -73, 75, 152};\n');
    fprintf(fid,'Curve Loop(22) = {89, -90};\n');
    fprintf(fid,'Plane Surface(16) = {21, 22};\n');
    fprintf(fid,'Curve Loop(23) = {82, 78, 79, -75, -163};\n');
    fprintf(fid,'Plane Surface(17) = {23};\n');
    fprintf(fid,'Curve Loop(24) = {96, 92, -153, 95};\n');
    fprintf(fid,'Plane Surface(18) = {24};\n');
    fprintf(fid,'Curve Loop(25) = {148, 146, 164, 165, 166};\n');
    fprintf(fid,'Plane Surface(19) = {25};\n');
    fprintf(fid,'Curve Loop(26) = {98, -97};\n');
    fprintf(fid,'Plane Surface(20) = {26};\n');
    fprintf(fid,'Curve Loop(27) = {101, 102, 94, 98, 99, 100};\n');
    fprintf(fid,'Plane Surface(21) = {27};\n');
    fprintf(fid,'Curve Loop(28) = {132, -130, 133, 131};\n');
    fprintf(fid,'Plane Surface(22) = {28};\n');
    fprintf(fid,'Curve Loop(29) = {110, 115, 116, -107, -106, 108, 109};\n');
    fprintf(fid,'Plane Surface(23) = {29};\n');
    fprintf(fid,'Curve Loop(30) = {138, -139, -136, -134, 135, 137};\n');
    fprintf(fid,'Plane Surface(24) = {30};\n');
    fprintf(fid,'Curve Loop(31) = {142, 144, 145, -128, -143, -138, -137, 140, 141};\n');
    fprintf(fid,'Plane Surface(25) = {31};\n');
    fprintf(fid,'Curve Loop(32) = {122, -123, -121, 124};\n');
    fprintf(fid,'Plane Surface(26) = {32};\n');
    fprintf(fid,'Curve Loop(39) = {105, -106, 108, 109, 110, 115, 116, 118, -143, -139, -136, -134, 135, 140, 141, 142, 144, 145, 125, 126, 121, 123, -122, -120, -119, -101, -100, -99, -97, -151, 96, -104};\n');
    fprintf(fid,'Plane Surface(27) = {39};\n');
    fprintf(fid,'Curve Loop(36) = {33, 34, 16, -38, 150};\n');
    fprintf(fid,'Plane Surface(28) = {36};\n');
    fprintf(fid,'// bodies\n');
    fprintf(fid,'For i In {1:28}\n');
    fprintf(fid,'vanorm~{i} = Extrude {1500, 0, 0}{ Surface{i}; };\n');
    fprintf(fid,'EndFor\n');
    fprintf(fid,'Coherence;\n');
    fprintf(fid,'// land Box\n');
    fprintf(fid,'Box ( 1001 ) = { -20000,-20000,400,40000,40000,19600 };\n');
    fprintf(fid,'Box ( 1004 ) = { -20000,-20000,100,40000,19250,300 };\n');
    fprintf(fid,'Box ( 1005 ) = { -20000,750,154,40000,19250,246 };\n');
    fprintf(fid,'Curve Loop(422) = {505, 511, 531, 743, 639, 741, 658, 739, 671, 669, 804, -661, -659, 797, -791, -548, 806, -753, -513, -508};\n');
    fprintf(fid,'Plane Surface(387) = {422};\n');
    fprintf(fid,'Curve Loop(423) = {506, 512, 532, 744, 640, 742, 687, 740, 700, 698, 805, -690, -688, 798, -792, -582, 807, -754, -519, -509};\n');
    fprintf(fid,'Plane Surface(388) = {423};\n');
    fprintf(fid,'vbulk1[] = Extrude {-19250, 0, 0}{ Surface{387};};\n');
    fprintf(fid,'vbulk2[] = Extrude {19250, 0, 0}{ Surface{388};};\n');
    fprintf(fid,'// air Box\n');
    fprintf(fid,'Box ( 1011 ) = { -20000,-20000,-20000,40000,40000,20000 };\n');
    fprintf(fid,'Box ( 1014 ) = { -20000,-20000,0,40000,19250,100 };\n');
    fprintf(fid,'Box ( 1015 ) = { -20000,750,0,40000,19250,154 };\n');
    fprintf(fid,'Curve Loop(485) = {735, -737, -671, -739, -658, -741, -639, -743, -531, -511, -505, -745};\n');
    fprintf(fid,'Plane Surface(449) = {485};\n');
    fprintf(fid,'Curve Loop(486) = {736, -738, -700, -740, -687, -742, -640, -744, -532, -512, -506, -746};\n');
    fprintf(fid,'Plane Surface(450) = {486};\n');
    fprintf(fid,'vbulk11[] = Extrude {-19250, 0, 0}{ Surface{449};};\n');
    fprintf(fid,'vbulk12[] = Extrude {19250, 0, 0}{ Surface{450};};\n');
    fprintf(fid,'// Coherence after all volumes, and physical entities\n');
    fprintf(fid,'Coherence;\n');
    fprintf(fid,'vv1() = BooleanUnion { Volume{3};Delete; } { Volume{1011,1014,1015,1016,1017};Delete; };\n');
    fprintf(fid,'vv2() = BooleanUnion { Volume{29};Delete; } { Volume{31,9,28,16,1001,1004,1005,1006,1007};Delete; };\n');
    fprintf(fid,'Coherence;\n');
    fprintf(fid,'Physical Volume("Air", 11) = {vv1(0)};\n');
    fprintf(fid,'Physical Volume("LeftAndRight_4400", 12) = {vv2(0)};\n');
    fprintf(fid,'Physical Volume("LowStatra_20", 13) = {4,21,23,25,26};\n');
    fprintf(fid,'Physical Volume("Anorm1_4", 14) = {1,8,10,11,12,5,6,15,19,22,24,20,17};\n');
    fprintf(fid,'Physical Volume("Anorm2_1175", 15) = {7,14};\n');
    fprintf(fid,'Physical Volume("Anorm3_385", 16) = {30,18};\n');
    fprintf(fid,'// source, source in air, and source field\n');
    fprintf(fid,'p1=newp;Point(p1)={'+string(coor(i,1))+','+string(coor(i,2))+','+string(coor(i,3))+'};\n');
    fprintf(fid,'l1=newreg; Circle(l1) = {'+string(coor(i,1))+','+string(coor(i,2))+','+string(coor(i,3))+',2};\n');
    fprintf(fid,'Physical Curve ( "lsource1",51 ) = {l1};\n');
    fprintf(fid,'Curve {l1} In Volume {1008};\n');
    fprintf(fid,'Field[1] = Distance;\n');
    fprintf(fid,'Field[1].NNodesByEdge = 60;\n');
    fprintf(fid,'Field[1].EdgesList = {l1};\n');
    fprintf(fid,'Field[2] = Threshold;\n');
    fprintf(fid,'Field[2].IField = 1;\n');
    fprintf(fid,'Field[2].LcMin = 0.2;\n');
    fprintf(fid,'Field[2].LcMax = 4000;\n');
    fprintf(fid,'Field[2].DistMin = 0;\n');
    fprintf(fid,'Field[2].DistMax = 20000;\n');
    fprintf(fid,'// receiver and receiver field\n');
    fprintf(fid,'Field[3] = Distance;\n');
    fprintf(fid,'Field[3].NodesList = {p1};\n');
    fprintf(fid,'Field[4] = Threshold;\n');
    fprintf(fid,'Field[4].IField = 3;\n');
    fprintf(fid,'Field[4].LcMin = 0.1;\n');
    fprintf(fid,'Field[4].LcMax = 4000;\n');
    fprintf(fid,'Field[4].DistMin = 0;\n');
    fprintf(fid,'Field[4].DistMax = 20000;\n');
    fprintf(fid,'Field[5] = Min;\n');
    fprintf(fid,'Field[5].FieldsList = {2,4};\n');
    fprintf(fid,'Background Field = 5;\n');
    fclose(fid);
    % meshing
%     [~,~] = system('./gmsh mod3_small.geo -3 -o mod3_small.msh');
    if exist('mod3_small.mesh','file')==2
        delete('mod3_small.mesh');
    end
    [~,~] = system('./gmsh mod3_small.geo -3 -o mod3_small.mesh');
    % simulation
    if exist('mod3_small.field','file')==2
        delete('mod3_small.field');
    end
    if exist('mod3_small.field2','file')==2
        delete('mod3_small.field2');
    end
    system(prog_mpi+" -n "+mpi_thread+" -genv OMP_NUM_THREADS="+omp_thread+" -genv I_MPI_PIN_DOMAIN=omp ./FETD -f mod3_small -v none");
    % EXTRACT
    tmp = importdata('mod3_small.field2');
    field(:,i) = tmp(:,3);
    save(sprintf('field%d.mat',i),'field');
    toc
end

[Y,X] = meshgrid(linspace(-750,750,16),linspace(-750,750,16));
field1 = reshape(field(1,:),16,16);
field2 = reshape(field(2,:),16,16);
field3 = reshape(field(3,:),16,16);
field4 = reshape(field(4,:),16,16);
field5 = reshape(field(5,:),16,16);
field6 = reshape(field(6,:),16,16);

f1=subplot(2,3,1);
contourf(X,Y,field1);
set(f1,'ColorScale','log');
title('t=5×10^{-6}s');
axis equal;
set(gca,'XTick',(-750:500:750));
set(gca,'YTick',(-750:500:750));
% ylabel('y(m)','position',[-1000,-1000]);
h1=colorbar('southoutside');
f2=subplot(2,3,2);
contourf(X,Y,field2);
set(f2,'ColorScale','log');
title('t=10^{-5}s');
axis equal;
set(gca,'XTick',(-750:500:750));
set(gca,'YTick',(-750:500:750));
h2=colorbar('southoutside');
f3=subplot(2,3,3);
contourf(X,Y,field3);
set(f3,'ColorScale','log');
title('t=5×10^{-5}s');
axis equal;
set(gca,'XTick',(-750:500:750));
set(gca,'YTick',(-750:500:750));
h3=colorbar('southoutside');
f4=subplot(2,3,4);
contourf(X,Y,field4);
set(f4,'ColorScale','log');
title('t=10^{-4}s');
axis equal;
set(gca,'XTick',(-750:500:750));
set(gca,'YTick',(-750:500:750));
h4=colorbar('southoutside');
f5=subplot(2,3,5);
contourf(X,Y,field5);
set(f5,'ColorScale','log');
title('t=5×10^{-4}s');
axis equal;
set(gca,'XTick',(-750:500:750));
set(gca,'YTick',(-750:500:750));
% xlabel('x(m)');
h5=colorbar('southoutside');
f6=subplot(2,3,6);
contourf(X,Y,field6);
set(f6,'ColorScale','log');
title('t=10^{-3}s');
axis equal;
set(gca,'XTick',(-750:500:750));
set(gca,'YTick',(-750:500:750));
h6=colorbar('southoutside');

set(f1,'position',[.05 .62 .32 .32]);
set(f2,'position',[.35 .62 .32 .32]);
set(f3,'position',[.65 .62 .32 .32]);
set(f4,'position',[.05 .13 .32 .32]);
set(f5,'position',[.35 .13 .32 .32]);
set(f6,'position',[.65 .13 .32 .32]);

set(h1,'Position',[.1 .56 .22 .012]);
h1.FontAngle = 'italic';
h1.FontWeight = 'bold';
set(h2,'Position',[.4 .56 .22 .012]);
h2.FontAngle = 'italic';
h2.FontWeight = 'bold';
set(h3,'Position',[.7 .56 .22 .012]);
h3.FontAngle = 'italic';
h3.FontWeight = 'bold';
set(h4,'Position',[.1 .07 .22 .012]);
h4.FontAngle = 'italic';
h4.FontWeight = 'bold';
set(h5,'Position',[.4 .07 .22 .012]);
h5.FontAngle = 'italic';
h5.FontWeight = 'bold';
set(h6,'Position',[.7 .07 .22 .012]);
h6.FontAngle = 'italic';
h6.FontWeight = 'bold';

exportgraphics(gcf,'figure23.png','Resolution',300);