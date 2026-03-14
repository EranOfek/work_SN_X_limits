%%

Args.SearchRadius = 0.1;  % [deg]
Args.AperRadius = 5;
Args.Annulus    = [30 80];
Args.ReturnCol  = {'TIME','X','Y','RAWX','RAWY','DETX','DETY',}
CritEdgeDist    = max(Args.Annulus);
PixScale        = 2.36;
PSF             = imUtil.kernel2.king([5.8./PixScale -1.55],[15 15]);

PWD = pwd;




% get List of SNe
[T, CsvFile, ZipFile]=VO.TNS.downloadAll;
% select SNe from 1990 till 2025


Nsn = size(T,1);

% for each SNe
for Isn=1:1:Nsn

    RA  = T.ra(Isn);
    Dec = T.dec(Isn);
    % query for Swift observations
    [Txrt,DU]  = VO.Swift.swiftXRT_obs(RAD, Dec, Args.SearchRadius);
    

    % select observations with large ExpTime
    Flag = Txrt.ExpTime>1000;
    Txrt = Txrt(Flag,:);
    % download all observations
    Txrt = VO.Swift.swiftXRT_wget(Txrt);



    %% find all dir of observations
    F=io.files.rdir('*xpcw*po_cl.evt*');
    
    Nobs = numel(F);
    for Iobs=1:1:Nobs
        cd(F(Iobs).folder)
        system('gzip -d *.gz');
        File = io.files.rdir('*xpcw*po_cl.evt');
        if numel(File)==1
            Info   = imUtil.swiftXRT.parseSwiftXRTFilenames(File.name);
            CCDSEC = imUtil.swiftXRT.getSwiftXRT_RawCorners(Info.Instrument, Info.Window)
            switch lower(Info.Window)


                
            % Read Evt file
            PL = PhotonsList.readPhotonsList1(File.name);
            PL.pi2energy;
    
            PL.selectEnergy([200 10000]);
    
            % analyze data
            PL.addSkyCoo;
            PL.markEventsNearBoundries(CCDSEC, CritEdgeDist);
    
    
    
            % Create Image;
            PL.Image;
    
            % Read ExpTime
            ExpTime = PL.HeaderData.getVal('EXPOSURE');
    
            % Select photons around source:
            % Inside AperRad, and with Annulus(1) to Annulus(2)
            
            % need to take into account image edges!
            Src = PL.getSrcPhotons(RA, Dec, 'InUnits','deg', 'Annulus',Args.Annulus, 'SearchRadius',Args.AperRadius, 'SearchRadiusUnits','pix',...
                                            'ReturnCol',PL.Events.ColNames);
            
            [X,Y]=PL.WCS.sky2xy(RA,Dec);
            % Match filter with PMF
            PMF = imUtil.poissNoise.poissonMatchedFilter(PSF, Src.BackPerPix, 3)
    
            FilteredImage = imUtil.filter.filter2_fast(PL.Image, PMF);
            FilteredImage(round(Y),round(X))
    
            % Read Flux
    
            % Estimate significance
    
            % Photometry
    
            % write data
    
        else
            % problem?
            numel(File)
        end
    
    % for each observation
       
end
