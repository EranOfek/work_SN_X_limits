%%

% TNS=VO.TNS.downloadAll;

% TNS.mat, include Z for host galaxies (excluding FRBs)
% generated in the work_SN_RadioPRS project on 2026-Mar-7
load TNS.mat


%%

Args.SearchRadius  = 0.2;  % [deg]
Args.MinExpTime    = 400;
Args.GetData       = true;
Args.AnalyzeData   = true;
Args.PauseDownload = 5;
Args.FileTemplate  = '*pcw*po_cl.evt';
Args.PixSize       = 2.36;  % Swift-XRT [arcsec/pix]
Args.FieldRadius   = 0.15; % removing near edge / 0.196;  % [deg]
Args.Annulus       = [15 120];  % arcsec
Args.AperRadius    = [5 10];    % arcsec

RAD = 180./pi;
ARCSEC_DEG = 3600;

AnnulusArea  = pi.*(Args.Annulus(2).^2 - Args.Annulus(1).^2);  % [arcsec^2]
AperArea     = pi.*Args.AperRadius(:).'.^2;

Ntns = numel(TNS);
FoundCounter = 0;
K            = 0;
tic;
for Itns=1:1:1000
    %Ntns
    RA   = TNS.ra(Itns);
    Dec  = TNS.declination(Itns);
    [T,DU]  = VO.Swift.swiftXRT_obs(RA, Dec, Args.SearchRadius);
    Nt = size(T,1);

    if Nt>0
        FlagET = T.xrt_exposure>Args.MinExpTime;
        T      = T(FlagET,:);

        Nt = size(T,1);
        if Nt>0
            FoundCounter = FoundCounter + 1;
            Result(FoundCounter).Itns     = Itns;
            Result(FoundCounter).TableXRT = T;  % Store the observation data

            % wget data
            if Args.GetData
                VO.Swift.swiftXRT_wget(Result(FoundCounter).TableXRT);

                % analyze data
                if Args.AnalyzeData
                    % List of obsid
                    ListObsID = Result(FoundCounter).TableXRT.obsid;
                    Nobsid    = numel(ListObsID);
                    % for each obsid
                    for Iobsid=1:1:Nobsid
                        % wait for download to finish
                        DownloadNotFinshed = true;
                        while DownloadNotFinshed
                            PathObsID = io.files.findDirBySubString(ListObsID{Iobsid});
                            if any(contains(PathObsID,'wget')) || numel(PathObsID)>1
                                % not done yet
                                pause(Args.PauseDownload);
                            else
                                DownloadNotFinshed = false;
                            end
                        end

                        PWD = pwd;
                        
                        cd(PathObsID{1});

                        Files = dir('sw*.ev*');
                        % uncompress
                        io.files.uncompress({Files.name});

                        % locate relevant file for analysis
                        FileForAnalysis = dir(Args.FileTemplate);
                        Nffa = numel(FileForAnalysis);

                        % extract photon list
                        for Iffa=1:1:Nffa
                            try
                                PL = PhotonsList(FileForAnalysis(Iffa).name);
                    
                                %--- analyze data ---
                                % Is the data clean?
                                %PL.populateBadTimes
                                %[a,b] = PL.nphotons;
                                % energy cut
                                PL = PL.pi2energy;
                                PL.selectEnergy([200 8000]);
                                [Nph,NphG] = PL.nphotons;
                                % Add RA/Dec to catalog
                                PL.addSkyCoo
                                % populate the image field, X/Y correspinds to image [default is sky coo]
                                [~,Image] = PL.constructImage; 
                                
    
                                
                                [X_Target, Y_Target] = PL.sky2xy(RA, Dec);
                                X_ph = PL.Events.getCol('x',false,false,'CaseSens',false);
                                Y_ph = PL.Events.getCol('y',false,false,'CaseSens',false);
    
                                RA_ph   = PL.Events.getCol('RA');
                                Dec_ph  = PL.Events.getCol('Dec');
                                RA_PNT  = PL.HeaderData.Key.RA_PNT;
                                Dec_PNT = PL.HeaderData.Key.DEC_PNT;
                                ExpTime = PL.HeaderData.Key.EXPOSURE;
                                DistCenter = celestial.coo.sphere_dist_fast(RA./RAD, Dec./RAD, RA_PNT./RAD, Dec_PNT./RAD);
    
                                DistTarget = celestial.coo.sphere_dist_fast(RA./RAD, Dec./RAD, RA_ph./RAD, Dec_ph./RAD);
                                %if sum(X_ph>X_Target)>5 && sum(X_ph<X_Target)>5 && sum(Y_ph>Y_Target)>5 && sum(Y_ph<Y_Target)>5
                                    
                                % source is not near edge.
    
                                    
                                % background
                                DistTargetAS = DistTarget.*RAD.*ARCSEC_DEG;
                                FlagBack = DistTargetAS>=Args.Annulus(1) & DistTargetAS<=Args.Annulus(2);
                                BackN    = sum(FlagBack);
                                BackRate = BackN./AnnulusArea;   % [back rate per exosure per arcsec^2]
                                    
                                % aperture photometry
                                FlagAper = DistTargetAS(:)<=Args.AperRadius(:).';
                                AperFlux = sum(FlagAper, 1);  % including background
                                AperExpBack  = BackRate.*AperArea;
                                
                                
                                % Poisson noise matched filter in position
                        

                                % store data
                                K = K + 1;
                                Meas(K).Itns         = Itns;
                                Meas(K).FoundCounter = FoundCounter;
                                Meas(K).Iffa         = Iffa;
                                Meas(K).RA           = RA;
                                Meas(K).Dec          = Dec;
                                Meas(K).RA_PNT       = RA_PNT;
                                Meas(K).Dec_PNT      = Dec_PNT;
                                Meas(K).ExpTime      = ExpTime;
                                Meas(K).ObsJD        = PL.HeaderData.julday;
                                Meas(K).DistCenter   = DistCenter.*RAD;  % [deg]
                                %
                                Meas(K).BackN        = BackN;
                                Meas(K).BackRate     = BackRate;
                                Meas(K).AperFlux     = AperFlux;
                                Meas(K).AperExpBack  = AperExpBack;

                              
                                % probability that the flux is from the
                                % background:
                                Meas(K).ProbFromBack = poisscdf(AperFlux, AperExpBack, 'upper');

                                % Transient
                                Meas(K).TranDiscJD  = juliandate(TNS.discoverydate(Itns));
                                Meas(K).TranDiscMag = TNS.discoverymag(Itns);
                                Meas(K).TranZ       = TNS.redshift(Itns);
                                Meas(K).TranHostZ   = TNS.Z(Itns);

                                if Meas(K).ProbFromBack<1e-5
                                    'a'
                                end



                            %end
                            catch ME
                                [Itns]
                            end


                        end

                        cd(PWD);
                    end
                end
            end
        end
    end

end
toc

%%

VO.Swift.swiftXRT_wget(210.75, 54.3, 'searchRadius', 1);
% VO.Swift.swiftXRT_wget(T);