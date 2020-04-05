%converter  converter model class
%
%If a thyristor, diode, or vsc converter is added to the system, this class
%is invoked. Each converter in the system, regardless of the type, is
%represented with an object of the converter type. The converter class
%contains all properties that are of importance for all converter types
%even if the property isn't universally relevant for each type. It is
%important to note that some properties are assumed such as circuit
%impedances. The user can overwrite these properties however by directly
%accessign the properties of the relvant converter object and modifying
%them (this is possible since all properties are set as public). There are
%two main functions that each converter has, that is the constructor and
%the solveConverter function. The constructor is called when the converter
%object is created and added to the powerconverters object converters
%array. It takes in the relevant inputs from the user and converts these
%inputs over to p.u. values. It also creates three current source
%objects (1 for each phase), spectrum, and monitor objects in the OpenDSS
%file that contains the
%system that is the object of this harmonic case study. The solveConverter
%function performs two purposes. Given the particular converter properties
%as well as the PCC voltage harmonic spectrum it will solve for the 3 phase
%currents of the particular converter for all harmonics of interest.
%Once these currents have been
%solved, it will then update the spectrum object for each of the three
%phases that correspond to the analysed converter. The three isource
%objects are then updated to read in the new updated spectrum objects as
%this is not done automatically by OpenDSS.

classdef converter < handle
    %Some properties are specified in the inputs.m file; others are
    %hard-coded in p.u. and then based on the user inputs are - in the
    %constructor - scaled accordingly to actual values.
    properties 
        %name that the particular converter object is referred to, type of
        %converter, bus # that the converter is shunt connected,
        %h = number of harmonics the user wants analysed
        name, type, bus, h
        Cdc, Idc, Q, Ldc, P, L, R %units: F, A, kvar, H, kW, H, Ohm
        
        %Vbase is in l-l,rms; Sbase is in 3 Phase.
        Vbase, Sbase, Zbase, Rdc, Vdc, Ibase %V, kVA, kVA, Ohm, Ohm, V, A
        mf, ma, sigma, mu, alpha, phaseAng, gamma
        L_pu, R_pu, Rdc_pu, Vdc_pu, Cdc_pu, Idc_pu, Q_pu, Ldc_pu, P_pu, err
        
        %Properties relevant to "thyristor6" and "thyristor12" converter
        %types:
        %Ldc, mu, L, R, Vbase, Sbase, Zbase, Rdc, Vdc, Ibase, P,alpha,
        %L_pu, R_pu, Rdc_pu, Vdc_pu, Ldc_pu, P_pu
        
        %properties relevant to "diode6" and "diode12" converter types:
        %L, R, Vbase, Sbase, Zbase, Rdc, Vdc, Ibase
        %Ldc, mu, phaseAng, gamma, L_pu, R_pu,Rdc_pu, Vdc_pu, Ldc_pu
        
        %Properties relevant to "vsc" converter types:
        %L, R, Vbase, Sbase, Zbase, Rdc, Vdc, Ibase
        %Cdc, Idc, Q, mf, ma, sigma, P, Cdc_pu, Idc_pu, Q_pu, P_pu, L_pu,
        %R_pu, Rdc_pu, Vdc_pu
    end
    methods
        %Constructor - Whenever a new converter object is added to the
        %system this function block is called.
        %This block takes in the inputs and corresondingly changes all the
        %properties that are relevant to that particular converter type. 
        function c_obj = converter(type,name,bus,Srated,Vdc,P,Q)
            global numHarmonics DSSText DSSCircuit
            c_obj.err = 0;
            c_obj.Sbase = Srated;
            c_obj.type = type;
            c_obj.name = name;
            c_obj.bus = bus;
            c_obj.h = numHarmonics;
            c_obj.P = P;
            
            %this block is used to extract the base voltage at the bus that
            %the converter is to be interfaced to---the base voltage is
            %specified in the original openDSS circuit file as kV l-n,rms. 
            DSSCircuit.SetActiveBus(c_obj.bus);
            DSSBus = DSSCircuit.ActiveBus;
            c_obj.Vbase = DSSBus.kVBase*1000;
            
            %calculating the per-phase base impedance from the Vbase and
            %Sbase values.
            c_obj.Zbase = 3*c_obj.Vbase^2/(c_obj.Sbase*1000); 
            %calculating the per-phase current base value.
            c_obj.Ibase = (c_obj.Sbase*1000)/(c_obj.Vbase)/sqrt(3);
            
            if strcmp(type,'vsc') %if one will explore this branch
                c_obj.Vdc = Vdc;
                c_obj.Rdc = 10^8*c_obj.Zbase; %Rdc set so that the p.u. is 10^8
                c_obj.Q = Q;
                c_obj.Cdc = abs(6*Srated/Vdc^2); %Set so that the dc-link cap will have at least 2 J of energy storage per kW rated. This is a general rule of thumb that limits the ripple on the dc-link to within a few percent.
                c_obj.mf = 15; %the switching frequency of the converter
                c_obj.Idc = P*1000/Vdc; %calculated based on the specified power and average dc-link voltage
            else
                c_obj.Ldc = 0.5*c_obj.Zbase/377; %Ldc set so that the p.u. is 0.5              
                c_obj.mu = -1;
                c_obj.Vdc = Vdc;
                if strcmp(type,'thyristor6') || strcmp(type,'thyristor12')
                    c_obj.Rdc = Vdc^2/(P*1000);
                else
                    c_obj.Vdc = 0;
                    if strcmp(type,'diode6')
                        c_obj.Rdc = 1.654^2/(P/c_obj.Sbase)*c_obj.Zbase;
                    else
                        c_obj.Rdc = 3.3080^2/(P/c_obj.Sbase)*c_obj.Zbase;
                    end
                end
            end
            c_obj.R = 0.02*c_obj.Zbase; %R set so that the p.u. is 0.02
            if strcmp(type,'diode12')
                c_obj.L = 0.2*c_obj.Zbase/377;
            elseif strcmp(type,'thyristor12')
                c_obj.L = 0.02*c_obj.Zbase/377;
            else
            	c_obj.L = 0.2*c_obj.Zbase/377; %L set so that the p.u. is 0.02
            end
            %c_obj.L = c_obj.L_pu*c_obj.Zbase/377;

            %the p.u. values are calculated based on the actual values and
            %the base values above. Keep in mind the p.u. values assume a
            %base frequency of 1 rad/s whereas the actual values assume a
            %base frequency of 60 Hz.
            c_obj.L_pu = 377*c_obj.L/c_obj.Zbase;
            c_obj.R_pu = c_obj.R/c_obj.Zbase;
            c_obj.Rdc_pu = c_obj.Rdc/c_obj.Zbase;
            c_obj.P_pu = c_obj.P/c_obj.Sbase;
            c_obj.Vdc_pu = c_obj.Vdc/(c_obj.Vbase*sqrt(2));

            if strcmp(type,'thyristor6') || strcmp(type,'thyristor12')
                c_obj.Rdc_pu = c_obj.Vdc_pu^2/c_obj.P_pu;
                c_obj.Rdc = c_obj.Rdc_pu*c_obj.Zbase;
            end
            if strcmp(type,'vsc')
                c_obj.Cdc_pu = 2*pi*60*c_obj.Cdc*c_obj.Zbase*sqrt(3);
                c_obj.Idc_pu = c_obj.Idc/1000/c_obj.Sbase*c_obj.Vbase*sqrt(2);
                c_obj.Q_pu = c_obj.Q/c_obj.Sbase;
                c_obj.Vdc_pu = c_obj.Vdc/(c_obj.Vbase*sqrt(2));
            else
                c_obj.Ldc_pu = 377*c_obj.Ldc/c_obj.Zbase;
            end
            
            %Adding 3 Isource objects to the OpenDSS circuit file, each
            %corresponding to a single phase. Assumed, at this point, to
            %contribute no current to the system.
            DSSText.Command = ['new Isource.',c_obj.name,'a Bus1=',c_obj.bus,'.1 basefreq = 60 Phases=1 amps=0'];
            DSSText.Command = ['new Isource.',c_obj.name,'b Bus1=',c_obj.bus,'.2 basefreq = 60 Phases=1 amps=0'];
            DSSText.Command = ['new Isource.',c_obj.name,'c Bus1=',c_obj.bus,'.3 basefreq = 60 Phases=1 amps=0'];
            %Defining new spectrum objects
            initSpectrum = zeros(numHarmonics,3);
            for i = 1:numHarmonics
                initSpectrum(i,1) = i;
            end
            
            %Adding three spectrum objects to the OpenDSS circuit file, one
            %for each of the three phases. Spectrum objects contain the harmonic magnitude and phase relative to the fundamental component. As far as I can tell, OpenDSS
            %does not allow you to directly update the spectrum object and
            %so first a csv file is created for each of the three phases.
            %The spectrum objects are then created for this converter. Each
            %new spectrum object is linked to the related csv file that was
            %just created. Thos files contain the harmonic spectra of that
            %particular converter. Finally, in order to monitor harmonics
            %at a particular node or line in the system, objects called
            %monitors must be placed at that particular node in the
            %network. Otherwise, the harmonic spectra of the converter's
            %PCC with the system will not be solved for which is the
            %desired outcome of each iteration. Each time the system is
            %solved, the csv's corresponding to the monitor objects are
            %updated. From these, we can extract the PCC voltage harmonic
            %spectra for the particular iteration.
            csvwrite([c_obj.name,'a_spectrum.csv'],initSpectrum);
            csvwrite([c_obj.name,'b_spectrum.csv'],initSpectrum);
            csvwrite([c_obj.name,'c_spectrum.csv'],initSpectrum);
            DSSText.Command = ['new Spectrum.',c_obj.name,'a numharm=',...
                num2str(numHarmonics),' csvfile=',c_obj.name,'a_spectrum.csv'];
            DSSText.Command = ['new Spectrum.',c_obj.name,'b numharm=',...
                num2str(numHarmonics),' csvfile=',c_obj.name,'b_spectrum.csv'];
            DSSText.Command = ['new Spectrum.',c_obj.name,'c numharm=',...
                num2str(numHarmonics),' csvfile=',c_obj.name,'c_spectrum.csv'];
            DSSText.Command = ['new monitor.',c_obj.name,...
                'a Isource.',c_obj.name,'a'];
            DSSText.Command = ['new monitor.',c_obj.name,...
                'b Isource.',c_obj.name,'b'];
            DSSText.Command = ['new monitor.',c_obj.name,...
                'c Isource.',c_obj.name,'c'];
        end
        %this function reads in the csv files containing the PCC voltage
        %harmonic spectra of the converter as input. It then solves the
        %converter model for the injection current harmonic spectra. This
        %harmonic spectra is used to update the spectrum and isource
        %objects of the openDSS file.
        function solveConverter(c_obj) %the only output we need is the current spectrum
            global OpenDSSFile invSymMtx DSSText numHarmonics
            
            if c_obj.err == 0
                lenHarm = 2*c_obj.h+1;
                
                %this reads in the abc phase voltages at the PCC from the csv
                %files from the most recent simulation.
                Va = csvread([OpenDSSFile,'_Mon_',c_obj.name,'a_1.csv'],1,2,...
                    [1,2,c_obj.h,3]);
                Vb = csvread([OpenDSSFile,'_Mon_',c_obj.name,'b_1.csv'],1,2,...
                    [1,2,c_obj.h,3]);
                Vc = csvread([OpenDSSFile,'_Mon_',c_obj.name,'c_1.csv'],1,2,...
                    [1,2,c_obj.h,3]);
                
                %this loop transforms the abc voltages into the positive and
                %negative sequence refererence frame. Then the real and
                %imaginary components of the positive and negative sequence
                %frame values are extracted and added to the voltage vector V.
                V = zeros(lenHarm*2,1);
                count = 0;
                for i = 1:c_obj.h
                    Vav = Va(i,1)*exp(1i*Va(i,2)*pi/180);
                    Vbv = Vb(i,1)*exp(1i*Vb(i,2)*pi/180);
                    Vcv = Vc(i,1)*exp(1i*Vc(i,2)*pi/180);
                    pnArr = invSymMtx*[Vav;Vbv;Vcv];
                    V(2*c_obj.h+3+count:2*c_obj.h+4+count) = [real(pnArr(2));imag(pnArr(2))];
                    V(2*c_obj.h-1-count:2*c_obj.h-count) = [real(pnArr(3));imag(pnArr(3))];
                    count = count + 2;
                end
                
                %the converter solver takes p.u. values as the input so this
                %operation pu-tizes the PCC voltages.
                V = V/c_obj.Vbase;
                
                outputFile = [c_obj.name,'V.mat']; %creates string called outputFile---the name of the file we want to create
                save(outputFile,'V');%saves the vector V into the .mat file with the name contained in the variable outputFile
                
                %In this "if" block, the object type is used to determine which
                %converter solver function to run. Each converter type operates
                %uniquely and thus a separate solving function is called for
                %each of the 5 different converter types. The output of each of
                %the solver functions is the injection current harmonic spectra
                %for the converter given the inputs and the PCC voltages.
                if strcmp(c_obj.type,'vsc')
                    [Ia Ib Ic fundCurrents c_obj.ma c_obj.sigma Rdc_new c_obj.err] = solveVSC(V,...
                        c_obj.L_pu,c_obj.R_pu,c_obj.mf,c_obj.Rdc_pu,c_obj.Cdc_pu,...
                        c_obj.Idc_pu,c_obj.h,c_obj.Vdc_pu,c_obj.Q_pu);
                elseif strcmp(c_obj.type,'thyristor6')
                    [Ia Ib Ic fundCurrents c_obj.alpha c_obj.mu] = solveThyristor6(V,c_obj.L_pu,...
                        c_obj.R_pu,c_obj.Rdc_pu,c_obj.Ldc_pu,...
                        c_obj.h,c_obj.Vdc_pu,c_obj.P_pu);
                elseif strcmp(c_obj.type,'diode6')
                    [Ia Ib Ic fundCurrents c_obj.phaseAng c_obj.gamma c_obj.mu c_obj.P_pu Rdc_new] = solveDiode6(V,...
                        c_obj.L_pu,c_obj.R_pu,c_obj.Rdc_pu,c_obj.Ldc_pu,...
                        c_obj.h,c_obj.Vdc_pu,c_obj.P_pu,c_obj.mu,c_obj.gamma,c_obj.phaseAng);
                elseif strcmp(c_obj.type,'thyristor12')
                    [Ia Ib Ic fundCurrents c_obj.alpha c_obj.mu] = solveThyristor12(V,c_obj.L_pu,...
                        c_obj.R_pu,c_obj.Rdc_pu,c_obj.Ldc_pu,...
                        c_obj.h,c_obj.Vdc_pu,c_obj.P_pu);
                elseif strcmp(c_obj.type,'diode12')
                    [Ia Ib Ic fundCurrents c_obj.phaseAng c_obj.gamma c_obj.mu c_obj.P_pu Rdc_new] = solveDiode12(V,...
                        c_obj.L_pu,c_obj.R_pu,c_obj.Rdc_pu,c_obj.Ldc_pu,...
                        c_obj.h,c_obj.Vdc_pu,c_obj.P_pu,c_obj.mu,c_obj.gamma,c_obj.phaseAng);
                end
                
                if ~strcmp(c_obj.type,'thyristor12') && ~strcmp(c_obj.type,'thyristor6')
                    c_obj.Rdc_pu = Rdc_new;
                end
                
                %To update the OpenDSS spectrum objects correspondign to the
                %converter, the csv files of the spectrum objects must first be
                %updated. Then the spectrum objects have to be updated to read
                %in the updated csv current harmonic spectrum values.
                csvwrite([c_obj.name,'a_spectrum.csv'],Ia);
                csvwrite([c_obj.name,'b_spectrum.csv'],Ib);
                csvwrite([c_obj.name,'c_spectrum.csv'],Ic);
                
                %Spectrums are being asked to read in the values of the updated
                %csv files
                DSSText.Command = ['edit Spectrum.',c_obj.name,'a numharm=',...
                    num2str(numHarmonics),' csvfile=',c_obj.name,'a_spectrum.csv'];
                DSSText.Command = ['edit Spectrum.',c_obj.name,'b numharm=',...
                    num2str(numHarmonics),' csvfile=',c_obj.name,'b_spectrum.csv'];
                DSSText.Command = ['edit Spectrum.',c_obj.name,'c numharm=',...
                    num2str(numHarmonics),' csvfile=',c_obj.name,'c_spectrum.csv'];
                
                %The spectrum objects contain relative values only. For our
                %case the values are relative to the fundamental frequency
                %component. As such, the fundamental current magnitude and
                %phase for the abc phases needs to be provided. This is done by
                %updating the isource objects directly. This block of code does
                %exactly this. The fundamental magnitude (in pu) and phase for
                %each of the three-phases is outputed from the solver
                %functions. So the first few lines deal with extracting these
                %values and transforming the pu magnitude value into an actual
                %current in Amps. The isource objects corresponding to the
                %converter are then updated with the fundamental magnitude and
                %phase as well as the updated spectrum objects.
                fundCurrents(1:end,1) = fundCurrents(1:end,1)*c_obj.Ibase*sqrt(2);
                ampVal = num2str(fundCurrents(1,1));
                angleVal = num2str(fundCurrents(1,2));
                DSSText.Command = ['Edit Isource.',c_obj.name,'a enabled=true amps = ',ampVal, ' angle = ',angleVal,' spectrum=',c_obj.name,'a'];
                ampVal = num2str(fundCurrents(2,1));
                angleVal = num2str(fundCurrents(2,2));
                DSSText.Command = ['Edit Isource.',c_obj.name,'b enabled=true amps = ',ampVal, ' angle = ',angleVal,' spectrum=',c_obj.name,'b'];
                ampVal = num2str(fundCurrents(3,1));
                angleVal = num2str(fundCurrents(3,2));
                DSSText.Command = ['Edit Isource.',c_obj.name,'c enabled=true amps = ',ampVal, ' angle = ',angleVal,' spectrum=',c_obj.name,'c'];
            end
        end
    end
end
