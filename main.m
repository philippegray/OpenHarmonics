%INPUTS OpenHarmonics User Interface.
%
%This file is where the user interfaces with OpenHarmonics. The 
%simulation parameters, converter and monitor objects are specified in
%this file. Once, all inputs have been added, the user must hit Run (F5)
%and OpenHarmonics will simulate the network. The outputs of the simulation 
%will appear in the Workspace block on the top right side of the MATLAB GUI.
%

%clear all;
tic;
global numHarmonics DSSText OpenDSSFile DSSCircuit %don't change

%numHarmonics - changes the number of harmonic orders that you would like to analyze.
%For example: numHarmonics = 27, would generate the harmonic spectrum up to and
%including, the 27th harmonic.
numHarmonics = 17;
%This line references the full path of the OpenDSS circuit file
%containing the network to be simulated. Make sure that the OpenDSS file is
%located in the \OpenDSSFiles subfolder.
OpenDSSFile = 'IEEE13Nodeckt';
OpenDSSFileLoc = ...
    '"Include Directory Path Here"\supporting_opendss_files\';
OpenDSSFile = [OpenDSSFileLoc,OpenDSSFile];

%Initializes the communication interface (COM) between MATLAB and OpenDSS
[DSSObj DSSText DSSCircuit DSSSolution DSSMonitors busBasekVs] = startOpenDSS(OpenDSSFile); %don`t change

%Creates the pwrconverters object for the simulated network
sys = pwrconverters(); %don't change line

%USER INPUT - Add Converter ObDjects...Here*********************************
%format = add(sys,'type','name','bus','Srated','Vdc','P','Q')
%Refer to the user manual further explanation of the input format
%vsc,diode6,diode12,thyristor6,thyristor12

sys = add(sys,'vsc','vsc1','634',100,862,55,-10);
%sys = add(sys,'vsc','Conv1','634',100,480/sqrt(3)*sqrt(2)*2.2,66,-70); %type,name,bus,Srated(kVA),Vdc(V),P(kW),Q(kvar)%sys = add(sys,'vsc','Conv2','670',1,2000,0.5,-0.1);
%sys = add(sys,'thyristor12','Conv4','634',100,480/sqrt(3)*sqrt(2)*0.6,40,'');
sys = add(sys,'diode6','Conv5','634',50,'',30,'');%
%sys = add(sys,'diode12','diode1','634',100,'',23,'');%
%sys = add(sys,'thyristor12','thy1','634',50,480/sqrt(3)*sqrt(2)*1.0,40,'');
%sys = add(sys,'thyristor12','thy2','634',10,0,6,'');%
%sys = add(sys,'thyristor6','Conv2','670',100,2400/sqrt(3)*sqrt(2)*1.6,88,0);
%sys = add(sys,'thyristor6','Conv3','634',100,480/sqrt(3)*sqrt(2)*0.3,99,'');
%sys = add(sys,'thyristor12','Conv4','634',100,480/sqrt(3)*sqrt(2)*2.2,100,'');
%**************************************************************************

%USER INPUT - List the bus names that you would like the voltage harmonic
%spectrums for into VSolve.
VSolve = {'670'};
ISolve = {}; %This isn't working currently.

%Calls and excecutes the OpenHarmonics solution engine.
ActualHarmonics(DSSText,DSSCircuit,DSSSolution,DSSMonitors,sys,ISolve,...
    VSolve);

%This loads the results of the simulation into the MATLAB Workspace
load([OpenDSSFileLoc,'outputs.mat']);
type([OpenDSSFileLoc,'errorFile.txt']);

toc;

 