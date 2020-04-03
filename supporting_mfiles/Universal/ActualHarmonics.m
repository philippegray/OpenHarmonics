%ActualHarmonics OpenHarmonics Solution Engine
%
%This is the solution engine that takes the inputs from the
%inputs.m file and drives towards a solution.
%

function ActualHarmonics(Text,Circuit,DSSSolution,DSSMonitors,sys,ISolve,VSolve)

global DSSText SymMtx invSymMtx DSSCircuit numHarmonics OpenDSSFile

DSSText = Text;
DSSCircuit = Circuit;


%"tol" refers to % tolerance between the results from one iteration to the
%next. More specifically, the desired result is a current harmonic spectra
%for the converter. After each iteration of solving for this current
%harmonic spectra, because a difference in load current harmonic spectra
%influences the PCC voltages, the current harmonic spectra will be
%different. This "tol" refers to the difference between each of the
%individual harmonic spectra terms from one iteration to the next. As an
%example, if tol = 1, we solve another iteration if the percent difference
%between the harmonic magnitudes for frequencies -h,...,+h is greater than
%1%.
tol = 5;
SymMtx = [1 1 1;1 exp(1i*240*pi/180) exp(1i*120*pi/180);...
    1 exp(1i*120*pi/180) exp(1i*240*pi/180)];
invSymMtx = inv(SymMtx);

%MonitorList is a string array that contains the names of ALL monitor
%objects in the system
MonitorList = DSSMonitors.AllNames;
%BusNames is a string array that contains the names of ALL buses in the system 
%BusNames = DSSCircuit.AllBusNames;

%%%%%%%%%%%%%%%%%%%%%
%The fundamental frequency solution is solved. This is a necessary first
%step before the harmonic frequencies are solved.
DSSSolution.Solve; 

%to test to see if placing monitors worked
%CircuitElementList = DSSCircuit.AllElementNames;

%harmVec is used to update the OpenDSS file. harmVec is used to tell
%OpenDSS to solve for the harmonic frequencies 1...numHarmonics.
harmVec = harmonicsVector(numHarmonics);

DSSText.Command = ['set harmonics = (',harmVec,')'];%tells OpenDSS which harmonics frequencies to solve the system at when the solve harmonics command is called.
DSSText.Command = 'solve mode = harmonics'; %tells OpenDSS to now simulate the system at those harmonics specified in the step before.

%Monitor objects at this point contain the harmonic spectra at their
%respective nodes. In order to manipulate and use this data though, the
%data must be exported into csv files that can be read into as
%manipulatable arrays. That is what the exportMonitors function does.

if length(MonitorList) ~= 1
    exportMonitors(MonitorList); %exports the results of the monitors into individual .csv files. There is one .csv file per Monitor object.
    
    Vprev = zeros(2*numHarmonics+1,length(MonitorList)/3);
    
    pass = 0;
    loopCount = 0;
    
    %This loop iterates until one of two conditions are met. First, if the
    %criteria condition is met, which is whether the PCC voltage harmonic
    %spectra between iterations for all buses connected to a converter remain
    %with the tolerance percentage specified above by tol. Second, if this loop
    %has already gone through 10 iterations.
    while(pass == 0 && loopCount < 10)
        solveConverters(sys); %the solveConverters fuction of the pwrconverters class object sys is called. For more information read the comments in the pwrconverters class file.
        
        DSSText.Command = 'set mode = snapshot';
        DSSText.Command = 'solve';
        DSSText.Command = 'solve mode = harmonics';
        exportMonitors(MonitorList);
        
        [pass Vcurr] = checkIfCriteriaMet(MonitorList,Vprev,tol); %returns 1
        Vprev = Vcurr;
        loopCount = loopCount + 1;
    end
end
%Setting up system to solve for the Bus voltages required by the user
for i = 1:length(VSolve)
    DSSText.Command = ['New Isource.V_bus',char(VSolve(i)),'a Bus1=',char(VSolve(i)),'.1 basefreq = 60 Phases=1 amps=0'];
    DSSText.Command = ['New Isource.V_bus',char(VSolve(i)),'b Bus1=',char(VSolve(i)),'.2 basefreq = 60 Phases=1 amps=0'];
    DSSText.Command = ['New Isource.V_bus',char(VSolve(i)),'c Bus1=',char(VSolve(i)),'.3 basefreq = 60 Phases=1 amps=0'];
    DSSText.Command = ['New Monitor.V_bus',char(VSolve(i)),'a Element=Isource.V_bus',char(VSolve(i)),'a Terminal=1'];
    DSSText.Command = ['New Monitor.V_bus',char(VSolve(i)),'b Element=Isource.V_bus',char(VSolve(i)),'b Terminal=1'];
    DSSText.Command = ['New Monitor.V_bus',char(VSolve(i)),'c Element=Isource.V_bus',char(VSolve(i)),'c Terminal=1'];
end

%Setting up system to solve for the Line currents required by the user

DSSText.Command = 'Solve Mode=snapshot';  % reset everything
DSSText.Command = 'solve';
DSSText.Command = 'solve mode = harmonics';
MonitorList = DSSMonitors.AllNames;
exportMonitors(MonitorList);

getOutputs(numHarmonics,MonitorList,sys,VSolve,ISolve);

