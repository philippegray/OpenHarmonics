%INPUTS startOpenDSS Sets up interface between OpenHarmonics and OpenDSS.
%
%This function initiates the communication between OpenHarmonics and
%OpenDSS. Calling this function returns a number of objects that directly
%reference certain properties of the OpenDSS circuit file. By calling
%functions on these objects we can simulate the circuit, change circuit
%properties, add components to the circuit, change simulation parameters.
%

function [DSSObj DSSText DSSCircuit DSSSolution DSSMonitors busBasekVs] = startOpenDSS(OpenDSSFile)

DSSObj = actxserver('OpenDSSEngine.DSS');

%Start = DSSObj.Start(0);
%DSSText is a pointer to the DSSObj.Text object. DSSObj.Text is used when
%you want to write new lines to the OpenDSS file.
DSSText = DSSObj.Text;

%DSSCircuit points to the DSSObj.ActiveCircuit object. Through this object
%you can access certain data from circuit objects of interest.
DSSCircuit=DSSObj.ActiveCircuit;

%By invoking the .solve function on DSSSolution you can simulate the
%circuit for the current set of simulation parameters.
DSSSolution=DSSCircuit.Solution;

%DSSMonitors points to an array of MonitorObjects. We strickly use
%DSSMonitors to grab all the names of the Monitor Objects in the network so
%that we can locate the associated .csv files that contain the relevant
%harmonic data.
DSSMonitors = DSSCircuit.Monitors;

%This command is necessary first step before any simulations of the system
%are performed
DSSText.command = ['Compile ', OpenDSSFile,'.dss'];

busBasekVs = zeros(DSSCircuit.NumBuses,1);

%BusNames contains the names of all buses in the system. A bus is not on a
%per phase basis. The bus definition from OpenDSS is the definition used in
%OpenHarmonics in defining a bus. That is a bus is defined as a collection
%of nodes connected - and this connection line or point is the bus.
BusNames = DSSCircuit.AllBusNames; 

%This loop goes through each bus in the circuit and stores the base voltage
%in units of kV (l-l, rms) into the array busBasekVs. This base voltage is
%used later to per-unit-ize the bus voltages.
for i = 1:DSSCircuit.NumBuses
    DSSCircuit.SetActiveBus(char(BusNames(i)));
    DSSBus = DSSCircuit.ActiveBus;
    busBasekVs(i) = DSSBus.kVBase;
end